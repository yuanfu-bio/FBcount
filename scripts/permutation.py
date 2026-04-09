#! /usr/bin/env python

import os
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
import numpy as np
import pandas as pd
from numba import njit
from typing import Dict, List
import logging
import argparse

from statsmodels.stats.multitest import fdrcorrection 
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="Barcode Validation.")
    parser.add_argument("-s", "--samples", required=True, help="sample name")
    parser.add_argument("-o", "--output", required=True, help="Path to the output path")
    parser.add_argument("-c", "--config", required=True, help="Path to the config json")
    args = parser.parse_args()
    return args

def load_single_df(sample, output_dir):
    file_path = f"{output_dir}/{sample}/05_rmMP/df_rmMP_WL.tsv.gz"
    df = pd.read_csv(file_path, sep="\t", compression="gzip")
    return sample, df

@njit(cache=True)
def compute_single_permutation_swap(M_binary: np.ndarray, edges_times: int = 5, seed: int = 0) -> np.ndarray:

    if seed != 0:
        np.random.seed(seed)
        
    M_rand = M_binary.copy()
    row_indices, col_indices = np.where(M_rand == 1)
    num_edges = len(row_indices)
    
    if num_edges < 2:
        return M_rand.T @ M_rand
        
    num_swaps = num_edges * edges_times 
    
    for _ in range(num_swaps):
        idx1 = np.random.randint(0, num_edges)
        idx2 = np.random.randint(0, num_edges)
        
        if idx1 == idx2:
            continue
            
        r1, c1 = row_indices[idx1], col_indices[idx1]
        r2, c2 = row_indices[idx2], col_indices[idx2]
        
        if (r1 != r2) and (c1 != c2) and (M_rand[r1, c2] == 0) and (M_rand[r2, c1] == 0):
            M_rand[r1, c1], M_rand[r2, c2] = 0, 0
            M_rand[r1, c2], M_rand[r2, c1] = 1, 1
            col_indices[idx1], col_indices[idx2] = c2, c1

    return M_rand.T @ M_rand

class CoOccurrencePermutationTest:
    def __init__(self, 
                 n_permutations: int = 1000, 
                 edges_times: int = 5, 
                 n_jobs: int = -1):
        self.n_permutations = n_permutations
        self.edges_times = edges_times
        self.n_jobs = n_jobs

    def preprocess(self, df: pd.DataFrame, white_list: List[str], black_list: List[str]) -> pd.DataFrame:
        df_clean = df.copy()
        if white_list:
            df_clean = df_clean[df_clean["Info"].isin(white_list)]
        if black_list:
            df_clean = df_clean[~df_clean["Info"].isin(black_list)]

        df_unique = df_clean[['PB', 'Info']].drop_duplicates()
        pb_counts = df_unique['PB'].value_counts()
        valid_pbs = pb_counts[pb_counts >= 2].index
        df_filtered = df_unique[df_unique['PB'].isin(valid_pbs)]

        df_matrix = pd.crosstab(df_filtered['PB'], df_filtered['Info'])
        df_matrix = df_matrix.clip(upper=1) 
        return df_matrix

    def fit_sample(self, df: pd.DataFrame, sample_name: str, white_list: List[str], black_list: List[str]) -> Dict:
        logging.info(f"[{sample_name}] Original rows: {len(df)}")
        
        df_matrix = self.preprocess(df, white_list, black_list)
        if df_matrix.empty:
            logging.warning(f"[{sample_name}] No valid data after filtering.")
            return {}

        cols = df_matrix.columns
        n_cols = len(cols)
        M_binary_real = df_matrix.values.astype(np.float32)
        real_vals = M_binary_real.T @ M_binary_real

        base_rng = np.random.default_rng(seed=42) 
        seeds = base_rng.integers(0, 1e8, size=self.n_permutations)

        # 1. 运行置换检验
        random_matrices_list = Parallel(n_jobs=self.n_jobs)(
            delayed(compute_single_permutation_swap)(M_binary_real, self.edges_times, seeds[i])
            for i in tqdm(range(self.n_permutations), desc=f"Permuting {sample_name}", leave=False)
        )
        
        random_matrices_array = np.array(random_matrices_list)
        
        # 2. 基础统计量
        mean_rand = np.mean(random_matrices_array, axis=0)
        std_rand = np.std(random_matrices_array, axis=0)
        diff_matrix = real_vals - mean_rand
        
        with np.errstate(divide='ignore', invalid='ignore'):
            z_score_matrix = np.divide(diff_matrix, std_rand, out=np.zeros_like(diff_matrix), where=std_rand != 0)

        # ==========================================
        # 核心优化：计算经验 P 值、FDR 和 Log2FC
        # ==========================================
        
        # 3. 经验 P 值 (Empirical P-value)
        # 计算富集(互作) P值: 随机打乱中，大于等于实际观测值的比例
        p_val_enrich = np.sum(random_matrices_array >= real_vals, axis=0) / self.n_permutations
        # 为了避免 P=0 导致下游对数计算报错，通常加上伪计数 1/N
        p_val_enrich = np.maximum(p_val_enrich, 1 / self.n_permutations)

        # 计算排斥(互斥) P值: 随机打乱中，小于等于实际观测值的比例
        p_val_deplete = np.sum(random_matrices_array <= real_vals, axis=0) / self.n_permutations
        p_val_deplete = np.maximum(p_val_deplete, 1 / self.n_permutations)

        # 4. 效应量 (Log2 Fold Change)
        # 加上伪计数 1 避免 log2(0)，反映观测共现是随机期望的多少倍
        log2fc_matrix = np.log2((real_vals + 1) / (mean_rand + 1))

        # 5. FDR 多重假设检验校正 (只针对上三角矩阵计算，避免对称矩阵导致惩罚过重)
        iu = np.triu_indices(n_cols, k=1) # 提取上三角索引（不包含对角线）
        
        q_val_enrich = np.ones((n_cols, n_cols))
        q_val_deplete = np.ones((n_cols, n_cols))
        
        if len(iu[0]) > 0:
            # 对富集和排斥分别进行 FDR 校正
            _, q_upper_enrich = fdrcorrection(p_val_enrich[iu])
            _, q_upper_deplete = fdrcorrection(p_val_deplete[iu])
            
            # 还原为对称矩阵
            q_val_enrich[iu] = q_upper_enrich
            q_val_enrich.T[iu] = q_upper_enrich
            q_val_deplete[iu] = q_upper_deplete
            q_val_deplete.T[iu] = q_upper_deplete

        # 清理对角线 (自身与自身比较无意义)
        for mat in [diff_matrix, z_score_matrix, std_rand, p_val_enrich, p_val_deplete, q_val_enrich, q_val_deplete, log2fc_matrix]:
            np.fill_diagonal(mat, np.nan)

        return {
            "real": pd.DataFrame(real_vals, index=cols, columns=cols), 
            "z_score": pd.DataFrame(z_score_matrix, index=cols, columns=cols),
            "diff": pd.DataFrame(diff_matrix, index=cols, columns=cols),
            "std_rand": pd.DataFrame(std_rand, index=cols, columns=cols),
            "log2fc": pd.DataFrame(log2fc_matrix, index=cols, columns=cols),
            "p_enrich": pd.DataFrame(p_val_enrich, index=cols, columns=cols),
            "q_enrich": pd.DataFrame(q_val_enrich, index=cols, columns=cols),
            "p_deplete": pd.DataFrame(p_val_deplete, index=cols, columns=cols),
            "q_deplete": pd.DataFrame(q_val_deplete, index=cols, columns=cols),
        }

if __name__ == "__main__":
    args = setup_and_parse_args()
    output_dir = args.output
    summary_dir = os.path.join(output_dir, "00_summary")
    summary_permutation_dir = os.path.join(summary_dir, "Permutation")
    os.makedirs(summary_permutation_dir, exist_ok=True)

    samples = args.samples.split()
    samples.sort()

    # 读取移除MP白名单中的数据
    results = Parallel(n_jobs=-1)(
        delayed(load_single_df)(sample, output_dir) for sample in tqdm(samples)
    )
    dfs = dict(results)

    N_PERMUTATIONS = 1000
    white_list = []
    black_list = []

    tester = CoOccurrencePermutationTest(n_permutations=N_PERMUTATIONS, n_jobs=-1)

    results = {}
    for sample in samples:
        results[sample] = tester.fit_sample(dfs[sample], sample, white_list, black_list)
        logging.info(f"[{sample}] Done.")

    metrics = ['real', 'z_score', 'diff', 'std_rand', 'log2fc', 'p_enrich', 'q_enrich', 'p_deplete', 'q_deplete']
    for metric in metrics:
        file_name = f"{summary_permutation_dir}/Permutation_{metric}.xlsx"
        with pd.ExcelWriter(file_name) as writer:
            for sample_name, data_dict in results.items():
                df = data_dict[metric]
                sheet_name = sample_name[:31] 
                df.to_excel(writer, sheet_name=sheet_name)
