#! /usr/bin/env python

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import numba
from scipy import sparse
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram 
import os
import argparse
# from tqdm_joblib import tqdm_joblib
# from tqdm import tqdm

def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="Barcode Validation.")
    parser.add_argument("-s", "--samples", required=True, help="sample name")
    parser.add_argument("-o", "--output", required=True, help="Path to the output path")
    parser.add_argument("-c", "--config", required=True, help="Path to the config json")
    args = parser.parse_args()
    return args

class LargeScaleProteinAnalysis:
    def __init__(self, matrix, sample):
        
        if isinstance(matrix, pd.DataFrame):
            self.matrix = matrix.values.astype(np.int8)
            self.protein_names = matrix.columns.tolist()
        else:
            self.matrix = matrix.astype(np.int8)
            self.protein_names = [f"Protein_{i}" for i in range(matrix.shape[1])]
        
        self.n_events, self.n_proteins = self.matrix.shape
        print(f"{sample} 数据规模: {self.n_events} 次检测, {self.n_proteins} 种蛋白质")
        
        # 使用稀疏矩阵存储以节省内存
        self.sample = sample
        self.sparse_matrix = sparse.csr_matrix(self.matrix)
        self.observed_cooccurrence = None
        self.null_distributions = None
        self.z_scores = None
        self.p_values = None
    
    @staticmethod
    @numba.jit(nopython=True, parallel=True)
    def fast_cooccurrence(matrix):
        """使用numba加速的共现矩阵计算"""
        n_proteins = matrix.shape[1]
        cooccurrence = np.zeros((n_proteins, n_proteins))
        
        for i in numba.prange(n_proteins):
            for j in range(i + 1, n_proteins):
                co_occur = 0
                for k in range(matrix.shape[0]):
                    if matrix[k, i] == 1 and matrix[k, j] == 1:
                        co_occur += 1
                cooccurrence[i, j] = co_occur
                cooccurrence[j, i] = co_occur
        
        return cooccurrence
    
    def calculate_cooccurrence_fast(self, matrix):
        """快速计算共现矩阵"""
        return self.fast_cooccurrence(matrix)

    def curveball_shuffle(self, matrix, n_trades_factor=5, seed=None):

        if seed is not None:
            np.random.seed(seed)

        n_rows, n_cols = matrix.shape

        # 将每行的 1 的列索引转换为集合
        row_sets = [set(np.where(matrix[i] == 1)[0]) for i in range(n_rows)]

        n_trades = n_trades_factor * n_rows

        for _ in range(n_trades):
            # 随机选两行
            i, j = np.random.choice(n_rows, 2, replace=False)

            A = row_sets[i]
            B = row_sets[j]

            # 行太类似就跳过
            if len(A) == 0 or len(B) == 0:
                continue

            # 公共列
            common = A & B

            # 各自独有的列
            Au = list(A - common)
            Bu = list(B - common)

            if len(Au) == 0 and len(Bu) == 0:
                continue

            # 合并 pool 并随机打乱
            pool = Au + Bu
            np.random.shuffle(pool)

            # 根据原本独有列的数量重新分配
            A_new = set(pool[:len(Au)]) | common
            B_new = set(pool[len(Au):]) | common

            row_sets[i] = A_new
            row_sets[j] = B_new

        # 转回矩阵
        shuffled = np.zeros_like(matrix)
        for i, cols in enumerate(row_sets):
            shuffled[i, list(cols)] = 1

        return shuffled
        
    def parallel_shuffle(self, n_permutations=100, n_jobs=-1):
        print("开始并行 Monte Carlo 模拟...")

        self.observed_cooccurrence = self.calculate_cooccurrence_fast(self.matrix)

        def single_permutation(perm_id):
            shuffled_matrix = self.curveball_shuffle(self.matrix)
            return self.calculate_cooccurrence_fast(shuffled_matrix)

        # # 使用 tqdm_joblib 包装
        # with tqdm_joblib(tqdm(total=n_permutations)) as progress_bar:
        #     results = Parallel(n_jobs=n_jobs)(
        #         delayed(single_permutation)(i)
        #         for i in range(n_permutations)
        #     )

        results = Parallel(n_jobs=n_jobs)(
                delayed(single_permutation)(i)
                for i in range(n_permutations))

        self.null_distributions = np.array(results)
        return self.null_distributions
    
    def calculate_p_values_fast(self):
        """快速计算p值"""
        if self.null_distributions is None:
            raise ValueError("请先运行parallel_shuffle方法")
        
        p_values = np.zeros((self.n_proteins, self.n_proteins))
        z_scores = np.zeros((self.n_proteins, self.n_proteins))
        
        for i in range(self.n_proteins):
            for j in range(i + 1, self.n_proteins):
                observed = self.observed_cooccurrence[i, j]
                null_dist = self.null_distributions[:, i, j]
                
                # 使用更高效的p值计算
                p_value = (np.sum(null_dist >= observed) + 1) / (len(null_dist) + 1)
                p_values[i, j] = p_value
                p_values[j, i] = p_value
                
                # 计算z-score
                mean_null = np.mean(null_dist)
                std_null = np.std(null_dist)
                if std_null > 0:
                    z_score = (observed - mean_null) / std_null
                else:
                    z_score = 0
                z_scores[i, j] = z_score
                z_scores[j, i] = z_score
            
        self.z_scores = z_scores
        self.p_values = p_values
        
        # return p_values, z_scores
    
    def get_significant_pairs_fast(self, alpha=0.05, min_cooccurrence=10):
        """快速获取显著共现的蛋白质对，可设置最小共现阈值"""
        self.calculate_p_values_fast()
        p_values = self.p_values
        z_scores = self.z_scores
        
        significant_pairs = []
        for i in range(self.n_proteins):
            for j in range(i + 1, self.n_proteins):
                if (p_values[i, j] < alpha and 
                    self.observed_cooccurrence[i, j] >= min_cooccurrence):
                    significant_pairs.append({
                        'protein1': self.protein_names[i],
                        'protein2': self.protein_names[j],
                        'observed_cooccurrence': self.observed_cooccurrence[i, j],
                        'expected_cooccurrence': np.mean(self.null_distributions[:, i, j]),
                        'std_null': np.std(self.null_distributions[:, i, j]),
                        'p_value': p_values[i, j],
                        'z_score': z_scores[i, j],
                        'fold_change': (self.observed_cooccurrence[i, j] + 1) / 
                                     (np.mean(self.null_distributions[:, i, j]) + 1)
                    })
        
        # 按p值排序
        significant_pairs.sort(key=lambda x: x['z_score'])
        return significant_pairs

def cluster_proteins(Z):
    model = AgglomerativeClustering(
        metric="euclidean",     # 新版参数
        linkage="average",
        n_clusters=None,
        distance_threshold=0     # 输出完整层次结构
    )
    model.fit(Z)

    linkage_matrix = np.column_stack([
        model.children_,
        model.distances_,
        np.zeros(model.children_.shape[0])
    ])

    dendro = dendrogram(linkage_matrix, no_plot=True)
    order = dendro["leaves"]
    return order


if __name__ == "__main__":
    args = setup_and_parse_args()
    summary_dir = os.path.join(args.output, "00_summary")
    summary_permutation_dir = os.path.join(summary_dir, "permutation")
    os.makedirs(summary_permutation_dir, exist_ok=True)
    
    writer_pair = pd.ExcelWriter(f"{summary_permutation_dir}/Permutation_top_pairs.xlsx")
    writer_z = pd.ExcelWriter(f"{summary_permutation_dir}/Permutation_z_scores.xlsx") 
    writer_p = pd.ExcelWriter(f"{summary_permutation_dir}/Permutation_p_scores.xlsx")

    for sample in args.samples.split():
        rmMP_dir = os.path.join(args.output, sample,"05_rmMP")
        permutation_dir = os.path.join(args.output, sample,"06_permutation")
        os.makedirs(permutation_dir, exist_ok=True)

        df_rmMP_WL_file = f"{rmMP_dir}/df_rmMP_WL.tsv.gz"
        df_rmMP_WL = pd.read_csv(df_rmMP_WL_file, sep="\t", compression="gzip")
        df_matrix = df_rmMP_WL.pivot_table(index='PB', columns='Info', values='UMI', aggfunc='nunique', fill_value=0)
        df_binary = df_matrix.copy()
        df_binary[df_binary > 0] = 1
        df_binary = df_binary[df_binary.sum(axis=1) > 1]

        print(f"Binary matrix shape: {df_binary.shape}")
        random_seed = 42
        np.random.seed(random_seed)
        n_permutations=100
        max_samples = 10000
        n_jobs=21
        data = df_binary.sample(n=min(max_samples, len(df_binary)), random_state=random_seed)
        
        analyzer = LargeScaleProteinAnalysis(data, sample)
        analyzer.parallel_shuffle(n_permutations=n_permutations, n_jobs=n_jobs)

        results = analyzer.get_significant_pairs_fast(alpha=1, min_cooccurrence=10)

        # 1. 保存 pairs
        results_df = pd.DataFrame(results)
        if not results_df.empty:
            results_df.sort_values(by='z_score', ascending=False, inplace=True)
            results_df.to_excel(writer_pair, sheet_name=sample, index=False)
            results_df.to_csv(f"{permutation_dir}/{sample}_top_pairs.csv", index=False)

            # 2. 获取原始数据
            Z_raw = analyzer.z_scores.copy()
            P_raw = analyzer.p_values.copy()
            names_raw = np.array(analyzer.protein_names)

            np.fill_diagonal(Z_raw, 0)

            # 3. 计算聚类顺序 (在 log 后的数据上聚类通常效果更好)
            Z_scaled = np.sign(Z_raw) * np.log1p(np.abs(Z_raw))
            order = cluster_proteins(Z_scaled)

            # 4. 根据聚类结果重排矩阵
            Z_ordered = Z_raw[np.ix_(order, order)]
            P_ordered = P_raw[np.ix_(order, order)]
            names_ordered = names_raw[order]

            # 5. 保存重排后的矩阵到 Excel
            pd.DataFrame(Z_ordered, index=names_ordered, columns=names_ordered).to_excel(writer_z, sheet_name=sample)
            pd.DataFrame(P_ordered, index=names_ordered, columns=names_ordered).to_excel(writer_p, sheet_name=sample)
            pd.DataFrame(Z_ordered, index=names_ordered, columns=names_ordered).to_csv(f"{permutation_dir}/{sample}_z_scores.csv")
            pd.DataFrame(P_ordered, index=names_ordered, columns=names_ordered).to_csv(f"{permutation_dir}/{sample}_p_values.csv")
    writer_pair.close()
    writer_z.close()
    writer_p.close()