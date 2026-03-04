#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from tqdm import tqdm
import matplotlib.colors as mcolors
import argparse
import os
import matplotlib.font_manager as fm
import matplotlib as mpl
from utils import custom_fonts

def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="Barcode Validation.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output path")
    args = parser.parse_args()
    return args

def plot_heatmap(z_scores, p_values, protein_names, sample_name,
                 custom_order=None, log_scale=True, p_cutoff=0.05, 
                 mask_nosig=True, cmap="coolwarm", output_file=None):
                 
    # 1. 确保数据类型为 numpy array，方便后续的矩阵高级索引操作
    Z = np.array(z_scores).copy()
    P = np.array(p_values).copy()
    names = np.array(protein_names)
    
    # 2. 对角线设为 0 (蛋白与自身的共现通常不参与可视化)
    np.fill_diagonal(Z, 0)
    
    # 3. 如果指定了自定义顺序，则进行矩阵重排
    if custom_order is not None:
        name_to_idx = {name: i for i, name in enumerate(names)}
        order = [name_to_idx[name] for name in custom_order if name in name_to_idx]
        
        Z = Z[np.ix_(order, order)]
        P = P[np.ix_(order, order)]
        names = names[order]

    # 4. Z-score 缩放 (对数转换，同时保留原有的正负号)
    if log_scale:
        Z_plot = np.sign(Z) * np.log1p(np.abs(Z))
    else:
        Z_plot = Z.copy()

    # 5. 遮罩不显著的点 (将 P >= p_cutoff 的位置强制设为 0，即显示为中心色)
    if mask_nosig:
        Z_plot[P >= p_cutoff] = 0


    num_proteins = len(names)
    
    # 设定每个格子的基础大小为 0.35 英寸
    cell_size = 0.2
    
    # 动态计算高和宽 (设置保底尺寸：高最小为8，宽最小为10)
    fig_height = max(8, num_proteins * cell_size)
    # 宽度额外加 2.0 是为了给右侧的 Colorbar 留出充足空间，防止挤压
    fig_width = max(10, num_proteins * cell_size + 2.0) 
    
    plt.figure(figsize=(fig_width, fig_height))
    
    sns.heatmap(
        Z_plot, 
        xticklabels=names, 
        yticklabels=names, 
        cmap=cmap, 
        center=0,
        square=True,
        cbar_kws={
            "shrink": 0.7,
            "aspect": 30
        }
    )
    
    plt.title(f"{sample_name} Z-score")
    plt.tight_layout()
    
    # 7. 保存文件并释放内存
    if not output_file:
        output_file = f"{sample_name}_heatmap.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":

    colors = ["#1f77b4", "#E6F9FF", "#F48FB1"]
    custom_cmap = mcolors.LinearSegmentedColormap.from_list("bluepink", colors)

    custom_fonts(default_font = "Arial")

    args = setup_and_parse_args()
    summary_dir = os.path.join(args.output, "00_summary")
    summary_permutation_dir = os.path.join(summary_dir, "Permutation")
    heatmap_dir = os.path.join(summary_permutation_dir, "Heatmaps")
    os.makedirs(heatmap_dir, exist_ok=True)

    z_file = pd.ExcelFile(f"{summary_permutation_dir}/Permutation_z_scores.xlsx")
    p_file = pd.ExcelFile(f"{summary_permutation_dir}/Permutation_p_scores.xlsx")

    my_custom_order = None 
    # my_custom_order = ["CD172a", "CD35", "CD16", "CD22", "CD32"]

    for sample in tqdm(z_file.sheet_names):
        df_z = pd.read_excel(z_file, sheet_name=sample, index_col=0)
        df_p = pd.read_excel(p_file, sheet_name=sample, index_col=0)
        
        plot_heatmap(
            z_scores=df_z.values,
            p_values=df_p.values,
            protein_names=df_z.columns.tolist(),
            sample_name=sample,
            custom_order=my_custom_order,
            log_scale=True, 
            p_cutoff=0.05, 
            mask_nosig=True, 
            cmap=custom_cmap,
            output_file=f"{heatmap_dir}/Heatmap_{sample}.pdf"
        )

        plot_heatmap(
            z_scores=df_z.values,
            p_values=df_p.values,
            protein_names=df_z.columns.tolist(),
            sample_name=sample,
            custom_order=my_custom_order,
            log_scale=True, 
            p_cutoff=1, 
            mask_nosig=True, 
            cmap=custom_cmap,
            output_file=f"{heatmap_dir}/Heatmap_all_{sample}.pdf"
        )