#! /usr/bin/env python

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from scipy.cluster import hierarchy
from utils import custom_fonts
import argparse

def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="Barcode Validation.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output path")
    args = parser.parse_args()
    return args

def get_cluster_order(df):
    df_for_clustering = df.fillna(0)

    Z = hierarchy.linkage(df_for_clustering.values, method='ward')
    dendro = hierarchy.dendrogram(Z, no_plot=True)

    order_indices = dendro['leaves']
    cluster_order = df.index[order_indices].tolist()

    return cluster_order


def plot_pyramid_heatmap(
        df,
        order,
        title=None,
        cmap='RdBu_r',
        vlim=None):

    num_proteins = len(df)
    cell_size = 0.2
    fig_height = max(8, num_proteins * cell_size)
    fig_width = max(10, num_proteins * cell_size + 2.0) 
    
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    df = df.reindex(index=order, columns=order)

    labels = df.columns.tolist()
    n = len(labels)
    data = df.values

    verts = []
    values = []

    for i in range(n):
        for j in range(i + 1, n):

            x = (i + j) / 2.0
            y = (j - i) / 2.0

            diamond = [
                (x, y - 0.5),
                (x + 0.5, y),
                (x, y + 0.5),
                (x - 0.5, y)
            ]

            verts.append(diamond)
            values.append(data[i, j])

    coll = PolyCollection(
        verts,
        cmap=cmap,
        edgecolors='white',
        linewidths=1
    )

    coll.set_array(np.array(values))

    if vlim is None:
        vmax = np.nanmax(np.abs(values))
        vlim = (-vmax, vmax)

    coll.set_clim(*vlim)

    ax.add_collection(coll)
    ax.autoscale()
    ax.set_aspect('equal')
    ax.axis('off')

    cbar = fig.colorbar(coll, ax=ax, shrink=0.4, aspect=20, pad=0.05)
    cbar.ax.tick_params(labelsize=8)
    cbar.set_label('Z-Score', fontsize=12)

    label_y_pos = -0.1

    label_fontsize = 8 
    for k, label in enumerate(labels):
        ax.text(
            k,
            label_y_pos,
            label,
            rotation=90,
            ha='center',
            va='top',
            fontsize=label_fontsize
        )
        
    if title:
        ax.set_title(title, fontsize=14, pad=20)

    return fig, ax, coll

if __name__ == "__main__":
    custom_fonts(default_font = "Arial")
    args = setup_and_parse_args()
    summary_dir = os.path.join(args.output, "00_summary")
    summary_permutation_dir = os.path.join(summary_dir, "Permutation")
    heatmap_dir = os.path.join(summary_permutation_dir, "Heatmaps")
    os.makedirs(heatmap_dir, exist_ok=True)
    
    zscore_dfs = pd.read_excel(f"{summary_permutation_dir}/Permutation_z_score.xlsx", sheet_name=None, index_col=0)

    df_first = zscore_dfs[list(zscore_dfs.keys())[0]]
    cluster_order = get_cluster_order(df_first)

    for i, (sample, df_zscore) in enumerate(zscore_dfs.items()):
        fig, ax, coll = plot_pyramid_heatmap(
            df_zscore,
            order=cluster_order,
            title=sample,
            cmap='RdBu_r',
        )
        
        fig.savefig(f"{heatmap_dir}/Heatmap_{sample}.pdf", dpi=300, bbox_inches='tight')
        plt.close(fig)