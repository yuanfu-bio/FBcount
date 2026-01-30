#! /usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse


def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="Saturation Correction.")
    parser.add_argument("-s", "--sample", required=True, help="sample name")
    parser.add_argument("-i", "--input_dir", required=True, help="Path to the input path")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the output path")
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = setup_and_parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    sample = args.sample

    df = pd.read_csv(f"{input_dir}/{sample}_Downsample.tsv", sep='\t')
    x = df["Downsample Ratio"].values
    y = df["UMI Types"].values

    # 定义Michaelis-Menten模型
    def michaelis_menten(x, Vmax, K):
        return (Vmax * x) / (K + x)

    popt, _ = curve_fit(michaelis_menten, x, y, bounds=(0, np.inf))
    Vmax, K = popt

    # 绘图
    x_fit = np.linspace(0, 5.0, 200)
    y_fit = michaelis_menten(x_fit, *popt)

    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, label="Observed Data", color="blue")
    plt.plot(x_fit, y_fit, label=f"Michaelis-Menten Fit\n$V_{{max}}$ = {Vmax:.0f}, K = {K:.4f}", color="red")
    plt.axhline(Vmax, color='gray', linestyle='--', label=f"$V_{{max}}$ ≈ {Vmax:.0f}")
    plt.xlabel("Downsample Ratio")
    plt.ylabel("UMI Types")
    plt.title("UMI Types Saturation Prediction (Michaelis-Menten Fit)")
    plt.legend()
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{sample}_Estimation.pdf", bbox_inches="tight")

    with open(f"{output_dir}/{sample}_Estimation.tsv", "w") as f:
        f.write(f"Vmax\tK\n{Vmax}\t{K}")


