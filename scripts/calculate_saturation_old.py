#!/usr/bin/env python3

import random
import json
import pandas as pd
from collections import Counter
import argparse
from utils import read_json_config

random.seed(42)

def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="Saturation Correction.")
    parser.add_argument("-s", "--sample", required=True, help="sample name")
    parser.add_argument("-i", "--input_dir", required=True, help="Path to the output path")
    parser.add_argument("-c", "--config", required=True, help="Path to the config json")
    args = parser.parse_args()
    return args

def generate_elements(bins, counts):
    n = 0
    elements = []
    rev_bins, rev_counts = bins[::-1], counts[::-1]
    for idx, bin in enumerate(rev_bins):
        count = rev_counts[idx]
        for _ in range(count):
            for _ in range(bin):
                elements.append(n)
            n += 1
    return elements

def calculate_bin_counts(elements):
    if not elements:
        return [], []
    element_counts = Counter(elements)
    value_counts = {}
    for value in element_counts.values():
        if value in value_counts:
            value_counts[value] += 1
        else:
            value_counts[value] = 1
    sorted_keys = sorted(value_counts.keys())
    return sorted_keys, [value_counts[key] for key in sorted_keys]

def compute_seq_saturation(bins, counts):
    if len(counts) == 0:
        return round(0, 2), 0, 0, round(0, 2)
    single = counts[0]
    n_duplicate_set = sum(counts)
    total, dup = 0, 0
    for i in range(0, len(bins)):
        total += bins[i] * counts[i]
        dup += (bins[i] - 1) * counts[i]
    PERCENT_DUPLICATION = (dup) * 100 / (total)
    seq_saturation = (1 - (single / n_duplicate_set)) * 100
    return round(seq_saturation, 2), single, n_duplicate_set, round(PERCENT_DUPLICATION, 2)

def downsample(bins, counts):
    # 生成元素
    elements = generate_elements(bins, counts)
    random.shuffle(elements)
    # 初始化数据收集
    total_elements = len(elements)
    sampling_results = [(0, 0, 0, 0, 0)] 
    for num in [0.0001, 0.001, 0.01, 0.1]:
        for i in range(1, 10):
            sample_ratio = round(num * i, 4)
            sample_size = int(total_elements * sample_ratio)
            sample_elements = elements[:sample_size]
            new_bins_transformed, new_counts_transformed = calculate_bin_counts(sorted(sample_elements))
            stats = compute_seq_saturation(new_bins_transformed, new_counts_transformed)
            sampling_results.append((sample_ratio, stats[0], stats[1], stats[2], stats[3]))

    full_stats = compute_seq_saturation(bins, counts)
    sampling_results.append((1.0, full_stats[0], full_stats[1], full_stats[2], full_stats[3]))
    
    df = pd.DataFrame(sampling_results, columns=['Downsample Ratio', 'Sequencing Saturation', "single", "n_duplicate_set", 'Duplication Rate'])
    return df

def calculate_bins(file, barcode):
    with open(file, "r") as f:
        dict_b = json.load(f)
    barcode_dict = {key: value for key, value in dict_b.items() if key.endswith(f"_{barcode}")}
    UMI_counts = []
    for subdict in barcode_dict.values():
        UMI_counts.extend(subdict.values())
    counts_matrix = pd.Series(UMI_counts).value_counts()
    counts_matrix = counts_matrix.sort_index()
    bins = counts_matrix.index.to_list()
    counts = counts_matrix.to_list()
    return bins, counts

if __name__ == "__main__":
    args = setup_and_parse_args()
    sample = args.sample
    input_dir = args.input_dir
    config = read_json_config(args.config)

    FB_info_file = config["feature_barcode_info"]
    if FB_info_file == "":
        print("feature_barcode_info is not specified, continue.")
    else:
        UMI_json_file = f"{input_dir}/{sample}_dic_B.json"
        FB_info = pd.read_csv(FB_info_file, sep = "\t", header = None)
        FB_info.columns = ["Num", "FB", "Info"]
        FB_dict = FB_info.set_index('FB')['Info'].to_dict()

        correct_UMI = {}
        saturation = {}
        df_dict = {}
        downsample_file_path = f"{input_dir}/{sample}_Downsample.xlsx"
        with pd.ExcelWriter(downsample_file_path, engine='xlsxwriter') as writer:
            for FB in FB_info["FB"]:
                bins, counts = calculate_bins(UMI_json_file, FB)
                df_downsample = downsample(bins, counts)
                # 存储抽样结果
                df_dict[FB_dict[FB]] = df_downsample
                # 存储测序饱和度最高时的UMI作为最终UMI
                correct_UMI[FB_dict[FB]] = df_downsample.loc[df_downsample['Sequencing Saturation'].idxmax()]["n_duplicate_set"].astype("int")
                # 存储最高测序饱和度
                saturation[FB_dict[FB]] = df_downsample['Sequencing Saturation'].max()
                df_downsample.to_excel(writer, sheet_name=str(FB_dict[FB]), index=False)

        df_result = pd.DataFrame({'UMI Counts': correct_UMI, 'Saturation': saturation})
        df_result.to_csv(f"{input_dir}/{sample}_Saturation.tsv", sep="\t")