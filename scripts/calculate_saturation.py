#!/usr/bin/env python3

import random
import json
import pandas as pd
from collections import Counter
import argparse
from utils import read_json_config, fa2dict, get_bc_umi_counts, write_dict_to_tsv

random.seed(42)

def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="Saturation Correction.")
    parser.add_argument("-s", "--sample", required=True, help="sample name")
    parser.add_argument("-i", "--input_dir", required=True, help="Path to the input path")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the output path")
    parser.add_argument("-c", "--config", required=True, help="Path to the config json")
    args = parser.parse_args()
    return args

def parse_json_config(config):
    barcode2_ref = config["feature_barcode"]

    barcode1_fq_types = config["barcode_struct"]["barcode1"]
    barcode_end = 0
    per_barcode1_len = []
    for fq_type in barcode1_fq_types:
        set_info = config["barcode"][fq_type]
        barcode1_len = set_info[2] - set_info[1]
        barcode_end += barcode1_len
        per_barcode1_len.append(barcode_end)

    return barcode2_ref, per_barcode1_len

def gen_sample_pool(dct):
    # Create a sample pool by expanding each barcode-umi according to its occurrence count
    sample_pool = []
    for barcode, umis in dct.items():
        for umi, count in umis.items():
            sample_pool.extend([f"{barcode}-{umi}"] * count)
    return sample_pool

def compute_seq_saturation(barcode_dict):
    # Compute sequencing saturation and duplication ratio
    UMI_counts = []
    for subdict in barcode_dict.values():
        UMI_counts.extend(subdict.values())
    counts_matrix = pd.Series(UMI_counts).value_counts()
    counts_matrix = counts_matrix.sort_index()
    bins = counts_matrix.index.to_list()
    counts = counts_matrix.to_list()

    if len(counts) == 0:
        return round(0, 2), 0, 0, round(0, 2)
    single = counts[0]
    n_duplicate_set = sum(counts)
    total, dup = 0, 0
    for i in range(0, len(bins)):
        total += bins[i] * counts[i]
        dup += (bins[i] - 1) * counts[i]
    duplication_ratio = (dup) * 100 / (total)
    seq_saturation = (1 - (single / n_duplicate_set)) * 100
    return round(seq_saturation, 2), single, n_duplicate_set, round(duplication_ratio, 2)

def downsample(dct, downsample_ratio):
    # Downsampling functions
    sample_pool = gen_sample_pool(dct)
    sample_size = int(len(sample_pool) * downsample_ratio)
    sampled_items = random.sample(sample_pool, sample_size)
    # Create a new dictionary to store the sampled results
    sampled_dict = {}
    # Group the sampled items by barcode and umi
    for item in sampled_items:
        barcode, umi = item.split('-')
        if barcode not in sampled_dict:
            sampled_dict[barcode] = {}
        if umi not in sampled_dict[barcode]:
            sampled_dict[barcode][umi] = 0
        sampled_dict[barcode][umi] += 1
    return sampled_dict

if __name__ == "__main__":
    args = setup_and_parse_args()
    sample = args.sample
    input_dir = args.input_dir
    output_dir = args.output_dir
    config = read_json_config(args.config)

    barcode2_ref, per_barcode1_len = parse_json_config(config)

    file = f"{input_dir}/{sample}_dic_B.json"
    with open(file, "r") as f:
            dict_b = json.load(f)

    # Calculate saturation after downsampling with different ratio
    sampling_results = [(0, 0, 0, 0, 0)] 
    for num in [0.0001, 0.001, 0.01, 0.1]:
        for i in range(1, 10):
            downsample_ratio = round(num * i, 4)
            downsample_result = downsample(dict_b, downsample_ratio)
            stats = compute_seq_saturation(downsample_result)
            sampling_results.append((downsample_ratio, stats[0], stats[1], stats[2], stats[3]))

    # Calculate saturation without downsampling
    stats = compute_seq_saturation(dict_b)
    sampling_results.append((1, stats[0], stats[1], stats[2], stats[3]))
    df = pd.DataFrame(sampling_results, columns=['Downsample Ratio', 'Sequencing Saturation', "UMI Types", "UMI Counts", 'Duplication Ratio'])
    df.to_csv(f"{output_dir}/{sample}_Downsample.tsv", sep="\t")

    # Get the optimal ratio of sampling results and generate the final result
    optimal_ratio = df.loc[df['Sequencing Saturation'].idxmax(), 'Downsample Ratio']
    downsample_result = downsample(dict_b, optimal_ratio)
    downsample_result_path = f"{output_dir}/{sample}_dic_after_downsample.json"
    with open(downsample_result_path, 'w') as file:
            json.dump(downsample_result, file, indent=4)

    barcode2_dict = fa2dict(barcode2_ref)
    per_bc_umi_count_after_downsample = get_bc_umi_counts(downsample_result)
    write_dict_to_tsv(per_bc_umi_count_after_downsample, f"{output_dir}/{sample}_per_bc_umi_count_after_downsample.map", per_barcode1_len, barcode2_dict)