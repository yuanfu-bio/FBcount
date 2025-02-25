#!/usr/bin/env python3

import sys
import os
import json
from copy import deepcopy
import gzip
import argparse
from utils import read_json_config, read_fq, fa2dict

def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="UMI Counting.")
    parser.add_argument("-i", "--input_dir", required=True, help="input fastq.gz of r2")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the input path")
    parser.add_argument("-s", "--sample", required=True, help="sample name")
    parser.add_argument("-c", "--config", required=True, help="Path to the config json")
    args = parser.parse_args()
    return args

def parse_json_config(config):
    barcode2_ref = config["feature_barcode"]

    barcode_start = 0

    barcode1_fq_types = config["barcode_struct"]["barcode1"]
    barcode_end = 0
    per_barcode1_len = []
    for fq_type in barcode1_fq_types:
        set_info = config["barcode"][fq_type]
        barcode1_len = set_info[2] - set_info[1]
        barcode_end += barcode1_len

        per_barcode1_len.append(barcode_end)

    umi_end = umi_start = barcode_end
    for _, value in config["umi"].items():
        umi_end += value[2] - value[1]

    return barcode2_ref, barcode_start, barcode_end, umi_start, umi_end, per_barcode1_len



def get_pibc_raw_umis(r1, r2, barcode_start, barcode_end, umi_start, umi_end):
    '''step1'''
    dic_A = {}

    total_reads = 0
    for (_, seq1, _), (_, seq2, _) in zip(read_fq(r1), read_fq(r2)):
        total_reads += 1

        barcode1 = seq1[barcode_start:barcode_end]
        umi = seq1[umi_start:umi_end]
        barcode2 = seq2

        barcode = f"{barcode1}_{barcode2}"
        if barcode in dic_A:
            if umi in dic_A[barcode]:
                dic_A[barcode][umi] += 1
            else:
                dic_A[barcode][umi] = 1
        else:
            dic_A[barcode] = {}
            dic_A[barcode][umi] = 1
    return total_reads, dic_A


def is_below_hamming_threshold(str1, str2, threshold):
    distance = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            distance += 1
            if distance > threshold:
                return 0
    return 1

def correct_umi(umi_count):
    MAXDIST_CORRECT_umi = 1
    umi_correct_mapping = {}
    keys_sorted, values_sorted = zip(*sorted(umi_count.items(), key=lambda item: item[1], reverse=False))
    keys_sorted_rev, values_sorted_rev = zip(*sorted(umi_count.items(), key=lambda item: item[1], reverse=True))
    umi_raw_kinds = len(keys_sorted)
    min_count_ten_mul = values_sorted[0]*10
    
    for idx in range(umi_raw_kinds):
        raw_umi = keys_sorted_rev[idx]
        mul_low = values_sorted_rev[idx]//10
        for i in range(umi_raw_kinds):
            less_umi = keys_sorted[i]
            if less_umi not in umi_correct_mapping:
                if values_sorted[i] <= mul_low:
                    if is_below_hamming_threshold(keys_sorted[i], raw_umi, MAXDIST_CORRECT_umi):
                        umi_correct_mapping[less_umi] = raw_umi
                else:
                    break
    return umi_correct_mapping


def get_pibc_new_umis(dic_A):
    '''step2'''
    dic_B = deepcopy(dic_A)
    # 定义一个新字典, 用来记录哪些umi得到了校正
    correct_list = {}

    n = 0
    for bc, umi_counts in dic_A.items():
        n += 1
        umi_count_new = deepcopy(umi_counts)
        umi_correct_mapping = correct_umi(umi_counts)

        correct_list[bc] = umi_correct_mapping

        for key, value in umi_correct_mapping.items():
            umi_count_new[umi_correct_mapping[key]] += umi_counts[key]
            umi_count_new[key] = 0
                
        umi_count_new_bak = deepcopy(umi_count_new)
        for key, value in umi_count_new_bak.items():
            if value == 0:
                del umi_count_new[key]
        
        dic_B[bc] = umi_count_new
    return dic_B, correct_list


def get_bc_umi_counts(dic_B):
    '''step 3'''
    per_bc_umi_count = {}
    for bc, umi_counts in dic_B.items():
        per_bc_umi_count[bc] = len(umi_counts)
    return per_bc_umi_count


def export_nested_dict_to_json(nested_dict, file_name):
    with open(file_name, 'w') as file:
        json.dump(nested_dict, file, indent=4)


def write_dict_to_tsv(data_dict, filename, per_barcode1_len, barcode2_dict):
    # Sort the dictionary by value in ascending order
    sorted_items = sorted(data_dict.items(), key=lambda item: item[1])

    # Open the file for writing
    with open(filename, 'w') as file:
        for key, value in sorted_items:
            barcode1, barcode2 = key.split("_")

            # 初始化变量，用于存放每段条形码
            barcode1_parts = []
            start = 0

            # 根据 per_barcode1_len 中的索引位置切割条形码
            for end in per_barcode1_len:
                barcode1_parts.append(barcode1[start:end])
                start = end  # 更新起始索引为上一个终止位置

            # 将各段条形码使用 "+" 连接
            barcode1 = "+".join(barcode1_parts)

            # 从 barcode2_dict 获取 barcode2 的值
            barcode2 = barcode2_dict[barcode2]

            # Write key and value separated by a tab, and end with a newline
            file.write(f"{barcode1}\t{barcode2}\t{value}\n")


def output_results(per_barcode1_len, barcode2_dict, dic_A, dic_B, correct_list, per_bc_umi_count_a_correct, per_bc_umi_count_b_correct, total_reads, dir, prefix):
    '''step 4'''
    dic_A_out = os.path.join(dir, f"{prefix}_dic_A.json")
    dic_B_out = os.path.join(dir, f"{prefix}_dic_B.json")
    per_bc_umi_count_a_correct_out = os.path.join(dir, f"{prefix}_per_bc_umi_count_after_correct.map")
    per_bc_umi_count_b_correct_out = os.path.join(dir, f"{prefix}_per_bc_umi_count_before_correct.map")
    log_out = os.path.join(dir, f"{prefix}_correct_umi.log")

    export_nested_dict_to_json(dic_A, dic_A_out)
    export_nested_dict_to_json(dic_B, dic_B_out)

    correct_list_new = deepcopy(correct_list)
    for key, value in correct_list_new.items():
        if len(value) == 0:
            del correct_list[key]
    log_dict = {"total_reads": total_reads, "correct_umi_stat": correct_list}
    export_nested_dict_to_json(log_dict, log_out)

    write_dict_to_tsv(per_bc_umi_count_a_correct, per_bc_umi_count_a_correct_out, per_barcode1_len, barcode2_dict)
    write_dict_to_tsv(per_bc_umi_count_b_correct, per_bc_umi_count_b_correct_out, per_barcode1_len, barcode2_dict)

    with open(os.path.join(dir, f"{prefix}.log"), "w") as log:
        log.write("Finished!\n")


if __name__ == "__main__":
    args = setup_and_parse_args()
    
    input_dir = args.input_dir
    out_dir = args.output_dir
    sample = args.sample
    config_file = args.config

    config = read_json_config(config_file)
    (barcode2_ref, 
     barcode_start, barcode_end, umi_start, umi_end, 
     per_barcode1_len) = parse_json_config(config)

    barcode2_dict = fa2dict(barcode2_ref)

    r1 = os.path.join(input_dir, f"{sample}_r1.fq.gz")
    r2 = os.path.join(input_dir, f"{sample}_r2.fq.gz")
    total_reads, dic_A = get_pibc_raw_umis(r1, r2, barcode_start, barcode_end, umi_start, umi_end)
    
    print(f"Start correct umi for {sample}.")
    dic_B, correct_list = get_pibc_new_umis(dic_A)
    print(f"Finish correct umi for {sample}.")
    per_bc_umi_count_a_correct = get_bc_umi_counts(dic_B)
    per_bc_umi_count_b_correct = get_bc_umi_counts(dic_A)
    output_results(per_barcode1_len, barcode2_dict, dic_A, dic_B, correct_list, per_bc_umi_count_a_correct, per_bc_umi_count_b_correct, total_reads, out_dir, sample)
    print(f"Finish count umi for {sample}.")