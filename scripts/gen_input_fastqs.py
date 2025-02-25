#!/usr/bin/env python3

import os
import gzip
import json
import argparse
import pickle
from utils import read_fq, read_json_config

def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="Generate input fastqs.")
    parser.add_argument("-s", "--sample", required=True, help="sample name")
    parser.add_argument("-r1", "--raw_r1", required=True, help="raw fastq.gz of r1")
    parser.add_argument("-r2", "--raw_r2", required=True, help="raw fastq.gz of r2")
    parser.add_argument("-l", "--logs", required=True, help="Path to the log file")
    parser.add_argument("-o", "--out", required=True, help="Path to the out dir")
    parser.add_argument("-c", "--config", required=True, help="Path to the json file")
    args = parser.parse_args()
    return args

def parse_json_config(config):
    # 提取 barcode 信息
    barcode_struct = config['barcode_struct']
    barcode1, barcode2 = barcode_struct["barcode1"], barcode_struct["barcode2"]
    umi_info = config["umi"]

    # 初始化两个列表来存储长度
    per_barcode1_len = []
    per_barcode2_len = []

    # 处理 barcode1 的长度
    for barcode in barcode1:
        start = config['barcode'][barcode][1]
        end = config['barcode'][barcode][2]
        per_barcode1_len.append(end - start)

    # 处理 barcode2 的长度
    for barcode in barcode2:
        start = config['barcode'][barcode][1]
        end = config['barcode'][barcode][2]
        per_barcode2_len.append(end - start)

    return barcode1, barcode2, umi_info, per_barcode1_len, per_barcode2_len

def get_umi(umi_info, r1_seq, r2_seq, r1_qual, r2_qual):
    umi_seq = []  # 用于存储提取的 UMI 序列片段
    umi_qual = []  # 用于存储提取的 UMI 质量值片段

    # 遍历 umi_info 中的每一段 UMI 信息
    for _, info in umi_info.items():
        read_type, start, end = info  # 提取 r1/r2 信息和 start/end 位置信息

        # 根据 umi 在 r1 或 r2 上进行选择并提取 UMI 片段
        if read_type == "r1":
            umi_seq.append(r1_seq[start:end])
            umi_qual.append(r1_qual[start:end])
        elif read_type == "r2":
            umi_seq.append(r2_seq[start:end])
            umi_qual.append(r2_qual[start:end])

    # 将多个 UMI 片段的序列和质量值分别拼接起来
    seq = ''.join(umi_seq)
    qual = ''.join(umi_qual)

    return seq, qual


if __name__=="__main__":

    args = setup_and_parse_args()

    sample = args.sample
    raw_r1 = args.raw_r1
    raw_r2 = args.raw_r2
    logs = args.logs
    out = args.out
    config_file = args.config

    config = read_json_config(config_file)
    (barcode1_fq_types,
      barcode2_fq_types, umi_info, 
      per_barcode1_len, per_barcode2_len) = parse_json_config(config)

    level_qual_map = {"A": "G", "B": "F", "C": "9", "D": "8"}

    all_barcode_dict = {"barcode1": [], "barcode2": []}
    for fq_type in barcode1_fq_types:
        pkl_f = f"{os.path.join(logs, sample)}_{fq_type}.barcode.pkl"
        with open(pkl_f, 'rb') as file:
            data = pickle.load(file)
        all_barcode_dict["barcode1"].append(data)
    for fq_type in barcode2_fq_types:
        pkl_f = f"{os.path.join(logs, sample)}_{fq_type}.barcode.pkl"
        with open(pkl_f, 'rb') as file:
            data = pickle.load(file)
        all_barcode_dict["barcode2"].append(data)

    out_r1, out_r2 = os.path.join(out, f"{sample}_r1.fq.gz"), os.path.join(out, f"{sample}_r2.fq.gz")
    with gzip.open(out_r1, 'wt') as out1, gzip.open(out_r2, "wt") as out2:
        count, valid = 0, 0
        for (name, seq_r1, qual_r1), (_, seq_r2, qual_r2) in zip(read_fq(raw_r1), read_fq(raw_r2)):
            count += 1
            if count % 1000000 == 0:
                print(f"processed {count}.")
            barcode1_lst, barcode1_level_lst = [], []
            for data in all_barcode_dict["barcode1"]:
                fq_type_info = data[name]
                barcode1_lst.append(fq_type_info[0])
                barcode1_level_lst.append(fq_type_info[1])

            barcode2_lst, barcode2_level_lst = [], []
            for data in all_barcode_dict["barcode2"]:
                fq_type_info = data[name]
                barcode2_lst.append(fq_type_info[0])
                barcode2_level_lst.append(fq_type_info[1])

            if "" in barcode1_lst + barcode2_lst:
                continue
            valid += 1

            barcode1_qual_lst = [level_qual_map[level] for level in barcode1_level_lst]
            barcode1_qual_lst = [qual*per_barcode1_len[idx] for idx, qual in enumerate(barcode1_qual_lst)]
            new_barcode1 = ''.join(barcode1_lst)
            new_barcode1_qual = ''.join(barcode1_qual_lst)

            umi_seq, umi_qual = get_umi(umi_info, seq_r1, seq_r2, qual_r1, qual_r2)

            barcode2_qual_lst = [level_qual_map[level] for level in barcode2_level_lst]
            barcode2_qual_lst = [qual*per_barcode2_len[idx] for idx, qual in enumerate(barcode2_qual_lst)]
            new_barcode2 = ''.join(barcode2_lst)
            new_barcode2_qual = ''.join(barcode2_qual_lst)

            # 写出到输出文件
            out1.write(f"@{name}\n{new_barcode1 + umi_seq}\n+\n{new_barcode1_qual + umi_qual}\n")
            out2.write(f"@{name}\n{new_barcode2}\n+\n{new_barcode2_qual}\n")

    with open(os.path.join(logs, f"{sample}_gen_input_fastqs.log"), "w") as file:
        file.write(f"total reads: {count}, valid reads: {valid}.\n")