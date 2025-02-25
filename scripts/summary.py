#!/usr/bin/env python3

import json
import os
import pandas as pd
import argparse
from utils import read_json_config

def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="Barcode Validation.")
    parser.add_argument("-s", "--samples", required=True, help="sample name")
    parser.add_argument("-o", "--output", required=True, help="Path to the output path")
    parser.add_argument("-c", "--config", required=True, help="Path to the config json")
    args = parser.parse_args()
    return args

def parse_json_config(config):
    barcode_struct = config['barcode_struct']
    barcode1, barcode2 = barcode_struct["barcode1"], barcode_struct["barcode2"]
    return barcode1, barcode2

def dict2df(dict):
    df = pd.DataFrame.from_dict(dict, orient='index')
    df = df.sort_index()
    df.index.name = 'sample'
    return df

if __name__ == "__main__":
    args = setup_and_parse_args()
    config = read_json_config(args.config)
    barcode1, barcode2 = parse_json_config(config)
    barcodes = barcode1 + barcode2
    summary_dir = os.path.join(args.output, "00_summary")
    os.makedirs(summary_dir, exist_ok=True)

    validation_dict = {}
    df_counts_dict = {}
    df_saturation_dict = {}
    for sample in args.samples.split():
        log_dir = os.path.join(args.output, sample, "00_logs")
        count_dir = os.path.join(args.output, sample, "03_counts")
        # 汇总barcode有效率
        validation_dict[sample] = {}
        for barcode in barcodes:
            barcode_info = os.path.join(log_dir, f"{sample}_{barcode}.barcode.info")
            with open(barcode_info, "r") as f:
                barcode_dict = json.load(f)
                validation_dict[sample][barcode] = barcode_dict['barcode_valid_percent']
                total_reads = int(barcode_dict['total_reads'])
        # 汇总counts与Saturation
        counts_saturation_file = os.path.join(count_dir, f"{sample}_Saturation.tsv")
        if os.path.exists(counts_saturation_file):
            df_counts_saturation = pd.read_csv(counts_saturation_file, sep = "\t", index_col=0)
            df_counts_dict[sample] = df_counts_saturation["UMI Counts"]
            df_saturation_dict[sample] = df_counts_saturation["Saturation"]
        else:
            print(f"{sample}_Saturation.tsv dosn't exist.")

    df_validation = dict2df(validation_dict)
    if os.path.exists(counts_saturation_file):
        df_counts = dict2df(df_counts_dict)
        df_saturation = dict2df(df_saturation_dict)

    file_path = os.path.join(summary_dir, f"summary.xlsx")
    with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
        df_validation.to_excel(writer, sheet_name="Validation")
        if os.path.exists(counts_saturation_file):
            df_counts.to_excel(writer, sheet_name="Counts")
            df_saturation.to_excel(writer, sheet_name="Saturation")