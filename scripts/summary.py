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
    FB_info_file = config["feature_barcode_info"]
    return barcode1, barcode2, FB_info_file

def dict2df(dict):
    df = pd.DataFrame.from_dict(dict, orient='index')
    df = df.sort_index()
    df.index.name = 'sample'
    return df

if __name__ == "__main__":
    args = setup_and_parse_args()
    config = read_json_config(args.config)
    barcode1, barcode2, FB_info_file = parse_json_config(config)
    barcodes = barcode1 + barcode2
    summary_dir = os.path.join(args.output, "00_summary")
    
    if FB_info_file == "":
        print("feature_barcode_info is not specified, continue.")
    else:
        FB_info = pd.read_csv(FB_info_file, sep = "\t", header = None)
        FB_info.columns = ["Code", "FB", "Info"]
        FB_dict = FB_info.set_index('Code')['Info'].to_dict()

    validation_dict = {}
    counts_dict = {}
    for sample in args.samples.split():
        log_dir = os.path.join(args.output, sample, "00_logs")
        saturation_dir = os.path.join(args.output, sample, "04_saturation")

        # 汇总barcode有效率
        validation_dict[sample] = {}
        for barcode in barcodes:
            barcode_info = os.path.join(log_dir, f"{sample}_{barcode}.barcode.info")
            with open(barcode_info, "r") as f:
                barcode_dict = json.load(f)
                validation_dict[sample][barcode] = barcode_dict['barcode_valid_percent']
                total_reads = int(barcode_dict['total_reads'])

        # 汇总counts
        final_map = os.path.join(saturation_dir, f"{sample}_per_bc_umi_count_after_downsample.map")
        df_map = pd.read_csv(final_map, sep="\t", header=None)
        df_map.columns = ["Barcode1", "Code", "Counts"]
        df_summary = df_map.groupby("Code")["Counts"].sum().reset_index()
        df_counts = df_summary.merge(FB_info, on="Code", how="right")
        df_counts.set_index("Info", inplace=True)
        counts_dict[sample] = df_counts["Counts"]

    df_validation = dict2df(validation_dict)
    df_counts = dict2df(counts_dict)

    file_path = os.path.join(summary_dir, f"summary.xlsx")
    with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
        df_validation.to_excel(writer, sheet_name="Validation")
        df_counts.to_excel(writer, sheet_name="Counts")