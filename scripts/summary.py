#! /usr/bin/env python

import json
import os
import pandas as pd
import argparse
import sys
from utils import read_json_config, draw_stacked_bar_plot
import matplotlib.pyplot as plt
import numpy as np

def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="Barcode Validation.")
    parser.add_argument("-s", "--samples", required=True, help="sample name")
    parser.add_argument("-o", "--output", required=True, help="Path to the output path")
    parser.add_argument("-c", "--config", required=True, help="Path to the config json")
    parser.add_argument("-mp", action="store_true", help="Enable multi-PI task")
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

def calc_pct(df, drop1=False, as_str=False):
    if as_str:
        df.index = df.index.map(lambda x: str(x) if x < 8 else '8+')
    if drop1:
        df = df.loc[df.index != ('1' if as_str else 1)]
    df = df.groupby(df.index).sum()
    return df.div(df.sum(axis=0), axis=1) * 100

if __name__ == "__main__":
    args = setup_and_parse_args()
    config = read_json_config(args.config)
    barcode1, barcode2, FB_info_file = parse_json_config(config)
    barcodes = barcode1 + barcode2
    summary_dir = os.path.join(args.output, "00_summary")

    # read FB info file
    if FB_info_file == "":
        print("feature_barcode_info is not specified, continue.")
        sys.exit(1)
    else:
        FB_info = pd.read_csv(FB_info_file, sep = "\t", header = None)
        FB_info.columns = ["Code", "FB", "Info"]
        FB_dict = FB_info.set_index('Code')['Info'].to_dict()

    meta_dict = {}
    counts_dict = {}
    df_counts_rmMP = pd.DataFrame()
    pb_umis = {}
    pb_fbs = {}
    df_mp_reports = []

    for sample in args.samples.split():
        log_dir = os.path.join(args.output, sample, "00_logs")
        saturation_dir = os.path.join(args.output, sample, "04_saturation")
        rmMP_dir = os.path.join(args.output, sample, "05_rmMP")

        # meta summary for barcode validation
        meta_dict[sample] = {}
        for barcode in barcodes:
            barcode_info = os.path.join(log_dir, f"{sample}_{barcode}.barcode.info")
            with open(barcode_info, "r") as f:
                barcode_dict = json.load(f)
                meta_dict[sample][barcode] = barcode_dict['barcode_valid_percent']
                total_reads = int(barcode_dict['total_reads'])
        meta_dict[sample]["Total reads"] = total_reads

        # meta summary for saturation
        downsample_file = os.path.join(saturation_dir, f"{sample}_Downsample.tsv")
        df_downsample = pd.read_csv(downsample_file, sep="\t")
        final_saturation = df_downsample["Sequencing Saturation"].max()
        final_duplication = df_downsample.loc[df_downsample["Sequencing Saturation"].idxmax(), "Duplication Ratio"]
        final_UMI_types = df_downsample.loc[df_downsample["Sequencing Saturation"].idxmax(), "UMI Types"]
        meta_dict[sample]["Saturation"] = f"{final_saturation}%"
        meta_dict[sample]["Duplication"] = f"{final_duplication}%"
        meta_dict[sample]["Observed"] = int(final_UMI_types)

        # estimation for FB counts
        estimation_file = os.path.join(saturation_dir, f"{sample}_Estimation.tsv")
        df_estimation = pd.read_csv(estimation_file, sep="\t")
        Vmax = df_estimation["Vmax"]
        meta_dict[sample]["Estimated"] = int(Vmax.iloc[0])

        if args.mp:
            MP_report = os.path.join(rmMP_dir, "MP_Report.tsv")
            df_mp_report = pd.read_csv(MP_report, sep="\t", header=0, index_col=0)
            df_mp_reports.append(df_mp_report)

        # summary for FB counts
        final_map = os.path.join(saturation_dir, f"{sample}_per_bc_umi_count_after_downsample.map")
        df_map = pd.read_csv(final_map, sep="\t", header=None)
        df_map.columns = ["Barcode1", "Code", "Counts"]
        df_summary = df_map.groupby("Code")["Counts"].sum().reset_index()
        ## merge the count matrix and FB information according to the input panel provided
        df_counts = df_summary.merge(FB_info, on="Code", how="right")
        df_counts.set_index("Info", inplace=True)
        counts_dict[sample] = df_counts["Counts"]

        if args.mp:
            # summary for FB counts afetr rmMP
            df_rmMP_WL_file = f"{rmMP_dir}/df_rmMP_WL.tsv.gz"
            df = pd.read_csv(df_rmMP_WL_file, sep="\t", compression="gzip")
            df_counts_rmMP[sample] = df["Info"].value_counts().fillna(0).astype(int)

            # calculate PB-UMI and PB-FB counts
            pb_umis[sample] = df.groupby('PB')['UMI'].nunique().value_counts()
            pb_fbs[sample] = df.groupby('PB')['FB_num'].nunique().value_counts()
        
    df_meta = dict2df(meta_dict)
    df_counts = dict2df(counts_dict)

    if args.mp:
        df_mp_summary = pd.concat(df_mp_reports, keys=args.samples.split(), names=['sample'])
        df_meta = df_meta.merge(df_mp_summary, on="sample", how="left")

        df_counts_rmMP = df_counts_rmMP.reindex(FB_info["Info"])
        df_counts_rmMP = df_counts_rmMP.T
        df_counts_rmMP.index.name = "sample"
        df_counts_rmMP = df_counts_rmMP.sort_index()

        df_pb_umis = pd.DataFrame(pb_umis).fillna(0)
        df_pb_umis = df_pb_umis[sorted(df_pb_umis.columns)]
        df_pct_umi = calc_pct(df_pb_umis.copy(), drop1=False, as_str=True)

        df_pb_fbs = pd.DataFrame(pb_fbs).fillna(0)
        df_pb_fbs = df_pb_fbs[sorted(df_pb_fbs.columns)]
        df_pct_fb = calc_pct(df_pb_fbs.copy(), drop1=False, as_str=True)
        
        draw_stacked_bar_plot(df_pct_umi, df_pct_fb, args.samples.split())
        plt.savefig(os.path.join(summary_dir, "UMI_FB_per_PB.pdf"), bbox_inches="tight")

    file_path = os.path.join(summary_dir, f"summary.xlsx")
    with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
        df_meta.to_excel(writer, sheet_name="Meta")
        df_counts.to_excel(writer, sheet_name="Counts")
        if args.mp:
            df_counts_rmMP.to_excel(writer, sheet_name="Counts_rmMP")
            df_pb_umis.to_excel(writer, sheet_name='UMI_per_PB_Count', index=True)
            df_pct_umi.to_excel(writer, sheet_name='UMI_per_PB_Pct', index=True)
            df_pb_fbs.to_excel(writer, sheet_name='FB_per_PB_Count', index=True)
            df_pct_fb.to_excel(writer, sheet_name='FB_per_PB_Pct', index=True)