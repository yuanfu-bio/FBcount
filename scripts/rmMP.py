#! /usr/bin/env python

import pandas as pd
from collections import defaultdict
import json
import argparse
import os
from utils import read_json_config, fa2df


def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="Saturation Correction.")
    parser.add_argument("-s", "--sample", required=True, help="sample name")
    parser.add_argument("-i", "--input_dir", required=True, help="Path to the input path")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the output path")
    parser.add_argument("-c", "--config", required=True, help="Path to the config json")
    args = parser.parse_args()
    return args

def rmMP(json_file, FB_info, FB_info_all, sample_name):

    with open(json_file, "r") as f:
        data = json.load(f)

    fbumi_pb_reads = defaultdict(dict)
    for pb_fb, umi_dict in data.items():
        pb, fb = pb_fb.split("_")
        for umi, umi_reads in umi_dict.items():
            fb_umi = f"{fb}_{umi}"
            fbumi_pb_reads[fb_umi][pb] = umi_reads

    df = pd.DataFrame([
        [fbumi, pb, reads]
        for fbumi, inner_dict in fbumi_pb_reads.items()
        for pb, reads in inner_dict.items()
    ], columns=['FB_UMI', 'PB', 'Reads'])

    # --- 2. 统计计算 ---
    threshold = 0.8
    summary = df.groupby('FB_UMI').agg(
        num_pbs=('PB', 'nunique'),
        total_reads=('Reads', 'sum'),
        max_reads=('Reads', 'max')
    ).reset_index()
    
    # 筛选 Multi-PB
    multi_pb_df = summary[summary['num_pbs'] > 1].copy()
    multi_pb_df['dominance_ratio'] = multi_pb_df['max_reads'] / multi_pb_df['total_reads']
    multi_pb_df['is_dominant'] = multi_pb_df['dominance_ratio'] > threshold
    
    # 标记 Top PB
    top_pb_info = df.sort_values(['FB_UMI', 'Reads'], ascending=[True, False]).drop_duplicates('FB_UMI')
    top_pb_info = top_pb_info[['FB_UMI', 'PB']].rename(columns={'PB': 'TOP_PB'})
    result_df = pd.merge(multi_pb_df, top_pb_info, on='FB_UMI', how='left')

    # --- 3. 移除 Multi-PI 逻辑 ---
    # 情况 A: 存在主导 PB
    dominant_keeps = result_df[result_df['is_dominant'] == True][['FB_UMI', 'TOP_PB', 'max_reads']]
    dominant_keeps.columns = ['FB_UMI', 'PB', 'Reads'] 
    
    # 情况 B: 单一 PB
    single_pb_keeps = summary[summary['num_pbs'] == 1][['FB_UMI', 'max_reads']]
    single_pb_keeps = pd.merge(single_pb_keeps, df, left_on=['FB_UMI', 'max_reads'], right_on=['FB_UMI', 'Reads'])
    single_pb_keeps = single_pb_keeps[['FB_UMI', 'PB', 'Reads']]
    
    df_cleaned = pd.concat([dominant_keeps, single_pb_keeps], ignore_index=True)

    # --- 4. 白名单处理 ---
    df_cleaned["FB"] = df_cleaned["FB_UMI"].apply(lambda x: x.split("_")[0])
    df_cleaned["UMI"] = df_cleaned["FB_UMI"].apply(lambda x: x.split("_")[1])
    df_cleaned_WL = df_cleaned.merge(FB_info, on="FB")

    # --- 5. Top 5 Not in WL 处理 ---
    df_not_in_WL = df_cleaned[~df_cleaned["FB"].isin(FB_info["FB"])]
    df_not_in_WL = df_not_in_WL.merge(FB_info_all, on="FB", how="left")
    FB_not_in_WL = df_not_in_WL["FB_num"].value_counts().reset_index()
    FB_not_in_WL.columns = ["FB_num", "count"]
    FB_not_in_WL = FB_not_in_WL.merge(FB_info_all, on="FB_num", how="left")

    # --- 6. 生成统计指标 DataFrame ---
    # 计算各项数值
    total_fb_umi = len(summary)
    multi_pb_cnt = len(result_df)
    dominant_cnt = result_df['is_dominant'].sum()
    final_cleaned_cnt = len(df_cleaned)
    final_wl_cnt = len(df_cleaned_WL)

    multi_pb_pct = multi_pb_cnt / total_fb_umi if total_fb_umi > 0 else 0
    dominant_in_multi_pct = dominant_cnt / multi_pb_cnt if multi_pb_cnt > 0 else 0
    clean_remove_pct = 1 - final_cleaned_cnt / total_fb_umi if total_fb_umi > 0 else 0
    wl_pct = final_wl_cnt / final_cleaned_cnt if final_cleaned_cnt > 0 else 0

    stats_data = {
        "Sample": sample_name,
        "Total_FB_UMI": total_fb_umi,
        "Multi_PB_Count": multi_pb_cnt,
        "Multi_PB_Pct": multi_pb_pct,
        "Dominant_PB_Count": dominant_cnt,
        "Dominant_in_Multi_Pct": dominant_in_multi_pct,
        "Cleaned_Count": final_cleaned_cnt,
        "Removed_Pct": clean_remove_pct,
        "Final_WL_Count": final_wl_cnt,
        "WL_Pct": wl_pct
    }
    stats_df = pd.DataFrame([stats_data])
    
    return df_cleaned, df_cleaned_WL, FB_not_in_WL, stats_df

if __name__ == "__main__":
    args = setup_and_parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    sample = args.sample
    config = read_json_config(args.config)
    print(sample)

    FB_info_file = config["feature_barcode_info"]
    FB_info_all_file = config["feature_barcode"]

    FB_info = pd.read_csv(FB_info_file, sep="\t", header=None)
    FB_info.columns = ["FB_num", "FB", "Info"]
    print(FB_info.head())

    FB_info_all = fa2df(FB_info_all_file, col_names=["FB_num", "FB"])

    json_file = os.path.join(input_dir, f"{sample}_dic_after_downsample.json")

    print(json_file)
    
    df_cleaned, df_cleaned_WL, FB_not_in_WL, stats_df = rmMP(json_file, FB_info, FB_info_all, sample)

    stats_df.to_csv(f"{output_dir}/MP_Report.tsv", index=False, sep="\t")
    FB_not_in_WL.to_csv(f"{output_dir}/FB_not_in_WL.tsv", index=False, sep="\t")

    df_cleaned.to_csv(f"{output_dir}/df_rmMP_all.tsv.gz", sep="\t", index=False, compression="gzip")
    df_cleaned_WL.to_csv(f"{output_dir}/df_rmMP_WL.tsv.gz", sep="\t", index=False, compression="gzip")