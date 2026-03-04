import json
import re
import gzip
from itertools import zip_longest
from datetime import datetime
import pandas as pd
import matplotlib.font_manager as fm
import matplotlib as mpl

def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def log_info(message):
    print(f"{timestamp()} [INFO] {message}")

# 部分代码参考10x gonomics
def read_json_config(file_path):
    # with open(file_path, 'r', encoding='utf-8') as f:
    with open(file_path, 'r') as f:
        config = json.load(f)
    return config
 
def open_maybe_gzip(file_path, mode, **kwargs):
    """支持自动处理 gzip 文件并传递编码参数"""
    if file_path.endswith('.gz'): 
        return gzip.open(file_path,  mode, **kwargs)
    else:
        return open(file_path, mode, **kwargs)

def read_fq(file_path):
    """读取 fastq 文件，兼容 gzip 压缩和 utf-8 编码"""
    with open_maybe_gzip(file_path, "rt", encoding="utf-8") as file:
        lines = (line.strip()  for line in file)
        for name_line, seq, _, qual in zip_longest(*[lines]*4):
            name = re.split(r'[  /]', name_line)[0][1:]
            yield name, seq, qual 

def fa2dict(file_path):
    info = {}
    with open_maybe_gzip(file_path, "rt", encoding="utf-8") as file:
        key = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                key = line[1:]
            else:
                info[line] = key
    return info    

def fa2df(file_path, col_names):
    records = []
    with open(file_path, "r") as f:
        header = ""
        seq = ""
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    records.append({col_names[0]: header, col_names[1]: seq})
                
                header = line[1:].strip()
                seq = ""
            else:
                seq += line
                
        if header:
            records.append({col_names[0]: header, col_names[1]: seq})
            
    # 将字典列表转换为 DataFrame
    return pd.DataFrame(records)


def load_barcode_whitelist(fn, ordered=False):
    """ Barcode whitelists are text files of valid barcodes, one per line.
    Lines containing the '#' character are ignored
    """
    with open(fn, 'r') as infile:
        if ordered:
            barcodes = [line.strip() for line in infile if '#' not in line]
        else:
            barcodes = {line.strip() for line in infile if '#' not in line}
    return barcodes

def read_generator_fastq(fastq_file, paired_end=False):
    """ Returns an interator over a fastq file tha produces (name, seq, qual)
    If paired_end, returns both reads (assuming interleaving fastq)
    """
    line_index = 0
    for line in fastq_file:
        if line_index == 0:
            name1 = line.strip()[1:]
        elif line_index == 1:
            seq1 = line.strip()
        elif line_index == 3:
            qual1 = line.strip()
        elif line_index == 4:
            name2 = line.strip()[1:]
        elif line_index == 5:
            seq2 = line.strip()
        elif line_index == 7:
            qual2 = line.strip()

        line_index += 1
        if not(paired_end) and line_index == 4:
            line_index = 0

        if line_index == 8:
            line_index = 0

        if line_index == 0:
            if paired_end:
                yield (name1, seq1, qual1, name2, seq2, qual2)
            else:
                yield (name1, seq1, qual1)

def get_bc_umi_counts(dic_B):
    per_bc_umi_count = {}
    for bc, umi_counts in dic_B.items():
        per_bc_umi_count[bc] = len(umi_counts)
    return per_bc_umi_count

def write_dict_to_tsv(data_dict, filename, per_barcode1_len, barcode2_dict):
    # Sort the dictionary by value in ascending order
    sorted_items = sorted(data_dict.items(), key=lambda item: item[1], reverse=True)

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

def save_dict_to_pkl(dictionary, filepath):
    with open(filepath, 'w') as f:
        for read_name, (barcode, quality) in dictionary.items():
            barcode = barcode if barcode else ""
            f.write(f"{read_name}\t{barcode}\t{quality}\n")

def load_dict_from_pkl(filepath):
    dictionary = {}
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 3:
                raise ValueError(f"Line format error: {line}")
            read_name, barcode, quality = parts
            dictionary[read_name] = [barcode, quality]
    return dictionary

def custom_fonts(default_font = "Arial", 
                 font_dir = "/work/xulab/xulab-seq/fonts"):
    font_files = fm.findSystemFonts(fontpaths=[font_dir])
    for font_path in font_files:
        fm.fontManager.addfont(font_path)
    mpl.font_manager._load_fontmanager(try_read_cache=False)

    mpl.rcParams['font.family'] = default_font
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42

def draw_stacked_bar_plot(df_pct_umi, df_pct_fb, samples):
    import matplotlib.pyplot as plt
    import numpy as np
        
    custom_fonts(default_font = "Arial")

    color_custom = [
        "#1f77b4", "#A3D3F5", "#aec7e8", "#A2E4F8",
        "#00c19f", "#8FCC91", "#4ECDC4", "#C1E8C4",
    ]

    fig_width = max(4, len(samples) * 0.5)
    # fig_width = len(samples) * 0.5
    _, ax = plt.subplots(figsize=(fig_width, 6))

    samples = df_pct_umi.columns
    x = np.arange(len(samples))
    width = 0.35

    # UMI stacked bar
    bottom = np.zeros(len(samples))
    for col, color in zip(df_pct_umi.index, color_custom):
        ax.bar(
            x - width/2, 
            df_pct_umi.loc[col],
            width*0.9, 
            bottom=bottom, 
            label=col if col not in ax.get_legend_handles_labels()[1] else "", 
            color=color
        )
        bottom += df_pct_umi.loc[col].values

    # FB stacked bar
    bottom = np.zeros(len(samples))
    for col, color in zip(df_pct_fb.index, color_custom):
        ax.bar(
            x + width/2, 
            df_pct_fb.loc[col], 
            width*0.9, 
            bottom=bottom, 
            label=col if col not in ax.get_legend_handles_labels()[1] else "", 
            color=color,
        )
        bottom += df_pct_fb.loc[col].values

    ax.set_ylabel("Percentage (%)")
    ax.set_title("nUMIs per PB vs FB per PB")
    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right")
    ax.legend(title=None, bbox_to_anchor=(1.02, 1), loc="upper left")