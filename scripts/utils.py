import json
import re
import gzip
from itertools import zip_longest

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
