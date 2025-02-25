#!/usr/bin/env python3

import os
import datetime
import multiprocessing
import pickle
import gzip
import itertools
import json
import argparse
import numpy as np
from utils import open_maybe_gzip, load_barcode_whitelist, read_generator_fastq

## 固定参数
DNA_ALPHABET = 'AGCT'
ALPHABET_MINUS = {char: {c for c in DNA_ALPHABET if c != char} for char in DNA_ALPHABET}
ALPHABET_MINUS['N'] = set(DNA_ALPHABET)
MAXDIST_CORRECT = 1
ILLUMINA_QUAL_OFFSET = 33
bc_confidence_threshold = 0.975
shiftCorrection = 2

def setup_and_parse_args():
    parser = argparse.ArgumentParser(description="Correct barcodes.")
    parser.add_argument("-f", "--fq_dir", required=True, help="file dir to _bc1.fq")
    parser.add_argument("-s", "--sample", required=True, help="sample name")
    parser.add_argument("-r1", "--raw_r1", required=True, help="raw fastq.gz of r1")
    parser.add_argument("-r2", "--raw_r2", required=True, help="raw fastq.gz of r2")
    parser.add_argument("-l", "--logs", required=True, help="Path to the log file")
    parser.add_argument("-c", "--config", required=True, help="Path to the config json")
    args = parser.parse_args()
    return args

def read_json_config(file_path):
    with open(file_path, 'r', encoding='utf-8') as f:
        config = json.load(f)
    return config

def parse_json_config(config, raw_r1, raw_r2):
    # 提取 barcode 信息
    barcode_struct = config['barcode_struct']
    barcode1, barcode2 = barcode_struct["barcode1"], barcode_struct["barcode2"]
    
    whitelists = []
    barcode_types = []
    starts = []
    ends = []
    raw_fq = []

    for barcode_type in barcode1 + barcode2:
        # 提取各个值
        # barcode_type = barcode
        config_barcode = config["barcode"][barcode_type]
        r_type = config_barcode[0]
        start = config_barcode[1]
        end = config_barcode[2]
        whitelist = config_barcode[5]

        # 填充到相应列表
        barcode_types.append(barcode_type)
        starts.append(start)
        ends.append(end)
        whitelists.append(whitelist)

        # 根据 r1 和 r2 判断 raw_fq
        if r_type == "r1":
            raw_fq.append(raw_r1)
        elif r_type == "r2":
            raw_fq.append(raw_r2)

    return whitelists, barcode_types, starts, ends, raw_fq

def save_dict_to_pkl(dictionary, filepath):
    with open(filepath, 'wb') as file:
        pickle.dump(dictionary, file)

def load_dict_from_pkl(filepath):
    with open(filepath, 'rb') as file:
        dictionary = pickle.load(file)
    return dictionary

def read_fastq_to_dict(filepath):
    fastq_dict = {}
    with open_maybe_gzip(filepath, 'rb') as file:
        while True:
            name = file.readline().strip()
            if not name:
                break
            seq = file.readline().strip()
            file.readline()  # Skip the '+' line
            qual = file.readline().strip()

            # 解码字节对象为字符串
            name = name.decode('utf-8')
            seq = seq.decode('utf-8')
            qual = qual.decode('utf-8')
            
            read_name = name[1:].split("/")[0].split(" ")[0]  # Remove the initial '@' from the read name
            fastq_dict[read_name] = (seq, qual)
    return fastq_dict

def get_bc_counts(fastq_dict, wl_idxs):
    bc_counts = [0] * len(wl_idxs)
    for _, (seq, _) in fastq_dict.items():
        idx = wl_idxs.get(seq)
        if idx is not None:
            bc_counts[idx] += 1
    return bc_counts       

def correct_barcode(bc_confidence_threshold, seq, qual, wl_idxs, wl_dist, maxdist=3):
    """Estimate the correct barcode given an input sequence, base quality scores, a barcode whitelist, and a prior
    distribution of barcodes.  Returns the corrected barcode if the posterior likelihood is above the confidence
    threshold, otherwise None.  Only considers corrected sequences out to a maximum Hamming distance of 3.
    """
    wl_cand = []
    likelihoods = []

    # qvs = np.fromstring(qual, dtype=np.byte) - ILLUMINA_QUAL_OFFSET
    qvs = np.frombuffer(qual.encode(), dtype=np.byte) - ILLUMINA_QUAL_OFFSET
    # Restrict QVs within a narrow range to prevent one base error from overly dominating others
    qvs[qvs < 3.0] = 3.0
    qvs[qvs > 40.0] = 40.0

    if seq in wl_idxs:
        if (qvs > 24).all():
            return seq, "uncorrected"

        wl_cand.append(seq)
        likelihoods.append(wl_dist[wl_idxs[seq]])

    for test_str, error_probs in gen_nearby_seqs(seq, qvs, wl_idxs, maxdist):
        idx = wl_idxs.get(test_str)
        p_bc = wl_dist[idx]
        log10p_edit = error_probs / 10.0
        likelihoods.append(p_bc * (10 ** -log10p_edit))
        wl_cand.append(test_str)

    posterior = np.array(likelihoods)
    posterior /= posterior.sum()

    if len(posterior) > 0:
        pmax = posterior.max()
        if pmax > bc_confidence_threshold:
            return wl_cand[np.argmax(posterior)], "corrected"
    return None, "failed"

def gen_nearby_seqs(seq, qvs, wl_idxs, maxdist=3):
    """Generate all sequences with at most maxdist changes from seq that are in a provided whitelist, along with the
    quality values of the bases at the changed positions.
    """
    allowed_indices = [i for i in range(len(seq)) if seq[i] != 'N']
    required_indices = tuple([i for i in range(len(seq)) if seq[i] == 'N'])
    mindist = len(required_indices)
    if mindist > maxdist:
        return

    for dist in range(mindist + 1, maxdist + 1):
        for modified_indices in itertools.combinations(allowed_indices, dist - mindist):
            indices = set(modified_indices + required_indices)
            error_probs = qvs[np.array(list(indices))]
            for substitutions in itertools.product(
                    *[ALPHABET_MINUS[base] if i in indices else base for i, base in enumerate(seq)]):
                new_seq = ''.join(substitutions)
                if new_seq in wl_idxs:
                    yield new_seq, error_probs.sum()

def get_barcodes_from_pos(seq, qual, seq_start, seq_end, shiftCorrection):
    lst_bc, lst_qual = [], []
    while 1:
        lst_bc.append(seq[seq_start:seq_end])
        lst_qual.append(qual[seq_start:seq_end])
        if seq_start == 0 or shiftCorrection == 0:
            break
        else:
            seq_start -= 1
            seq_end -= 1
            shiftCorrection -= 1
    return lst_bc, lst_qual
    

def correct_barcode_file(raw_fq_gz, seq_start, seq_end, wl_idxs, bc_dist, bc_confidence_threshold,
                         MAXDIST_CORRECT, fq_dict, shiftCorrection):
    log_dict = {"total_reads": 0, 
                "linker_right":{"bc_right_count": {"without_correct": 0, "need_correct": 0}, "bc_wrong_count": 0},
                "linker_wrong":{"bc_right_count": {}, "bc_wrong_count": 0}}
    read_name_bc_dict = {}
    for i in range(shiftCorrection + 1):
        log_dict["linker_wrong"]["bc_right_count"][f"shift_{i}_need_correct"] = 0 # 按照shift设定值给定日志字典的键值对
        log_dict["linker_wrong"]["bc_right_count"][f"shift_{i}_without_correct"] = 0 # 按照shift设定值给定日志字典的键值对

    fq = open_maybe_gzip(raw_fq_gz, 'rb')
    for read_group in read_generator_fastq(fq):
        level = 0 # 初始化校正等级为0
        log_dict["total_reads"] += 1
        if log_dict["total_reads"] % 500000 == 0:
            print(f"finish {log_dict['total_reads']}")
        name, seq, qual = read_group[0].decode("ASCII"), read_group[1].decode("ASCII"), read_group[2].decode("ASCII")
        name = name.split("/")[0].split(" ")[0]

        # 判断name是否在fq_dict, 若在, 则校正的过程不需要移位
        if name in fq_dict:
            barcode = fq_dict[name][0]
            barcode_qual = fq_dict[name][1]
            corrected_bc, correct_flag = correct_barcode(bc_confidence_threshold, 
                                                         barcode, barcode_qual, wl_idxs, bc_dist, MAXDIST_CORRECT)
            if not corrected_bc:
                log_dict["linker_right"]["bc_wrong_count"] += 1
                level = "E"
            elif correct_flag == "uncorrected":
                log_dict["linker_right"]["bc_right_count"]["without_correct"] += 1
                level = "A"
            else:
                level = "B"
                log_dict["linker_right"]["bc_right_count"]["need_correct"] += 1
        else:
            # 若不在, 则按位置取barcode并进行校正, 如果可能的话, 进行移位操作
            barcode_lst, barcode_qual_lst = get_barcodes_from_pos(seq, qual, 
                                                                  seq_start, seq_end, shiftCorrection)
            for idx, barcode in enumerate(barcode_lst):
                barcode_qual = barcode_qual_lst[idx]
                corrected_bc, correct_flag = correct_barcode(bc_confidence_threshold, 
                                                             barcode, barcode_qual, wl_idxs, bc_dist, MAXDIST_CORRECT)
                if corrected_bc:
                    if correct_flag == "uncorrected":
                        log_dict["linker_wrong"]["bc_right_count"][f"shift_{idx}_without_correct"] += 1
                        level = "C" if idx else "A"
                    else:
                        log_dict["linker_wrong"]["bc_right_count"][f"shift_{idx}_need_correct"] += 1
                        level = "D" if idx else "B"
                    break
            if not corrected_bc:
                log_dict["linker_wrong"]["bc_wrong_count"] += 1
                level = "E"
                
        # 执行校正后, 若corrected_bc有内容, 则校正成功/不需要校正, 若无内容, 则校正失败
        if corrected_bc:
            read_name_bc_dict[name] = [corrected_bc, level]
        else:
            read_name_bc_dict[name] = ['', level]
            
    fq.close()
    return log_dict, read_name_bc_dict

def write_nested_dict_to_json(file_path, data):
    with open(file_path, 'w', encoding='utf-8') as file:
        json.dump(data, file, indent=4, ensure_ascii=False)    

def main(whitelist, fq, raw_fq_gz, seq_start, seq_end, barcode_type, logs, sample,
         bc_confidence_threshold, MAXDIST_CORRECT, shiftCorrection):
    bc = load_barcode_whitelist(whitelist)
    wl_idxs = {bc: idx for (idx, bc) in enumerate(sorted(list(bc)))}
    fq_dict = read_fastq_to_dict(fq)
    bc_counts  = get_bc_counts(fq_dict, wl_idxs)

    bc_dist = np.array(bc_counts, dtype=float) + 1.0
    bc_dist = bc_dist / bc_dist.sum()

    log_dict, read_name_bc_dict = correct_barcode_file(raw_fq_gz, seq_start, seq_end, 
                                                       wl_idxs, bc_dist, bc_confidence_threshold,
                                                       MAXDIST_CORRECT, fq_dict, shiftCorrection)
    total_wrong_percent = (
        log_dict["linker_right"]["bc_wrong_count"] + log_dict["linker_wrong"]["bc_wrong_count"])*100/log_dict["total_reads"]
    log_dict["barcode_valid_percent"] = f"{(100-total_wrong_percent):.2f}%"
    log_json = os.path.join(logs, f"{sample}_{barcode_type}.barcode.info")
    write_nested_dict_to_json(log_json, log_dict)
    read_name_bc_dict_pkl = os.path.join(logs, f"{sample}_{barcode_type}.barcode.pkl")
    save_dict_to_pkl(read_name_bc_dict, read_name_bc_dict_pkl)
    return 1

def worker(args):
    result = main(*args)
    return result


if __name__=="__main__":
    args = setup_and_parse_args()

    fq_dir = args.fq_dir
    sample = args.sample
    raw_r1 = args.raw_r1
    raw_r2 = args.raw_r2
    logs = args.logs
    config_file = args.config

    config = read_json_config(config_file)
    (whitelists, barcode_types, starts, ends, raw_fq) = parse_json_config(config, raw_r1, raw_r2)

    # Required arguments
    fq = []
    for barcode_type in barcode_types:
        clipped_barcode_file = f"{os.path.join(fq_dir, sample)}_{barcode_type}.fq.gz"
        if not os.path.exists(clipped_barcode_file):
            with gzip.open(clipped_barcode_file, 'w') as f:
                pass
        fq.append(clipped_barcode_file)

    # 开始运行
    params = [(whitelists[i], fq[i], raw_fq[i], starts[i], ends[i], barcode_types[i],
               logs, sample, bc_confidence_threshold, MAXDIST_CORRECT, shiftCorrection) for i in range(len(whitelists))]
    print(params)
    
    multiprocessing.set_start_method('fork')
    with multiprocessing.Pool(processes=len(whitelists)) as pool:
        read_name_barcode_dict_lst = pool.map(worker, params)
    
    now = datetime.datetime.now()
    formatted_time = now.strftime("%H:%M:%S")
    print("Finish barcode correct. Current time:", formatted_time)

    with open(os.path.join(logs, f"{sample}_correct_attach.log"), "w") as file:
        file.write('Finished.\n')
