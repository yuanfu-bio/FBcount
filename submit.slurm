#! /bin/bash
#SBATCH -J FBcount
#SBATCH -p defq-ssd,defq,cpu003,gpu11,gpu12,gpu13,gpu14,gpu15,gpu16,gpu17,gpu18,gpu19,gpu20,gpu21
#SBATCH -c 28
#SBATCH -n 1
#SBATCH -N 1

input_dir="/work/xulab/xulab-seq/fy/pi-seq/01_result/AAY00021/01_raw"
output_dir="/work/xulab/xulab-seq/fy/pi-seq/01_result/AAY00021/02_output"
config="/work/xulab/xulab-seq/fy/pi-seq/00_submit/FBcount/config/STD-F15U10.json"

eval "$(conda shell.bash hook)"
conda activate FBcount

./FBcount \
-i "${input_dir}" \
-o "${output_dir}" \
-c "${config}"