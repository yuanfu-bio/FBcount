#! /bin/bash
#SBATCH -J FBcount
#SBATCH -p defq-ssd,defq,cpu003,gpu11,gpu12,gpu13,gpu14,gpu15,gpu16,gpu17,gpu18,gpu19,gpu20,gpu21
#SBATCH -c 28
#SBATCH -n 1
#SBATCH -N 1

# input_dir="/work/xulab/xulab-seq/fy/pi-seq/01_result/AAA00052/01_raw"
# output_dir="/work/xulab/xulab-seq/fy/pi-seq/01_result/AAA00052/02_output"
# config="/work/xulab/xulab-seq/fy/pi-seq/00_submit/FBcount/configs/STD-F10.json"

eval "$(conda shell.bash hook)"
conda activate FBcount

batch=AAY00029
raw_dir="/work/xulab/xulab-seq/data/${batch}/00.mergeRawFq"
input_dir="/work/xulab/xulab-seq/fy/pi-seq/01_result/${batch}/01_raw"
output_dir="/work/xulab/xulab-seq/fy/pi-seq/01_result/${batch}/02_output"
config="/work/xulab/xulab-seq/fy/pi-seq/00_submit/FBcount/configs/STD-F10.json"

if [[ ! -s ${input_dir} ]];then
    mkdir -p ${input_dir}
    ln -s ${raw_dir}/* ${input_dir}
fi

./FBcount \
-i "${input_dir}" \
-o "${output_dir}" \
-c "${config}"