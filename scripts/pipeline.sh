#!/bin/bash

SCRIPT_DIR=$(dirname "$(realpath "$0")")
root=$(dirname "$SCRIPT_DIR")

# 0. 参数设定
input_dir=$1
output_dir=$2
sample=$3
config=$4

# echo "input_dir: $input_dir"
# echo "output_dir: $output_dir"
# echo "sample: $sample"
# echo "config: $config"

# 1. 初始化变量, 切换到工作目录, 创建结果目录
echo "start step 1."
log_dir=${output_dir}/${sample}/00_logs
barcode_dir=${output_dir}/${sample}/01_barcodes
fastqs_dir=${output_dir}/${sample}/02_fastqs
counts_dir=${output_dir}/${sample}/03_counts
raw_r1=${input_dir}/${sample}/${sample}_raw_1.fq.gz
raw_r2=${input_dir}/${sample}/${sample}_raw_2.fq.gz

mkdir -p ${log_dir}
mkdir -p ${barcode_dir}
mkdir -p ${fastqs_dir}
mkdir -p ${counts_dir}

# 0. 获取涉及到的barcode名称
bcs=$(jq -r '.barcode | keys[]' ${config})

echo "start step 2. start cut barcodes of ${sample}"

# 根据linker, 使用cutadapt获取初始barcode
for bc in $bcs; do
    if [ ! -s ${barcode_dir}/${sample}_${bc}.fq.gz ]; then
        # 提取config文件中的barcode信息
        read_on_rn=$(jq -r ".barcode.${bc}[0]" ${config})
        if [ "$read_on_rn" = "r1" ]; then
            raw_r_file="$raw_r1"
        elif [ "$read_on_rn" = "r2" ]; then
            raw_r_file="$raw_r2"
        else
            echo "Invalid sample for read_on_rn"
        fi

        start=$(jq -r ".barcode.${bc}[1]" ${config})
        end=$(jq -r ".barcode.${bc}[2]" ${config})
        bc_len=$(awk -v a="$end" -v b="$start" 'BEGIN { print a - b }')

        adapt_5=$(jq -r ".barcode.${bc}[3]" ${config})
        adapt_3=$(jq -r ".barcode.${bc}[4]" ${config})
        
        if [ -n "$adapt_5" ] && [ -n "$adapt_3" ]; then
            cutadapt \
                -a "${adapt_5}...${adapt_3}" \
                -q 20 \
                -m "${bc_len}" \
                -M "${bc_len}" \
                -j 4 \
                -o "${barcode_dir}/${sample}_${bc}.fq.gz" \
                "${raw_r_file}" \
                > "${log_dir}/${sample}_${bc}.cut.summary" 2>&1 &
        else
            touch ${barcode_dir}/${sample}_${bc}.fq.gz
        fi
    fi
done
wait

echo "correct barcodes of ${sample}"
# 批量barcode校正
if [ ! -s ${log_dir}/${sample}_correct_attach.log ]; then
    ./scripts/correct_barcodes.py \
        -f ${barcode_dir} \
        -s ${sample} \
        -r1 ${raw_r1} \
        -r2 ${raw_r2} \
        -l ${log_dir} \
        -c ${config}
fi

# 3. 生成输入格式要求的r1和r2文件
echo "start step 3, generate r1 and r2 need for umi correct."
if [ ! -s ${log_dir}/${sample}_gen_input_fastqs.log ]
then
    ./scripts/gen_input_fastqs.py \
        -s ${sample} \
        -r1 ${raw_r1} \
        -r2 ${raw_r2} \
        -l ${log_dir} \
        -o ${fastqs_dir} \
        -c ${config}
else
    echo "${sample}样本已经生成用于umi校正的输入文件, 跳过..."
fi

# 4. Run UMI correct
echo "start step 5, UMI correct."
if [ ! -s ${counts_dir}/${sample}.log ]
then
    python ./scripts/count_UMI.py \
        -i ${fastqs_dir} \
        -o ${counts_dir} \
        -s ${sample} \
        -c ${config}
else
    echo "${sample}已经结束UMI校正, 跳过..."
fi

# 4. Run Saturation calculation
echo "start step 5, UMI correct."
if [ ! -s ${counts_dir}/${sample}_Saturation.tsv ]
then
    ./scripts/calculate_saturation.py \
        -i ${counts_dir} \
        -s ${sample} \
        -c ${config}
else
    echo "${sample}已经完成饱和度计算"
fi