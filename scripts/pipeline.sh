#!/bin/bash

input_dir=$1
output_dir=$2
sample=$3
config=$4

source ./scripts/utils.sh

# 1. Initialise variables, switch to working directory, and create result directory

log_dir=${output_dir}/${sample}/00_logs
barcode_dir=${output_dir}/${sample}/01_barcodes
fastqs_dir=${output_dir}/${sample}/02_fastqs
counts_dir=${output_dir}/${sample}/03_counts
saturation_dir=${output_dir}/${sample}/04_saturation
raw_r1=${input_dir}/${sample}/${sample}_raw_1.fq.gz
raw_r2=${input_dir}/${sample}/${sample}_raw_2.fq.gz

if [ ! -s ${output_dir}/${sample} ]; then
    log_info "Step 1. Creat a working directory for ${sample}"
    mkdir -p ${log_dir}
    mkdir -p ${barcode_dir}
    mkdir -p ${fastqs_dir}
    mkdir -p ${counts_dir}
    mkdir -p ${saturation_dir}
else
    log_info "Step 1. Working directory has been created for ${sample}"
fi

# 2. According to the linker, use cutadapt to get the initial barcode.

# Get barcode names in config files.
log_info "Step 2. Cut barcodes with cutadapt for ${sample}"
bcs=$(jq -r '.barcode | keys[]' ${config})

for bc in $bcs; do
    if [ ! -s ${log_dir}/${sample}_${bc}.cut.summary ]; then
        # Extract the barcodes information.
        read_on_rn=$(jq -r ".barcode.${bc}[0]" ${config})
        if [ "$read_on_rn" = "r1" ]; then
            raw_r_file="$raw_r1"
        elif [ "$read_on_rn" = "r2" ]; then
            raw_r_file="$raw_r2"
        else
            log_error "Invalid sample for read_on_rn"
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
            # If no linker is specified, touch an empty file.
            touch ${barcode_dir}/${sample}_${bc}.fq.gz
        fi
    fi
done
wait


# 3. Get barcodes with positions, and correct them.

if [ ! -s ${log_dir}/${sample}_correct_attach.log ]; then
    log_info "Step 3. Get barcodes with positions, and correct for ${sample}"
    ./scripts/correct_barcodes.py \
        -f ${barcode_dir} \
        -s ${sample} \
        -r1 ${raw_r1} \
        -r2 ${raw_r2} \
        -l ${log_dir} \
        -c ${config}
else
    log_info "Step 3. Barcodes Correcting completed for ${sample}"
fi


# 4. Generate input files needed for umi counting.

if [ ! -s ${log_dir}/${sample}_gen_input_fastqs.log ]; then
    log_info "Step 4. Generate r1 and r2 needed for umi counting for ${sample}"
    ./scripts/gen_input_fastqs.py \
        -s ${sample} \
        -r1 ${raw_r1} \
        -r2 ${raw_r2} \
        -l ${log_dir} \
        -o ${fastqs_dir} \
        -c ${config}
else
    log_info "Step 4. The inputs for umi counting of ${sample} has been generated"
fi

# 5. Run UMI counting and correcting.

if [ ! -s ${counts_dir}/${sample}.log ]; then
    log_info "Step 5. Run UMI counting for ${sample}"
    ./scripts/count_UMI.py \
        -i ${fastqs_dir} \
        -o ${counts_dir} \
        -s ${sample} \
        -c ${config}
else
    log_info "Step 5. UMI has been counted for ${sample}"
fi

# 6. Run Saturation calculation

if [ ! -s ${saturation_dir}/${sample}_Downsample.tsv ]; then
    log_info "Step 6. Downsample and calculate the sequence saturation for ${sample}"
    ./scripts/calculate_saturation.py \
        -i ${counts_dir} \
        -o ${saturation_dir} \
        -s ${sample} \
        -c ${config}
else
    log_info "Step 6. Saturation has been calculated for ${sample}"
fi