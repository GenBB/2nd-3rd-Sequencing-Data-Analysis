#!/bin/bash

base_dir="/rd2/wangwbx/project/xyz/260323/dataprocess"
fastp_res="${base_dir}/fastp_res"
multiqc_out="${base_dir}/multiqc_res"

mkdir -p "$fastp_res" "$multiqc_out"
cd "$base_dir"

for r in W1-C18-CD1 W1-C18-Ch2 W2-C18-SDS1 W1-C18-CD2 W2-C18-651 W2-C18-SDS2 W1-C18-Ch1 W2-C18-652
do
    fastq_1=${r}_1.fq.gz
    fastq_2=${r}_2.fq.gz

    fastp \
        -i "$fastq_1" -I "$fastq_2" \
        -j "${fastp_res}/${r}_fp.json" \
        -h "${fastp_res}/${r}_fp.html" \
        -w 16
done

multiqc "$fastp_res" -o "$multiqc_out"
