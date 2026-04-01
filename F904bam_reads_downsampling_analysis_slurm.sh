#!/bin/bash
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -J Downsample_Process
#SBATCH -o %j.out
#SBATCH -e %j.err

DATE="20260331"
TARGET_READS=66000000
SUFFIX="norm66M"
THREADS=32
WORK_DIR="/public5/home/t6s009028/project/20260331/dataprocess"
OUTPUT_DIR="${WORK_DIR}/results"
SCRIPT_DIR="/public5/home/t6s009028/workscript"
SAMPLES="W1-C18-Ch2 W1-C18-CD2 W2-C18-SDS2 W2-C18-652"

mkdir -p ${OUTPUT_DIR}/tmp_sort_dir
cd ${OUTPUT_DIR}

for f in coverage discordant splitmapping Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal total_base; do
    > ${OUTPUT_DIR}/${DATE}.${SUFFIX}.${f}.count
done

for r in ${SAMPLES}; do
    out_prefix="${OUTPUT_DIR}/${r}"
    
    total_reads=$(samtools idxstats ${out_prefix}.pe.F904.s.bam | awk '{sum += $3} END {print sum}')

    norm_F4="${out_prefix}.pe.F4.${SUFFIX}.s.bam"
    norm_F904="${out_prefix}.pe.F904.${SUFFIX}.s.bam"

    if [ "$total_reads" -le "$TARGET_READS" ]; then
        cp ${out_prefix}.pe.F4.s.bam $norm_F4
        cp ${out_prefix}.pe.F904.s.bam $norm_F904
    else
        fraction=$(awk -v t="$TARGET_READS" -v tot="$total_reads" 'BEGIN {printf "%.4f", t/tot}')
        seed_frac=$(awk -v f="$fraction" 'BEGIN { printf "42%s", substr(f, 2) }')
        
        samtools view -@ 16 -s "$seed_frac" -b ${out_prefix}.pe.F4.s.bam > $norm_F4
        samtools view -@ 16 -bS -F0x904 $norm_F4 > $norm_F904
    fi

    samtools index -@ 16 $norm_F4
    samtools index -@ 16 $norm_F904

    mosdepth_prefix="${out_prefix}.F904.${SUFFIX}"
    mosdepth -n -t 16 -b 100000 $mosdepth_prefix $norm_F904

    out_relative="${out_prefix}.pe.F904.${SUFFIX}.s.bam.100000.relative.depth"
    python3 ${SCRIPT_DIR}/mosdepth_relative_depth.py ${mosdepth_prefix}.regions.bed.gz $out_relative

    total_count=$(awk '$1=="total" {print $2}' ${mosdepth_prefix}.mosdepth.summary.txt)
    over0_count=$(awk -v tot="$total_count" '$1=="total" && $2==1 {printf "%d", $3*tot}' ${mosdepth_prefix}.mosdepth.global.dist.txt)
    echo -e "${r}\t${over0_count}\t${total_count}" >> ${OUTPUT_DIR}/${DATE}.${SUFFIX}.coverage.count

    samtools view -@ 16 $norm_F904 | awk -v id="$r" 'BEGIN{all_count=0; split_count=0} {if(int($2/64)%2 == 1){all_count++; if($7 != "=" || $9 > 1000 || $9 < -1000){split_count++}}} END{print id"\t"split_count"\t"all_count}' >> ${OUTPUT_DIR}/${DATE}.${SUFFIX}.discordant.count
    
    samtools view -@ 16 $norm_F904 | awk -v id="$r" 'BEGIN{total_count = 0; sa_count = 0} {total_count++; if($0 ~ /SA:Z:/){sa_count++}} END{print id"\t"sa_count"\t"total_count}' >> ${OUTPUT_DIR}/${DATE}.${SUFFIX}.splitmapping.count

    samtools view -@ 16 $norm_F904 | awk -v id="$r" 'BEGIN {r1=0; r2=0} {
        if(int($2/64)%2 == 1) {
            r1 += length($10)
        } else if(int($2/128)%2 == 1) {
            r2 += length($10)
        }
    } END {
        print id"\tread_1\t"r1
        print id"\tread_2\t"r2
    }' >> ${OUTPUT_DIR}/${DATE}.${SUFFIX}.total_base.count

    samtools sort -n -@ 16 -m 5G -T ${OUTPUT_DIR}/tmp_sort_dir/${r}_ns_tmp -o ${out_prefix}.F4.${SUFFIX}.name_s.bam $norm_F4
    inbam_ns=${out_prefix}.F4.${SUFFIX}.name_s.bam
    
    python3 ${SCRIPT_DIR}/get_chimera.Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal.pe.quick.py $inbam_ns $r >> ${OUTPUT_DIR}/${DATE}.${SUFFIX}.Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal.count &

done

wait

rm -rf ${OUTPUT_DIR}/tmp_sort_dir/
