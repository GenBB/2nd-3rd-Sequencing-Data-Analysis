#!/bin/bash

DATE="20260324"
TARGET_READS=66000000
SUFFIX="norm66M"
WORK_DIR="/rd2/wangwbx/project/xyz/260323/dataprocess"
SAMPLES="W1-C18-Ch2 W1-C18-CD2 W2-C18-SDS2 W2-C18-652"

cd ${WORK_DIR}

mkdir -p ./tmp_sort_dir

for f in coverage discordant splitmapping Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal total_base; do
    > ${DATE}.${SUFFIX}.${f}.count
done

for r in ${SAMPLES}; do
    total_reads=$(samtools idxstats ${r}.pe.F904.s.bam | awk '{sum += $3} END {print sum}')

    norm_F904="${r}.pe.F904.${SUFFIX}.s.bam"

    if [ "$total_reads" -le "$TARGET_READS" ]; then
        cp ${r}.pe.F904.s.bam $norm_F904
    else
        fraction=$(awk -v t="$TARGET_READS" -v tot="$total_reads" 'BEGIN {printf "%.4f", t/tot}')
        seed_frac=$(awk -v f="$fraction" 'BEGIN { printf "42%s", substr(f, 2) }')
        
        samtools view -@ 16 -s "$seed_frac" -b ${r}.pe.F904.s.bam > $norm_F904
    fi

    samtools index -@ 16 $norm_F904

    mosdepth_prefix="${r}.F904.${SUFFIX}"
    /rd2/wangwbx/software/mosdepth/mosdepth -n -t 16 -b 100000 $mosdepth_prefix $norm_F904

    out_relative="${r}.pe.F904.${SUFFIX}.s.bam.100000.relative.depth"
    python3 /rd2/wangwbx/project/mosdepth_relative_depth.py ${mosdepth_prefix}.regions.bed.gz $out_relative

    total_count=$(awk '$1=="total" {print $2}' ${mosdepth_prefix}.mosdepth.summary.txt)
    over0_count=$(awk -v tot="$total_count" '$1=="total" && $2==1 {printf "%d", $3*tot}' ${mosdepth_prefix}.mosdepth.global.dist.txt)
    echo -e "${r}\t${over0_count}\t${total_count}" >> ${DATE}.${SUFFIX}.coverage.count

    samtools view -@ 8 $norm_F904 | awk -v id="$r" 'BEGIN{all_count=0; split_count=0} {if(int($2/64)%2 == 1){all_count++; if($7 != "=" || $9 > 1000 || $9 < -1000){split_count++}}} END{print id"\t"split_count"\t"all_count}' >> ${DATE}.${SUFFIX}.discordant.count
    
    samtools view -@ 8 $norm_F904 | awk -v id="$r" 'BEGIN{total_count = 0; sa_count = 0} {total_count++; if($0 ~ /SA:Z:/){sa_count++}} END{print id"\t"sa_count"\t"total_count}' >> ${DATE}.${SUFFIX}.splitmapping.count

    samtools view -@ 8 $norm_F904 | awk -v id="$r" 'BEGIN {r1=0; r2=0} {
        if(int($2/64)%2 == 1) {
            r1 += length($10)
        } else if(int($2/128)%2 == 1) {
            r2 += length($10)
        }
    } END {
        print id"\tread_1\t"r1
        print id"\tread_2\t"r2
    }' >> ${DATE}.${SUFFIX}.total_base.count

    samtools sort -n -@ 8 -m 5G -T ./tmp_sort_dir/${r}_ns_tmp -o ${r}.F904.${SUFFIX}.name_s.bam $norm_F904
    inbam_ns=${r}.F904.${SUFFIX}.name_s.bam
    
    python3 /rd2/wangwbx/project/get_chimera.Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal.pe.quick.py $inbam_ns $r >> ${DATE}.${SUFFIX}.Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal.count &

done

wait

rm -rf ./tmp_sort_dir/
