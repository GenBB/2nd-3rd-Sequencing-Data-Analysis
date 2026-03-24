#!/bin/bash

cd /rd2/wangwbx/project/xyz/260323/dataprocess

mkdir -p ./tmp_sort_dir

### 2026.3.24
RAW_DIR="/nfs119/rd2/wangwbx/rawdata/xyz/20260323"

for file in ${RAW_DIR}/*/RawData/*/*.fq.gz
do
    filename=$(basename "$file")
    shortname="${filename#*L01_}"
    ln -sf "$file" "$shortname"
    echo "Linked: $shortname"
done

for f in coverage discordant splitmapping Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal total_base; do
    > MDA_20260323.${f}.count
done

dm6_ref=/rd/caiya/hg38.fa
for r in W1-C18-CD1 W1-C18-Ch2 W2-C18-SDS1 W1-C18-CD2 W2-C18-651 W2-C18-SDS2 W1-C18-Ch1 W2-C18-652
do
    fastq_1=${r}_1.fq.gz
    fastq_2=${r}_2.fq.gz

    bwa mem -t 36 $dm6_ref $fastq_1 $fastq_2 | samtools view -u -F0x4 -@ 4 - | samtools sort -@ 4 -m 3G -T ./tmp_sort_dir/${r}_tmp -o ${r}.pe.F4.s.bam -
        
    samtools index -@ 16 ${r}.pe.F4.s.bam

    samtools view -@ 16 -bS -F0x904 ${r}.pe.F4.s.bam > ${r}.pe.F904.s.bam
    samtools index -@ 16 ${r}.pe.F904.s.bam
    
    /rd2/wangwbx/software/mosdepth/mosdepth -n -t 16 -b 100000 ${r}.F904 ${r}.pe.F904.s.bam

    out_relative="${r}.pe.F904.s.bam.100000.relative.depth"
    python3 /rd2/wangwbx/project/mosdepth_relative_depth.py ${r}.F904.regions.bed.gz $out_relative

    total_count=$(awk '$1=="total" {print $2}' ${r}.F904.mosdepth.summary.txt)
    over0_count=$(awk -v tot="$total_count" '$1=="total" && $2==1 {printf "%d", $3*tot}' ${r}.F904.mosdepth.global.dist.txt)
    echo -e "${r}\t${over0_count}\t${total_count}" >> MDA_20260323.coverage.count

    samtools view -@ 8 ${r}.pe.F904.s.bam | awk -v id="$r" 'BEGIN{all_count=0; split_count=0} {if(int($2/64)%2 == 1){all_count++; if($7 != "=" || $9 > 1000 || $9 < -1000){split_count++}}} END{print id"\t"split_count"\t"all_count}' >> MDA_20260323.discordant.count
    
    samtools view -@ 8 ${r}.pe.F904.s.bam | awk -v id="$r" 'BEGIN{total_count = 0; sa_count = 0} {total_count++; if($0 ~ /SA:Z:/){sa_count++}} END{print id"\t"sa_count"\t"total_count}' >> MDA_20260323.splitmapping.count

    samtools sort -n -@ 4 -m 3G -T ./tmp_sort_dir/${r}_ns_tmp -o ${r}.F4.name_s.bam ${r}.pe.F4.s.bam
    inbam=${r}.F4.name_s.bam
    python3 /rd2/wangwbx/project/get_chimera.Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal.pe.quick.py $inbam $r >> MDA_20260323.Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal.count &

    pigz -dc $fastq_1 | awk 'BEGIN{all_base = 0}{if(NR % 4 == 2){all_base = all_base + length($0)}}END{print "'$r'""\tread_1\t"all_base}' >> MDA_20260323.total_base.count &
    p1=$!
    pigz -dc $fastq_2 | awk 'BEGIN{all_base = 0}{if(NR % 4 == 2){all_base = all_base + length($0)}}END{print "'$r'""\tread_2\t"all_base}' >> MDA_20260323.total_base.count &
    p2=$!
    wait "$p1" "$p2"

done

wait

rm -rf ./tmp_sort_dir/
# rm *.F4.s.bam
