#!/bin/bash
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH -n 50
#SBATCH -J scPTA_Process
#SBATCH -o %j.out
#SBATCH -e %j.err

DATE="20260416"
THREADS=50
WORK_DIR="/public5/home/t6s009028/project/20260416/dataprocess"
RAW_DIR="/public5/home/t6s009028/project/20260416/20260412_E260412002_U5107_FQ260400505/RawData"
OUTPUT_DIR="${WORK_DIR}/results"
TMP_SORT_DIR="${WORK_DIR}/tmp_sort_dir"
REF_GENOME="/public5/home/t6s009028/workscript/mm39/mm39.fa"
SCRIPT_DIR="/public5/home/t6s009028/workscript"
SAMPLES="scPTA-0-1 scPTA-50-1 scPTA-75-1 scPTA-0-2 scPTA-50-2 scPTA-75-2"

mkdir -p $OUTPUT_DIR
mkdir -p $TMP_SORT_DIR
cd $WORK_DIR

for file in ${RAW_DIR}/*/*.fq.gz ; do
    filename=$(basename "$file")
    shortname="${filename#*L01_}"
    [ ! -e "$shortname" ] && ln -sf "$file" "$shortname"
done

for f in coverage discordant splitmapping chimera total_base; do
    > ${OUTPUT_DIR}/${DATE}.${f}.count
done

for r in ${SAMPLES}; do
    fastq_1=${r}_1.fq.gz
    fastq_2=${r}_2.fq.gz
    out_prefix="${OUTPUT_DIR}/${r}"

    bwa mem -t $THREADS $REF_GENOME $fastq_1 $fastq_2 | \
    samtools view -u -F0x4 -@ 4 - | \
    samtools sort -@ 20 -T ${TMP_SORT_DIR}/${r}_tmp -o ${out_prefix}.pe.F4.s.bam -
    
    samtools index -@ $THREADS ${out_prefix}.pe.F4.s.bam

    samtools view -@ $THREADS -b -F0x904 ${out_prefix}.pe.F4.s.bam > ${out_prefix}.pe.F904.s.bam
    samtools index -@ $THREADS ${out_prefix}.pe.F904.s.bam
    
    mosdepth_prefix="${out_prefix}.F904"
    mosdepth -n -t $THREADS -b 100000 $mosdepth_prefix ${out_prefix}.pe.F904.s.bam

    out_relative="${out_prefix}.pe.F904.s.bam.100000.relative.depth"
    python3 ${SCRIPT_DIR}/mosdepth_relative_depth.py ${mosdepth_prefix}.regions.bed.gz $out_relative

    total_count=$(awk '$1=="total" {print $2}' ${mosdepth_prefix}.mosdepth.summary.txt)
    over0_count=$(awk -v tot="$total_count" '$1=="total" && $2==1 {printf "%d", $3*tot}' ${mosdepth_prefix}.mosdepth.global.dist.txt)
    echo -e "${r}\t${over0_count}\t${total_count}" >> ${OUTPUT_DIR}/${DATE}.coverage.count

    samtools view -@ 16 ${out_prefix}.pe.F904.s.bam | awk -v id="$r" 'BEGIN{all_count=0; split_count=0} {if(int($2/64)%2 == 1){all_count++; if($7 != "=" || $9 > 1000 || $9 < -1000){split_count++}}} END{print id"\t"split_count"\t"all_count}' >> ${OUTPUT_DIR}/${DATE}.discordant.count
    
    samtools view -@ 16 ${out_prefix}.pe.F904.s.bam | awk -v id="$r" 'BEGIN{total_count = 0; sa_count = 0} {total_count++; if($0 ~ /SA:Z:/){sa_count++}} END{print id"\t"sa_count"\t"total_count}' >> ${OUTPUT_DIR}/${DATE}.splitmapping.count

    samtools sort -n -@ 20 -T ${TMP_SORT_DIR}/${r}_ns_tmp -o ${out_prefix}.F4.name_s.bam ${out_prefix}.pe.F4.s.bam
    
    python3 ${SCRIPT_DIR}/get_chimera.Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal.pe.quick.py \
        ${out_prefix}.F4.name_s.bam $r >> ${OUTPUT_DIR}/${DATE}.chimera.count &
    chimera_pid=$!

    pigz -dc $fastq_1 | awk -v s="$r" 'BEGIN{b=0}{if(NR%4==2)b+=length($0)} END{print s"\tread_1\t"b}' >> ${OUTPUT_DIR}/${DATE}.total_base.count &
    p1=$!
    pigz -dc $fastq_2 | awk -v s="$r" 'BEGIN{b=0}{if(NR%4==2)b+=length($0)} END{print s"\tread_2\t"b}' >> ${OUTPUT_DIR}/${DATE}.total_base.count &
    p2=$!
    
    wait "$p1" "$p2" "$chimera_pid"
    rm -f "${out_prefix}.F4.name_s.bam" 2>/dev/null
done

# 清理统一的临时文件夹
wait
rm -rf $TMP_SORT_DIR
