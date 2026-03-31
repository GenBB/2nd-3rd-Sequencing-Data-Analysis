#!/bin/bash
#SBATCH -p v6_384
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -J MDA_Process
#SBATCH -o %j.out
#SBATCH -e %j.err

DATE="20260331"
THREADS=32
WORK_DIR="/public5/home/t6s009028/project/20260331/dataprocess"
RAW_DIR="/public5/home/t6s009028/project/20260331/20260327_E260327003_U5107_FQ260301311/RawData"
OUTPUT_DIR="${WORK_DIR}/results"
REF_GENOME="/public5/home/t6s009028/workscript/dm6_1/dm6.fa"
SCRIPT_DIR="/public5/home/sch5655/wmx"
SAMPLES="PTA-S3-1 PTA-S2-1 PTA-S0-1 A-B H-B A-BH H-BH A-SH H-SH PTA-S3-2 PTA-S2-2 PTA-S0-2"

mkdir -p $OUTPUT_DIR
cd $WORK_DIR

for file in ${RAW_DIR}/*/*.fq.gz ; do
    filename=$(basename "$file")
    shortname="${filename#*L01_}"
    [ ! -e "$shortname" ] && ln -s "$file" "$shortname"
done

for f in coverage discordant splitmapping chimera total_base; do
    > ${OUTPUT_DIR}/${DATE}.${f}.count
done

for r in ${SAMPLES}; do
    fastq_1=${r}_1.fq.gz
    fastq_2=${r}_2.fq.gz
    out_prefix="${OUTPUT_DIR}/${r}"

    bwa mem -t $THREADS $REF_GENOME $fastq_1 $fastq_2 | \
    samtools view -u -F0x4 - | \
    samtools sort -@ 8 -o ${out_prefix}.pe.F4.s.bam -
    
    samtools index -@ $THREADS ${out_prefix}.pe.F4.s.bam
    samtools view -@ $THREADS -b -F0x904 ${out_prefix}.pe.F4.s.bam > ${out_prefix}.pe.F904.s.bam
    samtools index -@ $THREADS ${out_prefix}.pe.F904.s.bam
    
    indepth=${out_prefix}.pe.F904.s.bam.depth
    samtools depth -a ${out_prefix}.pe.F904.s.bam > $indepth
    python3 ${SCRIPT_DIR}/relative_depth.py $indepth

    cat $indepth | awk -v s="$r" 'BEGIN{t=0; o=0}{t++; if($3>0)o++} END{print s"\t"o"\t"t}' >> ${OUTPUT_DIR}/${DATE}.coverage.count

    samtools view -@16 ${out_prefix}.pe.F904.s.bam | awk -v s="$r" 'BEGIN{a=0; sp=0}{if($1 in all){c=1} else{if($7!="="||$9>1000||$9<-1000){sp++}; a++}} END{print s"\t"sp"\t"a}' >> ${OUTPUT_DIR}/${DATE}.discordant.count
    
    samtools view -@16 ${out_prefix}.pe.F904.s.bam | awk -v s="$r" 'BEGIN{t=0; sa=0}{if($0~/SA:Z:/){sa++}; t++} END{print s"\t"sa"\t"t}' >> ${OUTPUT_DIR}/${DATE}.splitmapping.count

    python3 ${SCRIPT_DIR}/get_chimera.Interchromosomal_Inverted_Outward_Large_Insert_Unclassified_Normal.pe.quick.py \
        ${out_prefix}.pe.F4.s.bam $out_prefix >> ${OUTPUT_DIR}/${DATE}.chimera.count &
    chimera_pid=$!

    pigz -dc $fastq_1 | awk -v s="$r" 'BEGIN{b=0}{if(NR%4==2)b+=length($0)} END{print s"\tread_1\t"b}' >> ${OUTPUT_DIR}/${DATE}.total_base.count &
    p1=$!
    pigz -dc $fastq_2 | awk -v s="$r" 'BEGIN{b=0}{if(NR%4==2)b+=length($0)} END{print s"\tread_2\t"b}' >> ${OUTPUT_DIR}/${DATE}.total_base.count &
    p2=$!
    
    wait "$p1" "$p2" "$chimera_pid"
done
