#!/bin/bash
set -e

# --- 1. Configuration ---
SAMPLES=("B-20-8" "BH-20-8" "SH-20-8")
REF="/rd/caiya/dm6.fa"
BED="/rd2/wangwbx/ref/dm6_CDS_merged.bed"
TMP_DIR="/rd2/wangwbx/project/260311_HI/snp/snp_tmp"
THREADS=16
JAVA_OPTS="-Djava.io.tmpdir=$TMP_DIR -Xmx20G"
BG_THRESHOLD=${#SAMPLES[@]}

# Create numbered directories
mkdir -p "$TMP_DIR" 01_qc_metrics/fastp 01_qc_metrics/mosdepth 02_processed_bam 03_vcf_results 04_background 05_final_reports

# --- 2. Core Analysis Loop ---
for SAMPLE in "${SAMPLES[@]}"; do
    echo "Processing ${SAMPLE}"

    # [01] Fastp QC
    if [ ! -f "01_qc_metrics/fastp/${SAMPLE}_clean_1.fq.gz" ]; then
        fastp -i "${SAMPLE}_1.fq.gz" -I "${SAMPLE}_2.fq.gz" -o "01_qc_metrics/fastp/${SAMPLE}_clean_1.fq.gz" -O "01_qc_metrics/fastp/${SAMPLE}_clean_2.fq.gz" -h "01_qc_metrics/fastp/${SAMPLE}.html" -j "01_qc_metrics/fastp/${SAMPLE}.json" --thread 8
    fi

    # [02] Mapping & MarkDuplicates (GATK via Conda)
    if [ ! -f "02_processed_bam/${SAMPLE}.marked.bam" ]; then
        bwa mem -t 20 -M -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:MGI" \
            "$REF" "01_qc_metrics/fastp/${SAMPLE}_clean_1.fq.gz" "01_qc_metrics/fastp/${SAMPLE}_clean_2.fq.gz" | \
        samtools view -@ 8 -bS - | samtools sort -@ 8 -o "02_processed_bam/${SAMPLE}.sorted.bam" -

        gatk --java-options "$JAVA_OPTS" MarkDuplicates \
            -I "02_processed_bam/${SAMPLE}.sorted.bam" \
            -O "02_processed_bam/${SAMPLE}.marked.bam" \
            -M "01_qc_metrics/${SAMPLE}.markdup_metrics.txt" \
            --MAX_RECORDS_IN_RAM 500000

        samtools index -@ 8 "02_processed_bam/${SAMPLE}.marked.bam"
        rm -f "02_processed_bam/${SAMPLE}.sorted.bam"
    fi

    # [03] Coverage & Statistics
    samtools stats "02_processed_bam/${SAMPLE}.marked.bam" > "01_qc_metrics/${SAMPLE}_samstats.txt"
    /rd2/wangwbx/software/mosdepth/mosdepth -t "$THREADS" -b "$BED" -T 10 "01_qc_metrics/mosdepth/${SAMPLE}" "02_processed_bam/${SAMPLE}.marked.bam"

    # [04] Mutation Detection (HaplotypeCaller)
    if [ ! -f "03_vcf_results/${SAMPLE}_highqual.vcf.gz" ]; then
        gatk --java-options "$JAVA_OPTS" HaplotypeCaller -R "$REF" -I "02_processed_bam/${SAMPLE}.marked.bam" -O "03_vcf_results/${SAMPLE}_raw.vcf" -L "$BED"

        gatk --java-options "$JAVA_OPTS" VariantFiltration -R "$REF" -V "03_vcf_results/${SAMPLE}_raw.vcf" -O "03_vcf_results/${SAMPLE}_filtered.vcf" \
            --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0" --filter-name "Filter"

        gatk --java-options "$JAVA_OPTS" SelectVariants -V "03_vcf_results/${SAMPLE}_filtered.vcf" --exclude-filtered true -O "03_vcf_results/${SAMPLE}_highqual.vcf"

        bgzip -c "03_vcf_results/${SAMPLE}_highqual.vcf" > "03_vcf_results/${SAMPLE}_highqual.vcf.gz"
        tabix -p vcf "03_vcf_results/${SAMPLE}_highqual.vcf.gz"
        rm -f 03_vcf_results/${SAMPLE}_raw.vcf 03_vcf_results/${SAMPLE}_filtered.vcf 03_vcf_results/${SAMPLE}_highqual.vcf
    fi
done

# --- 3. Background Extraction (Coding Region) ---
VCF_LIST=""
for SAMPLE in "${SAMPLES[@]}"; do VCF_LIST+="03_vcf_results/${SAMPLE}_highqual.vcf.gz "; done

echo "Generating background mutation(Real SNP in coding region)"
bcftools isec -n=${BG_THRESHOLD} -p 04_background/isec_bg $VCF_LIST
bcftools view -R 04_background/isec_bg/sites.txt $(echo $VCF_LIST | awk '{print $1}') -O v -o 04_background/common_background.vcf
bgzip -f 04_background/common_background.vcf && tabix -f -p vcf 04_background/common_background.vcf.gz
rm -rf 04_background/isec_bg

# --- 4. MultiQC ---
multiqc 01_qc_metrics/ -o 01_qc_metrics/multiqc_report

# --- 5. Final Report ---
REPORT="05_final_reports/Final_Fidelity_Report.txt"
bg_snps=$(bcftools view -H 04_background/common_background.vcf.gz | wc -l)

echo -e "Background SNPs: $bg_snps\nSample\tGlobal_Err(%)\tUnique_SNPs\tCallable_Bases\tEPM\tDeamin(%)\tTs/Tv" > "$REPORT"

for SAMPLE in "${SAMPLES[@]}"; do
    bcftools isec -p "05_final_reports/isec_${SAMPLE}" "03_vcf_results/${SAMPLE}_highqual.vcf.gz" 04_background/common_background.vcf.gz
    mv "05_final_reports/isec_${SAMPLE}/0000.vcf" "05_final_reports/${SAMPLE}_pure_errors.vcf"
    rm -rf "05_final_reports/isec_${SAMPLE}"

    mismatches=$(grep "^SN" "01_qc_metrics/${SAMPLE}_samstats.txt" | grep "mismatches:" | cut -f 3)
    mapped_bases=$(grep "^SN" "01_qc_metrics/${SAMPLE}_samstats.txt" | grep "bases mapped (cigar):" | cut -f 3)
    snps=$(bcftools view -H "05_final_reports/${SAMPLE}_pure_errors.vcf" | grep -v "INDEL" | wc -l)
    callable_bases=$(zcat "01_qc_metrics/mosdepth/${SAMPLE}.thresholds.bed.gz" | awk '!/^#/ {sum += $5} END {print sum+0}')

    global_err=$(awk -v m="$mismatches" -v b="$mapped_bases" 'BEGIN{printf "%.4f", m/b*100}')
    epm=$(awk -v s="$snps" -v b="$callable_bases" 'BEGIN{printf "%.2f", (b>0?s/b*1000000:0)}')

    stats=$(bcftools query -i 'TYPE="snp"' -f '%REF>%ALT\n' "05_final_reports/${SAMPLE}_pure_errors.vcf" | awk -v s="$snps" '
    BEGIN {ts=0; tv=0; de=0} {
        split($0,a,">"); r=a[1]; t=a[2]
        if((r=="A"&&t=="G")||(r=="G"&&t=="A")||(r=="C"&&t=="T")||(r=="T"&&t=="C")) ts++; else tv++
        if((r=="C"&&t=="T")||(r=="G"&&t=="A")) de++
    } END {printf "%.1f\t%.2f", (s>0?de/s*100:0), (tv>0?ts/tv:0)}')

    echo -e "${SAMPLE}\t${global_err}\t${snps}\t${callable_bases}\t${epm}\t${stats}" >> "$REPORT"
done

column -t "$REPORT"
