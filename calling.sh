#!/bin/bash

# ==== [1] Input parsing ====
INPUT_BAM="$1"
gatk="java -jar /usr/bin/gatk.jar"
if [ -z "$INPUT_BAM" ]; then
  echo "Usage: $0 <input_bam>"
  exit 1
fi

SAMPLE=$(basename "$INPUT_BAM" _single.bam)
REF="/disk1/250616sra2/fastq/Homo_sapiens_assembly38.fasta"
DBSNP="/disk1/oneomics_analysis/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz"
INDELS="/disk1/oneomics_analysis/references/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

# ==== [2] MarkDuplicates ====
#${gatk} MarkDuplicates \
#    -I "$INPUT_BAM" \
#    -O "${SAMPLE}.markdup.bam" \
#    -M "${SAMPLE}.metrics.txt"

# ==== [3] Sort & Index ====
#samtools sort -o "${SAMPLE}.sorted.bam" "${SAMPLE}.markdup.bam"
#samtools index "${SAMPLE}.sorted.bam"

# ==== [4] BaseRecalibrator ====
#${gatk} BaseRecalibrator \
#    -R "$REF" \
#    -I "${SAMPLE}.sorted.bam" \
#    --known-sites "$DBSNP" \
#    --known-sites "$INDELS" \
#    -O "${SAMPLE}.recal.table"

# ==== [5] ApplyBQSR ====
#${gatk} ApplyBQSR \
#    -R "$REF" \
#    -I "${SAMPLE}.sorted.bam" \
#    --bqsr-recal-file "${SAMPLE}.recal.table" \
#    -O "${SAMPLE}.recal.bam"

# ==== [6] HaplotypeCaller ====
${gatk} HaplotypeCaller \
    -R "$REF" \
    -I "${INPUT_BAM}" \
    -O "${SAMPLE}.g.vcf.gz" \
    -ERC GVCF -L cleaned_intervals.bed
