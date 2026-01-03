# Benchmarking Variant Calling Performance at Different Depths of Coverage (HG002, chr22)

This repository documents a **reproducible benchmarking workflow** to evaluate how sequencing depth affects germline variant calling performance using **GIAB HG002**, **GATK Best Practices**, and **hap.py**.

All analyses are restricted to **chromosome 22** to reduce computational cost while maintaining benchmarking validity.

---

## Overview

- **Sample:** GIAB HG002
- **Sequencing type:** Whole-exome Illumina
- **Reference genome:** GRCh38
- **Capture kit:** Agilent SureSelect V5 (hg38)
- **Coverage levels evaluated:** 2×, 10×, 40×, 80×
- **Variant caller:** GATK HaplotypeCaller
- **Benchmarking tool:** hap.py (GA4GH)
- **Truth set:** GIAB HG002 v4.2.1 (VCF + confident regions BED), chr22

---

## Project Structure

```text
benchmark_project/
├── raw_data/        # FASTQ files
├── ref/             # Reference genome + BED files
├── bam/             # BAM files at each processing stage
├── vcf/             # Raw VCF outputs
├── truth/           # GIAB truth VCF + BED
└── metrics/happy/   # hap.py benchmarking results

#!/bin/bash
set -euo pipefail

########################################
# My paths (project structure)
########################################
PROJECT=~/benchmark_project
RAW=${PROJECT}/raw_data
REF=${PROJECT}/ref
BAM=${PROJECT}/bam
VCF=${PROJECT}/vcf
TRUTH=${PROJECT}/truth
METRICS=${PROJECT}/metrics/happy

mkdir -p "$BAM" "$VCF" "$METRICS"

########################################
# 0) From original exome BAM → FASTQ
########################################
cd "$BAM"

samtools fastq \
  -1 "${RAW}/HG002.exome.R1.ns.fastq.gz" \
  -2 "${RAW}/HG002.exome.R2.ns.fastq.gz" \
  -0 /dev/null \
  -s /dev/null \
  -n \
  HG002.exome.markdup.bam

########################################
# 1) FASTQ sanity check
########################################
cd "$RAW"
zcat HG002.exome.R1.ns.fastq.gz | wc -l
zcat HG002.exome.R2.ns.fastq.gz | wc -l

########################################
# 2) Align to GRCh38 and sort
########################################
cd "$BAM"

bwa mem -t 4 \
  -R '@RG\tID:HG002_exome\tSM:HG002\tPL:ILLUMINA\tLB:OsloExome' \
  "$REF/hg38.fa" \
  "$RAW/HG002.exome.R1.ns.fastq.gz" \
  "$RAW/HG002.exome.R2.ns.fastq.gz" \
  | samtools sort -o HG002.hg38.exome.sorted.bam

samtools index HG002.hg38.exome.sorted.bam

samtools quickcheck HG002.hg38.exome.sorted.bam && echo "BAM looks OK"

########################################
# 3) Mark duplicates (GATK)
########################################
gatk MarkDuplicates \
  -I HG002.hg38.exome.sorted.bam \
  -O HG002.hg38.exome.markdup.bam \
  -M HG002.hg38.exome.markdup.metrics.txt \
  --CREATE_INDEX true

########################################
# 4) Restrict to chr22 and estimate depth
########################################
samtools view -b HG002.hg38.exome.markdup.bam chr22 \
  > HG002.hg38.chr22.full.bam

samtools index HG002.hg38.chr22.full.bam

samtools depth HG002.hg38.chr22.full.bam \
  | awk '{sum+=$3; cnt++} END {print "Mean depth chr22 =", sum/cnt}'

########################################
# 5) Prepare SureSelect V5 chr22 BED
########################################
cd "$REF"

mv SureSelect_V5_hg38_targets_dir/S04380110_Regions.bed \
   SureSelect_V5_hg38_targets.bed

grep -w 'chr22' SureSelect_V5_hg38_targets.bed \
  > SureSelect_V5_hg38_chr22.bed

########################################
# 6) Restrict BAM to capture regions
########################################
cd "$BAM"

samtools view -b \
  -L ../ref/SureSelect_V5_hg38_chr22.bed \
  HG002.hg38.chr22.full.bam \
  > HG002.hg38.chr22.V5_ontarget.bam

samtools index HG002.hg38.chr22.V5_ontarget.bam

samtools depth \
  -b ../ref/SureSelect_V5_hg38_chr22.bed \
  HG002.hg38.chr22.V5_ontarget.bam \
  | awk '{sum+=$3; cnt++} END {print "Mean on-target depth =", sum/cnt}'

########################################
# 7) Subsample to target depths
########################################
samtools view -b -s 42.3617 HG002.hg38.chr22.V5_ontarget.bam > HG002.hg38.chr22.V5_ontarget_80x.bam
samtools view -b -s 42.1808 HG002.hg38.chr22.V5_ontarget.bam > HG002.hg38.chr22.V5_ontarget_40x.bam
samtools view -b -s 42.0452 HG002.hg38.chr22.V5_ontarget.bam > HG002.hg38.chr22.V5_ontarget_10x.bam
samtools view -b -s 42.0090 HG002.hg38.chr22.V5_ontarget.bam > HG002.hg38.chr22.V5_ontarget_2x.bam

samtools index HG002.hg38.chr22.V5_ontarget_{2x,10x,40x,80x}.bam

########################################
# 8) Validate depth for each subset
########################################
for depth in 2x 10x 40x 80x; do
  samtools depth \
    -b ../ref/SureSelect_V5_hg38_chr22.bed \
    HG002.hg38.chr22.V5_ontarget_${depth}.bam \
    | awk '{sum+=$3; cnt++} END {print "Mean depth", "'"$depth"'", "=", sum/cnt}'
done

########################################
# 9) GATK HaplotypeCaller
########################################
for depth in 2x 10x 40x 80x; do
  gatk HaplotypeCaller \
    -R "$REF/hg38.fa" \
    -I "$BAM/HG002.hg38.chr22.V5_ontarget_${depth}.bam" \
    -L chr22 \
    -O "$VCF/HG002_chr22_${depth}.raw.vcf.gz"
done

########################################
# 10) Index VCFs
########################################
for depth in 2x 10x 40x 80x; do
  gatk IndexFeatureFile \
    -I "$VCF/HG002_chr22_${depth}.raw.vcf.gz"
done

########################################
# 11) Benchmark with hap.py
########################################
for depth in 2x 10x 40x 80x; do
  hap.py \
    "$TRUTH/HG002_GRCh

