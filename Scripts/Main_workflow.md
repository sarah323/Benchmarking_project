
```markdown
# Benchmarking Variant Calling Performance at Different Depths (HG002, WGS)

This repository provides a **reproducible benchmarking workflow** to evaluate how **sequencing depth** affects **germline variant calling** performance using **GIAB HG002**, **GATK Best Practices**, and **hap.py**.

**Key difference between pipelines:**
- **Pipeline A:** Align FASTQ → whole GRCh38 → subsample BAM to chr22  
- **Pipeline B:** Align FASTQ directly → GRCh38_chr22

**All analyses are restricted to chromosome 22** to reduce computational cost while maintaining benchmarking validity.

---

## Overview

- **Sample:** GIAB HG002  
- **Sequencing type:** Whole-exome Illumina  
- **Reference genome:** GRCh38  
- **Capture kit:** Agilent SureSelect V5 (hg38)  
- **Coverage levels:** 2×, 10×, 40×, 80×  
- **Variant caller:** GATK HaplotypeCaller  
- **Benchmarking tool:** hap.py (GA4GH)  
- **Truth set:** GIAB HG002 v4.2.1 (VCF + confident regions BED), chr22  

---

## Project Structure

```

benchmark_project/
├── raw_data/        # FASTQ files
├── ref/             # Reference genome + BED files
├── bam/             # BAM files at each processing stage
├── vcf/             # Raw VCF outputs
├── truth/           # GIAB truth VCF + BED
└── metrics/happy/   # hap.py benchmarking results

````

---

## Setup

### Paths
```bash
PROJECT=~/benchmark_project
RAW=${PROJECT}/raw_data
REF=${PROJECT}/ref
BAM=${PROJECT}/bam
VCF=${PROJECT}/vcf
TRUTH=${PROJECT}/truth
METRICS=${PROJECT}/metrics/happy
mkdir -p "$BAM" "$VCF" "$METRICS" "$RAW" "$TRUTH" "$REF"
````

### Tools

```bash
samtools --version
bwa
gatk --version
bcftools --version

# Install hap.py in separate environment
conda create -n happenv -c bioconda -c conda-forge hap.py=0.3.12 python=3.7
conda activate happenv
which hap.py
```

---

## Data Preparation

### Download BAM & reference

```bash
wget <GIAB BAM URLs>
wget <GRCh38 reference URL>
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna hg38.fa
samtools faidx hg38.fa
bwa index hg38.fa
gatk CreateSequenceDictionary -R hg38.fa -O hg38.dict
```

### Download truth set

```bash
wget <GIAB truth VCF & BED URLs>
mv HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed HG002_GRCh38_1_22_v4.2.1_benchmark.bed
```

---

## Workflow Steps

### 0. Convert original BAM → FASTQ

```bash
samtools fastq -1 "${RAW}/HG002.exome.R1.ns.fastq.gz" \
               -2 "${RAW}/HG002.exome.R2.ns.fastq.gz" \
               -0 /dev/null -s /dev/null -n \
               HG002_posiSrt.markDup.namesort.bam
```

### 1. FASTQ sanity check

```bash
zcat HG002.exome.R1.ns.fastq.gz | wc -l
zcat HG002.exome.R2.ns.fastq.gz | wc -l
```

### 2. Align to GRCh38 and sort

```bash
bwa mem -t 4 -R '@RG\tID:HG002_exome\tSM:HG002\tPL:ILLUMINA\tLB:OsloExome' \
        "$REF/hg38.fa" \
        "$RAW/HG002.exome.R1.ns.fastq.gz" \
        "$RAW/HG002.exome.R2.ns.fastq.gz" \
  | samtools sort -o HG002.hg38.exome.sorted.bam
samtools index HG002.hg38.exome.sorted.bam
```

### 3. Mark duplicates (GATK)

```bash
gatk MarkDuplicates \
  -I HG002.hg38.exome.sorted.bam \
  -O HG002.hg38.exome.markdup.bam \
  -M HG002.hg38.exome.markdup.metrics.txt \
  --CREATE_INDEX true
```

### 4. Restrict to chr22 & depth estimation

```bash
samtools view -b HG002.hg38.exome.markdup.bam chr22 > HG002.hg38.chr22.full.bam
samtools index HG002.hg38.chr22.full.bam
```

### 5. Prepare SureSelect V5 chr22 BED

* Filter original BED to chr22:

```bash
grep -w 'chr22' SureSelect_V5_hg38_targets.bed > SureSelect_V5_hg38_chr22.bed
```

### 6. Restrict BAM to capture regions

```bash
samtools view -b -L ../ref/SureSelect_V5_hg38_chr22.bed \
  HG002.hg38.chr22.full.bam > HG002.hg38.chr22.V5_ontarget.bam
samtools index HG002.hg38.chr22.V5_ontarget.bam
```

### 7. Subsample to target depths

```bash
samtools view -b -s 42.3617 HG002.hg38.chr22.V5_ontarget.bam > HG002.hg38.chr22.V5_ontarget_80x.bam
# repeat for 40x, 10x, 2x
```

### 8. Validate depth

```bash
for depth in 2x 10x 40x 80x; do
  samtools depth -b ../ref/SureSelect_V5_hg38_chr22.bed \
    HG002.hg38.chr22.V5_ontarget_${depth}.bam \
    | awk '{sum+=$3; cnt++} END {print "Mean depth", "'"$depth"'", "=", sum/cnt}'
done
```

### 9. GATK HaplotypeCaller

```bash
for depth in 2x 10x 40x 80x; do
  gatk HaplotypeCaller -R "$REF/hg38.fa" \
                       -I "$BAM/HG002.hg38.chr22.V5_ontarget_${depth}.bam" \
                       -L chr22 \
                       -O "$VCF/HG002_chr22_${depth}.raw.vcf.gz"
done
```

### 10. Index VCFs

```bash
for depth in 2x 10x 40x 80x; do
  gatk IndexFeatureFile -I "$VCF/HG002_chr22_${depth}.raw.vcf.gz"
done
```

### 11. Benchmark with hap.py

```bash
for depth in 2x 10x 40x 80x; do
  hap.py truth/HG002_GRCh38_chr22_v4.2.1_benchmark.vcf.gz \
         vcf/HG002.hg38.chr22_V5_ontarget_${depth}.raw.vcf.gz \
         -f truth/HG002_GRCh38_chr22_v4.2.1_benchmark.bed \
         -r ref/hg38.fa \
         -o metrics/happy/HG002.chr22_${depth}
done
```

---

## Notes

* PipelineB aligns **directly to chr22**, PipelineA aligns **whole genome then subsamples**.
* Analyses are restricted to **chr22** for computational efficiency.
* Coverage evaluation: 2×, 10×, 40×, 80×.
* Benchmarking performed using **hap.py** with GIAB truth set.

```

---


