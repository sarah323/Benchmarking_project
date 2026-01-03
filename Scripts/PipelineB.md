

# HG002 WES chr22 Subsampling and Variant Calling Workflow (GRCh38) (Direct chr22 alignment)

This repository documents a complete workflow for processing **HG002 whole-exome sequencing (WES)** data, restricting analysis to **chromosome 22**, subsampling to multiple coverage depths, performing **GATK variant calling**, and benchmarking results using **hap.py** against GIAB truth sets. (The main difference between pipeline B and pipeline A is that in pipeline B, we aligned the FASTQ R1 and R2 files directly to the GRCh38_chr22 reference from the start. In contrast, in pipeline A, we first aligned the FASTQ files to the entire GRCh38 reference and then subsampled the resulting BAM file to chromosome 22.).

---

## 1. Project Setup

### Create project directory

```bash
mkdir -p workdir/ngs2-project
cd workdir/ngs2-project
```

### Verify required tools

```bash
samtools --version
bwa
gatk --version
bcftools --version
```

---

## 2. hap.py Environment Setup

`hap.py` was installed in a separate Conda environment due to dependency constraints.

```bash
conda create -n hap_py_env -c bioconda -c conda-forge hap.py -y
conda activate hap_py_env
hap.py --help
```

Reactivate the main analysis environment:

```bash
conda activate ngs1
```

---

## 3. Download HG002 WES BAM (GRCh37)

The HG002 WES dataset was available as a **pre-aligned, mark-duplicated BAM** (GRCh37).

```bash
wget ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
```
```bash
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai
```
---

## 4. Inspect Original BAM

```bash
samtools view -c -f 1 151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
samtools flagstat 151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
```

---

## 5. Name-Sort BAM and Convert to FASTQ

```bash
samtools sort -n \
  -o HG002_posiSrt.markDup.namesort.bam \
  151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
```

### Convert BAM → FASTQ using Picard

```bash
picard SamToFastq \
  I=HG002_posiSrt.markDup.namesort.bam \
  F=HG002_R1.clean.fastq.gz \
  F2=HG002_R2.clean.fastq.gz \
  INCLUDE_NON_PF_READS=true \
  VALIDATION_STRINGENCY=SILENT
```

---

## 6. Reference Genome Preparation (GRCh38)

```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
mv Homo_sapiens_assembly38.fasta GRCh38.fa
```

```bash
samtools faidx GRCh38.fa
gatk CreateSequenceDictionary -R GRCh38.fa -O GRCh38.dict
bwa index GRCh38.fa
```

---

## 7. Restrict Reference to Chromosome 22

```bash
samtools faidx GRCh38.fa chr22 > GRCh38_chr22.fa
samtools faidx GRCh38_chr22.fa
bwa index GRCh38_chr22.fa
```

---

## 8. Align FASTQs to chr22 Only

```bash
bwa mem \
  -t 2 \
  -K 10000000 \
  -R '@RG\tID:HG002\tSM:HG002\tPL:ILLUMINA' \
  GRCh38_chr22.fa \
  HG002_R1.clean.fastq.gz \
  HG002_R2.clean.fastq.gz | \
samtools view -b - | \
samtools sort -@ 2 -m 512M -o HG002_chr22.sorted.bam
```

```bash
samtools index HG002_chr22.sorted.bam
samtools flagstat HG002_chr22.sorted.bam
```

---

## 9. Coverage Calculation (Whole chr22)

```bash
samtools depth -a HG002_chr22.sorted.bam | \
awk '{sum+=$3; cnt++} END {print "Average depth =", sum/cnt}'
```

---

## 10. Identify Exome Capture Kit

BAM header inspection:

```bash
samtools view -H 151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
```

Metadata indicated:

```
excap51 → Agilent SureSelect Human All Exon V5 (~51 Mb)
```

---

## 11. Prepare Capture BED (Agilent V5)

Downloaded from **Agilent SureDesign** (hg38).

```bash
unzip S04380110_hg38.zip
grep -w "^chr22" S04380110_Regions.bed > S04380110_chr22_Regions.bed
```

---

## 12. Mark Duplicates (Picard)

```bash
picard MarkDuplicates \
  I=HG002_chr22.sorted.bam \
  O=HG002_chr22.sorted.markdup.bam \
  M=HG002_chr22.markdup.metrics.txt \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=STRICT
```

---

## 13. Restrict BAM to On-Target Regions

```bash
samtools view -b \
  -L S04380110_chr22_Regions.bed \
  HG002_chr22.sorted.markdup.bam \
  > HG002_chr22_V5_ontarget.markdup.bam

samtools index HG002_chr22_V5_ontarget.markdup.bam
```

### On-target coverage

```bash
samtools depth \
  -b S04380110_chr22_Regions.bed \
  HG002_chr22_V5_ontarget.markdup.bam | \
awk '{sum+=$3; cnt++} END {print "Mean on-target depth =", sum/cnt}'
```

---

## 14. Subsample BAM to Target Depths

Using fixed seed (`42`) for reproducibility.

```bash
samtools view -b -s 42.336 HG002_chr22_V5_ontarget.markdup.bam > HG002_chr22_V5_ontarget.markdup_80x.bam
samtools view -b -s 42.168 HG002_chr22_V5_ontarget.markdup.bam > HG002_chr22_V5_ontarget.markdup_40x.bam
samtools view -b -s 42.042 HG002_chr22_V5_ontarget.markdup.bam > HG002_chr22_V5_ontarget.markdup_10x.bam
samtools view -b -s 42.0084 HG002_chr22_V5_ontarget.markdup.bam > HG002_chr22_V5_ontarget.markdup_2x.bam
```

```bash
samtools index HG002_chr22_V5_ontarget.markdup_*.bam
```

---

## 15. Validate Subsampled Coverage

```bash
samtools depth -b S04380110_chr22_Regions.bed HG002_chr22_V5_ontarget.markdup_40x.bam | awk '{sum+=$3; cnt++} END {print sum/cnt}'
samtools depth -b S04380110_chr22_Regions.bed HG002_chr22_V5_ontarget.markdup_80x.bam | awk '{sum+=$3; cnt++} END {print sum/cnt}'
samtools depth -b S04380110_chr22_Regions.bed HG002_chr22_V5_ontarget.markdup_10x.bam | awk '{sum+=$3; cnt++} END {print sum/cnt}'
samtools depth -b S04380110_chr22_Regions.bed HG002_chr22_V5_ontarget.markdup_2x.bam  | awk '{sum+=$3; cnt++} END {print sum/cnt}'
```

---

## 16. Variant Calling (GATK HaplotypeCaller)

```bash
for depth in 2x 10x 40x 80x
do
  gatk HaplotypeCaller \
    -R ref/GRCh38.fa \
    -I HG002_chr22_V5_ontarget.markdup_${depth}.bam \
    -L chr22 \
    -O HG002_chr22_${depth}.raw.vcf.gz
done
```

---

## 17. GIAB Truth Set Preparation

```bash
mkdir -p truth
cd truth
```

```bash
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
```

```bash
bcftools view -r chr22 HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz -Oz -o HG002_chr22_truth.vcf.gz
tabix -p vcf HG002_chr22_truth.vcf.gz
grep '^chr22' HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed > HG002_chr22_truth.bed
```

---

## 18. hap.py Benchmarking

```bash
conda activate hap_py_env
conda install -c bioconda rtg-tools
```

```bash
export HGREF=/full/path/to/ref/GRCh38_chr22.fa
```

```bash
for depth in 2x 10x 40x 80x
do
  hap.py \
    truth/HG002_chr22_truth.vcf.gz \
    HG002_chr22_${depth}.raw.vcf.gz \
    -f truth/HG002_chr22_truth.bed \
    -r ref/GRCh38_chr22.fa \
    -o hap_py_chr22_${depth} \
    --engine=vcfeval
done
```

---

## Notes

* All steps are reproducible and chromosome-restricted for performance.
* Fixed random seeds ensure deterministic subsampling.
* Workflow is suitable for benchmarking variant calling sensitivity across coverage depths.

---

