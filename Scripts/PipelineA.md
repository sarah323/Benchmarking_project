# Benchmarking Variant Calling Performance at Different Depths of Coverage (HG002) (Whole genome alignment)

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

mkdir -p "$BAM" "$VCF" "$METRICS" "$raw_data" "$truth" "$ref"

########################################
##Make sure that all the necessary packages and tools are installed:
samtools --version
bwa
gatk --version
bcftools –version

##Installing hap.py in a separate environment because of compatibility issues 
conda create -n happenv -c bioconda -c conda-forge hap.py=0.3.12 python=3.7
conda activate happenv
which hap.py

########################################
##Activating ngs1 environment again
conda activate ngs1
########################################
##Downloading and preparing files:
#######################################
#downloading the GIAB reference (WES)
Wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai

#Downloading GRCh38 reference and unzipping it
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

#Renaming the file
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna hg38.fa

#indexing
samtools faidx hg38.fa
bwa index hg38.fa

#creating a sequence dictionary (describing chromosomes)
gatk CreateSequenceDictionary \
   -R hg38.fa \
   -O hg38.dict

#downloading the truthset and its index
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi
wget https://ftp.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

#renaming the file
mv HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
>    HG002_GRCh38_1_22_v4.2.1_benchmark.bed


#######################################
##exploring the original bam file and here is the code along with its output:
samtools view -c -f 1 151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
##150538907

samtools flagstat 151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
##The output
150538907 + 0 in total (QC-passed reads + QC-failed reads)
152131 + 0 secondary
0 + 0 supplementary
7305282 + 0 duplicates
150113954 + 0 mapped (99.72% : N/A)
150386776 + 0 paired in sequencing
75193388 + 0 read1
75193388 + 0 read2
148839498 + 0 properly paired (98.97% : N/A)
149715400 + 0 with itself and mate mapped
246423 + 0 singletons (0.16% : N/A)
307022 + 0 with mate mapped to a different chr
247662 + 0 with mate mapped to a different chr (mapQ>=5)
#######################################
##Sorting the original bam file
samtools sort -n -o HG002_posiSrt.markDup.namesort.bam   151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam
#This resulted in a 15GB sorted file

#######################################
# 0) From original exome BAM → FASTQ
########################################
cd "$BAM"

samtools fastq \
  -1 "${RAW}/HG002.exome.R1.ns.fastq.gz" \
  -2 "${RAW}/HG002.exome.R2.ns.fastq.gz" \
  -0 /dev/null \
  -s /dev/null \
  -n \
  HG002_posiSrt.markDup.namesort.bam

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
  ##Mean depth chr22 =17

########################################
# 5) Prepare SureSelect V5 chr22 BED
########################################
###The next step is to restrict our bam file to WES capture kit used for the HG002 sample
##First, we viewed the head of the original bam file to see the type of sequencing kit used, using this command
samtools view -H 151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam

##We saw this multiple times in the read group and sample metadata:
#Sample_Diag-excap51-HG002-EEogPU
#Project_Diag-excap51-2015-09-23
##excap51 is a very common internal shorthand used by diagnostic labs for:
#Exome capture ~51 Mb → Agilent SureSelect Human All Exon V5

##Then, an account on the Agilent SureDesign was created then we searched for and downloaded the matching bed file “SureSelect Human All Exon V5” ==> from Suredesign ==> sureselect DNA ==> Agilent Catalog ==> SureSelect Human All Exon V5 (hg38) ==> Agilent SureDesign

##Copying it to my work directory
cp ../../Downloads/S04380110_hg38.zip .

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
  ##Mean on-target depth =221

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
conda activate happenv

#1) Run hap.py once on the 80x BAM (as a test)
---------------------------------
hap.py \
  truth/HG002_GRCh38_chr22_v4.2.1_benchmark.vcf.gz \
  vcf/HG002.hg38.chr22_V5_ontarget_80x.raw.vcf.gz \
  -f truth/HG002_GRCh38_chr22_v4.2.1_benchmark.bed \
  -r ref/hg38.fa \
  -o metrics/happy/HG002.chr22_80x

#2) Run hap.py for all four depths (2x, 10x, 40x, 80x)
-----------------------------------------------------
for depth in 2x 10x 40x 80x
do
  echo "Running hap.py for depth ${depth} ..."
  hap.py \
    truth/HG002_GRCh38_chr22_v4.2.1_benchmark.vcf.gz \
    vcf/HG002.hg38.chr22_V5_ontarget_${depth}.raw.vcf.gz \
    -f truth/HG002_GRCh38_chr22_v4.2.1_benchmark.bed \
    -r ref/hg38.fa \
    -o metrics/happy/HG002.chr22_${depth}
done

#3) Extract SNP/INDEL PASS lines from each summary file
------------------------------------------------------
for depth in 2x 10x 40x 80x
do
  echo "===== ${depth} SNP PASS ====="
  grep "SNP, PASS" metrics/happy/HG002.chr22_${depth}.summary.csv
  echo
  echo "===== ${depth} INDEL PASS ====="
  grep "INDEL, PASS" metrics/happy/HG002.chr22_${depth}.summary.csv
  echo
done




