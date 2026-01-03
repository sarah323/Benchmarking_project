

# NGS Data Analysis: Troubleshooting and Solutions

## 1. Project Setup

* Created project directory and activated environment:

```bash
mkdir workdir/ngs2-project && cd workdir/ngs2-project
conda activate ngs1
```

* Verified tools installed: `samtools`, `bwa`, `gatk`, `bcftools`.
* Installed `hap.py` in a separate environment due to compatibility issues.

---

## 2. Input Data

* Downloaded HG002 WES BAM (aligned to GRCh37, mark-duplicated):

```bash
wget ftp://ftp-trace.ncbi.nlm.nih.gov/.../HG002.posiSrt.markDup.bam
```

* Explored BAM:

```bash
samtools flagstat HG002.posiSrt.markDup.bam
```

* Observed ~150 million reads, 99.7% mapped, 7.3M duplicates.

---

## 3. Problem: Converting BAM → FASTQ for realignment to GRCh38

* Initial attempt using `samtools fastq`:

```bash
samtools fastq -1 HG002_R1.fastq.gz -2 HG002_R2.fastq.gz -s HG002_singletons.fastq.gz -0 /dev/null -n HG002.posiSrt.markDup.bam
```

**Issue:**

* ~135M reads went to the singleton file.
* Resulting R1/R2 FASTQ files were too small (~690MB vs 9.7GB BAM).

**Reason:**

* BAM was coordinate-sorted; `samtools fastq -n` expects name-sorted BAM → paired reads were not adjacent → discarded as singletons.

---

## 4. Solution: Name-sorting the BAM

```bash
samtools sort -n -o HG002.namesort.bam HG002.posiSrt.markDup.bam
samtools fastq -F 1024 -1 HG002_R1.fastq.gz -2 HG002_R2.fastq.gz -s /dev/null -n HG002.namesort.bam
```

✅ Outcome: Minimal singletons discarded, resulting FASTQ files ~6.2–6.4 GB.

**Alternative:** Used `Picard SamToFastq` for extra validation:

```bash
picard SamToFastq I=HG002.namesort.bam F=HG002_R1.clean.fastq.gz F2=HG002_R2.clean.fastq.gz INCLUDE_NON_PF_READS=true VALIDATION_STRINGENCY=SILENT
```

---

## 5. Reference Preparation

* Downloaded GRCh38 reference:

```bash
wget https://storage.googleapis.com/.../Homo_sapiens_assembly38.fasta
mv Homo_sapiens_assembly38.fasta GRCh38.fa
samtools faidx GRCh38.fa
gatk CreateSequenceDictionary -R GRCh38.fa -O GRCh38.dict
bwa index GRCh38.fa
```

---

## 6. Problem: Whole-genome Alignment Failures

* Initial BWA MEM alignment to GRCh38 failed repeatedly:

```text
[E::sam_parse1] SEQ and QUAL are of different length
[W::sam_read1] Parse error ... samtools sort: truncated file
```

* Integrity checks showed FASTQ files were valid.
* Likely causes: large genome, excessive memory, ambiguous pairings, or supplementary alignments.

---

## 7. Solution: Aligning to Chromosome 22 Only (2 workflows were followed to accomodate for different computational powers)

* Created chr22 reference:

```bash
samtools faidx GRCh38.fa chr22 > GRCh38_chr22.fa
bwa index GRCh38_chr22.fa
```

* Aligned reads to chr22:

```bash
bwa mem -t 2 -K 10000000 -R '@RG\tID:HG002\tSM:HG002\tPL:ILLUMINA' GRCh38_chr22.fa HG002_R1.clean.fastq.gz HG002_R2.clean.fastq.gz | samtools view -b - | samtools sort -@2 -m512M -o HG002_chr22.sorted.bam
```

✅ Outcome: BAM 12GB, all reads preserved, reasonable singleton percentage (~1%).

**Observation:** BAM larger than original due to:

* No duplicate marking
* Retained unmapped and supplementary reads

---

## 8. Depth of Coverage

* Whole chr22:

```bash
samtools depth -a HG002_chr22.sorted.bam | awk '{sum+=$3; cnt++} END {print sum/cnt}'
```

* Average coverage: 27.17x

---

## 9. Using Agilent SureSelect Human All Exon V5 BED File

* Determined capture kit from BAM header: `excap51` → Agilent SureSelect V5 (~51 Mb).
* Downloaded BED from Agilent SureDesign and organized project:

```bash
mkdir agilentbedfiles
unzip S04380110_hg38.zip -d agilentbedfiles/
grep -w "^chr22" agilentbedfiles/S04380110_Regions.bed > S04380110_chr22.bed
```

---

## 10. Restrict BAM to Capture Regions

```bash
samtools view -b -L S04380110_chr22.bed HG002_chr22.sorted.bam > HG002_chr22_V5_ontarget.bam
samtools index HG002_chr22_V5_ontarget.bam
```

* Average coverage (all targets including 0-depth): 251.79x
* Average coverage (covered bases only): 160.75x

**Note:** These values are inflated for mapped reads only; including 0-depth gives true capture performance.

---

## 11. Marking Duplicates After Alignment

```bash
picard MarkDuplicates I=HG002_chr22.sorted.bam O=HG002_chr22.sorted.markdup.bam M=HG002_chr22.markdup.metrics.txt CREATE_INDEX=true VALIDATION_STRINGENCY=STRICT
```

* Then restricted duplicate-marked BAM to Agilent V5 regions:

```bash
samtools view -b -L S04380110_chr22_Regions.bed HG002_chr22.sorted.markdup.bam > HG002_chr22_V5_ontarget.markdup.bam
samtools index HG002_chr22_V5_ontarget.markdup.bam
```

✅ Outcome: BAM ready for variant calling; coverage and alignment validated.

---

### Key Lessons / Troubleshooting Highlights

1. **BAM → FASTQ conversion:** Name-sort BAM to preserve pairing, or use Picard.
2. **Alignment failures:** Large genome → memory limits; alignment to single chromosome resolves parsing errors.
3. **BAM size inflation:** Includes duplicates, unmapped, and supplementary reads.
4. **On-target analysis:** Use Agilent SureSelect BED file to restrict coverage calculations accurately.
5. **Duplicate marking:** Should be done **before restricting to capture regions** to ensure accurate downstream analyses.

---
