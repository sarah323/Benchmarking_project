# Benchmarking Project: Evaluating Variant Calling Performance at Different Depths of Coverage

## Overview
This project benchmarks variant calling pipelines by evaluating the effect of sequencing coverage on variant detection accuracy. Using whole-exome Illumina reads from the HG002 reference sample, the analysis compares variant calling performance at multiple depths of coverage (2×, 10×, 40×, and 80×), restricted to **chromosome 22** to reduce computational cost.

The Genome in a Bottle (GIAB) consortium provides well-characterized human reference samples with high-confidence variant calls, serving as a gold standard for benchmarking.  
[GIAB Project Page](https://www.nist.gov/programs-projects/genome-bottle)

## Objectives
- Assess how depth of sequencing coverage affects variant calling accuracy.
- Generate benchmarking metrics for SNPs and indels.
- Explore the trade-off between sequencing cost and variant calling performance.
- Optionally, compare an alternative variant caller (e.g., DeepVariant or FreeBayes) at 40× coverage.

## Data
- **Reference genome:** GRCh38 human reference genome
- **Gold truth dataset:** GIAB HG002 VCF + BED (confident regions) restricted to chromosome 22
- **Input reads:** Whole-exome Illumina reads from HG002

## Methods
1. **Variant Calling:** Apply GATK Best Practices pipeline for germline variant calling.
2. **Benchmarking Metrics:** Use `hap.py` from the GA4GH benchmarking toolkit to compute:
   - Precision (Positive Predictive Value)
   - Recall (Sensitivity)
   - F1-score
   - False Positive Rate
   - Genotype Concordance
   - Separate metrics for SNPs and indels  

Reference for benchmarking methodology:  
[Best practices for benchmarking germline small-variant calls in human genomes](https://www.nature.com/articles/s41587-019-0054-x)

3. **Analysis and Visualization:**  
   - Summarize metrics across coverage depths.
   - Generate plots/tables for Precision–Recall–F1 and SNP vs. indel performance.
   - Discuss trends and insights regarding coverage, accuracy, and sequencing cost.

