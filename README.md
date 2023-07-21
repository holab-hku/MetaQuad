# MetaQuad
MetaQuad: Shared Informative Variants Discovery in Metagenomic Samples

MetaQuad is a tool designed to detect **shared informative variants** in a population of shotgun metagenomic data. It uses a density clustering model to assess the informative variants among background noise. MetaQuad is an extension of [MQuad](https://github.com/single-cell-genetics/MQuad), specifically tailored to support shotgun metagenomic data.

We recommend the following pipeline for using MetaQuad:

  1. Filter and align reads from shotgun metagenomic sequencing to a reference database using alignment tools such as bwa or bowtie2. The recommended suffix for the final files is ".filtered.sorted.bam". <br/>

  2. Pile up variants using [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite). <br/> with the following recommended parameters: minMAF = 0.02 and minCOUNT = 100. A higher value is recommended to reduce computation time for subsequent steps. 

  3. Detect informative mutations using MetaQuad, with the following recommended parameters: minSample = 5% of the total sample size and num_of_clusters >= 2. 




## Installation

MetaQuad is available from the GitHub repository for latest version by the following command line: 

```
pip install -U git+https://github.com/holab-hku/MetaQuad.git
```

## Quick start

Once installed, the version and input parameters can be checked with metaquad -h

MetaQuad takes the output folder of cellSNP-lite: 

```
metaquad  -i cellSNP -o Metaquad --minSample=5
```

## Manual

### Expected output

**Column description for BIC_params.csv:**

- variant_names: Name of variants.

- num_of_clusters: Predicted number of clusters for a variant. 
  
- mean_DP, mean_AD: Average AD and DP (respectively) across all samples.

- number_DP: Number of non-zero DP samples.




## Calculation of nucleotide diversity

Nucleotide diversity is a measure of the degree of polymorphism within a population. To calculate nucleotide diversity, you can use the nucleotide_diversity_code.R script, which requires the R package "vcfR". The script should be run as follows: 

```
Rscript nucleotide_diversity_code.R cellSNP.cells.vcf.gz BIC_params.csv genome_length.csv
```
The expected output from the nucleotide_diversity_code.R script is a CSV file called "Nucleotide_diversity.csv". This script requires two input files: cellSNP.cells_ARG.vcf.gz, which is one of the output files from cellsnp-lite, and genome_length.csv, which is used to normalize the raw data. The values in genome_length.csv can be either genome lengths or covbases for each gene and sample, depending on the specific situation. Genome lengths and covbases can be calculated with samtools and R script Genome_coverage.R: 

```
for i in  $(cat idlist)
do 
    samtools coverage ${i}.sorted.bam > ${i}.tab
    grep -Ff genomes.txt ${i}.tab > ${i}_coverage.tab
done

Rscript Genome_coverage.R
```
genomes.txt contains genomes that are interested, each lines for a genome or a gene. Genome_coverage.R will provide two csv files, genome_length.csv and cov_bases.csv. 






