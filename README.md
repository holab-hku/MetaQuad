# MetaQuad
MetaQuad: Shared Informative Variants Dis-covery in Metagenomic Samples

MetaQuad is a tool that can detect **shared informative variants** in a population of shotgun metagenomic data. It uses a binomial mixture model to assess the informative variants among background noise. False positive mutations can be significanly reduced via MetaQuad. MetaQuad is incorporated from [MQuad](https://github.com/single-cell-genetics/MQuad), and extended to support shotgun metagenomic data.

A recommended pipeline for MetaQuad:

  1. filter and align reads from shotgun metaganomic sequencing to a reference database with some alignment tools, including bwa and bowtie2. Recommended suffix for the final files are ".filtered.sorted.bam". <br/>

  2. pileup variant with [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite). <br/> Recommended parameter: minMAF = 0.02 and minCOUNT = 100. A higher value is recommended as the computation time of subsequent steps can be significantly reduced. 

  3. detect informative mutations with MetaQuad. Recommended parameter: minSample = 5% of total sample size and delta BIC cutoff (-t) = 10, as it is a very strong evidence. 




## Installation

MetaQuad is available from the GitHub repository for latest version by the following command line: 

```
pip install -U git+https://github.com/holab-hku/MetaQuad.git
```

## Quick start

Once installed, the version and input parameters can be checked with metaquad -h

MetaQuad takes the output folder of cellSNP-lite: 

```
metaquad  -i cellSNP -o metaquad_test -p 20 -t 10 --minDP=3 --batchSize 512
```

## Manual

### Expected output

**Column description for BIC_params.csv:**

- deltaBIC: score of informativeness, higher is better. In MetaQuad, a variant can be called as an informative mutation only if the delta BIC is larger than the cutoff value. 

- params1, params2, model1BIC, model2BIC: fitted parameteres for the binomial models. 
  
- num_samples_minor_cpt: no. of samples in the minor component, can be used to filtering variants that only happens in few samples. It can be a tradeoff between true positive and false positive. 

- PASS_cutoff: T/F, whether a mutation passed the delta BIC cutoff.

- PASS_MINSAMPLES: T/F, whether a mutation passed the minSample.


## Calculation of nucleotide diversity

Nucleotide diversity can be used to measure the degree of polymorphism within a population. Nucleotide diversity can be calculated via nucleotide_diversity_code.R, with R package "vcfR": 

```
Rscript nucleotide_diversity_code.R cellSNP.cells.vcf.gz BIC_params.csv genome_length.csv
```
The expected output is a csv file, "Nucleotide_diversity.csv". cellSNP.cells_ARG.vcf.gz is one of the output file from cellsnp-lite. genome_length.csv is used to normalised the raw data, and it can also be covbases for each gene and sample, according to the actual situation. Genome lengths and covbases can be calculated with samtools and R script Genome_coverage.R: 

```
for i in  $(cat idlist)
do 
    samtools coverage ${i}.sorted.bam > ${i}.tab
    grep -Ff genomes.txt ${i}.tab > ${i}_coverage.tab
done

Rscript Genome_coverage.R
```
genomes.txt contains genomes that are interested, each lines for a genome or a gene. Genome_coverage.R will provide two csv files, genome_length.csv and cov_bases.csv. 






