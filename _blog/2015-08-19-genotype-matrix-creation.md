---
layout: posts
title: "Population structure analysis input file generation"
date: 2015-08-19
excerpt: "Generating input files for NGSadmix and Entropy."
---

### NGSadmix (Skotte et al. 2013)
NGSadmix genotype matrices include a header line and two beginning columns (with headers) with the marker ID (scaffold and position) and the reference and alternative allele (all sites must be biallelic). Three genotype likelihoods are given for each sample and marker in a standardized format (sum to 1.0) and correspond to the likelihood of increasingless less reference alleles (homozygous reference, heterozygous, homozygous alternative). All values are space-delimited and missing data is coded as 0.000 across all three allele combinations. Here is an example with three samples at two markers:
```
Marker Ref. Alt. Sample1 Sample1 Sample1 Sample2 Sample2 Sample2 Sample3 Sample3 Sample3
scaffold1_100 A C 1.000 0.000 0.000 0.333 0.333 0.333 0.250 0.750 0.000
scaffold2_1000 G T 0.000 0.000 0.000 0.500 0.500 0.000 0.010 0.990 0.000
```
The following RADpipe command will create this as output from a filtered VCF:
```
python genotypes_from_VCF.py --samplesheet <samplesheet.txt> --filvcf <filteredVCF.vcf> --genotype 3 --locinfo --refalt --headers 2,4 --delimit 1 --prefix <output>
```

### Entropy (Gompert et al. 2014)
Entropy genotype matrices include three header lines: (1) one with the matrix dimensions (with a trailing 1 that is required), (2) one with the sample IDs, and (3) one with the population assignments. They also include a leading column (without a header) with the marker ID (scaffold and position). Three genotype likelihoods are given for each sample and marker in the raw PHRED likelihood format (0-255), with higher numbers indicating more support for a given allelic state. The three values are ordered from homozygous reference to homozygous alternative. All values in the matrix are space-delimited and missing data is coded as 0 across all three allele combinations. Here is an example with three samples at two markers.
```
3 2 1
Sample1 Sample1 Sample1 Sample2 Sample2 Sample2 Sample3 Sample3 Sample3
scaffold1_100 0 0 0 255 0 0 255 255 0
scaffold2_1000 0 12 125 175 110 0 255 10 0
```
The following RADpipe command will create this as output from a filtered VCF:
```
python genotypes_from_VCF.py --samplesheet <samplesheet.txt> --filvcf <filteredVCF.vcf> --genotype 1 --locinfo --headers 1,2,3 --delimit 1 --prefix <output>
```

### Entropy Chain Initialization via DAPC
Gompert et al. 2014 (see supplementary materials) initialize their Entropy MCMC chains with approximate values of admixture proportions derived from DAPC (Jombert et al. 2010). The RADpipe function entropyStart.R will perform the DAPC analysis and output initialization files for Entropy for desired values of K. entropyStart.R requires a genotype matrix with three header rows and a leading marker ID column. A single genotype uncertainty value (0-2) is given for each sample, and basically provides the fuzzy estimation of the number of non-reference alleles (i.e., rather than being discrete integers, the value can vary continuously between 0 and 2 such that it provides a measure of uncertainty). For example, a genotype uncertainty value of 1.90 represents confidence that there are two non-reference alleles (homozygous alternative allele state), but also demonstrates that there is a bit of uncertainty in this genotype inference (it would be 2.00 if we were totally sure about the call). Here is an example matrix with three samples at two markers:
```
3 2 1
Marker Sample1 Sample2 Sample3
scaffold1_100 0.01 1.90 1.00
scaffold2_1000 0.00 1.05 0.45
```
The following RADpipe command will create this as output from a filtered VCF:
```
python genotypes_from_VCF.py --samplesheet <samplesheet.txt> --filvcf <filteredVCF.vcf> --genotype 1 --locinfo --headers 1,2,3,4 --delimit 1 --prefix <output>
```
Note: Missing data for genotype uncertainty values is currently represented as '-9.00'. This may change to just be 'NA' so that the matrix can easily be manipulated in R. The entropyStart.R function has not yet been updated to deal with these missing data values. This note will be removed once these scripts have been adjusted accordingly.

### Citations
Gompert, et al. 2014. Admixture and the organization of genetic diversity in a butterfly species complex revealed through common and rare genetic variants. Molecular Ecology 23 (18): 4555-4573. [doi: 10.1111/mec.12811](http://doi.org/10.1111/mec.12811).

Jombart, et al. 2010. Discriminant analysis of principal components: a new method for the analysis of genetically structured populations. BMC Genetics 11 (94). [doi: 10.1186/1471-2156-11-94](http://doi.org/10.1186/1471-2156-11-94).

Skotte, et al. 2013. Estimating individual admixture proportions from next generation sequencing data. Genetics 195 (3): 693-702. [doi: 10.1534/genetics.113.154138](http://doi.org/10.1534/genetics.113.154138).
