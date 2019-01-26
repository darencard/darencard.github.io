---
layout: posts
title: "Useful One-liners for Calculating Population Genetic Statistics from VCF files"
date: 2017-07-18
excerpt: "Various Bash one-liners for summarizing variant data."
---

Useful one-liners that I have come up with to summarize variant data from VCF files.

The following commands require non-standard software like [BCFtools](http://www.htslib.org/) and [VCFtools](https://vcftools.github.io/index.html).

**Description:** Thin variants to prevent linkage biases and output the number of sampled alleles and the allele frequency for the reference allele.
```bash
vcftools --thin 10000 --recode --recode-INFO-all --stdout --gzvcf <my_variants.vcf.gz> | \
  bcftools query -f '%CHROM\t%POS[\t%GT]\n' - | \
  awk -v OFS="\t" '{ miss=0; hom_ref=0; hom_alt=0; het=0; \
    for (i=3; i<=NF; i++) \
    if ($i == "./.") miss += 1; \
    else if ($i == "0/0") hom_ref += 1; \
    else if ($i == "1/1") hom_alt += 1; \
    else if ($i == "0/1" || $i == "1/0") het += 1; \
    pres=hom_ref+hom_alt+het; \
    print $1, $2, pres*2, ((2*hom_ref)+het)/(2*pres) }'
```
