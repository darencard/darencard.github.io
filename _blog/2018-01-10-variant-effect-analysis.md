---
layout: posts
title: "Variant Effect Analysis using VEP"
date: 2018-01-10
excerpt: "Understand potential impact of variants based on VEP analysis."
---

When calling variants it can be very useful to know what impact a polymorphism might have on the biology of an organism. It is often difficult to extend a simple list of variants to specific phenotypes, but it is possible to broadly classify a variant based on its possible impact on protein coding genes (e.g., mis-sense mutations, etc.). A very useful tool for performing such an analysis is the [Varient Effect Predictor (`VEP`)](http://ensembl.org/info/docs/tools/vep/index.html) tool, which is produced by the folks at [Ensembl](http://ensembl.org/). This tutorial will describe how to perform such an analysis on your organism of interest using variants stored in a VCF file and annotations from a GFF file.

## Software and Data

1. `VEP` must be installed. Detailed information on this tool and its installation can be found [here](http://ensembl.org/info/docs/tools/vep/script/index.html).
2. `Boa_constrictor_SGA_7C_scaffolds.fa`: A genome FASTA file for your organism is required to run `VEP`.
3. `Boa_40x_remapping.variants.hard_filters.vcf`: Genomic variants inferred using one of many approaches not described here must be provided as an unzipped VCF file.
4. `Bcon_rnd3.all.maker.noseq.format.gff.gz`: Gene annotations from an existing source of produced using annotation software like MAKER must be provided as a bg-zipped, tabix-indexed GFF file.

Importantly, the GFF file must contain `contig` records that define the coordinates of scaffolds or we will get errors or warnings. In my analysis, I began with a raw GFF file produced from MAKER, with the trailing sequence FASTA lines removed. This file contains lots of information about how protein and transcript evidence map to the genome, which isn't needed here. We are interested in the MAKER gene annotations, which are the lines that have `maker` in the 2nd column. We can process this raw GFF file to produce what we would like.

```bash
# keeps only lines with "maker" in 2nd column and "contig" in 3rd column, sorts output, and bg-zips it
grep -v "^#" Bcon_rnd3.all.maker.noseq.gff | \
awk '{ if ($2 == "maker" || $3 == "contig") print $0 }' | \
sort -k1,1 -k4,4n -k5,5n -t $'\t' | \
bgzip -c \
> Bcon_rnd3.all.maker.noseq.format.gff.gz
# we also need to tabix index this file
tabix -p gff Bcon_rnd3.all.maker.noseq.format.gff.gz
```

We will now run `VEP` in a specific way, which is outlined [here](http://ensembl.org/info/docs/tools/vep/script/vep_cache.html#gff). Note that this is simply scratching the surface of the capability of `VEP`, and see the [full documentation](http://ensembl.org/info/docs/tools/vep/script/index.html) for more information on running VEP in many other advanced ways.

Here is our specific command, which should be pretty self explanatory. This will take a while to run (hours) on larger variant datasets, but we can parallelize this process by passing the number of available cores to the `--fork` option.

```bash
path/to/ensembl-vep/vep --fork 1 \
-i Boa_40x_remapping.variants.hard_filters.vcf \
-gff Bcon_rnd3.all.maker.noseq.format.gff.gz \
-fasta Boa_constrictor_SGA_7C_scaffolds.fa \
-o Boa_variant_effects
```

Two output files are now produced:

`Boa_variant_effects`: The full data table for all variants. Note that there is some redudancy here as certain variants are processed a couple times in certain circumstances.

`Boa_variant_effects_summary.html`: A summary of the table. This did not render for me and I don't know why, but the full data table is what is important. It is pretty easy to get a rough idea of the amounts of certain types of variant effects by issuing a command like this:

```bash
grep -v "^#" Boa_variant_effects | \
cut -f 7 | sort | uniq -c
```

Details on what each of these categories (referred to as Sequence Ontologies [SOs]) mean is available [here](https://ensembl.org/info/genome/variation/predicted_data.html#consequences).

Let's say we want to get counts of the sequence ontologies for SNPs and InDels separately. We can do so using the following commands. Note that some variants will have multiple SO terms, which we aren't keeping track of here. Also, frameshift variants show up under both SNP and InDel output, though I would classify them as a type of InDel (even if they are only 1bp features).

```bash
# SNPs
cat Boa_variant_effects | grep -v "^#" | \
awk -F ":|\t" '{ if ($3 !~ "-") print $0 }' | \
cut -f 7 | tr ',' '\n' | sort | uniq -c
# InDels
cat Boa_variant_effects | grep -v "^#" | \
awk -F ":|\t" '{ if ($3 ~ "-") print $0 }' | \
cut -f 7 | tr ',' '\n' | sort | uniq -c
```

We may also wish to know how InDel lengths are distributed within coding regions, as InDel lengths that aren't a multiple of 3 are more likely to cause protein-coding changes. We can do that with the following command, which outputs the +/- InDel length for coding-region InDels.

```bash
cat Boa_variant_effects.txt | grep -e "insertion" -e "deletion" -e "frameshift" | \
cut -f 12 | tr -dc '[:upper:]|\/|\n' | \
awk '{ if ($1 ~ "^/") print length($1)-1; else if ($1 ~ "/$") print "-"length($1)-1 }'
```

We can also summarize what I think are the more important fields in a slightly better way by issuing a command like this:

```bash
cat Boa_variant_effects.txt | grep -v "^#" | \
awk -F "\t|:|;|IMPACT=|STRAND=" -v OFS="\t" \
'{ if ($4 == "-" || length($4) > 1) \
print $2, $3-1, $3+0, $18, "INDEL", $4, $5, $8, $9, $10, $11, $12, $13, $16; \
else print $2, $3-1, $3+0, $18, "SNP", $4, $5, $8, $9, $10, $11, $12, $13, $16 }' | \
awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "." }; 1'
```
