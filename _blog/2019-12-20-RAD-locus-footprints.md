---
layout: posts
title: "Identifying genome-wide RAD-seq loci"
date: 2019-10-26
excerpt: "A useful script for identifying the footprints of RAD-seq loci from read mappings"
---

I have analyzed a lot of RAD-seq data largely using the pipelines built around traditional "best practices" whole genome resequencing data. The typical approach is to quality filter the read data, map to a reference genome, and call variants. However, the traditional variant file format (VCF) assumes that reads have been sampled from the entire genome, not necessarily from small portions of the genomes like what we get with RAD-seq data. This isn't a problem in most cases, but it is still worthwhile to determine the footprints/coordinates of RAD loci that were sequenced as part of an experiment. This helps when one wants to know how many independent loci they have good sequence sampling for and also helps the user know how many monomorphic loci are sampled alongside polymorphic ones. VCF only usually records variant (i.e., polymorphic) sites, leaving you naive to the number of RAD loci with no variation.

To overcome this, I wrote a script that attempts to infer the footprints/coordinates of RAD loci based on straight-forward (and easy to manipulate) criteria. It takes mapping data (i.e., BAM files) and a sample metadata sheet as input. The script will then look for sites that meet certain quality and coverage criteria and that are present in some proportion of a given population/experiment (i.e., missing data filter). From there, it merges adjacent sites together into loci, allowing a user-defined merging distance to span short lower coverage/quality stretches that can sometimes appear with stochastic read data; this prevents us from splitting one effective RAD loci into two just because there is a gap of, say, 2 bp that are lower coverage, for example. It then filters away shorter loci to only keep ones that are of a certain length. All these parameters are adjustable by the user, so it is easy to iterate over different parameters and see how it impacts the number and distribution of RAD loci.

There are a couple software dependencies that should be standard for anyone doing genomics: [samtools](http://www.htslib.org/) and [bedtools](https://bedtools.readthedocs.io/en/latest/). Output of this script are locus coordinates in BED format.

Otherwise, the input requirements and arguments (positional) are well-documented in the beginning of the script (see below). Hopefully, this is useful to others working with RAD-seq (or perhaps other reduced representation) data.

```bash
# script that identifies RAD loci for a given population sample based on
# several quality control metrics
#
# requires input metadata table (tab separated)
# with columns for sample IDs and for
# population/experiment assignment
#
# example:
# sample01  population01
# sample02  population01
# sample03  population02
# sample04  population02
#
# 9 positional arguments:
#
# ${1} = sample metadata table
# ${2} = table column with population/experiment assignments (2 in above example)
# ${3} = table column with sample IDs (1 in above example)
# ${4} = BAM extension (text after the sample ID in BAM files)
# ${5} = minimum mapping quality score (30 should be good)
# ${6} = minimum mapping depth (iterate as necessary)
# ${7} = maximum amount of missing data (as percentage; e.g., 50)
# ${8} = merge distance for spanning low-depth sites (10 should be good)
# ${9} = minimum locus length (iterate as necessary)

cat ${1} | \
grep -v "^#" | cut -f ${2} | sort | uniq | \
while read experiment; \
do \
samples=`cat ${1} | grep -v "^#" | grep "${experiment}" | \
cut -f ${3} | awk -v ext="${4}" '{ print $1ext }' | paste -s -d " "`; \
cmd="samtools depth -Q 30 ${samples}"; \
eval ${cmd} | \
awk -v OFS="\t" -v depth="${6}" -v miss="${7}" '{ count=0; for (i=3;i<=NF;i++) \
{if ($i > depth) count+=1} \
{ if (count / (NF-2) >= (1 - (miss/100)))  print $0 } }' | \
awk -v OFS="\t" '{ $2=$2"\t"$2; print $0 }' | \
bedtools merge -d ${8} | awk -v loclen="${9}" '{ if ($3 - $2 >= loclen) print $0 }' \
> ${experiment}_rad_loci.qual${5}.dp${6}.miss${7}.merge${8}.len${9}.bed; \
done

```
