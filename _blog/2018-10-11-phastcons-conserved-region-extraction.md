---
layout: posts
title: "Extracting conserved regions from PhastCons"
date: 2018-10-11
excerpt: "Extracting conserved regions sequences from PhastCons conservation tracks."
---

Commands for extracting conserved regions from 100-way PhastCons conservation tracks. Then BLAST can be used to extract orthologous regions from genome of interest for as many regions as possible.

```bash
# download the human gene annotations
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz

# convert human gene annotations to GTF file format
zcat refGene.txt.gz | cut -f 2- | genePredToGtf -utr file stdin stdout > refGene.gtf

# extract exon features from GTF file
cat refGene.gtf | awk '{ if ($3 == "exon") print $0 }' > refGene.exons.gtf

# download the human genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# download chromosome sizes file
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

# download MAF alignment files
mkdir maf
for i in `seq 8 22` M X Y; \
do \
rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/multiz100way/maf/chr${i}.maf.gz ./maf/; \
done

# download phastCons conservation tracks
mkdir phastCons100way
for i in `seq 8 22` M X Y; \
do \
rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/phastCons100way/hg38.100way.phastCons/chr${i}.phastCons100way.wigFix.gz ./phastCons100way/; \
done

# convert phastCons tracks into BedGraph
cd phas5Cons100way
for file in *.gz; \
do \
base=`basename ${file} .wigFix.gz`; \
echo ${base}; \
gunzip -c ${file} > ${base}.wigFix; \
wigToBigWig ${base}.wigFix ../hg38.chrom.sizes ${base}.bw; \
bigWigToBedGraph ${base}.bw ${base}.bedGraph; \
done

# download the annotation table for 100 vertebrate conserved elements
# 1. Go to http://genome.ucsc.edu/cgi-bin/hgTables
#
# 2. Settings:
# 	clade=Mammal
# 	genome=Human
# 	assembly=hg38
# 	group=Comparative Genomics
# 	track=Conservation
# 	table=100 Vert. El (phastConsElements100way
# 	output format=all fields from selected table
# 	output file=phastConsElements100way.txt
# 	**others**=defaults
# This downloads a BED-like table of conserved elements based on the 100-way vertebrate alignment.
# Here are some stats:
# 	count=10,350,729
# 	bases=162,179,256 (5.32%)
# 	smallest item=1
# 	average item=16
# 	biggest item=3,732
# 	smallest score=186
# 	average score=333
# 	biggest score=1,000

# extract conserved regions that fall outside coding exons (25% or less of the conserved region falls inside an exon)
cat phastConsElements100way.txt | \
tail -n +2 | cut -f 2- | \
awk '{ if ($3 - $2 >= 100) print $0 }' | \
bedtools intersect -wao -a stdin -b refGene.exons.gtf | \
awk -v OFS="\t" -F "\t" '{ print $0, $15 / ($3 - $2) }' | \
awk -F "\t" '{ if ($16 <= 0.25) print $0 }' | cut -f 1-5 \
> phastConsElements100way.noncoding.bed

# Here are some stats:
# 	count=138,744
# 	bases=25,180,979

# extract sequences for conserved regions from human genome
bedtools getfasta -fi hg38.fa -bed phastConsElements100way.noncoding.bed \
> phastConsElements100way.noncoding.fasta
```
