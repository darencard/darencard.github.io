---
layout: posts
title: "A tutorial on whole genome alignment"
date: 2019-11-01
excerpt: "A thorough tutorial on performing whole genome alignment with the Lastz/MultiZ pipeline"
---

# Note: This blog post/tutorial is in progress and will be updated as data analysis proceeds. At this point, very little will be useful to anyone. Please check back later for a hopefully complete post.

## links to integrate
lastz documentation: http://www.bx.psu.edu/~rsharris/lastz/README.lastz-1.04.00.html
Michael Hiller's software page (lots of useful tools for WGA): https://www.mpi-cbg.de/research-groups/current-groups/michael-hiller/software-and-data/
Hiller lab Github for WGA: https://github.com/hillerlab/GenomeAlignmentTools
144 vertebrate alignment paper: https://academic.oup.com/nar/article/45/14/8369/3875570
Human to zebrafish alignment paper: https://academic.oup.com/nar/article/41/15/e151/2411468
chainCleaner paper: https://academic.oup.com/bioinformatics/article/33/11/1596/2929344
*new* RepeatFiller preprint: https://www.biorxiv.org/content/10.1101/696922v1.full
Example of WGA in a reptile: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177708
bx-python software for MSAs: https://github.com/bxlab/bx-python
50 bird genomes alignment scripts: https://github.com/gigascience/paper-zhang2014/blob/master/Whole_genome_alignment/pairwise/bin/lastz_CNM.pl
UCSC WGA wiki: http://genomewiki.ucsc.edu/index.php/Whole_genome_alignment_howto
Bioconductor WGA: https://bioconductor.org/packages/release/bioc/vignettes/CNEr/inst/doc/PairwiseWholeGenomeAlignment.html
TBA guide (MultiZ): https://www.bx.psu.edu/miller_lab/dist/tba_howto.pdf
Toast/Roast guide: http://www.bx.psu.edu/~cathy/toast-roast.tmp/README.toast-roast.html

## preparation
#### SOFT mask reference genomes for repeats

#### create .2bit and .sizes files
faToTwoBit target_genome.fasta target_genome.2bit
faToTwoBit query_genome.fasta query_genome.2bit
faSize target_genome.fasta -detailed target_genome.sizes
faSize query_genome.fasta -detailed query_genome.sizes

#### split query into reasonable size number of chunks (i.e., 20)
split_fasta.pl query_genome.fasta query_genome_chunks/ 20 avg yes

## alignment
#### run the lastz alignment on each query chunk against the target genome
lastz lastz_parameters target_genome.2bit query_genome_chunks/query_genome_chunk_N.2bit > chunk_N_output.axt
K=2400
L=3000
Y=9400
H=2000
default scoring matrix

## chaining and processing
#### chain each lastz alignment
axtChain ?chain_parameters? chunk_N_output.axt target_genome.2bit query_genome_chunks/query_genome_chunk_N.2bit > chunk_N_output.chain

#### merge and sort all the chains
chainMergeSort chunk_1_output.axt ... chunk_N_output.axt > all_output.chain

#### patchchain for highly sensitive local alignments
patchChain.perl all_output.chain target_genome.2bit query_genome.2bit target_genome.sizes query_genome.sizes [?parameters?]

#### repeatfiller (may need to do this on chunks before merging; but given it relies on colinear chains, probably best to run on whole genome)
RepeatFiller.py -c all_output.chain -T2 target_genome.2bit -Q2 query_genome.2bit

#### chainCleaner to remove alignments that obscure the evolutionary history of genomic rearrangements
chainCleaner all_output.chain target_genome.2bit query_genome.2bit -tSizes target_genome.sizes -qSizes query_genome.sizes all_output_clean.chain removedSuspects.bed -linearGap=loose

## netting
#### filter the chains to remove chains that do not have a chance of being netted
chainPreNet all_output_clean.chain target_genome.sizes query_genome.sizes all_output_sort.chain

#### UCSC chainNet to form nets
chainNet all_output_sort.chain target_genome.sizes query_genome.sizes target_genome.prenet query_genome.prenet

#### Hiller chainNet to form nets (will probably use this one)
chainNet -rescore all_output_sort.chain target_genome.sizes query_genome.sizes target_genome.prenet query_genome.prenet

#### add synteny information
netSyntenic target_genome.prenet target_genome.net
netSyntenic query_genome.prenet query_genome.net

#### convert to axt format and sort (not sure why; perhaps only way to get to maf from net?)
netToAxt target_genome.net all_output_sort.chain target_genome.2bit query_genome.2bit all_output.axt
axtSort all_output.axt all_output_sort.axt

## multiz alignment
#### convert axt to maf format, which is required for MultiZ alignment
axtToMaf -tPrefix=target_name -qPrefix=query_name all_output_sort.axt target_genome.sizes query_genome.sizes all.maf

#### create the multiple alignment of all maf files
tba "newick_tree" *.maf tba_out.maf

#### projecting the multiple alignment onto a reference sequence (for visualization)
maf_project tba_out.maf "species"
