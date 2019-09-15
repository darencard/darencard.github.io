---
layout: posts
title: "Introductory RAD-seq Activity"
date: 2019-09-15
excerpt: "Experimental design activity to accompany an introductory lesson on RAD-seq"
---

## Introduction

Today, let's go through a brief experimental design activity to help demonstrate some of the features of RAD-seq and how they can be harnessed to explore numerous questions in ecology, evolutionary biology, and genetics. To perform this activity, we will use Harvard's [Odyssey compute cluster](https://www.rc.fas.harvard.edu/odyssey-3-the-next-generation/).

## Setting up Our Analysis Environment

We are going to use the very useful environments feature of [Anaconda Python](https://www.anaconda.com/distribution/) to manage the software dependencies we have for this activity. Anaconda is already installed on Odyssey and can be loaded as follows:

```bash
module load Anaconda3/5.0.1-fasrc02
```

Now we can create an analysis environment that will always be available for us to use for this activity, which will run using Python version 3.7.

```bash
conda create -n radseq python=3.7
```

This will take a few minutes to run and along the way, we will answer Yes (y) to allow installation. When it is done, we can activate our newly created environment.

```bash
conda activate radseq
```

We must also install some other important software in our environment, namely [Biopython](https://biopython.org/). For whatever reason, installing `biopython` through `conda` is not allowing later software to run properly for me, so we will install it using `pip` instead, which is Python's default package manager.

```bash
pip install --user biopython
```

We must also install [`bedtools`](https://bedtools.readthedocs.io/en/latest/), which is required for some of our software. This is available on Odyssey already, but let's install it into our `radseq` environment to demonstrate how to install using `conda`.

```bash
conda install -c bioconda bedtools
```

Now let's create a directory for our activity called `radseq_activity` and move into it so we can do our work there.

```bash
mkdir radseq_activity
cd radseq_activity
```

We are also going to be using some software that I wrote for performing *in silico* RAD-seq called [RADis](https://github.com/darencard/RADis). I have, unfortunately, not packaged this into easily-installable software, but we can still make use of it fairly easily. Let's download the GitHub repository for this purpose.

```bash
git clone https://github.com/darencard/RADis
```

One more thing we will need to do is add the software in this repository to our software `$PATH`, which is where the computer looks for software to run by default. Note that this is a temporary addition that will not persist beyond this session or if this directory is moved. There are ways to permanently add these software to your `$PATH` that I can show if anyone is interested.

```bash
cd RADis
export PATH=$PATH:`pwd`
cd ..
```

## *In Silico* RAD-seq Activity

In order to perform our *in silico* RAD-seq, we need a sequence for our digital restriction enzyme digestion. For time and ease, we will use chromosome 22 from the human genome, which is about 50 Mb in length. We can download this from the [UCSC Genome Browser](https://genome.ucsc.edu/) database.

```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz
gunzip chr22.fa.gz
```

We need to create a restriction enzymes file to run the *in silico* restriction digest. This software performs a double digest, so two enzymes need to be used. We are going to work in 4 groups, each with a different combination of enzymes as follows.

```bash
# group 1
sbfI	CCTGCAGG	6
sau3AI	GATC	0

# group 2
pstI	CTGCAG	1
sau3AI	GATC	0

# group 3
sbfI	CCTGCAGG	6
mspI	CCGG	1

# group 4
pstI	CTGCAG	1
mspI	CCGG	1
```

Each group should take their restriction enzymes and paste them into a new text file called `enzymes.txt`. This can be done using the terminal-based text editor `nano`. Note that there are `TAB` characters between each of the columns.

Once we have our reference genome (human chromosome 22) and our enzymes, we can run an *in silico* digestion with each enzyme using the following command. Sit back and be patient, as I'm not a professional software engineer and my program is slower than I would like (should take about 10 minutes at most).

```bash
restriction_digest_insilico.py --input chr22.fa --output chr22 --enzymes enzymes.txt
```

Once the command ends, you will see a couple new files in the working directory when we run `ls`.

```bash
ls
# group 1 output
# chr22_sbfI_CCTGCAGG_8.bed
# chr22_sau3AI_GATC_4.bed
```

Now that we have the individual restriction digest files for single enzymes, we can run the double digest. As part of this command, we need to also set the lower and upper limit on the fragment lengths we will select for. Given there are 4-5 people in a group, let's have each person select a different range of fragment sizes.

```bash
# person 1
200 - 400
# person 2
300 - 500
# person 3
400 - 600
# person 4
500 - 700
# person 5
600 - 800
```

Now we can run our double digest as follows. Be sure to use the output from the longer enzyme as the `-r` input. This will run much faster than the other command.

```bash
# example for group 1, person 1
extract_fragments.sh -g chr22.fa -r chr22_sbfI_CCTGCAGG_8.bed -c chr22_sau3AI_GATC_4.bed -l 200 -u 400 -o chr22_sbfI_sau3AI_200_400_output
```

You will see two outputs of this command.

```bash
ls
# example for group 1, person 1
# chr22_sbfI_sau3AI_200-400_output.bed.gz
# chr22_sbfI_sau3AI_200-400_output.fasta.gz
```

Both are gzipped, but we can easily unzip them.

```bash
gunzip *.gz
```

And look at their contents. One is a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file.

```bash
head chr22_sbfI_sau3AI_200_400_output.bed
# chr22	10520281	10520562	sbfI-sau3AI	281	+
# chr22	10527904	10528214	sbfI-sau3AI	310	-
# chr22	10528218	10528432	sbfI-sau3AI	214	+
# chr22	11054755	11054977	sbfI-sau3AI	222	-
# chr22	11265186	11265460	sbfI-sau3AI	274	+
# chr22	11438712	11439070	sbfI-sau3AI	358	+
# chr22	11553161	11553482	sbfI-sau3AI	321	+
# chr22	11716993	11717223	sbfI-sau3AI	230	+
# chr22	11953435	11953797	sbfI-sau3AI	362	-
# chr22	11956673	11957032	sbfI-sau3AI	359	-
```

The other is a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file.

```bash
head chr22_sbfI_sau3AI_200_400_output.fasta
# >chr22:10520281-10520562(+)
# ggtggggcttcactggggacccatcaccttccacccaggagcctgtctgcctcccacaaccatccatggcacccaggctgctggcatcaaggggcacttgcaggccagtgccaagccaccctcgtaccccctcatcttcccctcccatgctcctgctcctcagtgtccaaagtccagaaggggctgaggtggcaggggactgacatgtcagcactgcttccaatgtgtgcacacctggctgggcagtgacagcaccctgctgggtcccaaccccactctga
# >chr22:10527904-10528214(-)
# GGAAACCCCTTGAGGTCAAGACCCCACAATCAGACAAGGATGGAGTGGCTCACCCTCAGTAAACAGGCCAGACTCAAGGTGGTATAATGTCTTAACCAAGAGTGTGGGCCTCCAGGTCTGACTCCCATCTCAGTTCTCCTTTAATAACCACACTTTGATAATTCTCCTTAACAGGGGTTCCTggcaagtcagttcttcctcaggccttcggtttcctcacctacaagatgagagggctggaccagatgGAAATTCAGGGTGGAAGGGGATGTCCGCGCACAGCCCACCACCCCCCCCCAACAGGACCCTG
# >chr22:10528218-10528432(+)
# GGAAGAAGCAAATGGGTTTGCTTTCCTAGCTCTGTCCAGTATCTTAGGGACCCTGAGGACTGAAGAGATTCTTGTAGAGCCATCTGGTGTATGTCATGGGTGGGCCTTTTTTGAATGTCAGTCTGCCCAGTGAGCTGGCTCAGCCTGAATGAACTGTCTTGAATCTTTGGAGTTGTCTGTGTACTTTTAAGGGTTTCTCATCCTTGCACCAAAA
# >chr22:11054755-11054977(-)
# GGAGACCTTCTCCTGCCCGGACCAGAGAGGCTGCAGGGGTGCTGTGGGGAGAGGGCCAGATGGCCACAACACAGCTGGGAATGGCTCCAGAGGGGAGGGTGGAGAGTGCGCTGGCCCATGGGAGGACTCAGAGAGGAGTTGCGTGTGTTTGGCGCTTGTGGACATGGTGTGTGGCAGTGTGAGCCAACAGGCCAACGGGCAGCAAGCATAGGAGGCTTTATG
# >chr22:11265186-11265460(+)
# GGGCGGAGGTTCCCTAGGCAACGAGGGAGAGAGGGAGGGGTCTCCAGAAGGGAGAGACAGAAACCGGTTGCCCCAGGTTCGGTGAAGTCGGCCAGACCTCTCCCCGTGTCACCTCGACTTTCAATAACAGTGGCTGCTAGGTGATGCCCAAAGACAACCGATGCCTGCAAGTGTCAGTCAGCAGGGAAAAGAATGCatttatttatttatttatttatttatttagagacagagtctcactcaatctacagcccaggctgtagtgcagtggcgc
```

Now the critical piece of information we need for experimental design considerations is the number of fragments, as marker number/density enables different types of questions to be addressed. We can easily determine this by counting the number of lines in our BED file.

```bash
wc -l chr22_sbfI_sau3AI_200_400_output.bed
# 952
```

Or the number of sequences in our FASTA file.

```bash
grep -c ">" chr22_sbfI_sau3AI_200_400_output.fasta
# 952
```

Within your group, work together to determine how fragment number varies across different fragment size selections. What pattern do you see and why might this be important?

As a class, what do you notice about the number of fragments sampled based on the restriction enzymes and the fragment size windows used? What aspect of the restriction enzyme appears to impact fragment number?

Based on your research interests, think of an interesting question (or questions) you might want to ask and how genome-wide data could be useful. Assuming you are restricted to using RAD-seq only, what experimental design considerations must you consider to answer your research question?
