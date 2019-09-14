---
layout: posts
title: "Determining whether RNA-seq data is stranded"
date: 2019-09-03
excerpt: "How-to guide on quickly determining whether RNA-seq data is stranded"
---

### Determining Whether RNA-seq Data is Stranded or Unstranded

I have a mixture of RNA-seq data including newly generated data that I produced and some older data from other researcher. When sending my RNA for library preparation, I elected to have stranded RNA-seq libraries produced. However, you never know what can happen when you make intentions known to a sequencing core and I have no idea what these inherited libraries are. [Here](https://www.ecseq.com/support/ngs/how-do-strand-specific-sequencing-protocols-work) is a nice overview of stranded vs. unstranded RNA-seq libraries and data.

So I took a look around for an answer and found this [issue](https://github.com/alexdobin/STAR/issues/258) on the [STAR RNA-seq mapper GitHub repository](https://github.com/alexdobin/STAR). If you have stranded data, all (or almost all) of you reads will map to one of the DNA strands versus mapping at about 50/50 ratios to each strand when you have unstranded data. This will look something like this:

![Stranded/Unstranded DNA](https://cloud.githubusercontent.com/assets/631218/25322008/39407a7c-28f6-11e7-9434-1852dcbacae0.png)

So I decided to take a look at the two different groups of libraries I had. I quickly subsetted both my BAM mapping file and my reference down to the first 1 Mb of the first scaffold, to make it easier. I then visualized the mappings in [IGV](https://software.broadinstitute.org/software/igv/) and turned on the option `Colour Alignments By -> First-of-Pair Strand` by right clicking on the data track.

My genome is not yet annotated, so I can't focus on regions with genes specifically, so I just looked for regions with decent pileups of reads, like the following:

![My Data](https://github.com/darencard/darencard.github.io/raw/master/assets/images/blog/igv_stranded_unstranded.png)

It is not quite as clean as the example data, but there is an obvious bias in the top sample towards one strand or the other. In the bottom sample, it is a pretty equal mixture of strands. This indicates that my new data are stranded (confirming what I expected to see) and that the older data are not stranded.

### Splitting BAM into Sense and Antisense Strand BAMs

It is also sometimes necessary to take an existing mapping BAM from stranded RNA-seq data and split it into two BAM files for the sense ("plus") and antisense ("minus") strands of DNA. That is pretty easily accomplished using `samtools`. To do this, I ended up using the following bash script:

```bash
#!/usr/bin/env bash

# Get the bam file (1st argument) and output directory (2nd argument) from the command line
BAM=$1
TARGET_D=$2

FILE=$(basename $BAM)
NAME=$(basename $BAM)
BAMF1=${TARGET_D}/${NAME}.fwd1.bam
BAMF2=${TARGET_D}/${NAME}.fwd2.bam
BAMF=${TARGET_D}/${NAME}.fwd.bam
BAMR1=${TARGET_D}/${NAME}.rev1.bam
BAMR2=${TARGET_D}/${NAME}.rev2.bam
BAMR=${TARGET_D}/${NAME}.rev.bam

# Forward strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse strand
#
# 0x1 - paired
# 0x2 - properly paired
# 0x20 - partner on reverse strand
# 0x40 - read one
# FLAGs 0x1 + 0x2 + 0x20 + 0x40 = 0x63 = 99 in decimal
samtools view -bh -f 99 $BAM > $BAMF1
samtools index $BAMF1
# 0x1 - paired
# 0x2 - properly paired
# 0x10 - on reverse strand
# 0x80 - read two
# FLAGs 0x1 + 0x2 + 0x10 + 0x80 = 0x93 = 147 in decimal
samtools view -bh -f 147 $BAM > $BAMF2
samtools index $BAMF2

#
# Combine alignments that originate on the forward strand.
#
samtools merge -f $BAMF $BAMF1 $BAMF2
samtools index $BAMF

# Reverse strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#

# 0x1 - paired
# 0x2 - properly paired
# 0x10 - reverse strand
# 0x40 - read one
# FLAGs 0x1 + 0x2 + 0x10 + 0x40 = 0x53 = 83 in decimal
samtools view -bh -f 83 $BAM > $BAMR1
samtools index $BAMR1
# 0x1 - paired
# 0x2 - properly paired
# 0x30 - partner on reverse strand
# 0x80 - read two
# FLAGs 0x1 + 0x2 + 0x20 + 0x80 = 0xA3 = 163 in decimal
samtools view -bh -f 163 $BAM > $BAMR2
samtools index $BAMR2

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -f $BAMR $BAMR1 $BAMR2
samtools index $BAMR
```

I modified this slightly from the script presented in [this](https://josephcckuo.wordpress.com/2016/11/18/splitting-reads-in-a-bam-file-by-strands/) excellent blog post. But I also found some other useful information on this process [here](https://www.biostars.org/p/348134/#348350), [here](https://www.biostars.org/p/92935/), and [here](http://seqanswers.com/forums/showthread.php?t=29399#7). The trickiest part is understanding what stranded libraries are/mean (see above) and the flags used in the SAM format to indicate mapping type - [this tool](https://broadinstitute.github.io/picard/explain-flags.html) can be used to help with the latter.
