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
