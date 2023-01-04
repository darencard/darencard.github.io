---
layout: posts
title: "Performing the Bonferonni and Benjamini–Hochberg Procedures on a Large Dataset"
date: 2022-10-22
excerpt: "How to control the false discovery rate without relying on R."
---

## Introduction

A common procedure in statistics when performing many statistical tests is to control for false positives using one of many procedures devised for doing so. I commonly use the Bonferonni or Benjamini–Hochberg procedures, depending on the situation. The easiest way to take advantage of these statistical approaches is to use R, which has several options for "correcting" p-values, transforming the raw p-value to a "corrected" form that minimizes false positives. In this circumstance, one considers tests with a corrected p-value below some threshold (usually 0.05) to be statistically significant.

While this procedure is straightforward, it relies on being able to load your dataset, containing raw p-values, into R. However, sometimes datasets are very large (i.e., 10s / 100s of megabytes or even gigabytes) and R does not play nice with very large files because of the memory requirements, which your computer may not support. So, how does one determine which statistical tests are significant without using R? Below, I outline how one can perform either the Bonferonni or Benjamini-Hochberg procedure to control for false positives using the Unix shell.

Dependencies: Following these procedures should not require any software installations, as it uses standard Unix utilities that should be available on your system. The one tool that may not be available is [datamash](https://www.gnu.org/software/datamash/), which should be easy to install on your system. However, you will need to install a tool for handling the genomics file formats we are working with, which is described when needed below.

## Background Information 

I found [a great overview](https://stats.libretexts.org/Bookshelves/Applied_Statistics/Book%3A_Biological_Statistics_(McDonald)/06%3A_Multiple_Tests/6.01%3A_Multiple_Comparisons) of both the Bonferonni and Benjamini-Hochberg procedures, and controlling for false positives overall, in a section of the free, online textbook [*Biological Statistics*](https://stats.libretexts.org/Bookshelves/Applied_Statistics/Book%3A_Biological_Statistics_(McDonald)) by [John H. McDonald (University of Deleware)](https://udel.edu/~mcdonald/).

In reading through the background information, I found that one does not need to strictly "correct" p-values. Rather, one can perform a procedure that will identify statistically significant tests based on the raw p-values. Moreover, the operations performed in these procedures are relatively straightforward, which makes Unix a great candidate for performing this task on a dataset that is too large for R.

## Preparing the Dataset

Let's start by gathering a dataset of reasonable scale that justifies why I am writing this blog post. We will go with this very large dataset of PhyloP scores produced from a [241-way alignment of mammalian genomes](https://cglgenomics.ucsc.edu/data/cactus/) produced as part of the [Zoonomia Project](https://zoonomiaproject.org/). [PhyloP scores](http://compgen.cshl.edu/phast/phyloP-tutorial.php) measure evolutionary conservation or acceleration of genomic sites based on an alignment of genomics regions and a phylogenetic tree. The file we will be using was produced using the procedure outlined in Example 1 [here](http://compgen.cshl.edu/phast/phyloP-tutorial.php). Because this test is performed across all sites in a 241-way alignment of mammals, a giant file is produced in [Wiggle format](https://genome.ucsc.edu/goldenPath/help/wiggle.html). This file was converted to [bigWig format](https://genome.ucsc.edu/goldenPath/help/bigWig.html) to compress the data for a lower disk footprint. We can retrieve this file with positions based on the human reference genome from the Zoonomia data repository.

```bash
# big file! this may take a while to download
wget https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/Homo_sapiens/241-mammalian-2020v2.bigWig
```

This bigWig file is ~21 GB. However, given it is compressed / binary, we cannot actually look at it or perform operations. We first need to convert from bigWig to a human readable format like Wiggle. Let's make this conversion using the tool [`bigWigToWig`](https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig) (warning: link directly downloads program) from UCSC (distributed as part of the [Kent utilities](https://hgdownload.soe.ucsc.edu/admin/exe/)). Let's look at the beginning of the file to get a sense of the format.

```bash
bigWigToWig 241-mammalian-2020v2.bigWig /dev/stdout | head
```

Here is what the output looks like.

```
#bedGraph section chr1:10074-11098
chr1	10074	10075	0.053
chr1	10075	10076	0.064
chr1	10076	10077	0.064
chr1	10077	10078	0.064
chr1	10078	10079	-2.109
chr1	10079	10080	0.053
chr1	10080	10081	0.053
chr1	10081	10082	0.064
chr1	10082	10083	0.064
```

Notice that there are lines beginning with `#`, which are comment lines that we can ignore. Then we have a pretty standard genomic coordinate format for each line indicating the chromosome / scaffold and start and end positions of each site in the human genome where a PhyloP score was calculated. The 4th column is the PhyloP score. This score requires some interpretation. First, it can be either positive or negative, where positive scores are considered conserved sites and negative scores are considered accelerated sites evolving faster than expected under neutral evolution. Scores of 0 encode sites evolving neutrally. Beyond that, the numeric score is not a p-value, but rather a -log10(p-value) from the likelihood ratio test performed in PhyloP. We can use the formula `10^-(|PhyloP|)` to convert to a raw p-value. The following command will classify sites as conserved, accelerated, or neutral and will perform the conversion to a p-value for all sites within the bigWig alignment file. We will write this to a tabular, plain text file.

```bash
bigWigToWig 241-mammalian-2020v2.bigWig /dev/stdout | grep -v "^#" | awk -v OFS="\t" '{ if ($4 > 0) print $0, "conserved", 10^-$4; else if ($4 < 0) print $0, "accelerated", 10^$4; else print $0, "neutral", 10^$4 }' > 241-mammalian-2020v2.parse.txt
```

We only added a couple additional columns of information, but between that and uncompressing the data from the bigWig format, we now have a file that is ~134 GB. Good luck loading that into R!

## Bonferroni Correction in Unix

We will use Unix to follow procedures for identifying sites with p-values that are statistically significant after controlling for false positives. First, let's perform the Bonferonni correction. Bonferroni correction is a simple procedure: the alpha value used in the statistical comparison (usually 0.05) is divided by the number of tests performed and this value is the new threshold for statistical significance. This is a quite conservative correction.

```bash
# this will take a while to run!
total=`cat 241-mammalian-2020v2.parse.txt | wc -l` && cat 241-mammalian-2020v2.parse.txt | awk -v OFS="\t" -v total="${total}" '{ if ($6 <= (0.05 / total)) print $0, "TRUE"; else print $0, "FALSE" }' > 241-mammalian-2020v2.parse.bonferroni.txt
```

The output was the same as the input with an additional binary column where sites that remain significant after Bonferroni correction are encoded as `TRUE` and non-signficant sites are encoded as `FALSE`. It is also possible to view the new Bonferroni threshold for statistical significance.

```bash
cat 241-mammalian-2020v2.parse.txt | wc -l | awk -v total="${total}" -v OFS="\t" '{ print total, 0.05 / total }'
# 2852623265	1.75277e-11
```

Pretty simple! That's a ton of statistical tests / sites and thus the Bonferonni threshold is quite low. Here it is with commas so it is more readable: 2,852,623,265. 2.85 billion!

## Benjamini-Hochberg Procedure in Unix

Now let's turn to the Benjamini-Hochberg (BH) procedure, which is more complex than the simple Bonferroni correction. But the BH procedure is less conservative and may detect statistical differences in comparisons that the Bonferroni does not. You can read more details about the procedure in the [Background Information](#background-information) above, but let's walk through the BH procedure simply as it is laid out there.

1. In the BH procedure, raw p-values are sorted from lowest to highest (a step that is not necessary in the Bonferroni correction).
2. Based on the sorted p-values, a rank is assigned to each test (i.e., site) from 1 to N where N is the number of tests.
3. A BH critical value is then calculated based on the ranks using the formula `(i/m)Q`, where `i` is the rank, `m` is the total number of tests performed (as is used in the Bonferroni correction), and `Q` is the false discovery rate you choose (usually 0.05).
4. Finally, the BH critical value is compared to the raw p-value. Specifically, we are looking for circumstances where the raw p-value is below the BH critical value. Many tests could meet this criterion but we are most interested in the *largest* p-value that is below its BH critical value. Tests with a raw p-value at or below this *largest* p-value are all considered statistically significant.

Again, if you need more information and a simple example with a far smaller dataset, visit the 
Information above.

We can use this procedure in our giant dataset by applying the following code. We are tabulating the total number of tests (i.e., sites) in the original file, sorting the raw p-values using 4 cores (`--parallel=4`) and the current working directory for temporary files (`-T $(pwd)`), calculating the ranks and BH critical values for each test (i.e., site), and identifying all sites where the raw p-value is below the BH critical value.

```bash
# this will take a while to run!
total=`cat 241-mammalian-2020v2.parse.txt | wc -l` && cat 241-mammalian-2020v2.parse.txt | sort -k6,6g --parallel=4 -T $(pwd) | awk -v OFS="\t" '{ print $0, NR }' |
awk -v OFS="\t" -v total="${total}" '{ print $0, ($7 / total)*0.05 }' |
awk -v OFS="\t" '{ if ($6 < $8) print $0, "TRUE"; else print $0, "FALSE" }' > 241-mammalian-2020v2.parse.fdr.txt
```

Unlike the Bonferroni procedure, this one will require multiple steps / commands. This first one produced a useful intermediate with all the information we need. Now we need to determine the largest raw p-value where the raw p-value is below the BH critical value.

```bash
# this will take a while to run!
cat 241-mammalian-2020v2.parse.fdr.txt | awk '{ if ($9 == "TRUE") print $0 }' | tail -n 1 | awk '{ print $6 }'
# 0.00229615
```

Now we can gather those tests (i.e., sites) where the raw p-value is lesser than or equal to 0.00229615. We will also re-sort the data by chromosome / scaffold and coordinates, which is more logical (using the same tricks as before to run with 4 cores).

```bash
cat 241-mammalian-2020v2.parse.fdr.txt | awk '{ if ($6 <= 0.00229615) print $0 }' |
sort -k1,1 -k2,2n --parallel=4 -T $(pwd) > 241-mammalian-2020v2.parse.fdr.significant.txt
```

## Summarizing the Results

Remember, the Bonferroni correction is far more conservative than the BH FDR procedure, so we expect far fewer sites to be identified as statistically significant. Let's see how that shakes out.

```bash
# number of statistically significant sites after Bonferroni correction
wc -l 241-mammalian-2020v2.parse.fdr.significant.txt
# 131124451 241-mammalian-2020v2.parse.fdr.significant.txt
# number of statistically significant sites after performing the BH procedure
cat 241-mammalian-2020v2.parse.bonferroni.txt | awk '{ if ($7 == "TRUE") print $0 }' > 241-mammalian-2020v2.parse.bonferroni.significant.txt
cat 241-mammalian-2020v2.parse.bonferroni.significant.txt | wc -l
# 175862
```

Just as we expected: 131 million sites vs. only 176 thousand sites. A very large difference! Let's dig further into the dataset of sites that are statistically significant after running the Benjamini-Hochberg procedure. How many sites are conserved versus accelerated? We can use `datamash` to investigate this in Unix without needing to use R, which would be difficult anyway (`241-mammalian-2020v2.parse.fdr.significant.txt` is 9.6 GB).

```bash
cat 241-mammalian-2020v2.parse.fdr.significant.txt | datamash --sort groupby 5 count 5
# accelerated	50558458
# conserved	80565993
```

We can even go further and summarize based on chromosome as well.

```bash
cat 241-mammalian-2020v2.parse.fdr.significant.txt | datamash --sort groupby 1,5 count 1
# chr1	accelerated	4215264
# chr1	conserved	7427904
# chr10	accelerated	2466657
# chr10	conserved	3769684
# chr11	accelerated	2243972
# chr11	conserved	4151872
# chr12	accelerated	2211997
# chr12	conserved	3893533
# chr13	accelerated	1560650
# chr13	conserved	2161264
# chr14	accelerated	1569982
# chr14	conserved	2936207
# chr15	accelerated	1464326
# chr15	conserved	2706121
# chr16	accelerated	1929270
# chr16	conserved	2821175
# chr17	accelerated	1791503
# chr17	conserved	3377866
# chr18	accelerated	1334192
# chr18	conserved	1944436
# chr19	accelerated	1393935
# chr19	conserved	1924627
# chr2	accelerated	3946703
# chr2	conserved	7261033
# chr20	accelerated	1385740
# chr20	conserved	1902183
# chr21	accelerated	921060
# chr21	conserved	703765
# chr22	accelerated	1145298
# chr22	conserved	1114539
# chr3	accelerated	2798629
# chr3	conserved	5486505
# chr4	accelerated	2850552
# chr4	conserved	3928253
# chr5	accelerated	2804820
# chr5	conserved	4893091
# chr6	accelerated	2753938
# chr6	conserved	4439519
# chr7	accelerated	2744352
# chr7	conserved	4098364
# chr8	accelerated	2454247
# chr8	conserved	3322716
# chr9	accelerated	2027966
# chr9	conserved	3732695
# chrX	accelerated	2253501
# chrX	conserved	2322307
# chrY	accelerated	289904
# chrY	conserved	246334
```

Pretty easy and cool! And that's essentially it! These commands could easily be put together into a little pipeline for files like this but I will not do that here. The amazing thing is that this analysis can be done in an hour or two of interactive coding on the terminal even though we are dealing with files that are 100s of gigabytes in size! The Unix shell does not really blink but we would be helpless with R. I hope this tutorial is useful to others out there working with genomics or other extremely large statistical datasets!
