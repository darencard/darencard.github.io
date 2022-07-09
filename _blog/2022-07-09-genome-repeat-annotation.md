---
layout: posts
title: "Repeat Annotation using RepeatModeler and RepeatMasker"
date: 2022-07-09
excerpt: "Thorough guide for annotating and masking repeats using RepeatModeler/RepeatMasker."
---

## Introduction

Repeat annotation and masking is a critical step in genome annotation, as it illuminates an important driver of genome evolution and cuts down on the number of spurious gene models in subsequent gene annotation steps of genome assembly, annotation, and curation. However, I have noticed anecdotally that many researchers put less effort into annotating repeats relative to gene models and often only perform very basic annotations based on existing repeat libraries. However, there are numerous ways to annotate repeats in genomes and tools seem to be multiplying quickly. I have previously outlined some of the steps I take for annotating repeats as part of my guide on [genome annotation using MAKER](https://darencard.net/blog/2017-05-16-maker-genome-annotation/), but I have made some changes to my approach more recently and thought a stand-alone guide on repeat annotation might be useful to others. This guide will outline *de novo* repeat identification using RepeatModeler with subsequent annotation and masking using RepeatMasker. These two pieces of software are probably the most widely used in repeat annotation and have a long (10+ year) track record in genomics research. That said, the same general approach could be followed using alternative repeat identification/annotation tools and new tools may eventually emerge that outperform RepeatModeler and RepeatMasker, so it is always good to keep an eye on popular software in this area.

## Software & Data

#### Software prerequisites:
1. [RepeatModeler2](http://www.repeatmasker.org/RepeatModeler/): RepeatModeler was recently updated to version 2 and a formal [publication](https://www.pnas.org/doi/10.1073/pnas.1921046117) on the software now exists (it never did before). This guide will work with the various version 1s of RepeatModeler, but using the most updated software is likely best. This guide is based on RepeatModeler v. 2.0.2a.
2. [RepeatMasker](http://www.repeatmasker.org/RMDownload.html): While a paper has never been published about the software, RepeatMasker has become very popular for identifying and masking repeats. It is basically a wrapper around BLAST and other tools which use sequence similarity with known repeats (provided as a "library") to identify repeats in a genome. This guide is based on RepeatMasker v. 4.1.2-p1.
3. Numerous dependencies required by RepeatModeler and RepeatMasker: See the installation guides for these software for more details. I used the most up-to-date versions of these dependencies. Details: [RECON](http://eddylab.org/software/recon/) v. 1.08, RepeatScout v. 1.0.6, [TRF](https://tandem.bu.edu/trf/trf.html) 4.0.9, [RMBlast](http://www.repeatmasker.org/RMBlast.html) v. 2.11, [CD-HIT](http://weizhong-lab.ucsd.edu/cd-hit/) v. 4.8.1, [MAFFT](https://mafft.cbrc.jp/alignment/software/) v. 7.481, [GenomeTools LtrHarvest](http://genometools.org/pub/) v. 1.5.9, [Ltr_retriever](https://github.com/oushujun/LTR_retriever) v. 2.9.0 and [NINJA](https://github.com/TravisWheelerLab/NINJA) v. 0.95-cluster_only. I use the BLAST module for identifying repeats in RepeatMasker but this guide should work if you decide to use alternatives like CrossMatch or HMMER instead. 

Both RepeatModeler and RepeatMasker can be a little complex to install, but the instructions for doing so are pretty clear. Generally, each dependency must be pre-installed and then when RepeatModeler/Masker is configured, you provide the paths to each dependency so everything can be properly setup. See the detailed instructions for all software for doing this. I am not involved with developing or maintaining any of this software and, therefore, cannot provide support.

#### Raw data/resources:
1. `reference-genome.fasta`: A reference genome for the species of interest in FASTA format.
2. Existing, curated repeat libraries sourced from an appropriate database. A good solution fot his are the database releases from [RepBase](http://www.girinst.org/repbase/), and I used database release version 3.3 (date: 2020-11-09) of this resource for this guide. Note that this database is now unfortunately behind a paywall, so it may take some more effort to source what you need.

#### Running commands:
Our campus HPC system uses SLURM for managing jobs and gathering STDOUT and STDERR streams from any of the commands I am running. But if you do not have a well-managed HPC system, I still have some recommendations. I have become accustomed to running most programs (especially those that take hours/days) using `screen`. I would create new screens for each command below. I also like to use `tee` to keep track of run logs. I will demo these `tee` commands below for those not on an HPC system and those who are can leave off these portions of the commands to have SLURM take care of everything.

## *De Novo* Repeat Identification

The first, and very important, step to genome annotation is identifying repetitive content. Existing libraries from Repbase or from internal efforts are great, but it is also important to identify repeats *de novo* from your reference genome using `RepeatModeler`. This is pretty easy to do and normally only takes a couple days using 8-12 cores.

```bash
# make a directory for storing logs
mkdir -p logs
# build new RepeatModeler BLAST database with a name that includes an ID (e.g., a species code, specimen ID, etc.) 
# and genus/species. Modify accordingly.
BuildDatabase -name ID_Genus-species -engine ncbi reference-genome.fasta
# now run RepeatModeler with 16 cores and send results from STDOUT and STDERR streams to 1_repeatmodeler.log
# in my experience, this command takes 1-3 days with vertebrate genomes
RepeatModeler -pa 16 -engine ncbi -database ID_Genus-species 2>&1 | tee 00_repeatmodeler.log
```

All RepeatModeler output is normally placed into a directory `RM_<dateinfo>` where `<datainfo>` includes the date/time of the run (a way of creating a unique directory for the results). You will see subdirectories for each round of RepeatModeler that is run as well as several files with the prefix `reference-genome` (varies depending on your file name). The most relevant outputs are the final, classified repeat consensus libraries in FASTA and Stockholm format: `reference-genome-families.fa` and `reference-genome-families.stk`, respectively. I do not know much about the Stockholm file and the rest of the guide will make use of the library in FASTA format. The FASTA library contains important repeat classification information in the header lines, which look something like this: `>rnd-1_family-174#LINE/L1`. A unique round and family identifier are included from the internal RepeatModeler run and after the `#`, repeat classification information is provided. This example element is an L1 LINE element. You will see consensus sequences for other repeat types as probably a fair number of unknown repeat elements that have the form `>rnd-2_family-408#Unknown`. Consensi from the first round of RepeatModeler, which uses RepeatScout for identification, will also have additional information in the form `( RepeatScout Family Size = 445, Final Multiple Alignment Size = 100, Localized to 501 out of 508 contigs )` while those in rounds 2 and later, which are based on RECON, instead have the form ` ( Recon Family Size = 19, Final Multiple Alignment Size = 19 )`.

I typically like to do a couple of things with the raw FASTA repeat library from RepeatModeler. First, it is best to add some sort of species/genome identifier to each consensus sequence ID so we can keep things straight and include some minor curation. This is especially useful if you are assessing repeats across multiple reference genomes in a comparative manner, which is becoming increasingly common as genomes quickly proliferate. I normally add some sort of species code. One that has become a bit of a standard in the comparative genomics community is to use a 6-letter species code with a number for the version of the reference genome. The 6-letter species code is based on the first 3 letters of the genus and species and for a new reference genome, the genome version will be 1. For human, for example, such a prefix may be `homSap1`. Here is a command that will add this prefix to each FASTA header line using the useful software `seqkit` and `awk`.

```bash
cat reference-genome-families.fa | seqkit fx2tab | awk '{ print "abcDef1_"$0 }' | seqkit tab2fx > reference-genome-families.prefix.fa
```

The above L1 element will now have a header line like this: `>abcDef1_rnd-1_family-174#LINE/L1`.

I also like to split my library into elements that were successfully classified and those that remain as unclassified or unknown elements. Here are some commands for doing so.

```bash
cat reference-genome-families.prefix.fa | seqkit fx2tab | grep -v "Unknown" | seqkit tab2fx > reference-genome-families.prefix.fa.known
cat reference-genome-families.prefix.fa | seqkit fx2tab | grep "Unknown" | seqkit tab2fx > reference-genome-families.prefix.fa.unknown
```

Then it is possible to quantify the number of known vs. unknown repeat consensi. RepeatModeler does its best to classify repeats using existing info from Repbase and other features of the repeat consensi, but you may still encounter several unknown elements. They can be especially prevalent when working with species from lineages not well represented among existing reference genomes. Therefore, the squamate reptiles I most commonly study normally have 50% or more unknown repeat consensi, but if you work on mammals and birds, you may get better results. Here is a command for quantifying known and unknown repeat elements.

```bash
# quantify number of classified elements
grep -c ">" reference-genome-families.prefix.fa.known
# quantify number of unknown elements
grep -c ">" reference-genome-families.prefix.fa.unknown
```

Further steps can be taken to annotate unknowns in the resulting library. This could be its own guide and generally devolves into a fair amount of manual curation, so I will only briefly provide some suggestions. First, a recent [paper](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-021-00259-7) lays out a lot of great advice for curating repeat elements and I recommend reading this paper and working to implement the suggestions (something I am still doing myself). 

I have built one tool that could also be useful, though I have not published on it and should do some thorough validation before anyone relies on it. I will describe it briefly here but please note that your mileage may vary and you should thorougly evaluate your results before relying on them for downstream inference.

My approach leverages RepeatMasker as an engine for sequence similarity searches. Given RepeatMasker is good at identifying repeats in a genome, I made the assumption that it will also do a decent job identifying unknown element consensus sequences if supplied with a library of known elements. These known elements could come from an existing database like Repbase, which is nicely incorporated into RepeatMasker. However, RepeatModeler is likely already using these same resources to initially identify elements, so this is not likely to make a large difference. Classified elements could also come from other curated libraries, including the known/classified elements identified by RepeatModeler in the same species. Or combined libraries of known elements from a clade of organisms, which are becoming more common.

With this in mind, my approach is iterative and uses one library of "known" elements to classify another library of "unknown" elements. I have implemented it in the program `repclassifier`, which is available [here](https://github.com/darencard/GenomeAnnotation/blob/master/repclassifier). See the beginning of the script for details about dependencies and how to run the software, but it basically has two modes. One mode uses the Repbase database and RepeatMasker to try to identify unknown elements with sequence similarity to curated repeat elements in Repbase. The other mode does the same thing but with a custom library, supplied as a FASTA file, instead of the Repbase database. I normally use both modes to perform an initial classification of unknowns.

```bash
# classifying unknowns (-u): run with 3 threads/cores (-t) and using the Tetrapoda elements (-d) from Repbase 
# and known elements (-k) from the same reference genome; append newly identified elements to the existing known 
# element library (-a) and write results to an output directory (-o)
repclassifier -t 3 -d Tetrapoda -u reference-genome-families.prefix.fa.unknown \
-k reference-genome-families.prefix.fa.known -a reference-genome-families.prefix.fa.unknown \
-o round-1_RepbaseTetrapoda-Self
```

The main results in `round-1_RepbaseTetrapoda-Self` are the files `round-1_RepbaseTetrapoda-Self.unknown`, which are those elements that remain unknowns, and `round-1_RepbaseTetrapoda-Self.known`, which are a combination of the initial known repeat library provided for classification (`reference-genome-families.prefix.fa.known`) and the newly classified elements from this first round of classification. But we can also continue running `repclassifier` iteratively using the newly annotated known elements to try to classify more unknown elements. I find that doing this for a few rounds allows many unknown elements to be identified early on leading to a pattern where the number of unkowns decreases quickly at first before plateauing with no additional gains after a few rounds. Elements annotated across multiple species could also be combined together as a library of known elements, which might also lead to gains in the number of known elements. We can run an artibrary number of additional rounds of classification from 2 to N using a command like the following.

```bash
# classifying unknowns (-u): run with 3 threads/cores (-t) and using only the known elements (-k) from the 
# same reference genome; append newly identified elements to the existing known element library (-a) and 
# write results to an output directory (-o). No Repbase classification is used here.
repclassifier -t 3 -u round-1_RepbaseTetrapoda-Self/round-1_RepbaseTetrapoda-Self.unknown \
-k round-1_RepbaseTetrapoda-Self/round-1_RepbaseTetrapoda-Self.known \
-a round-1_RepbaseTetrapoda-Self/round-1_RepbaseTetrapoda-Self.known -o round-2_Self
```

No matter whether/how you decide to classify unknown elements, I still recommend separating the known and unknown elements and annotating using these libraries separately (see below). This gives us the advantage of classifying any elements with sequence similarity to a known element before we introduce unknown elements that may have more similarity to certain repeats but include no functional information about the type of repeat. While it is best to identify/mask unknown elements, these annotations are not useful when characterizing repeat landscapes or making comparative assessments, so prioritizing known annotations is worthwhile. There are probably some reasonable arguments for annotating knowns and unknowns together, which cuts down on the rounds of RepeatMasker, but I have had good success annotating knowns and unknowns separately.

## Full Repeat Annotation and Masking

RepeatMasker, and repeat element identification/masking in general, is great because it can be run serially using several libraries with different levels of priority to annotate repeats in a flexible manner. By "run serially", I mean that the results of one round of repeat annotation/masking with RepeatMasker, which results in a genome with those elements masked as Ns, can be fed into a sebsequent round to annotate distinct elements that were not already identified. I like to take advantage of this when annotating repeats rather than just combining a bunch of repeat element consensi together into a single library for one round of annotation/masking. I will lay out this approach below to demonstrate the utility of serially annotating repeats.

I typically serially annotate/mask genomes with RepeatMasker in the following order:

1. Simple repeats
2. Well-curated repeats from existing databases like Repbase for the target organism
3. Known species-specific repeats from RepeatModeler
4. Unknown species-specific repeats from RepeatModeler

Note that more rounds of masking where additional libraries are used could be performed. For example, if one is annotating a bunch of in-group genomes, it may make sense to produce a combined library of consensus repeat elements from across the species and perform an additional round of annotation/masking based on this clade-level library. Alternatively, if one knows of existing artifacts in a repeat library/database that could be problematic or wishes to take great care in annotating certain elements, they could run an initial round of annotation/masking that addresses these issues before performing broader repeat annotations. For example, previous findings highlighted chimeric CR1 and BovB LINE elements in Repbase that may lead to erroneous classifications of repeats in squamate reptiles. By curating our own CR1/BovB LINE library and running an initial round of RepeatMasker with this library, we can overcome this issue by correctly identifying elements with our new library and masking them so that they are not erroneously annotated in subsequent rounds of repeat annotation/masking.

It is also worth briefly describing how repeat masking is performed. There are two types of sequence masking, which are commonly used to indicate stretches of the genome that are repetitive. The first is soft masking, which is where the features in question (repeats in this, and most, cases) are converted from uppercase nucleotides to lowercase nucleotides. Hard masking, on the other hand, is where these stetches of nucleotides are instead converted to an `N`, which stands for any, or totally ambiguous, nucleotide. Less information is lost during soft masking since the nucleotide identities are still present, but the lower case letters tell software that the sequence is distinct from other regions. Usually, soft masking is used for simple repeats, since these repeats can be pretty common in protein-coding regions of the genome and it makes sense to retain nucleotide identities while identifying a region as repetitive. Software that is designed to properly interpret this masking, such as MAKER, will treat these regions slighly differently. For example, when annotating genes, MAKER will allow gene models to be started only in non-masked regions of the genome (capital nucleotides) but not in regions that are masked in anyway. However, MAKER will allow gene models to be extended into regions with soft masking but not hard masking. This example hopefully helps you to understand why soft and hard masking are different and important.

With this all in mind, let's run through this serial annotation/masking pipeline. We will be running RepeatMasker multiple times but the commands are quite similar and it is worth outlining the details all at once. Here is a generic RepeatMasker command.

```bash
RepeatMasker -pa 16 -a -e ncbi -dir <out_dir> < ... other arguments... > reference-genome.fasta
```

The above aspects of the commands are the same for all runs. We are running RepeatMasker with 16 cores/threads using NCBI BLAST as the search engine for identifying similar elements. The output will go into a designated output directory and we are using the `-a` flag to write repeat alignments as a `.align` output file, which allows for some interesting downstream analyses (see below). For our first round of RepeatMasker, we will use the `reference-genome.fasta` file that we used to identify repeats using RepeatModeler, but keep in mind that in subsequent rounds we will tecnically be feeding in a different file where repeats from previous rounds are masked. The major difference between each run of RepeatMasker is included in the `< ... other arguments ... >` portion of the command, where one of three major options can be used. Runtime can vary quite a bit based on the reference genome and repeat library used, but the commands I outline below typically do not run more than a couple of days (sometimes they finish in minutes).

1. `-noint -xsmall`: We can annotate/mask only simple repeats using these arguments. These repeats will be soft masked and in subsequent rounds of RepeatMasker, the software may incorporate some of these simple repeats into broader annotations of complex, interspersed repeats and hard mask them. Simple repeat masking is performed whenever RepeatMasker is run but I have begun to do this all in one up-front step to save on time in subsequent rounds of annotation. This also helps keep simple and complex, interspersed repeats (e.g., transposable elements) separate, which helps because these repeat elements are often handled/masked differently.
2. `-nolow -species tetrapoda`: We can annotate/mask only complex, interspersed repeats (e.g., transposons) by setting the `-nolow` argument. In this case, we are using repeat consensus sequences from Tetrapoda included in the Repbase database to annotate these repeats. This taxonomic group could be modified to encompass other broader or narrower lineages.
3. `-nolow -lib round-1_RepbaseTetrapoda-Self/round-1_RepbaseTetrapoda-Self.known`: Again, we use the `-nolow` argument to annotate/mask only complex, interspersed repeats. However, in this case we are avoiding Repbase and instead providing our own custom repeat library in FASTA format. In this example I am supplying the known elements after one round of `repclassifier` classification but this file can be any library you would like to use.

Now, let's see how we use the above commands and information to run our serial repeat annotation/masking pipeline. First, let's create output directories for the results of each round of annotation/masking and a directory for log files.

```bash
mkdir -p logs 01_simple_out 02_tetrapoda_out 03_known_out 04_unknown_out
```

Now we can run our annotation/masking of the simple repeats in our reference genome.

```bash
# round 1: annotate/mask simple repeats
RepeatMasker -pa 16 -a -e ncbi -dir 01_simple_out -noint -xsmall reference-genome.fasta 2>&1 | tee logs/01_simplemask.log
```

RepeatMasker will initially create outputs and process them in a directory named `RM_<dateinfo>` where `<datainfo>` includes the date/time of the run (a way of creating a unique directory for the results), so do not mess with this! Upon finishing, the results will automatically be copied to the designated output directory `01_simple_out`. Five output files are created.

1. `reference-genome.fasta.align`: Resulting repeat alignments
2. `reference-genome.fasta.cat.gz`: Full RepeatMasker results (file can get large)
3. `reference-genome.fasta.masked`: The masked reference genome. So, a copy of the original reference genomes with repeat elements soft masked (for this round)
4. `reference-genome.fasta.out`: A tabular list of all repeat regions/elements annotated in this round of annotation/masking
5. `reference-genome.fasta.tbl`: A summary table of the overall composition of the genome based on major repeat element groupings

I like to rename the outputs of each round of RepeatMasker so that it is more clear what I have when looking at the file names.

```bash
# round 1: rename outputs
rename fasta simple_mask 01_simple_out/reference-genome*
rename .masked .masked.fasta 01_simple_out/reference-genome*
```

We could look at the summary table and begin to understand our repeat landscape, but remember, we have only masked the simple repeats. Let's continue our serial repeat annotation/masking pipeline and summarize everything at the end.

We will next mask Tetrapoda repeats from Repbase and then species specific known and unknown repeats, respectively, using three more rounds of RepeatMasker. Below, I will simply provide the commands for this example for each of these rounds.

```bash
# round 2: annotate/mask Tetrapoda elements sourced from Repbase using output from 1st round of RepeatMasker
RepeatMasker -pa 16 -a -e ncbi -dir 02_tetrapoda_out -nolow \
-species tetrapoda 01_simple_out/reference-genome.simple_mask.masked.fasta 2>&1 | tee logs/02_tetrapodamask.log
# round 2: rename outputs
rename simple_mask.masked.fasta tetrapoda_mask 02_tetrapoda_out/reference-genome*
rename .masked .masked.fasta 02_tetrapoda_out/reference-genome*

# round 3: annotate/mask known elements sourced from species-specific de novo repeat library using output froom 2nd round of RepeatMasker
RepeatMasker -pa 16 -a -e ncbi -dir 03_known_out -nolow \
-lib round-1_RepbaseTetrapoda-Self/round-1_RepbaseTetrapoda-Self.known \
02_tetrapoda_out/reference-genome.tetrapoda_mask.masked.fasta 2>&1 | tee logs/03_knownmask.log
# round 3: rename outputs
rename tetrapoda_mask.masked.fasta known_mask 03_known_out/reference-genome*
rename .masked .masked.fasta 03_known_out/reference-genome*

# round 4: annotate/mask unknown elements sourced from species-specific de novo repeat library using output froom 3nd round of RepeatMasker
RepeatMasker -pa 16 -a -e ncbi -dir 04_unknown_out -nolow \
-lib round-1_RepbaseTetrapoda-Self/round-1_RepbaseTetrapoda-Self.unknown \
03_known_out/reference-genome.known_mask.masked.fasta 2>&1 | tee logs/04_unknownmask.log
# round 4: rename outputs
rename known_mask.masked.fasta unknown_mask 04_unknown_out/reference-genome*
rename .masked .masked.fasta 04_unknown_out/reference-genome*
```

Great, now we can combine everything together to do build a full summary of the repeat element landscape of our reference genome.

```bash
# create directory for full results
mkdir -p 05_full_out

# combine full RepeatMasker result files - .cat.gz
cat 01_simple_out/reference-genome.simple_mask.cat.gz \
02_tetrapoda_out/reference-genome.tetrapoda_mask.cat.gz \
03_known_out/reference-genome.known_mask.cat.gz \
04_unknown_out/reference-genome.unknown_mask.cat.gz \
> 05_full_out/reference-genome.full_mask.cat.gz

# combine RepeatMasker tabular files for all repeats - .out
cat 01_simple_out/reference-genome.simple_mask.out \
<(cat 02_tetrapoda_out/reference-genome.tetrapoda_mask.out | tail -n +4) \
<(cat 03_known_out/reference-genome.known_mask.out | tail -n +4) \
<(cat 04_unknown_out/reference-genome.unknown_mask.out | tail -n +4) \
> 05_full_out/reference-genome.full_mask.out

# copy RepeatMasker tabular files for simple repeats - .out
cat 01_simple_out/reference-genome.simple_mask.out > 05_full_out/reference-genome.simple_mask.out

# combine RepeatMasker tabular files for complex, interspersed repeats - .out
cat 02_tetrapoda_out/reference-genome.tetrapoda_mask.out \
<(cat 03_known_out/reference-genome.known_mask.out | tail -n +4) \
<(cat 04_unknown_out/reference-genome.unknown_mask.out | tail -n +4) \
> 05_full_out/reference-genome.complex_mask.out

# combine RepeatMasker repeat alignments for all repeats - .align
cat 01_simple_outreference-genome.simple_mask.align \
02_tetrapoda_out/reference-genome.tetrapoda_mask.align \
03_known_out/reference-genome.known_mask.align \
04_unknown_out/reference-genome.unknown_mask.align \
> 05_full_out/reference-genome.full_mask.align
```

That takes care of our `.cat.gz`, `.out`, and `.align` files. To generate a new, combined `.tbl` summary table of repeat composition, we can run the command `ProcessRepeats`, distributed as part of RepeatModeler, with the combined `.cat.gz` full RepeatMasker output. This can take some time to run.

```bash
# resummarize repeat compositions from combined analysis of all RepeatMasker rounds
ProcessRepeats -a -species tetrapoda 05_full_out/reference-genome.full_mask.cat.gz 2>&1 | tee logs/05_fullmask.log
```

However, I have found that the summary of repeat element clades can be somewhat arbitrary and that it is sometimes better to summarize repeats in our own manner. Doing so is essentially an exercise in parsing the tabular list of all annotated/masked repeat regions/elements contained in the combined `.out` file. Here are some commands for doing this that I have put together, which result in distinct, but similar, values to those reported in the combined `.tbl` file produced by `ProcessRepeats`. I have not done systematic comparisons so I would avoid relying on these commands, but they may be useful.

I have combined everything into a single 'one-liner' command, but here is what is going on. Let's first summarize the nucleotide content of our reference genome, since this is needed to calculate percentages (i.e., this is the denominator). We can use [`seqtk`](https://github.com/lh3/seqtk) and [`datamash`](https://www.gnu.org/software/datamash/) to do so. I have quantified the total number of nucleotides in the genome and the number that are Ns (i.e., gaps). The first two lines of the command are doing this.

Then we can summarize each record in the combined `.out` tabular output of annotated repeats, which includes the major repeat group (e.g., LINE, SINE, DNA, etc.), subgroup (if applicable; e.g., CR1, hAT, etc.), and the total basepairs masked for each record. This information is then summarized per subfamily and fed through an `awk` command to calculate the proportion of the genome that each repeat subfamily composes as an additional column in a new `.tabulate` file. 

```bash
# calculate the length of the genome sequence in the FASTA
allLen=`seqtk comp reference-genome.fasta | datamash sum 2`; 
# calculate the length of the N sequence in the FASTA
nLen=`seqtk comp reference-genome.fasta | datamash sum 9`; 
# tabulate repeats per subfamily with total bp and proportion of genome masked
cat 05_full_out/reference-genome.full_mask.out | tail -n +4 | awk -v OFS="\t" '{ print $6, $7, $11 }' | 
awk -F '[\t/]' -v OFS="\t" '{ if (NF == 3) print $3, "NA", $2 - $1 +1; else print $3, $4, $2 - $1 +1 }' | 
datamash -sg 1,2,3 sum 4 | grep -v "\?" | 
awk -v OFS="\t" -v genomeLen="${allLen}" '{ print $0, $4 / genomeLen }' > 05_full_out/reference-genome.full_mask.tabulate
```

These proportions are based on the full genome length but one could easily calculate them based on the number of non-N bases in the genome instead by swapping in the following command: `awk -v OFS="\t" -v genomeLen="${allLen}" -v nLen="${nLen}" '{ print $0, $4 / (genomeLen - nLen) }'`.

By summing across the families or subfamilies of repeats in different ways, using, for example, `datamash` in the Unix shell or the [tidyverse](https://www.tidyverse.org/) in R, one can easily summarize repeat composition in a way that is more transparent than what is going on in RepeatMasker/ProcessRepeats (which has never been detailed or published). For example, I explicitly filter away the small number of repeat records that have a `?` in the repeat family annotation (see the `grep` command), as it is not clear what this means and if/how it is different from a record without the mark.

We should also create GFF files of our repeat annotations, since GFF has become a standard for genome annotation information. RepeatMasker has a distributed script called `rmOutToGFF3.pl` that will convert the combined `.out` file to GFF version 3. I did not like how this script names repeats in the 9th column of the GFF, so I created my [own script](https://github.com/darencard/GenomeAnnotation/blob/master/rmOutToGFF3custom) that produces the same file in a similar way, but with improved (in my opinion) annotation information in the 9th column. Use whichever you would prefer, but I will demo mine below.

```bash
# use Daren's custom script to convert .out to .gff3 for all repeats, simple repeats only, and complex repeats only
rmOutToGFF3custom -o 05_full_out/reference-genome.full_mask.out > 05_full_out/reference-genome.full_mask.gff3
rmOutToGFF3custom -o 05_full_out/reference-genome.simple_mask.out > 05_full_out/reference-genome.simple_mask.gff3
rmOutToGFF3custom -o 05_full_out/reference-genome.complex_mask.out > 05_full_out/reference-genome.complex_mask.gff3
```

Finally, we can produce the final masked genome with simple repeats soft masked and complex, interspersed repeats hard masked. This file is compatible with downstream annotation software that interprets masking, like MAKER (see above). In fact, you could feed this masked genome into the MAKER pipeline and turnoff the repeat masking portion of the pipeline (see `-RM_off` option) to annotate protein-coding genes. See [my guide for annotating protein-coding regions using MAKER](https://darencard.net/blog/2017-05-16-maker-genome-annotation/) for more information. I use [bedtools](https://bedtools.readthedocs.io/en/latest/) and my GFF files to perform the masking.

```bash
# create masked genome FASTA files
# create simple repeat soft-masked genome
bedtools maskfasta -soft -fi reference-genome.fasta -bed 05_full_out/reference-genome.simple_mask.gff3 \
-fo 05_full_out/reference-genome.fasta.simple_mask.soft.fasta
# create complex repeat hard-masked genome
bedtools maskfasta -fi 05_full_out/reference-genome.simple_mask.soft.fasta \
-bed 05_full_out/reference-genome.complex_mask.gff3 \
-fo 05_full_out/reference-genome.simple_mask.soft.complex_mask.hard.fasta
```

Finally, let's parse the output one last way that allows us to understand patterns of ancestral repeat element proliferation by using variation across the elements in the repeat alignments to date the relative timing. RepeatMasker has a supplementary script called `createRepeatLandscape.pl` that summarizes the data in this way but I have preferred to use one that does something similar created by [Aurelie Kapusta](https://scholar.google.com/citations?user=AW7uxjIAAAAJ&hl=en) called [`parseRM.pl`](https://github.com/4ureliek/Parsing-RepeatMasker-Outputs).

```bash
allLen=`seqtk comp reference-genome.fasta | datamash sum 2`;
parseRM.pl -v -i 05_full_out/reference-genome.full_mask.align -p -g ${allLen} -l 50,1 2>&1 | tee logs/06_parserm.log
```

This script creates five output files that use the input file name as a prefix, and are thus saved alongside that file. The suffixes are `.landscape.Div.Rclass.tab`, `.landscape.Div.Rfam.tab`, `.landscape.Div.Rname.tab`, `.parseRM.all-repeats.tab`, and `.parseRM.summary.tab`. Their contents are summarized in Aurelie's documentation.

Finally, let's do some file compression, as many of these files can get quite large and we have several copies of each. The file extensions for those files that get quite large are as follows: `.align`, `.fasta`, `.gff3`, and `.out`. Each can be used with the following `find` command to identify and compress all files for permanent storage.

```bash
find . -type f -name "*.align" | sort | while read file; do gzip ${file}; done
```

Armed with these outputs, one can summarize the repeat compositions/landscapes of one or more genomes in a variety of ways. I'll leave that to a future update at this point. I pieced this guide together from some scripts that were doing more complex processes, so I apologize in advance if I made any small mistakes in generalizing/simplifying everything. Hopefully, everything is clear enough to easily overcome any instances of that if you are doing a lot of copying and pasting from the guide. You can always submit a pull request against this blog post markdown file in my [website repository](https://github.com/darencard/darencard.github.io) to suggest edits.
