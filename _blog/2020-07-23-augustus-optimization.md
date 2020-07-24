---
layout: posts
title: "Optimizing Augustus"
date: 2020-07-23
excerpt: "Thorough guide for optimizing Augustus on well-suited BUSCO genes."
---

[Augustus](http://bioinf.uni-greifswald.de/augustus/) is a very popular piece of software for gene prediction and is commonly used for whole-genome gene annotation. One critical step in using Augustus - especially for new genomic resources for divergent taxa - is the optimization process. Fortunately, there are some ways of automating the process using tools, like [autoAugTrain](http://augustus.gobics.de/binaries/scripts/autoAugTrain.pl) (packaged with Augustus) and [BRAKER2](https://github.com/Gaius-Augustus/BRAKER). While both will work, I have found the former to be difficult to use and slow (can only use 1 core) and the latter is more of a black box (and is pretty new).

One method of training Augustus that I was drawn to and that [I have written about before](http://darencard.net/blog/2017-05-16-maker-genome-annotation/) is to use the popular software [BUSCO](https://busco.ezlab.org/). BUSCO's original intent was as a quality-control assessment of genome assemblies and annotations. The approach is to use a set of highly conserved genes across certain taxonomic groups (e.g., Insecta, Vertebrata, Mammalia) that are single-copy in most species examined. These genes have been identified and curated as part of the [OrthoDB project](https://www.orthodb.org/). This results in upwards of thousands of genes to profile and their presence/absence, completeness, and copy number are a great way of evaluating the quality of new genomic resources. BUSCO utilizes three key pieces of software: [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/), [HMMER](http://hmmer.org/), and Augustus. Importantly, in utilizing Augustus, the authors of BUSCO have included an option (off by default) that automates Augustus training using the complete BUSCOs identified. BUSCO can, therefore, be co-opted as a tool for automated Augustus training in new genomic resources using hundreds to thousands of conserved, single-copy (usually) genes and it also utilizes multiple cores.

While BUSCO is great for this purpose, one big disadvantage is that BUSCO will attempt to train Augustus with *all* "complete" genes, which with certain taxa sets results in thousands of genes. Even with multiple cores, this process can take a long time because of the many genes. If one is patient, that is fine, but the Augustus authors point out that beyond about 1000 genes, improvements in gene prediction training are very diminished, so it ends up being inefficient to have BUSCO use thousands of genes. Moreover, BUSCO can be somewhat of a black box because the user is locked into the genes included within the BUSCO taxon sets with no real ability to make adjustments to prioritize featuers the [Augustus authors note as important for training purposes](https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html): longer, non-redundant genes with numerous exons and introns. Perhaps most importantly, BUSCO does not set aside a subset of genes that can be used to test the Augustus HMMs after training, which is standard practice and is useful for understanding specificity and sensitivity.

In thinking about these shortcomings, I devised an approach that could be taken to speed up the Augustus training process while prioritizing the complete genes that have the most beneficial features for training. So while we are still limited to the genes in the BUSCO taxon sets, at least we are prioritizing those that are best for training and not overloading BUSCO/Augustus with many times more genes than are needed to properly train. The following sections outline the approach I took.

#### Software prerequisites:
1. [Augustus](http://bioinf.uni-greifswald.de/augustus/) version 3.3.2.
2. [BUSCO](http://busco.ezlab.org/) version 3.1.0.
3. [jq](https://stedolan.github.io/jq/) version jq-1.5.
4. [seqkit](https://bioinf.shenwei.me/seqkit/) version 0.10.2.
5. [BEDtools](https://bedtools.readthedocs.io/en/latest/) version 2.29.0.
6. [cd-hit](http://weizhongli-lab.org/cd-hit/) version 4.6.4.

#### Raw data/resources:
1. An appropriate [BUSCO taxon dataset](https://busco.ezlab.org/busco_v4_data.html). This tutorial assumes use of ODB10 and I used the [Vertebrata](https://busco-data.ezlab.org/v4/data/lineages/vertebrata_odb10.2019-11-20.tar.gz) dataset. Other taxon datasets would probably work but YMMV.
2. A reference genome that you want to optimize Augustus for or, if following my [MAKER tutorial](http://darencard.net/blog/2017-05-16-maker-genome-annotation/), a set of contigs for the regions surrounding annotated genes from a preliminary round of gene prediction (MAKER, Augustus, etc.). For the purposes of this tutorial, I'll just refer to a generic `genome.fasta` that could be either.

#### Running commands:
I've become accustomed to running most programs (especially those that take hours/days) using `screen`. I would create new screens for each command below. I also like to use `tee` to keep track of run logs, as you will see in the commands below. Most of these commands can be run interactively if you are using a shared computer cluster, with the exception of the gene prediction and optimization commands near the end, which should instead be submitted as a job.

## Ranking BUSCO genes

While the BUSCO genes are very useful, it is currently not possible to prioritize certain sets of genes for Augustus optimization. My main goal was to gather metadata about the BUSCO genes that would allow me to rank the genes based on features that the Augustus authors indicate as most important:

1. Genes with larger numbers of exons and introns are best for training Augustus to identify these features.

2. Related to #1, longer genes are better for training purposes.

3. The Augustus authors don't mention this specifically, but more conserved genes are easier to identify completely for training purposes.

While none of these features are overly-complicated, they are not included by default as metadata alongside the BUSCO genes included in certain taxa sets. This makes it impossible to rank the BUSCO genes to meet my overall goals. Fortunately, the OrthoDB website does include this type of information and is directly linked to the identifiers that are used in BUSCO ([example](https://www.orthodb.org/?query=73826at7742)). However, with thousands of genes, it is time-consuming and error-prone to extract these by hand. Luckily, OrthoDB has a nice API that we can query to extract this metadata. To query the database using the API, all one needs to do is construct an appropriate URL. When visiting this URL, the user will see the raw data underlying OrthoDB in a text-based JSON format. Parsing this output appropriately allows one to gather the relevent metadata we are after. See the following script for this purpose, which is based on the Vertebrata set and OrthoDB release 10 (odb10).

```bash
cat \
<(echo -e "#group\tgene_name\torthodb_url\tevolutionary_rate\tmedian_exon_count\tstdev_exon_count\tmedian_protein_length\tstdev_protein_length") \
<(cat vertebrata_odb10/links_to_ODB10.txt | cut -f 1 | \
while read id; \
do \
wget -qO - https://dev.orthodb.org/group?id=${id} | \
jq -r '. | [.data.id, .data.name, .url, .data.evolutionary_rate, .data.gene_architecture.exon_median_counts, .data.gene_architecture.exon_stdev_counts, .data.gene_architecture.protein_median_length, .data.gene_architecture.protein_stdev_length] | @tsv'; \
sleep 1s; \
done) > vertebrata_odb10.info.txt
```

The output includes several fields, but the ones we are most interested in is an evolutionary rate calculated by OrthoDB and the median and standard deviation of exon number and protein length across all species included for a given BUSCO. By extracting these fields for all BUSCO proteins, we can pretty easily sort to rank the BUSCO genes by some objective criteria. I choose to prioritize exon number (and by extension, intron number; more is better), followed by protein length (longer is better) and evolutionary rate (lower is better).

```bash
cat <(cat vertebrata_odb10.info.txt | head -1) \
<(cat vertebrata_odb10.info.txt | grep -v "^#" | sort -t $'\t' -k5,5nr -k7,7nr -k4,4n) > vertebrata_odb10.info.ranked.txt
```

Now we have our BUSCO gene set ranked based on features that are helpful for optimizing Augustus.

## Extracting complete, single-copy, non-redundant BUSCOs

Now that we have rankings of the Vertebrata BUSCOs, we can extract those that are complete and non-redundant for Augustus training purposes. To identify complete BUSCOs, we need to do a standard BUSCO run, but without the length internal Augustus training. These runs typically take a few hours (or maybe a day at most), so it is pretty quick.

```bash
# remember we're using genome.fasta as generic file
run_BUSCO.py -i genome.fasta  -o genome_busco_out -l vertebrata_odb10/ \
-m genome -c 24 -sp human --augustus_parameters='--progress=true'
```

As indicated in the above command, our BUSCO output will be within `run_genome_busco_out`. There is a lot of output that is organized hierarchically and you can find details on all these data on the [BUSCO website](https://busco.ezlab.org/busco_userguide.html#outputs). Several of these BUSCO outputs will be used below. Importantly, the output of BUSCO gives us a couple pieces of important information.

1. Overall, the amount of redundancy across BUSCOs (i.e., two or more BUSCOs with high sequence identity) is likely pretty low since duplicated genes are exluded by definition (i.e., no gene families to worry about). That said, there can still be some sequence redundancy across different BUSCOs and the authors of Augustus recommend collapsing redundant genes that are 80% similar to one. Apparently similar genes can lead Augustus to overfit parameters and results can be lower quality.

2. For all BUSCOs in a given taxon set, BUSCO will identify whether they are complete, missing, fragmented, or duplicated. Therefore, it is easy to cross-reference our ranked list of BUSCOs with those that are complete. With most genomic resources of reasonable quality, this basic BUSCO run reliably identifies at least 80% of BUSCOs as complete and single copy (in my experience). So while we might toss out some BUSCOs, we are still left with thousands (in the case of Vertebrata) and we get the complete, single copy genes we want most.

To remove redundancy, we can use a clustering tool like `cd-hit` to cluster complete BUSCO proteins that share 80% identity. This is quite well facilitated by the output of BUSCO, which includes protein FASTA sequences located within `run_genome_busco_out/augustus_output/extracted_proteins`. However, we only want to consider complete BUSCOs, but another output makes these very easy to identify: `run_genome_busco_out/full_table_genome.tsv`. So we can pretty easily gather the BUSCO IDs of complete, single-copy genes out of this table and use this to concatenate the relevant protein sequence FASTAs together into a single file.

Let's first extract the BUSCO IDs for genes that were complete and that contain only a single sequence in the protein FASTA file. You will see some errors related to files not being found, which is because "Missing" genes do not have an amino acid sequence FASTA that we can use to count the number of sequences, so this can safely be ignored.

```bash
cat vertebrata_odb10.info.ranked.txt | grep -v "^#" | cut -f 1 | \
while read id; \
do \
status=`cat run_genome_busco_out/full_table_genome.tsv | awk -v id="${id}" '{ if ($1 == id) print $2 }' | sort | uniq`; \
file="run_genome_busco_out/augustus_output/extracted_proteins/"${id}".faa.1"; \
count=`cat ${file} | grep -c ">"`; \
echo -e "${status}\t${count}\t${file}" | awk '{ if ($1 == "Complete" && $2 == 1) print $3 }'; \
done > complete_singlecopy.aafiles.txt
```

This resulted in 1829 BUSCOs of the original 3354 in the dataset, so quite a few genes were not complete and single copy. However, we still have plenty to work with. If you find that you are getting far fewer than this, you may have to adjust your approach. Having over 1000 retained BUSCOs would be ideal. From there, we can concatenate all of these files together and then use this to identify redundancy using `cd-hit`. While doing this, we're going to rename the headers in the FASTA so that they are the BUSCO IDs instead of coordinate data that is more obscure and less useful.

```bash
cat complete_singlecopy.aafiles.txt | \
while read file; \
do \
id=`basename ${file} .faa.1`; \
cat ${file} | seqkit fx2tab | awk -v id="${id}" -v OFS="\t" '{ print id, $2 }' | seqkit tab2fx; \
done > complete_singlecopy.seqs.rename.fasta
```

Now we can run `cd-hit` on this file to cluster by identity when identity exceeds 80%, which will remove redundancy. See `cd-hit` usage for more details. This command ran nearly instantaneously for me.

```bash
cd-hit -o complete_singlecopy.seqs.rename.cdhit -c 0.8 -i complete_singlecopy.seqs.rename.fasta -p 1 -d 0 -b 3 -T 0 -M 16000
# let's rename the default FASTA output to provide an extension
mv complete_singlecopy.seqs.rename.cdhit complete_singlecopy.seqs.rename.cdhit.fasta
```

We can then see how much redundancy was weeded out by counting the number of sequences in the two FASTA files.

```bash
grep -c ">" complete_singlecopy.seqs.rename.fasta
# pre-cd-hit = 1829
grep -c ">" complete_singlecopy.seqs.rename.cdhit
# post-cd-hit = 1818
```

So in my case, 9 proteins were identified as redundant with other proteins. This confirms that the BUSCO dataset is largely non-redundant to begin with, but we did filter out some genes that we would not want to use for Augustus training. You can take a look at the other `cd-hit` output file (`.clstr`) to see how similar these clustered sequences were. Here is what I observed with my dataset.

```bash
# ugly command to visualize redundancy, but it works
cat complete_singlecopy.seqs.rename.cdhit.clstr | grep -B1 -A1 "^1"
# Ignore the >Cluster lines
# 0	1201aa, >164538at7742... *
# 1	1201aa, >103284at7742... at 1:1201:1:1201/100.00%
# >Cluster 117
# --
# 0	1168aa, >38882at7742... *
# 1	1168aa, >325125at7742... at 1:1168:1:1168/100.00%
# >Cluster 126
# --
# 0	770aa, >321778at7742... *
# 1	545aa, >153051at7742... at 1:545:1:545/100.00%
# >Cluster 362
# --
# 0	589aa, >225169at7742... *
# 1	589aa, >276337at7742... at 1:589:1:589/100.00%
# 2	589aa, >261612at7742... at 1:589:1:589/100.00%
# --
# 0	471aa, >290048at7742... *
# 1	471aa, >215706at7742... at 1:471:1:471/100.00%
# >Cluster 859
# --
# 0	430aa, >280833at7742... at 1:430:5:434/100.00%
# 1	434aa, >239989at7742... *
# >Cluster 937
# --
# 0	391aa, >216073at7742... *
# 1	391aa, >401756at7742... at 1:391:1:391/100.00%
# >Cluster 1056
# --
# 0	149aa, >303510at7742... *
# 1	149aa, >386595at7742... at 1:149:1:149/100.00%
# >Cluster 1741
```

So we see that 8 proteins were redundant with one other protein and there was one case where three proteins were redundant with each other. Interestingly, these sequences are 100% similar and looking at a few, are the same exact length. To me, this seems problematic considering the goal of BUSCO, but given it is rare, I don't think it detracts from the utility of BUSCO.

## Extracting annotation GFF files for optimization

Now we need to cross-reference our ranked set of BUSCOs with these non-redundant `cd-hit` outputs. From there, we can extract the GFF annotations for these complete, non-redundant outputs, which will be important for doing the Augustus optimization. This is also a critical step, as it is where we trim the number of genes down to help the Augustus training run more efficiently. You can set this arbitrarily, but I'm going to subset down to 1500 genes and later we will use 1000 of them for training and 500 for testing. Feel free to adjust this accordingly, but the Augustus authors indicate that 200 genes for training is the absolute minimum (better to shoot for 500) and remember that over 1000 genes, gains in optimization are very diminished. I'm going to outline doing this a couple different ways depending on what type of input you used for this process.

If you ran BUSCO on your entire reference genome, this is easier and you can run the following. Note that all genes have the ID "g1" in the Augustus outputs, which is non-unique and would be an issue later on. So we need to modify this within the loop so gene IDs are unique. These names are arbitrary, so I did this by using a `counter` variable that increases each time through the loop.

```bash
counter=1; cat vertebrata_odb10.info.ranked.txt | grep -v "^#" | cut -f 1 | \
while read id; \
do \
grep -w "${id}" complete_singlecopy.seqs.rename.cdhit.fasta; \
done | \
tr -d ">" | head -1500 | \
while read busco; \
do \
cat run_genome_busco_out/augustus_output/predicted_genes/${busco}.out.1 | sed "s/g1/g${counter}/g"; \
let counter+=1; \
done | \
grep -v "^#" | sort -k1,1 -k4,4n > vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gff
```

You may instead be following the [approach I have outlined](http://darencard.net/blog/2017-05-16-maker-genome-annotation/) where you feed in the genomic regions surrounding annotated genes from a preliminary annotation from MAKER or other pipelines. In this case, the "genome" includes shorter contigs and the sequence headers are the coordinates of that putative gene-containing region instead of `chrN` or `scaffold-N`. However, ideally we go back to the reference genome to extract the training sets around annotated genes. Unfortunately, by feeding in subsets of the genome, you somewhat lose the broader genomic context and have to re-establish that. In other words, we need to convert the coordinates of our small chunks of the genome back to the full genome coordinates. Below I do this assuming the short contigs have the header name `chromosome:start-end`, which is the default output for several tools that will extract sequences based on coordinates (e.g., `bedtools`).

What we're doing is adding the start coordinate for the genomic range described by `chromosome:start-end` in the first column of the GFF to the start and end coordinates of the feature found in column 4 and 5, respectively, to convert back to the full genome context. Here is an example.

```bash
scaffold_154:449818-521019	AUGUSTUS	CDS	1282	1752	.	+	0	transcript_id "g9.t1"; gene_id "g9";
# becomes
scaffold_154	AUGUSTUS	CDS	451100	451570	.	+	0	transcript_id "g9.t1"; gene_id "g9";
# start=449818+1282=451100; end=449818+1752=451570
```

Hopefully, that makes sense and here is the command that will do it.

```bash
counter=1; cat vertebrata_odb10.info.ranked.txt | grep -v "^#" | cut -f 1 | \
while read id; \
do \
grep -w "${id}" complete_singlecopy.seqs.rename.cdhit.fasta; \
done | \
tr -d ">" | head -1500 | \
while read busco; \
do \
cat run_genome_busco_out/augustus_output/predicted_genes/${busco}.out.1 | sed "s/g1/g${counter}/g"; let counter+=1; \
done | \
grep -v "^#" | \
awk -F "\t" -v OFS="\t" '{split($1,a,/:/); print a[1], a[2], $2, $3, $4, $5, $6, $7, $8, $9 }' | \
awk -F "\t" -v OFS="\t" '{split($2,a,/-/); print $1, $3, $4, a[1]+$5, a[1]+$6, $7, $8, $9, $10 }' | \
sort -k1,1 -k4,4n > vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gff
```

> Tangent: Importantly, I checked carefully that this conversion is correct using the converted GFF coordinates and the predicted CDS sequence that Augustus outputs as part of BUSCO. The way to do this is to select a single, complete BUSCO and extract the coordinate-modified CDS lines and their corresponding sequences based on the full genome. From there, we can compare it to the CDS sequences in the Augustus GFF file produced in BUSCO. One useful hint with this is to pick a gene that is on the `+` strand of the reference, as when it is on the `-` strand you need to reverse complement the sequences to check that it all makes sense (I didn't think enough and learned this the hard way, but it was a good sanity check). Here is an example command for extracting those CDS lines and gathering the sequences in FASTA. Note that I just picked an arbitrary complete gene for this purpose (`109081at7742`; you can pick anything you'd like) and that `full_genome.fasta` represents the original genome sequences that the gene region contigs were extracted from.


> ```bash
> cat run_genome_busco_out/augustus_output/predicted_genes/109081at7742.out.1 | awk '{ if ($3 == "CDS") print $0 }' | \
> awk -F "\t" -v OFS="\t" '{split($1,a,/:/); print a[1], a[2], $2, $3, $4, $5, $6, $7, $8, $9 }' | \
> awk -F "\t" -v OFS="\t" '{split($2,a,/-/); print $1, $3, $4, a[1]+$5, a[1]+$6, $7, $8, $9, $10 }' | \
> bedtools getfasta -fi full_genome.fasta -bed - -fo -
> ```

> Then we can manually compare this combined CDS sequence to that generated by Augustus, which is commented with "#" at the beginning of the line in the file `run_genome_busco_out/augustus_output/predicted_genes/109081at7742.out.1`. You are looking for the sequence labeled `coding sequence = [seq]`, which you can manually compare with the sequences you generated above. When I did this, I found that the CDS sequence from my coordinate-adjusted GFF matched that generated by Augustus within BUSCO, which confirms that my genomic coordinate adjustments are correct.

Great, so now we have a set of gene models in GFF format. These are generated by Augustus as it is run in BUSCO. These genes are complete, single-copy, and non-redundant. And we also used rankings based on exon number, protein length, and overall conservation to prioritize those gene models that are better for Augustus optimization. So now we can proceed with that training with the expectation that it will be higher quality and potentially faster than if we had run Augustus optimization within BUSCO.

## Optimizing Augustus using our filtered gene set

The following steps follow the standard Augustus method of extracting training and testing sets, training Augustus, and evaluating how well training has done. I'll report the steps here briefly to complete this tutorial, but there are lots of other resources out there that outline this process and training in general, so I encourage everyone to be familiar with those ([here](https://vcru.wisc.edu/simonlab/bioinformatics/programs/augustus/docs/tutorial2015/training.html) and [here](https://currentprotocols.onlinelibrary.wiley.com/doi/abs/10.1002/cpbi.57) are some resources).

First we'll extract the genomic regions containing these genes, plus up to 1000 bp on either side.

```bash
# note we're using the full reference genome here, as no matter what your input was, we now have coordinates based on the full genome
gff2gbSmallDNA.pl vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gff full_genome.fasta 1000 vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gb
```

Let's confirm how many gene models we have.

```bash
grep -c "LOCUS" vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gb
# 1497
```

Looks like we lost a few somehow (probably because they overlap or are very nearby), but good enough for our purposes. Next we can randomly split this into two dataset: a training dataset and a testing dataset. Let's keep the 497 for testing and use 1000 for training.

```bash
randomSplit.pl vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gb 497
```

And we can confirm these counts also.

```bash
grep -c "LOCUS" vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gb*
# vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gb.test:497
# vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gb.train:1000
```

Then we can create a new species to hold our training data. This creates a new directory with the user specified species name at $AUGUSTUS_CONFIG_PATH/species.

```bash
new_species.pl --species=myspecies
```

Now we can do some initial training.

```bash
etraining --species=myspecies vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gb.train
```

Let's first see how this initial training works with our test dataset. So, first, we'll run Augustus. If you are on a computer cluster, you will want to submit this (and the optimize_augustus.pl below) command as a job rather than working interactive (which is otherwise fine in this tutorial).

```bash
augustus --species=myspecies vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gb.test | tee augustus_untrained.out
```

Now we can look at the evaluation and see how it did.

```bash
grep -A 22 Evaluation augustus_untrained.out
#
# *******      Evaluation of gene prediction     *******
#
# ---------------------------------------------\
#                  | sensitivity | specificity |
# ---------------------------------------------|
# nucleotide level |       0.908 |       0.933 |
# ---------------------------------------------/
#
# ----------------------------------------------------------------------------------------------------------\
#            |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
#            | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
#            | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
# ----------------------------------------------------------------------------------------------------------|
#            |        |        |      |                866 |               1346 |             |             |
# exon level |   5516 |   5996 | 4650 | ------------------ | ------------------ |       0.776 |       0.843 |
#            |   5516 |   5996 |      |  471 |   20 |  375 |  477 |   25 |  844 |             |             |
# ----------------------------------------------------------------------------------------------------------/
#
# ----------------------------------------------------------------------------\
# transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
# ----------------------------------------------------------------------------|
# gene level |   517 |   497 |   60 |  457 |  437 |       0.121 |       0.116 |
# ----------------------------------------------------------------------------/
```

Now we can do the optimization of the parameters for our gene set, which is the primary pupose of this tutorial. This is pretty easy to do with Augustus, but it could take quite a while to run (days). Where possible, it is best to provide multiple compute cores for this process. While this may be shorter than running the optimization automatically with BUSCO, it still takes some time, so be patient.

```bash
optimize_augustus.pl --cpus=24 --kfold=24 --species=myspecies vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gb.train
```

Now we can do the training one more time.

```bash
etraining --species=myspecies vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gb.train
```

And we can run Augustus on our test set again ...

```bash
augustus --species=myspecies vertebrata_odb10.info.ranked.complete.nonredundant.top1500.gb.test | tee augustus_trained.out
```

... and see how well it does relative to our intial, un-optimized run.

```bash
grep -A 22 Evaluation augustus_trained.out
# *******      Evaluation of gene prediction     *******
#
# ---------------------------------------------\
#                  | sensitivity | specificity |
# ---------------------------------------------|
# nucleotide level |       0.915 |       0.934 |
# ---------------------------------------------/
#
# ----------------------------------------------------------------------------------------------------------\
#            |  #pred |  #anno |      |    FP = false pos. |    FN = false neg. |             |             |
#            | total/ | total/ |   TP |--------------------|--------------------| sensitivity | specificity |
#            | unique | unique |      | part | ovlp | wrng | part | ovlp | wrng |             |             |
# ----------------------------------------------------------------------------------------------------------|
#            |        |        |      |                885 |               1322 |             |             |
# exon level |   5559 |   5996 | 4674 | ------------------ | ------------------ |        0.78 |       0.841 |
#            |   5559 |   5996 |      |  488 |   19 |  378 |  496 |   26 |  800 |             |             |
# ----------------------------------------------------------------------------------------------------------/
#
# ----------------------------------------------------------------------------\
# transcript | #pred | #anno |   TP |   FP |   FN | sensitivity | specificity |
# ----------------------------------------------------------------------------|
# gene level |   514 |   497 |   52 |  462 |  445 |       0.105 |       0.101 |
# ----------------------------------------------------------------------------/
```

My optimization took under 48 hours using the 24 cores I provided in the above command. For comparison, when I used this same dataset and ran BUSCO with the optimization step activated, which would do the training with all complete BUSCOs, it had run over 10 days on 32 cores and still had not completed. Given this took about 3 days of collective work, this is quite a time savings. Though in looking at the results, the optimization led to very modest differences in sensitivity and specificity, which I'll have to explore more (maybe that will be a topic for a later post).
