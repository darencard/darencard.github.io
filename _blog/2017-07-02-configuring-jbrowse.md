---
layout: posts
title: "Configuring JBrowse"
date: 2017-07-02
excerpt: "Instructions for configuring JBrowse for viewing genome annotations"
---

JBrowse is a handy genome browser and is especially useful for viewing the results of iterative rounds of [MAKER](http://gmod.org/wiki/MAKER). The documentation is decent, but for those not used to creating a data server, it can be difficult to understand. I struggled a bit at first.

## Software & Data

#### Software
1. [JBrowse](http://gmod.org/wiki/JBrowse) version 1.12.3 (though other versions should work just fine)

#### Data
1. GFF output files produced after each round of MAKER. These are created using the `gff3_merge` script supplied with MAKER. You can also pass the `-n` flag to avoid having the FASTA sequence written to the end of the GFF.

## Initializing JBrowse
The [Quick Start Tutorial](http://jbrowse.org/code/JBrowse-1.12.3/docs/tutorial/) and [Configuration Guide](http://gmod.org/wiki/JBrowse_Configuration_Guide) provided for JBrowse are pretty good at explaining how to get going. I'll reiterate the important steps here.

1. JBrowse relies on some dependencies that you'll need to install, mostly related to data compression. Here is how you will install these on Ubuntu, and you may need to compile from source on other systems like OSX.

```bash
 sudo apt-get install zlib1g-dev libpng-dev libgd2-noxpm-dev build-essential libexpat-dev libxml2-dev libdb-dev
```

2. Once dependencies are installed properly, you can download JBrowse and move it into the directory that is served by your web browser (i.e., apache2). On Linux machines, this will be either `/var/www/` or `/var/www/html`. The documentation doesn't provide info on OSX, but the web browser location for Macs is `/Library/WebServer/Documents`. I'm running mine on a Mac, but after this step everything should be generally the same.

```bash
cd /Library/WebServer/Documents
curl -O http://jbrowse.org/releases/JBrowse-1.12.3.zip
unzip JBrowse-1.12.3.zip
```

3. Now to initialize JBrowse and view the sample volvox genome data.

```bash
cd JBrowse-1.12.3
./setup.sh
```

4. Finally, to view the volvox data being served from your local machine, you simply have to open an internet browser and go the the following URL: [http://localhost/JBrowse-1.12.3/index.html?data=sample_data/json/volvox](http://localhost/JBrowse-1.12.3/index.html?data=sample_data/json/volvox). Note that this may be slightly different if you are working with a different version of JBrowse.

## Integrating Your Own Data

It is great to have the initialized JBrowse, but the real goal for most users is to integrate their own genome and tracks. For example, JBrowse is useful for evaluating gene models produced using MAKER. Here are the steps to follow to do this. With all commands, simply pass the `-h` flag to view help information.

#### 1. Subsetting Data

Loading genome tracks for an entire genome can take quite a while, so in many cases it is better just to load a subset. For this example, I'm going to use data from the first 25 scaffolds (based on IDs), which provide a decent mix of long and short scaffolds. In other cases, it might make sense to use the longest scaffolds and/or those with the most gene density to evaluate gene models. The following commands will subset the FASTA and GFF files.

```bash
# subset full genome FASTA file (requires seqtk)
seqtk seq Boa_constrictor_SGA_7C_scaffolds.fa | grep -A 1 -w -E 'scaffold-[0-9]|scaffold-1[0-9]|scaffold-2[0-4]' > Boa_constrictor_SGA_7C_scaffolds.first25.fa
# subset annotation GFF file
cat Bcon_rnd1.all.maker.noseq.gff | grep -w -E '^scaffold-[0-9]|^scaffold-1[0-9]|^scaffold-2[0-4]' > Bcon_rnd1.all.maker.noseq.first25.gff
```

#### 2. Preparing the Reference Genome

Prepare the reference sequence using `prepare-refseqs.pl` script. I like to provide an alternative output path so you can keep distinct browsers separate. In this example, I have all my FASTA and GFF files for the *Boa constrictor* genome I am annotating on my `~/Desktop` in a directory called `boa_annotation`. Let's prepare the reference genome I have stored in that location.

```bash
cd /Library/WebServer/Documents/JBrowse-1.12.3
bin/prepare-refseqs.pl --out data/json/boa_first25 --fasta ~/Desktop/boa_annotation/Boa_constrictor_SGA_7C_scaffolds.first25.fa
# notice I'm storing the output data in data/json/boa
```

If all went well and you received no errors, you can point your browser towards the location for these data by varying the ending to include the supplied path in the command: [http://localhost/JBrowse-1.12.3/index.html?data=data/json/boa_first25](http://localhost/JBrowse-1.12.3/index.html?data=data/json/boa_first25). You should see the genome viewer with the Reference sequence track available. You can change between the scaffolds or chromosomes and zoom in and out, as you can on normal genome browsers.

#### 3. Adding Genomic Tracks

A reference sequence is the first step, but obviously the real goal is to overlay different tracks of genomic features. We'll do this for two subsequent rounds of a MAKER annotation of the *Boa* genome. The documentation outlines how to do this, but I think it is much easier to use a Perl script created by Yannick Wurm's research group: [gff2jbrowse.pl](https://github.com/wurmlab/afra/blob/master/bin/gff2jbrowse.pl). I produced a [modified script](https://github.com/darencard/Genome_annotation/blob/master/gff2jbrowse.pl) to allow users to specify a prefix for a given run, as the output of two different MAKER runs will have similar tracks and it can be difficult to keep everything straight with the default naming used by this script.

```bash
perl gff2jbrowse.pl --out data/json/boa_first25 --prefix rnd1_ ~/Desktop/boa_annotation/Bcon_rnd1.all.maker.noseq.first25.gff
```

While this runs, you will see the different tracks, with the passed prefix added. This script might not know how to handle some of the tracks, but those are the ones that are less important in most cases. If you know some Perl you can use the existing code to overcome this.

Once this has run it is also useful to run an additional script that allows users to search for features with autocompletion.

```bash
bin/generate-names.pl --verbose --out data/json/boa_first25
```

Once this is complete you can now refresh the browser page you visited before and you should see the tracks from the MAKER GFF. It is easy to toggle them on and off and to get a better idea of what these features look like against your genome. In the context of genome assembly, the first round of MAKER usually uses empirical RNAseq or protein data to produce gene models, so it is useful to look at the degree of overlap between the `est2genome`, `protein2genome`, and `MAKER` tracks. You may also be interested in looking at the repeat element annotations.

#### 4. Downstream Considerations

Now that the first round of MAKER has been added, one can go forward and add subsequent rounds using the same series of commands outlined in Step 3. If evaluating gene models from subsequent rounds of MAKER, it is best to turn on the `Augustus`, `SNAP`, and `MAKER` tracks and evaluate the amount of overlap. If gene models from gene prediction software are nicely overlapping the MAKER models that also evaluate empirical evidence from RNA and proteins, then your annotation is probably pretty decent. There are other [scripts](https://github.com/mscampbell/Genome_annotation) available from Michael Campbell to evaluate genome annotation quality.
