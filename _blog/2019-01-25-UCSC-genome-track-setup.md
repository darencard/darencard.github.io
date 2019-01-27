---
layout: posts
title: "Visualizing Genome Annotations"
date: 2019-01-25
excerpt: "Creating a UCSC Genome Track for Viewing Genome Annotations"
---

# Creating a UCSC Genome Track for Viewing Genome Annotations

When creating a custom, in-house genome annotation, there is no straight-forward way to share it upon publication. NCBI does not accept and archive these, so most users just end up depositing the text files in an online repository. The [UCSC Genome Browser team](http://genome.ucsc.edu/) has fortunately come up with a decent way to share an annotation in a graphical manner, so readers can browse the assembly and annotation at some level without a lot of work.

There are several webpages that document how to prepare annotation files to do this properly, but I found that the gloss over the details a bit, which led me to struggle. Here are some of those resources.

* [https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html](https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html)
* [http://genomewiki.ucsc.edu/index.php/Assembly_Hubs](http://genomewiki.ucsc.edu/index.php/Assembly_Hubs)

I created this Gist to remind me of how to do this again later and to perhaps help others who want to do the same thing. This guide walks through creating the necessary data and configuration files for the genome of *Solenodon paradoxus*. The genome was assembled by folks at the Broad independent of me, but I was asked to help with the project by creating a genome assembly. I did this using the MAKER pipeline, following essentially the steps laid out in [this gist](). More details are available from the published paper (PENDING).

Commands issued below are from the following software dependencies, which will need to be installed:

* [UCSC kentUtils](https://github.com/ucscGenomeBrowser/kent)
* [bioawk](https://github.com/lh3/bioawk)
* [bedops](https://bedops.readthedocs.io/en/latest/)


## Creating Necessary Data Files

The first step is to create the necessary data files that UCSC requires to display the genome browser. These are analogous to FASTA, GFF, and BED files that most users will be familiar with, but are distinct.

First, we will create a `.2bit` file of our genome using the FASTA file of the genome assembly. This ends up being quite easy.

```
faToTwoBit SolPar.formated_assembly.fasta SolPar.formated_assembly.2bit
```

For the following steps, we also need to create a `.chrom.sizes` file, which is pretty straight-forward

```
cat SolPar.formated_assembly.fasta | bioawk -c fastx '{ print $name, length($seq) }' > SolPar.formated_assembly.chrom.sizes
```

Now that we have the backbone genome, we need to create the necessary data tracks for displaying annotation features. We will be creating `bigBed` files for this purpose, which are a special kind of BED file that can be more randomly accessed for displaying information in a genome browser. Let's begin by creating a `bigBed` file from an existing GFF file that contains only the repeat annotations outputted by MAKER. Put simply, we are converting the GFF into a BED file, cleaning up some of the name encoding to make it more readable, and then converting that into the `bigBed` file.

```
# unzip the repeats-only GFF file from maker
gunzip -c SolPar_full.all.maker.noseq.full_repeats.gff.gz > SolPar_full.all.maker.noseq.full_repeats.gff

# convert maker GFF3 repeats file into bed4
# formatting htmlencoding to make more readable
cat SolPar_full.all.maker.noseq.full_repeats.gff | gff2bed | awk -F "\t|;" -v OFS="\t" '{ print $1, $2, $3, $11 }' | sed -e 's/%2528/\(/g' -e 's/%2529/\)/g' -e 's/%252F/\//g' -e 's/Name=//g' | sort-bed - > SolPar_full.all.maker.noseq.full_repeats.bed4

# convert the bed4 into bigBed format
bedToBigBed -type=bed4 SolPar_full.all.maker.noseq.full_repeats.bed4 SolPar.formated_assembly.chrom.sizes SolPar_full.all.maker.noseq.full_repeats.bigBed
```

We can do something analogous to create a `bigBed` file for the gene annotations. The big difference is that we instead convert from GFF into a `genePred` format before converting to BED and `bigBed`, which is because of the way genes are formatted versus repeats (i.e., exons/introns).

```
# unzip the no repeat, no sequence GFF file from maker
gunzip -c SolPar_full.all.maker.noseq.no_repeats.gff.gz > SolPar_full.all.maker.noseq.no_repeats.gff

# great new proper GFF3 by extracting all maker gene annotation lines from GFF and
# replacing the _AED, _eAED, and _QI tags with aed, eaed, and qi (per standard)
cat <(echo "##gff-version 3") <(cat SolPar_full.all.maker.noseq.no_repeats.gff | awk '{ if ($2 == "maker") print $0 }' | sed -e 's/_AED/aed/g' -e 's/_eAED/eaed/g' -e 's/_QI/qi/g') > SolPar_full.all.maker.noseq.no_repeats.maker_only.gff

# convert new gene annotation GFF3 into a genePred file
gff3ToGenePred SolPar_full.all.maker.noseq.no_repeats.maker_only.gff SolPar_full.all.maker.noseq.no_repeats.maker_only.genePred

# convert the genePred file into a bed12 file
genePredToBed SolPar_full.all.maker.noseq.no_repeats.maker_only.genePred SolPar_full.all.maker.noseq.no_repeats.maker_only.bed12

# sort the bed12 file
sort-bed SolPar_full.all.maker.noseq.no_repeats.maker_only.bed12 > SolPar_full.all.maker.noseq.no_repeats.maker_only.sorted.bed12

# convert the bed12 into bigBed format
bedToBigBed SolPar_full.all.maker.noseq.no_repeats.maker_only.sorted.bed12 SolPar.formated_assembly.chrom.sizes SolPar_full.all.maker.noseq.no_repeats.maker_only.bigBed
```

The `2bit` and two `bigBed` files are the data files we need for our genome browser.

## Creating Genome Track Configuration Files

The other part of the genome browser is creating the necessary configuration files for rendering the data. UCSC keeps this quite simple overall, but some are still required. Fortunately, there is a nice Python package that makes this pretty straight-forward called [trackhub](https://daler.github.io/trackhub/index.html). See that webpage for installation and usage instructions.

Let's begin by creating symlinks that take our big ugly file names and make them more simple and analogous to what we will want our final data files to look like.

```
ln -s SolPar.formated_assembly.2bit solPar1.2bit
ln -s SolPar_full.all.maker.noseq.full_repeats.bigBed solPar1_repeats.bigBed
ln -s SolPar_full.all.maker.noseq.no_repeats.maker_only.bigBed solPar1_genes.bigBed
```

Now we can create a Python script that will use these files to create the configuration files. Here is the Python script (called `format_hub.py`) that I used, which is very similar to the one on the `trackhub` instructions webpage.

```
import trackhub
import re
import sys
import os
import glob

# In contrast to the example in the README, we do not use the
# `trackhub.default_hub` function but instead build up the hub from its
# component pieces.
hub = trackhub.Hub(
    "assembly_hub",
    short_label="solPar1 hub",
    long_label="Assembly hub for Solenodon paradoxus",
    email="daren.card@gmail.com")

# The major difference from a regular track hub is this object, which needs
# to be added to the genomes_file object:
genome = trackhub.Assembly(
    genome="solPar1",
    twobit_file="solPar1.2bit",
    organism="Hispaniolan solenodon",
    defaultPos="scaffold-1:0-1000000",
    scientificName="Solenodon paradoxus",
    description="Solenodon paradoxus V1",
    html_string="Solenodon paradoxus V1 INFO\n",
    orderKey=4800
)

genomes_file = trackhub.GenomesFile()
hub.add_genomes_file(genomes_file)

# we also need to create a trackDb and add it to the genome
trackdb = trackhub.TrackDb()
genome.add_trackdb(trackdb)

# add the genome to the genomes file here:
genomes_file.add_genome(genome)

# Find all bigBeds
for bb in glob.glob("solPar*.bigBed"):
    os.path.basename(bb).split(".")
    name, _ = os.path.basename(bb).split(".")
    track = trackhub.Track(
        name=trackhub.helpers.sanitize(name),
        source=bb,
        tracktype='bigBed')
    trackdb.add_tracks(track)

# Assembly hubs also need to have a Group specified. Here's how to do that:
main_group = trackhub.groups.GroupDefinition(
    "solPar1_tracks",
    label="solPar1 Tracks",
    priority=1,
    default_is_closed=False)

groups_file = trackhub.groups.GroupsFile([main_group])
genome.add_groups(groups_file)

# We can now add the "group" parameter to all the children of the trackDb
for track in trackdb.children:
    track.add_params(group="solPar1_tracks")

# render/stage the hub for upload
trackhub.upload.stage_hub(hub, staging="solPar1-staging")
```

The file names used in this Python script are the symlinks we just setup. This script must be run from the directory that contains these symlinks to work properly. Note that the first couple Python commands require information specific to the genome being worked with, so these will vary from instance to instance. We can now run this script to create the necessary files, which we are storing in a directory called `solPar1-staging` (see last Python command).

```
python ./format_hub.py
```

You will now see all of the necessary configuration file alongside new symlinks within the `solPar1-staging` directory.

## Hosting the Genome Hub Files

As mentioned above, UCSC will render these files into a genome browser, but they need to be stored somewhere accessible for this to happen. As documented by UCSC, (free) optins include [GitHub](https://github.com), [figshare](https://figshare.com/), and [Cyverse](https://www.cyverse.org/). GitHub has file size limits that are too restrictive and I could not understand how to export proper links to individual files for figshare, so I went with Cyverse. In order to use Cyverse, you must create an account, and the [Cyverse wiki](https://wiki.cyverse.org/wiki/dashboard.action) provides clear instructions on how to do this.

Once the user account is established, you can login through the Cyverse Discovery Environment at [https://de.cyverse.org/](https://de.cyverse.org/). You will be brought to a GUI interface similar to a computer desktop. You should create a new subfolder to store the genome hub files and then you can upload all files within the `solPar1-staging` directory to this subfolder. We are following the directions at [https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html#Hosting](https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html#Hosting). Notice that (at least on a Mac) when you click on the symlink to the data files the webpage is smart and follows those symlinks to the data file they represent, uploading these files with the names of the symlink, which is quite useful.

We also need to make some modifications to the configuration files, as we are now putting everything together on a remote server instead of your local computer. This means we need to change any pointer to a given file to the full URL where we can find that file online. This is also laid out ain the directions at [https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html#Hosting](https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html#Hosting). All URLs for files should follow the format `https://de.cyverse.org/anon-files/iplant/home/<user>/<subfolder>/<file.ext>`. You must be sure to make all necessary changes and the following table outlines how many changes were made for which files in the case of the Solenodon genome.

| File Name | Number of changes |
| ------------- | ------------- |
| assembly_hub.genome.txt | 4  |
| assembly_hub.hub.txt | 1  |
| trackDb.txt | 1 per track  |

Finally, as I figured out after having some trouble and finding [this post](https://groups.google.com/a/soe.ucsc.edu/forum/#!searchin/genome/cyverse%7Csort:date/genome/W9gqz8uLsRE/mAuc8sAMBQAJ) on the UCSC Genome Browser Public Support group, the permissions also need to be set to allow these files to be accessed by UCSC. This must be done with iCommands, which is part of iRODS. Cyverse has good documentation on setting up and running iCommands [here](https://wiki.cyverse.org/wiki/display/DS/Using+iCommands). Once setup and connected to Cyverse, you can navigate with `ils` and `icd` commands into the subfolder that was created for this project. Then for each of the files, issue the following command to make the file publically accessible: `ichmod read anonymous <file>`.

## Viewing the Genome Hub Browser

Finally, everything is setup properly for rendering your genome browser. To view it, go to the [Track Hubs Page](https://genome.ucsc.edu/cgi-bin/hgHubConnect) and click on the 'My Hubs' tab. You can then enter in the full URL to the `hubs.txt` file to load it, which will follow the following structure `https://de.cyverse.org/anon-files/iplant/home/<username>/<subfoder>/assembly_hub.hub.txt`. Once loaded, you can click the link under 'Assembles' to view the Genome Hub Browser. Navigating this genome browser should be pretty intuitive, but see the UCSC Genome Browser [webpage](https://genome.ucsc.edu/) for more information.
