---
layout: posts
title: "Batch Ensembl data download"
date: 2014-03-24
excerpt: "How-to on doing large downloads of Ensembl data."
---

Recently, I began work on a project that required downloading the genomes of currently sequenced vertebrates. Most of these genomes are available through Ensembl, an online repository of genomes and a website where you can browse genome annotations. In the past, I've only had to download a genome or two and always used an ftp connection through the terminal. Here are the commands you can use to do this:

```bash
cd /path/to/store/data
ftp
open ftp.ensembl.org
# Username: anonymous
# Password: <email>
# can then move around hierarchy using cd and ls commands
# to download your target file, once you navigate to it, type:
mget <file>
# answer y(es) when the server prompts.
```

If instead you want to download many files at once, you need to take a different approach using a program called rsync (which should be installed on Linux/OSX). Here are the commands you can adapt for your purpose:

```bash
cd /path/to/store/data
rsync -av rsync://ftp.ensembl.org/ensembl/pub/path/to/folder/targetfile.ext ./
# downloads targetfile.ext to your current directory (which is the \\
# path where you want to store data)
rsync -av rsync://ftp.ensembl.org/ensembl/pub/path/to/*/*.ext ./
# downloads all files with a .ext extension in all directories \\
# within the path/to/ directory
```

The command I used to download the most current releases of all complete Ensembl genomes was:

```bash
rsync -av rsync://ftp.ensembl.org/ensembl/pub/current_fasta/*/dna/*.dna.toplevel.fa.gz ./
```

You can modify this command to download all soft-masked assemblies, all cDNA assemblies, and many other file types.

One should also note the root directory structure for the Ensembl Genomes database, which includes Metazoans, Protozoans, etc. For Metazoans, it is `ftp.ensemblgenomes.org/all/pub/metazoa/current/` and this can be modified accordingly for the target genome files the user is seeking.

In all cases, these batch rsync downloads are much better than performing the above ftp steps 30+ times and will work in the background. Download rates will vary based on connection and data type, but I would expect to wait between 5 and 10 minutes per assembly file. Hopefully this helps save you some time as well.
