---
layout: posts
title: "SRA Batch download"
date: 2015-10-06
excerpt: "Tutorial on batch downloading NCBI SRA files using Bash."
---

It has been some time since I've posted anything, and I'm trying to start blogging regularly again. I've already described how to do a batch submission of data to the [NCBI Sequence Read Archive](http://www.ncbi.nlm.nih.gov/sra), but today I was trying to do a batch download of a set of SRA sequence data for a project. Turns out, it can be a bit difficult to setup and use the [SRA Toolkit](http://www.ncbi.nlm.nih.gov/books/NBK158900/), at least in my opinion, but it is certainly easier than uploading data. Fortunately, I found a nice Biostars thread that solved this issue for me, so I figured I'd reiterate it here for my own reference, and that of others.

You can find the thread by visiting [https://www.biostars.org/p/111040](https://www.biostars.org/p/111040/).

The command is a long pipe, as follows:
```bash
esearch -db sra -query | efetch --format runinfo | cut -d ',' -f 1 | grep SRR | xargs fastq-dump --split-files --bzip2
```

I'll break this all down so it is apparent what this command is doing.
1. First, `esearch`, part of the NCBI Entrez Direct utilities, queries the database chosen (here, the SRA) with an accession (in the case of a batch, a Bioproject is most appropriate).
2. The STDOUT from #1 is directed into `efetch`, which uses this metadata to format a report in the 'runinfo' format, which is a comma-separated table of information about the accession.
3. The STDOUT from #2 is then subsetted, such it splits columns by commas and takes only the first column, which corresponds to the SRR accession numbers.
4. The STOUT list of SRR accession (#3), plus the header you don't want, are then sent through grep so that only SRR accession numbers are passed along.
5. Finally, xargs is used to take the STDOUT from #4 and run fastq-dump, from the NCBI SRA Toolkit, on it. The --split-files argument splits the paired reads into two files, instead of interleaving them, and the --bzip2 flag allows you to compress the output fastq files (you could use --gzip instead).

It is worth pointing out that by correctly setting up your SRA Toolkit configuation, you'll notice some intermediate files being written to you `/path/to/ncbi/sra directory`, where the /path/to usually equals $HOME. The ultimate fastq(.gz/bz2) files are written to the directory you called this command from, or to an output directory you can set using fastq-dump (see documentation). This will take some time to actually run on a larger set of files.
