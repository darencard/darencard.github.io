---
layout: posts
title: "Extracting spliced sequences from GFF files"
date: 2018-01-18
excerpt: "Extract CDS and exons sequences from GFF files."
---

GFF is a common format for storing genetic feature annotations. In the case of gene annotations, subsets of elements are split over multiple lines, as things like exons and CDS features will have gaps based on the full genome sequence. Therefore, while it is easy to extract exon and CDS lines, it can be difficult to associate them together based on a parent (e.g., transcript) ID and perform downstream operations. Even extracting the full CDS sequence using a GFF file can be tricky for this reason, even though it seems trivial.

Here we'll overcome this difficulty using the [`gffread`](https://github.com/gpertea/gffread) tool. Installation is pretty easy and is documented in the GitHub README. `gffread` has a lot of options, but here we'll just document one that extracts the spliced CDS for each GFF transcript (`-x` option). Note that you can do the same thing for exons (`-w` option) and can also produce the protein sequence (`-y` option).

Let's extract the CDS sequences for each transcript using a genome sequence and a GFF annotation file.

```bash
gffread -x <out.fasta> -g <genome.fasta> <annotation.gff>
```

Pretty straight-forward. We can also set the command to output to STDOUT instead of a file by seeing the option `-x -`. Let's demo this by going a step further with the output and extracting every 3rd codon position (relies on [`bioawk`](https://github.com/lh3/bioawk)). Key is an `awk` command documented [here](https://stackoverflow.com/questions/22354082/print-every-nth-column-of-a-file).

```bash
gffread -x - -g <genome.fasta> <annotation.gff> | \
bioawk -c fastx '{ print $name, $seq }' | \
while read line; \
do \
name=$(echo $line | cut -f 1); \
echo $line | cut -f 2 | \
awk -F "" '{ for (i = 3; i <= NF; i += 3) \
printf "%s%s", $i, (i+3>NF?"\n":FS) }' | \
awk -v name="$name" '{ print ">"name; print $1 }'; \
done \
> <out.fasta>
```

Now it is easy to calculate something like GC3, the GC content of 3rd codon positions, using a tool like [`seqtk`](https://github.com/lh3/seqtk). We'll make the output in BED format.

```bash
gffread -x - -g <genome.fasta> <annotation.gff> | \
bioawk -c fastx '{ print $name, $seq }' | \
while read line; \
do \
name=$(echo $line | cut -f 1); \
echo $line | cut -f 2 | \
awk -F "" '{ for (i = 3; i <= NF; i += 3) \
printf "%s%s", $i, (i+3>NF?"\n":FS) }' | \
awk -v name="$name" '{ print ">"name; print $1 }'; \
done | \
seqtk comp - | \
awk -v OFS="\t" '{ print $1, "0", $2, ($4 + $5) / $2 }'
```

It is even easier if you want to look at GC across the entire CDS.


```bash
gffread -x - -g <genome.fasta> <annotation.gff> | \
seqtk comp - | \
awk -v OFS="\t" '{ print $1, "0", $2, ($4 + $5) / $2 }'
```


[![Analytics](https://ga-beacon.appspot.com/UA-112753663-1/gist/CDS_extract?pixel)](https://gist.github.com/darencard/9497e151882c3ff366335040e20b6714)
