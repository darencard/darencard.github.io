---
layout: posts
title: "Extracting paired FASTQ read data from a BAM mapping file"
date: 2018-09-07
excerpt: "Details on acquiring reads from BAM file."
---

Sometimes FASTQ data is aligned to a reference and stored as a BAM file, instead of the normal FASTQ read files. This is okay, because it is possible to recreate raw FASTQ files based on the BAM file. The following outlines this process. The useful software `samtools` and `bedtools` are both required.

From each bam, we need to extract:
1. reads that mapped properly as pairs
2. reads that didn’t map properly as pairs (both didn’t map, or one didn’t map)

For #1, the following command will work. This was taken from this [webpage](https://crazyhottommy.blogspot.com/2013/06/count-how-many-mapped-reads-in-bam-file.html).

```bash
samtools view -u -f 1 -F 12 lib_002.sorted.md.bam > lib_002_map_map.bam
```

The `-f` and `-F` filter using flags in column 2 of the BAM file. These aren't always intuitive, and I won't describe them more here, but you can use this [handy tool](https://broadinstitute.github.io/picard/explain-flags.html) to better understand. Also note that the `-u` flag creates uncompressed BAM output rather than default compressed BAM output, so the files will be larger. This helps with quicker reading in later steps, but it isn't necessary to include this if you want to save disk space. `samtools` is super fast either way.

Resolving #2 is more complicated, as there are three ways a read might not have mapped as a proper pair. **A.** The first read mapped but the paired read did not. **B.** The first read did not map but the paired read did. **C.** Neither paired read mapped at all. Again, flags will be used to filter the original BAM file. This information was found at this [webpage](http://www.novocraft.com/documentation/novoalign-2/novoalign-ngs-quick-start-tutorial/1040-2/).

```bash
# R1 unmapped, R2 mapped
samtools view -u -f 4 -F 264 lib_002.sorted.md.bam > lib_002_unmap_map.bam
# R1 mapped, R2 unmapped
samtools view -u -f 8 -F 260 lib_002.sorted.md.bam > lib_002_map_unmap.bam
# R1 & R2 unmapped
samtools view -u -f 12 -F 256 lib_002.sorted.md.bam > lib_002_unmap_unmap.bam
```

As you might expect, you have to then merge the three files that contain at least one unmapped pair.

```bash
samtools merge -u lib_002_unmapped.bam lib_002_unmap_map.bam lib_002_map_unmap.bam lib_002_unmap_unmap.bam
```

Next, these BAM files must be resorted so that they are ordered by read ID instead of location in the reference.

```bash
samtools sort -n lib_002_map_map.bam lib_002_mapped.sort
samtools sort -n lib_002_unmapped.bam lib_002_unmapped.sort
```

At this time, it is a good idea to check that you have the correct number of reads and no redundancy. You can summarize the original BAM file to get an idea of where you started.

```bash
samtools flagstat lib_002.sorted.md.bam
## output
# 160673944 + 0 in total (QC-passed reads + QC-failed reads)
# 44648692 + 0 duplicates
# 136861465 + 0 mapped (85.18%:-nan%)
# 160673944 + 0 paired in sequencing
# 80336972 + 0 read1
# 80336972 + 0 read2
# 123173914 + 0 properly paired (76.66%:-nan%)
# 123173914 + 0 with itself and mate mapped
# 13687551 + 0 singletons (8.52%:-nan%)
# 37581548 + 0 with mate mapped to a different chr
# 28422566 + 0 with mate mapped to a different chr (mapQ>=5)
```
Notice the toal number of input reads that is found on the first line. You want to be sure that the number of unmapped and mapped reads total this number. It is easy to check using the following commands.

```bash
samtools view -c lib_002_mapped.sort.bam
## output
# 123173914
samtools view -c lib_002_unmapped.sort.bam
## output
# 37500030
```

Note that one paired read is counted as two reads here. If you sum these two numbers, they should equal the number you noted above, as they do here.

If all is good, you can now extract the FASTQ reads into two paired read files, as follows.

```bash
bamToFastq -i lib_002_mapped.sort.bam -fq lib_002_mapped.1.fastq -fq2 lib_002_mapped.2.fastq
bamToFastq -i lib_002_unmapped.sort.bam -fq lib_002_unmapped.1.fastq -fq2 lib_002_unmapped.2.fastq
```

And then it also makes sense to combine both the first and paired reads together from the mapped and unmapped files.

```bash
cat lib_002_mapped.1.fastq lib_002_unmapped.1.fastq > lib_002.1.fastq
cat lib_002_mapped.2.fastq lib_002_unmapped.2.fastq > lib_002.2.fastq
```

These two files should now have the same number of reads that are exactly as you would have received them if they had come directly from the sequencer as FASTQ.

Please also note that all of the commands above can be piped together in bash using `|`, which will save on disk space and time. So it is best to combine commands where possible.
