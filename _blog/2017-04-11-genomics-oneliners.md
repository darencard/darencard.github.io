---
layout: posts
title: "Useful genomics one-liners"
date: 2017-04-11
excerpt: "Running list of shell one-liners that can be quite useful in genomics."
---

The following commands sometimes require non-standard software like [bioawk](https://github.com/lh3/bioawk) and [seqtk](https://github.com/lh3/seqtk).

rename scaffold headers with sequential numbers and lengths ("scaffold-N <len>")
```bash
bioawk -c fastx '{ print ">scaffold-" ++i" "length($seq)"\n"$seq }' < genome.fasta > new_genome.fasta
```

make association table of old and renamed scaffold names after above renaming command
```bash
paste <(grep ">" old_genome.fasta) <(grep ">" new_genome.fasta) | sed 's/>//g'
```

wrap fasta sequence lines to desired length (usually 80bp)
```bash
seqtk seq -l 80 seqs.fasta > wrapped_seqs.fasta
```

calculate sequence N50 length
```bash
seqtk comp genome.fasta | cut -f 2 | sort -rn | \
awk '{ sum += $0; print $0, sum }' | tac | \
awk 'NR==1 { halftot=$2/2 } lastsize>halftot && $2<halftot { print $1 } { lastsize=$2 }'
```

keep lines in a file if a certain field repeats at least N times.
for example: if you wanted to calculate scaffold-wide stats only from scaffolds with 10 or greater samples,
you'd use this to filter away lines from scaffolds with less than 10 samples.
Just vary each $1 to the desired column number and change the 10 in the term `if (c[i] >= 10)` to desired theshold.
```bash
cat file.txt | awk 'BEGIN { FS="\t" } { c[$1]++; l[$1,c[$1]]=$0 } END { for (i in c) { if (c[i] >= 10) for (j = 1; j <= c[i]; j++) print l[i,j] } }'
```

remove FASTA sequence lines from a GFF file (sometimes they are included)
```bash
sed '/^##FASTA$/,$d' <genome.fasta> > <genome.noseq.fasta>
```

rename FASTA headers using a 2-column lookup table
output format: >New|Old
lookup table: tab-delimited 2-column text file with Old ID in first column and New ID in second
```bash
awk 'FNR==NR { a[">"$1]=$2; next } $1 in a { sub(/>/,">"a[$1]"|",$1)}1' lookup.txt seqs.fasta
```

run a blast and automatically create GFF format output
must output blast in format 6
```bash
blast -db <db> -query <query> -outfmt 6 | \
  awk -v OFS="\t" '{ if ($10 > $9) print $2, "tblastn", "match", $9, $10, $12, "+", ".", "ID="$1; \
    else print $2, "tblastn", "match", $10, $9, $12, "-", ".", "ID="$1 }'
```

average of a set of sequential column per row
adjust 'X' in two places below if not all columns are needed
```bash
awk -v OFS="\t" '{ sum = 0; for(i=X; i<=NF; i++) sum += $i; print $0, sum, sum / (NF-(X-1)) }'
```

produce a lookup table of protein symbols from a FASTA sequence
this relies on `bioawk` for parsing FASTA and `jq` for parsing json
the command basically queries [mygene](mygene.info) using NCBI accession IDs and parses the output
can then use this lookup table to rename the FASTA sequence with `seqkit`

[bioawk](https://github.com/lh3/bioawk)
[jq](https://stedolan.github.io/jq/)
[seqkit](http://bioinf.shenwei.me/seqkit/)
```bash
# create lookup table
while read id; \
do paste \
<(echo ${id}) \
<(wget -qO- "http://mygene.info/v3/query?q=${id}&species=all&fields=name,symbol,taxid" | jq '.hits[0] | "\(._id)_\(.symbol)"' | sed 's/"//g'); \
done < \
<(bioawk -c fastx '{ print $name }' GCF_000090745.1_AnoCar2.0_protein.faa)

# create renamed fasta
bioawk -c fastx '{ print ">"$name; print $seq }' GCF_000090745.1_AnoCar2.0_protein.faa | \
seqkit replace -p '^(.+)$' -r "\${1}_{kv}" -k GCF_000090745.1_AnoCar2.0_protein_lookup.txt \
> GCF_000090745.1_AnoCar2.0_protein_rename.faa
```
