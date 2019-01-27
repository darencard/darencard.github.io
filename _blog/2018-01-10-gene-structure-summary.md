---
layout: posts
title: "Summarizing the structure of gene annotations"
date: 2018-01-10
excerpt: "Summarizing various measures about the structure of gene annotations."
---

When annotating genomes it is often desireable to know the overall structure of genes, including information like exon and intron lengths among other metrics. Here is a program `genestats` that will calculate such measures for a user.

```bash
#!/usr/bin/env bash

usage()
{
cat << EOF
genestats
Version 1.0 (10 January, 2018)
License: GNU GPLv2
To report bugs or errors, please contact Daren Card (dcard@uta.edu).
This script is provided as-is, with no support and no guarantee of proper or
desirable functioning.

This script calculates several gene structure measures based on a GFF annotation file.
The bgzip and tabix programs (http://www.htslib.org/doc/tabix.html) and bedtools
(https://bedtools.readthedocs.io/en/latest/) must be installed and in the \$PATH.
The user supplies the GFF3 file as a single argument to the command and the output is a
12 column tab-delimited text written to STDOUT with the following columns:
1. transcript ID
2. transcript sequence length
3. number of exons
4. total exon sequence length
5. number of introns
6. total intron sequence length
7. number of CDS chunks
8. total CDS sequence length
9. number of 5' UTR sequences
10. total 5' UTR sequence length
11. number of 3' UTR sequences
12. total 3' UTR sequence length

To work optimally and produce reliable results, the following tags in column 3 must be
included within the GFF file:
mRNA = transcript
exon = exons
CDS = coding sequence
five_prime_UTR = 5' UTR sequence
three_prime_UTR = 3' UTR sequence

The first three tags should be present in any GFF file but the latter two may not be
depending on the source of the file. In those cases, estimates of UTR features will be
erroneous, but the other features should be reported correctly. Note that temporary files
are created in the working directory as the program runs.

USAGE:
genestats <file.gff>

EOF
}

if [[ -z $1 ]] || [[ $1 == "-h" ]] || [[ $1 == "-help" ]]
then
	usage
	exit 1
fi

# this script basically pastes together 6 queries of the GFF file for each mRNA annotation
# for the mRNA, exons, introns, CDS, and UTRs (x2)
# for each, tabix is used to rapidly pull out the feature from the GFF file and
# the appropriate feature lines are then used to produce counts and sequence lengths
# for the introns, bedtools is used to essentially find the features not matching exons
# (i.e. inverse match)

cat $1 | sort -k1,1 -k4,4n -k5,5n | bgzip -c > $1.gz && \
tabix -p gff $1.gz && \
cat $1 | \
awk '{ if ($3 == "mRNA") print $0 }' | \
awk -F "ID=|;" '{ print $1, $2 }' | \
while read line; do \
parent=`echo $line | awk '{ print $9 }'`; \
chrom=`echo $line | awk '{ print $1 }'`; \
start=`echo $line | awk '{ print $4 }'`; \
end=`echo $line | awk '{ print $5 }'`; \
paste <(echo $parent) \
<(tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "mRNA") print $0 }' | \
	awk -v sum=0 -v count=0 -v OFS="\t" '{ count += 1; sum += $5 - $4 } END { print sum }') \
<(tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "exon") print $0 }' | \
	awk -v sum=0 -v count=0 -v OFS="\t" '{ count += 1; sum += $5 - $4 } END { print count, sum }') \
<(tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "mRNA") print $0 }' > mRNA.txt && \
	tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "exon") print $0 }' > exons.txt && \
	bedtools subtract -a mRNA.txt -b exons.txt | \
	awk -v sum=0 -v count=0 -v OFS="\t" '{ count += 1; sum += $5 - $4 } END { print count, sum }' && \
	rm mRNA.txt exons.txt) \
<(tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "CDS") print $0 }' | \
	awk -v sum=0 -v count=0 -v OFS="\t" '{ count += 1; sum += $5 - $4 } END { print count, sum }') \
<(tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "five_prime_UTR") print $0 }' | \
	awk -v sum=0 -v count=0 -v OFS="\t" '{ count += 1; sum += $5 - $4 } END { print count, sum }') \
<(tabix $1.gz $chrom:$start-$end | grep "$parent" | awk '{ if ($3 == "three_prime_UTR") print $0 }' | \
	awk -v sum=0 -v count=0 -v OFS="\t" '{ count += 1; sum += $5 - $4 } END { print count, sum }'); \
done
```

Let's run this on a test GFF file to demonstrate it. Here is the GFF file we will use, called `test.gff`.

```bash
##gff-version 3
scaffold-162	.	contig	1	2399525	.	.	.	ID=scaffold-162;Name=scaffold-162
scaffold-162	maker	gene	1481900	1501151	.	+	.	ID=BoaCon13386;Name=BoaCon13386;Alias=maker-scaffold-162-augustus-gene-5.4;
scaffold-162	maker	mRNA	1481900	1501151	.	+	.	ID=BoaCon13386-RA;Parent=BoaCon13386;Name=BoaCon13386-RA;Alias=maker-scaffold-162-augustus-gene-5.4-mRNA-1;_AED=0.05;_QI=0|0|0|0.8|0.75|0.6|5|0|297;_eAED=0.25;
scaffold-162	maker	exon	1481900	1481924	.	+	.	ID=BoaCon13386-RA:exon:218;Parent=BoaCon13386-RA;
scaffold-162	maker	exon	1485901	1486158	.	+	.	ID=BoaCon13386-RA:exon:219;Parent=BoaCon13386-RA;
scaffold-162	maker	exon	1490203	1490387	.	+	.	ID=BoaCon13386-RA:exon:220;Parent=BoaCon13386-RA;
scaffold-162	maker	exon	1495786	1495929	.	+	.	ID=BoaCon13386-RA:exon:221;Parent=BoaCon13386-RA;
scaffold-162	maker	exon	1500870	1501151	.	+	.	ID=BoaCon13386-RA:exon:222;Parent=BoaCon13386-RA;
scaffold-162	maker	CDS	1481900	1481924	.	+	0	ID=BoaCon13386-RA:cds;Parent=BoaCon13386-RA;
scaffold-162	maker	CDS	1485901	1486158	.	+	2	ID=BoaCon13386-RA:cds;Parent=BoaCon13386-RA;
scaffold-162	maker	CDS	1490203	1490387	.	+	2	ID=BoaCon13386-RA:cds;Parent=BoaCon13386-RA;
scaffold-162	maker	CDS	1495786	1495929	.	+	0	ID=BoaCon13386-RA:cds;Parent=BoaCon13386-RA;
scaffold-162	maker	CDS	1500870	1501151	.	+	0	ID=BoaCon13386-RA:cds;Parent=BoaCon13386-RA;
scaffold-162	maker	gene	1641672	1650129	.	-	.	ID=BoaCon13387;Name=BoaCon13387;Alias=maker-scaffold-162-augustus-gene-5.5;
scaffold-162	maker	mRNA	1641672	1650129	.	-	.	ID=BoaCon13387-RA;Parent=BoaCon13387;Name=BoaCon13387-RA;Alias=maker-scaffold-162-augustus-gene-5.5-mRNA-1;_AED=0.21;_QI=0|0|0|0.5|0.66|0.75|4|0|130;_eAED=0.06;
scaffold-162	maker	exon	1650055	1650129	.	-	.	ID=BoaCon13387-RA:exon:226;Parent=BoaCon13387-RA;
scaffold-162	maker	exon	1647269	1647309	.	-	.	ID=BoaCon13387-RA:exon:225;Parent=BoaCon13387-RA;
scaffold-162	maker	exon	1644412	1644586	.	-	.	ID=BoaCon13387-RA:exon:224;Parent=BoaCon13387-RA;
scaffold-162	maker	exon	1641672	1641773	.	-	.	ID=BoaCon13387-RA:exon:223;Parent=BoaCon13387-RA;
scaffold-162	maker	CDS	1650055	1650129	.	-	0	ID=BoaCon13387-RA:cds;Parent=BoaCon13387-RA;
scaffold-162	maker	CDS	1647269	1647309	.	-	0	ID=BoaCon13387-RA:cds;Parent=BoaCon13387-RA;
scaffold-162	maker	CDS	1644412	1644586	.	-	1	ID=BoaCon13387-RA:cds;Parent=BoaCon13387-RA;
scaffold-162	maker	CDS	1641672	1641773	.	-	0	ID=BoaCon13387-RA:cds;Parent=BoaCon13387-RA;
scaffold-162	maker	gene	2299416	2309756	.	+	.	ID=BoaCon13393;Name=BoaCon13393;Alias=maker-scaffold-162-augustus-gene-6.11;
scaffold-162	maker	mRNA	2299416	2309756	.	+	.	ID=BoaCon13393-RA;Parent=BoaCon13393;Name=BoaCon13393-RA;Alias=maker-scaffold-162-augustus-gene-6.11-mRNA-1;_AED=0.02;_QI=171|1|1|1|1|1|2|4907|158;_eAED=0.02;
scaffold-162	maker	exon	2299416	2299884	.	+	.	ID=BoaCon13393-RA:exon:227;Parent=BoaCon13393-RA;
scaffold-162	maker	exon	2304671	2309756	.	+	.	ID=BoaCon13393-RA:exon:228;Parent=BoaCon13393-RA;
scaffold-162	maker	five_prime_UTR	2299416	2299586	.	+	.	ID=BoaCon13393-RA:five_prime_utr;Parent=BoaCon13393-RA;
scaffold-162	maker	CDS	2299587	2299884	.	+	0	ID=BoaCon13393-RA:cds;Parent=BoaCon13393-RA;
scaffold-162	maker	CDS	2304671	2304849	.	+	2	ID=BoaCon13393-RA:cds;Parent=BoaCon13393-RA;
scaffold-162	maker	three_prime_UTR	2304850	2309756	.	+	.	ID=BoaCon13393-RA:three_prime_utr;Parent=BoaCon13393-RA;
scaffold-162	maker	gene	2319095	2328392	.	+	.	ID=BoaCon13395;Name=BoaCon13395;Alias=maker-scaffold-162-augustus-gene-6.12;
scaffold-162	maker	mRNA	2319095	2328392	.	+	.	ID=BoaCon13395-RA;Parent=BoaCon13395;Name=BoaCon13395-RA;Alias=maker-scaffold-162-augustus-gene-6.12-mRNA-1;_AED=0.08;_QI=0|1|0.8|1|1|1|5|789|478;_eAED=0.08;
scaffold-162	maker	exon	2319095	2319800	.	+	.	ID=BoaCon13395-RA:exon:229;Parent=BoaCon13395-RA;
scaffold-162	maker	exon	2322272	2322404	.	+	.	ID=BoaCon13395-RA:exon:230;Parent=BoaCon13395-RA;
scaffold-162	maker	exon	2323841	2323973	.	+	.	ID=BoaCon13395-RA:exon:231;Parent=BoaCon13395-RA;
scaffold-162	maker	exon	2325881	2326087	.	+	.	ID=BoaCon13395-RA:exon:232;Parent=BoaCon13395-RA;
scaffold-162	maker	exon	2327346	2328392	.	+	.	ID=BoaCon13395-RA:exon:233;Parent=BoaCon13395-RA;
scaffold-162	maker	CDS	2319095	2319800	.	+	0	ID=BoaCon13395-RA:cds;Parent=BoaCon13395-RA;
scaffold-162	maker	CDS	2322272	2322404	.	+	2	ID=BoaCon13395-RA:cds;Parent=BoaCon13395-RA;
scaffold-162	maker	CDS	2323841	2323973	.	+	1	ID=BoaCon13395-RA:cds;Parent=BoaCon13395-RA;
scaffold-162	maker	CDS	2325881	2326087	.	+	0	ID=BoaCon13395-RA:cds;Parent=BoaCon13395-RA;
scaffold-162	maker	CDS	2327346	2327603	.	+	0	ID=BoaCon13395-RA:cds;Parent=BoaCon13395-RA;
scaffold-162	maker	three_prime_UTR	2327604	2328392	.	+	.	ID=BoaCon13395-RA:three_prime_utr;Parent=BoaCon13395-RA;
scaffold-162	maker	gene	2049228	2051664	.	+	.	ID=BoaCon13389;Name=BoaCon13389;Alias=augustus_masked-scaffold-162-processed-gene-6.2;
scaffold-162	maker	mRNA	2049228	2051664	.	+	.	ID=BoaCon13389-RA;Parent=BoaCon13389;Name=BoaCon13389-RA;Alias=augustus_masked-scaffold-162-processed-gene-6.2-mRNA-1;_AED=0.17;_QI=0|0|0|0.33|1|1|3|0|99;_eAED=0.16;
scaffold-162	maker	exon	2049228	2049258	.	+	.	ID=BoaCon13389-RA:exon:234;Parent=BoaCon13389-RA;
scaffold-162	maker	exon	2049470	2049538	.	+	.	ID=BoaCon13389-RA:exon:235;Parent=BoaCon13389-RA;
scaffold-162	maker	exon	2051465	2051664	.	+	.	ID=BoaCon13389-RA:exon:236;Parent=BoaCon13389-RA;
scaffold-162	maker	CDS	2049228	2049258	.	+	0	ID=BoaCon13389-RA:cds;Parent=BoaCon13389-RA;
scaffold-162	maker	CDS	2049470	2049538	.	+	2	ID=BoaCon13389-RA:cds;Parent=BoaCon13389-RA;
scaffold-162	maker	CDS	2051465	2051664	.	+	2	ID=BoaCon13389-RA:cds;Parent=BoaCon13389-RA;
```

Now it is very easy to run `genestats`.

```bash
genestats test.gff
```

Which will write the following output to STDOUT.

```bash
BoaCon13386-RA	19251	5	889	4	18354	5	889	0	0	0	0
BoaCon13387-RA	8457	4	389	3	8062	4	389	0	0	0	0
BoaCon13393-RA	10340	2	5553	1	4785	2	475	1	170	1	4906
BoaCon13395-RA	9297	5	2221	4	7068	5	1432	0	0	1	788
BoaCon13389-RA	2436	3	297	2	2135	3	297	0	0	0	0
```

It is possible to redirect this to an output file, but the user may also wish to calculate some basic statistics like mean exon and intron lengths, like this.

```bash
genestats test.gff | \
awk -v OFS="\t" '{ exon_sum += $3; exon_len += $4; intron_sum += $5; intron_len += $6 } END { print exon_sum, exon_len / exon_sum, intron_sum, intron_len / intron_sum }'
```

Which provides the number of exons, mean exon length, the number of introns, and the mean intron length.
