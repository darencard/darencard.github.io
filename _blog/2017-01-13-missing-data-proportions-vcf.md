---
layout: posts
title: "Quantify missing data in VCF file"
date: 2017-01-13
excerpt: "One-liner to calculate proportion missing data per sample in VCF/BCF file."
---

Simply replace `<<FILE>>` with your properly formated VCF/BCF file name (2 places).
Required bcftools v. 1.2+.

```bash
paste \
<(bcftools query -f '[%SAMPLE\t]\n' <<FILE>> | head -1 | tr '\t' '\n') \
<(bcftools query -f '[%GT\t]\n' <<FILE>> | awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | sort -k1,1n | cut -f 2)
```

We can also filter based on these proportions. Simply replace `<<INPUT>>`, `<<OUTPUT>>`, and `<<PROP>>` with the input file name, output file name, and proportion missing data at which points samples begin to get excluded, repectively. For example, 0.75 means that samples with greater than 75% missing data are filtered away. Requires bcftools v. 1.2+.

```bash
bcftools view -S ^<(paste <(bcftools query -f '[%SAMPLE\t]\n' <<INPUT>> | head -1 | tr '\t' '\n') <(bcftools query -f '[%GT\t]\n' <<INPUT>> | awk -v OFS="\t" '{for (i=1;i<=NF;i++) if ($i == "./.") sum[i]+=1 } END {for (i in sum) print i, sum[i] / NR }' | sort -k1,1n | cut -f 2) | awk '{ if ($2 > <<PROP>>) print $1 }') <<INPUT>> | bgzip > <<OUTPUT>>
```
