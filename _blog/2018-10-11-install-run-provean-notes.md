---
layout: posts
title: "Notes from running Provean"
date: 2018-10-11
excerpt: "Notes from work installing and running Provean to predict protein impact of variants."
---

Provean input files were produced based on VEP output using commands below. Some trial runs were completed using a computer to understand how quickly Provean can be run in parallel to work through all annotated genes.

```bash
# running PROVEAN

# installation & dependencies
# 1. checked that blast was installed and also reinstalled cd-hit to avoid issue with certain version
# 2. installed the NCBI nr database
sudo mkdir /opt/ncbi_blast_nr_db_2018-01-29
sudo chmod 775 /opt/ncbi_blast_nr_db_2018-01-29
sudo chown administrator:administrator /opt/ncbi_blast_nr_db_2018-01-29
cd /opt/ncbi_blast_nr_db_2018-01-29
update_blastdb --passive --decompress --timeout 300 --force --verbose nr
# 3. install provean
tar xvf provean-1.1.5.tar.gz
cd provean-1.1.5
./configure BLAST_DB=/opt/ncbi_blast_nr_db_2018-01-29/nr
make
sudo make install

# prepare PROVEAN inputs using protein sequence fasta and the results of a VEP analysis
# create bed summary of output that we can easily do random searches from
cat Boco_wgs_BZ_HN_joint_variants_VEP_effects | \
awk '{ if ($11 ~ "/") print $0 }' | cut -f 5,10,11 | \
awk -v OFS="\t" '{ print $1, $2+0, $2+1, $3 }' | \
sort -k1,1 -k2,2n | bgzip -c \
> Boco_wgs_BZ_HN_joint_variants_VEP_effects.coding.bed.gz

# index this bed and the protein fasta
tabix -p bed Boco_wgs_BZ_HN_joint_variants_VEP_effects.coding.bed.gz
samtools faidx Bcon_rnd3.all.maker.proteins.fasta

# pull out protein variants and format properly, and gather protein fasta sequences
# variants with an X (any amino acid) are excluded
cat Boco_wgs_BZ_HN_joint_variants_VEP_effects | \
awk '{ if ($11 ~ "/") print $0 }' | cut -f5 | sort | uniq | \
while read protein; \
do \
echo ${protein}; \
tabix Boco_wgs_BZ_HN_joint_variants_VEP_effects.coding.bed.gz ${protein} | \
awk -v OFS="\t" '{ if ($4 !~ "X") print $0 }' | \
awk -F "\t|/" -v OFS="," '{ print $2, $4, $5 }' | sed 's/-/./g' \
> input_files_per_protein/${protein}.var.txt; \
samtools faidx Bcon_rnd3.all.maker.proteins.fasta ${protein} \
> input_files_per_protein/${protein}.seq.fasta; \
done

# I am able to run PROVEAN just fine, but it takes N minutes on 1 core for a protein
# of average length, mostly due to the BLAST. Given we need to run almost 15k proteins,
# we need to parallelize this a lot. I plan to run 4 sets of jobs in parallel on
# moonunits 0-2 and javelina (2000 proteins each, running on 10 cores, total run time
# of 8-10 days). I will also run 1 large set of jobs using jetstream (8000 proteins,
# running on 43 cores, total run time of 8-9 days). Math assumes 1 core per job and that
# jobs last approximately N minutes on average.
```
