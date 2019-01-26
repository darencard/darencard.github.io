---
layout: posts
title: "Tutorial for running OrthoMCL"
date: 2018-01-12
excerpt: "Step-by-step guide to running OrthoMCL ortholog detection software."
---

[`OrthoMCL`](http://orthomcl.org/orthomcl/) is the leading piece of software for inferring orthologs across several organisms. In this tutorial I will provide detailed instructions for running a set of protein annotations through `OrthoMCL`.

## Software and Data

1. `OrthoMCL`, and it's dependencies, must be installed. Detailed information on this tool and its installation can be found [here](http://orthomcl.org/common/downloads/software/v2.0/UserGuide.txt). I actually used a slightly modified version of `OrthoMCL` that was made available by the author of the `orthomcl-pipeline` (see below). There isn't much details on the ways this is different from the existing `OrthoMCL`, but this is available [here](https://github.com/apetkau/orthomclsoftware-custom).
2. `orthmcl-pipeline` must also be installed, as this is how we will automate the `OrthoMCL` process. Detailed information on this tool and its installation can be found [here](https://github.com/apetkau/orthomcl-pipeline). There is also a nice overview of the processes we'll be essentially going through below located [here](https://github.com/apetkau/microbial-informatics-2014/tree/master/labs/orthomcl).
3. Protein sequences from desired organisms must also be assembled. In this example, I am using a *de novo* annotation that I produced for one species and am combining these data with protein annotations retrieved from [NCBI](https://www.ncbi.nlm.nih.gov/) for 4 other species, as follows:

- *Boa constrictor* - My own annotation produced from MAKER. If wanting to recreate this analysis, you can leave this sample off and run with the others.
- [*Anolis carolinensis*](https://www.ncbi.nlm.nih.gov/genome/?term=txid28377[orgn])
- [*Python molurus*](https://www.ncbi.nlm.nih.gov/genome/?term=txid176946[orgn])
- [*Thamnophis sirtalis*](https://www.ncbi.nlm.nih.gov/genome/?term=txid35019[orgn])
- [*Protobothrops mucrosquamatus*](https://www.ncbi.nlm.nih.gov/genome/?term=txid103944[orgn])

We will be specifically downloading the protein sequences in FASTA format (see banner at top of the page) and will also follow the link to see the tabular annotation data, which we can download on the next page (you'll see why in a minute).

## Setting up `MySQL`

We will be using a MySQL database to store data when running `OrthoMCL`, so we need to set this up for the first time. First, we must change the default storage engine used by `MySQL` to be `myisam` and not `innodb`. This is because of an error that may appear later that is documented [here](https://www.biostars.org/p/163888/). Note that I tried increasing the buffer size for `innodb`, as documented on that page, but that this didn't work. The details of this aren't important, but `myisam` works based on my limited experience. To do so, open the `my.cnf` file (stored in `/etc/mysql/` on Linux systems) that stores settings for `MySQL` and append the following lines beneath the `[mysqld]` tag.

```bash
#
# * Set database engine
#
default-storage-engine  = myisam
```

It is easy to make this change by running `nano` with root access.

```bash
nano /etc/mysql/my.cnf`
```

This has the effect of making `myisam` the default storage engine for any new database. We must restart `MySQL` so that this change can take effect, which can be done by issuing a command like this:

```bash
sudo service mysql restart
```

Now we should create a new `MySQL` user so that `OrthoMCL` can create the desired database. To do this, we must log into `MySQL` as the `root` user.

```bash
mysql -u root
# a password usually isn't necessary if you've never used mysql
```

Now we can create a generic user by issuing the following command at the `MySQL` prompt. I just used these simple (but unsecure) settings, but feel free to customize the username (`orthomcl`) and password (`password`) as desired.

```mysql
CREATE USER 'orthomcl'@'localhost' IDENTIFIED BY 'password';
```

You can issue the command `exit` to get out of `MySQL`. Now let's create the database we would like to write our `OrthoMCL` data to. First make sure you are in the desired working directory for the downstream `MySQL` run. Now we can use a script from `orthomcl-pipeline` to automatically create the database. In this command we are supplying the username and password we created above and also naming the database (I'm using the generic name `orthomcl` but customize as desired). We just also write out settings to a file that will be used again later using `--outfile`.

```
orthomcl-setup-database.pl --user orthomcl --password password --host localhost --database orthomcl --outfile orthomcl.conf
```

Now let's do one more step to make sure that `OrthoMCL` will have the proper permissions it needs to write the data it generates to our `orthomcl` database.

```bash
mysql -u root
```

You can now make sure our `orthomcl` user has full privileges to this database.

```mysql
use orthomcl;
GRANT ALL PRIVILEGES ON orthomcl.* TO 'orthomcl'@'localhost' IDENTIFIED BY 'password';
```

It is also worth checking to make sure our database has the correct storage engine set (`myisam`).

```mysql
SELECT TABLE_NAME, ENGINE FROM information_schema.TABLES where TABLE_SCHEMA = 'database';
```

Our database is now ready for us to run `OrthoMCL`.

## Preparing our data

While the `orthomcl-pipeline` makes running `OrthoMCL` pretty straight-forward, we first need to prepare our data a little bit. First, the headers for the FASTA protein sequences currently aren't conducive to running `OrthoMCL`. We instead want them in the format `organism_ID|protein_ID`, where `organism_ID` is a short ID that is unique to each FASTA file and `protein_ID` is the protein ID set by NCBI (usually has a `NP_` or `XP_` prefix) or by MAKER. We can easily do this using [`bioawk`](https://github.com/lh3/bioawk). Here is how we would do this using our *Anolis* proteins.

```bash
bioawk -c fastx '{ print ">Acar|"$name; print $seq }' GCF_000090745.1_AnoCar2.0_protein.faa
```

Now or `organism_ID` is always a short, easy to discern `Acar` and the `protein_ID` is the first part of the existing FASTA header before the first space (i.e., the `NP_` portion).

Beyond formatting the FASTA headers, we also want to filter the protein data further. Let's only work with longer proteins of 50 amino acids or longer. We can use [`seqkit`](http://bioinf.shenwei.me/seqkit/) to easily filter our protein FASTA files, which excludes some proteins.

```bash
seqkit seq -m 50 GCF_000090745.1_AnoCar2.0_protein.faa
```

Finally, in the case of NCBI proteins (and likely proteins from other sources), we will have some redundancy in the form of different isoforms for the same gene. We would like to collapse this down and only keep a single isoform for each gene for the `OrthoMCL` analysis. This is where the tabular annotation data from NCBI will come into play. The following command will extract the longest protein product (i.e., isoform) for each distinct Gene ID to produce a new, somewhat smaller protein FASTA file with no redundancy.

```
tail -n +2 GCF_000090745.1_AnoCar2.0_protein_table.txt | sort -k9,9nr | sort -ut $'\t' -k6,6 | cut -f 8 | seqkit grep -f - GCF_000090745.1_AnoCar2.0_protein.faa
```

We can easily pipe all of these commands together to produced the processed FASTA protein files that we need to run `OrthoMCL`. Let's say we have our protein FASTA files (i.e., `_protein.faa`) and our annotation data (stored with a `_protein_table.txt` suffix) stored in a directory called `raw` in our working directory. `OrthoMCL` will run on all FASTA files in a specified directory, so let's write our processed protein FASTAs to a new directory called `processed` with the extension `.fasta` (required by `OrthoMCL`).

```
# directory structure
--/path/to/working/directory
|
-----raw/
|
-----processed/
```

```bash
# example with Anolis
tail -n +2 raw/GCF_000090745.1_AnoCar2.0_protein_table.txt | \
sort -k9,9nr | sort -ut $'\t' -k6,6 | cut -f 8 | \
seqkit grep -f - raw/GCF_000090745.1_AnoCar2.0_protein.faa | \
seqkit seq -m 50 | bioawk -c fastx '{ print ">Acar|"$name; print $seq }' \
> /processed/Acar.fasta
```

And here are the commands I used for all 5 samples (note that the processing for the *Boa* genome is simpler because we don't have issues with isoforms like we do with the NCBI genomes). In each case I'm using 4 letter `organism_ID` prefixes.

```bash
tail -n +2 raw/GCF_000090745.1_AnoCar2.0_protein_table.txt | \
sort -k9,9nr | sort -ut $'\t' -k6,6 | cut -f 8 | \
seqkit grep -f - raw/GCF_000090745.1_AnoCar2.0_protein.faa | \
seqkit seq -m 50 | bioawk -c fastx '{ print ">Acar|"$name; print $seq }' \
> processed/Acar.fasta
tail -n +2 raw/GCF_000186305.1_Python_molurus_bivittatus-5.0.2_protein_table.txt | \
sort -k9,9nr | sort -ut $'\t' -k6,6 | cut -f 8 | \
seqkit grep -f - raw/GCF_000186305.1_Python_molurus_bivittatus-5.0.2_protein.faa | \
seqkit seq -m 50 | bioawk -c fastx '{ print ">Pmol|"$name; print $seq }' \
> processed/Pmol.fasta
tail -n +2 raw/GCF_001077635.1_Thamnophis_sirtalis-6.0_protein_table.txt | \
sort -k9,9nr | sort -ut $'\t' -k6,6 | cut -f 8 | \
seqkit grep -f - raw/GCF_001077635.1_Thamnophis_sirtalis-6.0_protein.faa | \
seqkit seq -m 50 | bioawk -c fastx '{ print ">Tsir|"$name; print $seq }' \
> processed/Tsir.fasta
tail -n +2 raw/GCF_001527695.2_P.Mucros_1.0_protein_table.txt | \
sort -k9,9nr | sort -ut $'\t' -k6,6 | cut -f 8 | \
seqkit grep -f - raw/GCF_001527695.2_P.Mucros_1.0_protein.faa | \
seqkit seq -m 50 | bioawk -c fastx '{ print ">Pmuc|"$name; print $seq }' \
> processed/Pmuc.fasta
seqkit seq -m 50 raw/Bcon_rnd3.all.maker.proteins.fasta | \
bioawk -c fastx '{ print ">Bcon|"$name; print $seq }' \
> processed/Bcon.fasta
```

## Running `OrthoMCL`

We're finally ready to run `OrthoMCL` using the `orthomcl-pipeline`. This is the easiest portion of this tutorial. Let's first create a directory for the output in our working directory.

```bash
mkdir orthomcl_out
```

Now we can construct our command and put it into a shell script. We'll use settings that tell `orthomcl-pipeline` that our data is compliant with the required formatting (see above), that keeps all intermediate data (in case we want it), that splits our data into 10 chunks for running `BLAST` (I have 12 cores available, so I set this to 10), and we'll pass in the `orthomcl.conf` file we created earlier with the names of the input and output directories.

```bash
cat orthomcl_run.sh
# output:
# orthomcl-pipeline --no-cleanup --compliant -s 10 -m orthomcl.conf -i processed -o orthomcl_out
```

We can now run the pipeline, using `screen` to manage the terminal and writing a log file.

```bash
screen -S orthomcl
bash ./orthomcl_run.sh 2>&1 | tee orthomcl_run.log
```

That's it! With these data and my system, this only took about a day to run. Adding more data (i.e., species) will increase the run time.

## Summarizing the Results

It turns out there isn't much out there in the way of tutorial and tools for summarizing the results of `OrthoMCL`. So, here are some commands that will help when `OrthoMCL` is finished.

Proteins that are not clustered into othologous groups (i.e., present in at least 2 species) are exluced from the resulting `groups.txt` output file. Here is a command that will output a single line for all of the unclustered genes in all species, with `unclustered:` in the first column and the gene IDs for all species in subsequent columns (space delimited). You will have to adjust for the organism IDs you are using and the location of the input compliant FASTA sequences. You can put this in its own file or append it to the existing `groups.txt` file. I will do the latter and create a new file called `groups.all.txt`.

```bash
for species in Acar Bcon Pmol Pmuc Tsir; \
do cat \
<(cat groups/groups.txt | \
awk -v species="$species" '{ for(i=1; i<=NF; i++) if ($i ~ species) print $i }') \
<(bioawk -c fastx '{ print $name }' ../processed/$species.fasta) | \
sort | uniq -u; \
done | \
paste -s -d ' ' | \
paste -d ' ' <(echo "unclustered:") - | \
cat groups/groups.txt - > groups/groups.all.txt
```

It is also useful to have counts for each cluster for each species, rather than a long list of the proteins. That's pretty easy to do as well.

```bash
cat groups/groups.all.txt | \
awk '{ acar=bcon=pmol=pmuc=tsir=0; for(i=2; i<=NF; i++) \
if ($i ~ "Acar") acar += 1; \
else if ($i ~ "Bcon") bcon += 1; \
else if ($i ~ "Pmol") pmol += 1; \
else if ($i ~ "Pmuc") pmuc += 1; \
else if ($i ~ "Tsir") tsir += 1; \
print $1, acar, bcon, pmol, pmuc, tsir }' > groups/groups.all.counts.txt
```

Let's go further and tabulate some basic information about orthologs.

1. How many ortholog clusters have a single gene in each species (full 1:1 orthologs)?

```bash
cat groups/groups.all.counts.txt | \
awk '{ if ($2 == 1 && $3 == 1 && $4 == 1 && $5 == 1 && $6 == 1) print $0 }' | \
wc -l
```

2. How many have a single gene in at least two of the species (1:1 orthologs)?

```bash
cat groups/groups.all.counts.txt | \
awk '{ if (($2 == 1 || $2 == 0) && ($3 == 1 || $3 == 0) && ($4 == 1 || $4 == 0) && \
($5 == 1 || $5 == 0) && ($6 == 1 || $6 == 0)) print $0 }' | \
wc -l
```

3. How about how many of the above 1:1 orthologs are found in each genome? (The 1st column is the column number and the 2nd is the count.)

```bash
cat groups/groups.all.counts.txt | \
awk '{ if (($2 == 1 || $2 == 0) && ($3 == 1 || $3 == 0) && ($4 == 1 || $4 == 0) && \
($5 == 1 || $5 == 0) && ($6 == 1 || $6 == 0)) print $0 }' | \
awk '{ for (i=2; i<=NF; i++) sum[i] += $i } END { for (i in sum) print i, sum[i] }' | \
sort -k1,1
```

4. Now let's figure out how many groups contain at least two paralogs that are found in only 1 of the species. (The 1st column is the column number and the 2nd is the count.)

```bash
cat groups/groups.all.counts.txt | \
awk '{ count=0; for (i=2; i<=NF; i++) if ($i == 0) count += 1; if (count == NF - 2) print $0 }' |  \
awk '{ for (i=2; i<=NF; i++) sum[i] += $i } END {for (i in sum) print i, sum[i] }' | \
sort -k1,1
```

5. How many groups are expanded (size > 1) in each species? (Note we have to factor out the last row of unclustered genes. The 1st column is the column number and the 2nd is the count.)

```bash
cat groups/groups.all.counts.txt | \
awk '{ count=0; for (i=2; i<=NF; i++) if ($i == 0) count += 1; if (count < NF - 2) print $0 }' |  \
awk '{ for (i=2; i<=NF; i++) if ($i > 1) sum[i] += 1 } END { for (i in sum) print i, sum[i] - 1 }' | \
sort -k1,1
```

6. How many groups are expanded (size = 1) in each species? This encompasses 1:1 orthologs from above. (The 1st column is the column number and the 2nd is the count.)

```bash
cat groups/groups.all.counts.txt | \
awk '{ count=0; for (i=2; i<=NF; i++) if ($i == 0) count += 1; if (count < NF - 2) print $0 }' |  \
awk '{ for (i=2; i<=NF; i++) if ($i == 1) sum[i] += 1 } END { for (i in sum) print i, sum[i] }' | \
sort -k1,1
```

These types of data are usually plotted as stacked barplots. Also see [here](https://github.com/halexand/Ehux_HD/blob/master/orthoMCL_output/Ehux_Groups.ipynb) for another interesting way of displaying these data.
