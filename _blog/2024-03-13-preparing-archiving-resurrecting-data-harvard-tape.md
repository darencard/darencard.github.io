---
layout: posts
title: "Preparing, Archiving, and Resurrecting Data To and From Harvard Tape Storage"
date: 2024-03-13
excerpt: "How to use Harvard's tape storage with your data."
---

Note: This post is currently a work in progress so please bear with me as I work to complete it. In the meantime, please proceed with caution if you are referring to any information in the post.

Currently, I have outlined the post and written some of the introductory sections. I also copied in two scripts that form the basis of the computational workflow that is described in this post. For each, I have begun providing detailed descriptions of what is happening so that these can be used by others.

Eventually, I will fill in more details about moving the data to and from the tape storage (with screenshots). I will also outline a simple log system that I devised for when data is written to the tape storage. I also plan to include a working example with a small dataset that will allow others to reproduce what is happening here before they work with their own data.

## Preamble

Given the date on which I am beginning to compose this post, March 13, 2024, I find it fitting to comment briefly about what happened four short years ago when the COVID-19 pandemic became a national emergency. I was about 18 months into my postdoc at Harvard University and I remember the date so distinctly because it was Friday the 13th in 2020. The following 18-24 months are now honestly a bit of a blur but I'm grateful I made it through that trying time. Others were not so lucky and that fact continues to weigh on me to this day. This comment is otherwise not relevant to the topic of this blog post but I still felt compelled to include it because (A) it's my blog and I can do what I want and (B) in reaching "middle age", I increasingly find it helpful and fun to reflect on lived experiences in a public way that others (or me, in the future) may appreciate, learn from, or relate to. Perhaps this passing comment resonates with others and you can "pay it forward" in your own way.

## Introduction

Speaking of "paying it forward", the topic of this blog post is my effort to "pay it forward" on some work and resources I put together during my postdoc (and the pandemic). In reading this post you will learn how to "pay it forward" (I promise I won't write this anymore) to your future self by preparing, archiving, and resurrecting research data from Harvard University's tape storage.

First, some general context. Early on in my postdoc and before, PIs/labs at Harvard were granted generous storage quotas at low/no cost on the university high-performance computing cluster (current generation: Cannon). Logically, the financial burden on the university of this policy was becoming too much in the era of big data and during the pandemic, the decision was made to begin charging research labs for their longterm storage allocations. Given I was managing aspects of cluster maintenance for the Edwards lab, I ended up needing to implement a solution to the now-costly storage the lab was using. We devised a new data management plan that included (1) getting rid of a bunch of old data from former lab members after coordinating with these individuals, (2) reapproaching how we work with data on HPC to keep our data storage footprint smaller, and (3) taking advantage of tiered storage options for storing different types of data most cost-effectively, including utilizing the new tape storage option that was rolled out at the same time as the new pricing policies. This blog post will focus on #3 by outlining the mechanics of preparing data for tape storage, archiving the data to tape storage, and resurrecting data from tape storage. Perhaps I will write more about #1 and #2 in the future.

Paragraph on tape storage in general and specifically at Harvard

Disclaimer that I'm not a software engineer and have not rigorously tested everything. But given that no HPC staff-supported solution exists yet, this is better than nothing and is probably useful to others. Certainly to my fellow colleagues in the Edwards lab, but perhaps also to others at Harvard or other places.

## Preparing Data for Tape Storage

```bash
#!/bin/bash
#SBATCH -n 4                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH --constraint=intel  # only use faster intel cores
#SBATCH -t 2-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p edwards   # Partition to submit to
#SBATCH --mem=12000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH	-J curateData4Tape            # Job name for job
#SBATCH -o curateData4Tape_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e curateData4Tape_%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=FAIL,END        # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=dcard@fas.harvard.edu # Email to which notifications will be sent

# useful information
# SLURM environmental variables
# ${SLURM_CPUS_ON_NODE} = Number of cores/node
# ${SLURM_CPUS_PER_TASK} = Number of cores per task
# ${SLURM_JOB_ID} = Job ID
# ${SLURM_JOB_NAME} = Job Name
# ${SLURM_NTASKS} = Total number of cores for job
# ${SLURM_SUBMIT_DIR} = Submit Directory
# ${SLURM_SUBMIT_HOST} = Host submitted from

# version 0.2

# load modules and environments
conda activate zip

# commands to run
# create a full file manifest with sizes in bytes and megabytes
# and then identify files > 500 MB to be saved "as is" and write the relative subdirectory as column 10
# and then identify files < 500 MB to be 'tared' and write a placeholder tar file based on top-level directory as column 10
find . -mindepth 1 -type f ! -name "curateData4Tape*" ! -name "MANIFEST.txt" -exec ls -l --time-style=long-iso {} \; |
awk -v OFS="\t" -v OFMT='%f' '{ print $1, $2, $3, $4, $5, $5/1024/1024, $6, $7, $8 }' |
sort --parallel=${SLURM_NTASKS} -k9,9 |
while read line;
do
parent=`echo -e "${line}" | awk '{ print $9 }' | awk -F "/" '{ print $2 }'`;
subdir=`dirname $(echo -e "${line}" | awk '{ print $9 }')`;
echo -e "${line}" |
awk -v OFS="\t" -v parent="${parent}" -v subdir="${subdir}" '{ if ($6 > 500) print $0, "asis/"subdir; else if ($6 <= 500 && gsub(/\//,"/",$9) == 1) print $0, "tar/base.tar.gz"; else print $0, "tar/"parent".tar.gz" }'
done > MANIFEST.txt

# create a staging locations for the tape archive transfer
mkdir -p transfer_stage/tar transfer_stage/asis

# read the file manifest and for files that should be 'tared', create a tar.gz for each top-level directory with those files < 500 MB
cat MANIFEST.txt |
awk '{ if ($10 ~ "^tar") print $10 }' |
sort | uniq |
while read mytar;
do
echo $mytar;
cat MANIFEST.txt |
awk -v mytar="${mytar}" '{ if ($10 == mytar) print $9 }' |
tar --use-compress-program="pigz" -cvf transfer_stage/${mytar} --files-from -;
done

# read the file manifest and for files that can be stored as is, create a subdirectory under "transfer_stage/asis" based on the relative
# directory structure that the file is stored in (just in case of duplicates)
# then copy each of these files to this staging area for archival
cat MANIFEST.txt |
awk '{ if ($10 !~ "^tar") print $0 }' |
while read file;
do
source=`echo -e ${file} | awk '{ print $9 }'`;
destination=`echo -e ${file} | awk '{ print $10 }'`;
mkdir -p transfer_stage/${destination};
cp ${source} transfer_stage/${destination};
done

# generate MD5 checksums for all files in transfer_stage and write file alongside transfer_stage
cd transfer_stage
find . -type f -exec md5sum {} + > ../CHECKSUMS.md5.txt
cd ..
```

I will go through essentially each command/line and outline what is happening in detail. I will also provide examples of the inputs and outputs so you know what the data look like at different steps.

```bash
conda activate zip
```

This command activates a conda environment called "zip". This environment was provisioned with the utility pigz, which allows gzipping with more than one core. Therefore, ensure that pigz is properly installed and in your system `$PATH`.

```bash
find . -mindepth 1 -type f ! -name "curateData4Tape*" ! -name "MANIFEST.txt" -exec ls -l --time-style=long-iso {} \;
```

This command finds all files in the current working directory, except those that begin with `"curateData4Tape*"` (the name of the script) and `MANIFEST.txt` (the output of the script). Each of these files is listed using `ls` with the long listing format and a more detailed/parsed ISO string for the date/time. Here is a look at the output.

```
### Insert example
```

```bash
awk -v OFS="\t" -v OFMT='%f' '{ print $1, $2, $3, $4, $5, $5/1024/1024, $6, $7, $8 }'
```

This command adds a field after the 5th field (size of the file in bytes) recording the size of the file in MB.

```bash
sort --parallel=${SLURM_NTASKS} -k9,9
```

This command sorts by the last (9th) column, which is the relative path of each file. The sort is performed using multiple cores since the list of files can be quite large. The number of cores is inherited from the SLURM job parameters.

```bash
while read line;
do
...
done
```

This while loop reads each line and performs the operations within (`...`) on each line.

```bash
parent=`echo -e "${line}" | awk '{ print $9 }' | awk -F "/" '{ print $2 }'`;
```

This command writes an environmental variable recording the highest level parent directory of a given file (beneath the current working directory) after parsing the 9th field (relative file path).

```bash
subdir=`dirname $(echo -e "${line}" | awk '{ print $9 }')`;
```

This command writes an environmental variable recording the lowest level parent directory of a given file after parsing the 9th field (relative file path).

```bash
echo -e "${line}" |
```

The individual line (representing a single file) is now echoed to the terminal so it can be operated on. The output is piped forward.

```bash
awk -v OFS="\t" -v parent="${parent}" -v subdir="${subdir}" '{ if ($6 > 500) print $0, "asis/"subdir; else if ($6 <= 500 && gsub(/\//,"/",$9) == 1) print $0, "tar/base.tar.gz"; else print $0, "tar/"parent".tar.gz" }'
```

Each line is then operated on by an `awk` command. This awk command inherits the highest- and lowest-level directories of a file. It then evaluates whether the size of the file in MB exceeds 500. If it does, the file will be written as-is since the file size can be written natively to the tape. In this circumstance, an extra (10th) field is added with the prefix directory `asis/` and the lowest-level directory of the given file. If the size does not exceed 500 MB, this file and others like it within the same parent directory (i.e., highest-level directory beneath the current working directory) will be gathered together into a `tar` archive for storage, since the individual file size is too small to be written to tape. In this circumstance, two things can happen. If the file in question is not within a highest-level parent directory, a 10th field is added with the string for the downstream `.tar.gz` file called `base.tar.gz` for all base-level files. Note that this should be adjusted if there is a subdirectory called `base`. If the file in question is within a parent directory (this is ideal - try to place all files beneath at least one parent directory), add a 10th field with the string for the downstream `.tar.gz` file called `<parent>.tar.gz` where `<parent>` is the highest-level parent directory of the file. This 10th column is important for staging the files in the following commands.

```bash
mkdir -p transfer_stage/tar transfer_stage/asis
```

Makes a new directory called `transfer_stage/` where the files will be written during the preparation (i.e., "staging") process. Within this directory are two subdirectories: (1) `tar` for the `.tar.gz` files holding small (<500 MB) files and (2) `asis` for holding large (>500 MB) files as-is.

2nd command

```bash
cat MANIFEST.txt |
```

```bash
awk '{ if ($10 ~ "^tar") print $10 }' |
```

```bash
sort | uniq |
```

```bash
while read mytar;
do
...
done
```

```bash
echo $mytar;
```

```bash
cat MANIFEST.txt |
```

```bash
awk -v mytar="${mytar}" '{ if ($10 == mytar) print $9 }' |
```

```bash
tar --use-compress-program="pigz" -cvf transfer_stage/${mytar} --files-from -;
```

3rd command

```bash
cat MANIFEST.txt |
```

```bash
awk '{ if ($10 !~ "^tar") print $0 }' |
```

```bash
while read file;
do
...
done
```

```bash
source=`echo -e ${file} | awk '{ print $9 }'`;
```

```bash
destination=`echo -e ${file} | awk '{ print $10 }'`;
```

```bash
mkdir -p transfer_stage/${destination};
```

```bash
cp ${source} transfer_stage/${destination};
```

4th command

```bash
cd transfer_stage
```

```bash
find . -type f -exec md5sum {} + > ../CHECKSUMS.md5.txt
```

```bash
cd ..
```

## Archiving Data to Tape Storage

```bash
#!/bin/bash
#SBATCH -n 4                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH --constraint=intel  # only use faster intel cores
#SBATCH -t 2-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p edwards   # Partition to submit to
#SBATCH --mem=12000           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH	-J resurrectDataFromTape            # Job name for job
#SBATCH -o resurrectDataFromTape_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e resurrectDataFromTape_%j.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=FAIL,END        # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=dcard@fas.harvard.edu # Email to which notifications will be sent

# useful information
# SLURM environmental variables
# ${SLURM_CPUS_ON_NODE} = Number of cores/node
# ${SLURM_CPUS_PER_TASK} = Number of cores per task
# ${SLURM_JOB_ID} = Job ID
# ${SLURM_JOB_NAME} = Job Name
# ${SLURM_NTASKS} = Total number of cores for job
# ${SLURM_SUBMIT_DIR} = Submit Directory
# ${SLURM_SUBMIT_HOST} = Host submitted from

# load modules and environments
conda activate zip

# commands to run
# extract all tars in tar/
find tar/ -type f | sort |
while read tarfile;
do
tar -I pigz -xvf ${tarfile};
done

# copy over "as is" files to proper location after creating location (if necessary)
cat MANIFEST.txt |
awk '{ if ($10 !~ "^tar") print $0 }' |
while read line;
do
source=`echo -e "${line}" | awk '{ print $10 }'`;
destination=`dirname $(echo -e "${line}" | awk '{ print $9 }')`;
file=`basename $(echo -e "${line}" | awk '{ print $9 }')`;
mkdir -p ${destination}
cmd="cp ${source}/${file} ${destination}/${file}";
echo ${cmd};
eval ${cmd};
done
```

## Resurrecting Data from Tape Storage
