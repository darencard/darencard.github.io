---
layout: posts
title: "BaseSpace CLI Quickstart"
date: 2018-10-04
excerpt: "Installing, authenticating, and downloading using BaseSpace CLI."
---

## Installation on Mac computer

1. Install BaseSpace CLI to `$HOME/bin` directory and make executable.

```bash
wget "https://api.bintray.com/content/basespace/BaseSpaceCLI-EarlyAccess-BIN/latest/\$latest/amd64-osx/bs?bt_package=latest" \
-O $HOME/bin/bs
chmod u+x $HOME/bin/bs
```

2. Install the BaseSpace copy application to the `$HOME/bin` directory and make executable.

```bash
wget https://api.bintray.com/content/basespace/BaseSpace-Copy-BIN/\$latest/osx/bscp?bt_package=develop \
-O $HOME/bin/bs-cp
chmod u+x $HOME/bin/bs-cp
```

3. In your internet browser, log into BaseSpace through either your personal or a private domain (i.e., organization) login.

4. Once logged in in your browser, authenticate `bs`.

```bash
bs authenticate
```

Copy the URL supplied by the command and paste into your web browser, and since you are already signed, BaseSpace should prompt to accept terms and then properly authenticate. The command prompt in your shell will recycle when this is complete.

5. List available datasets in BaseSpace.

```bash
bs list datasets
```
For the sake of this tutorial, let's focus on the project `CVOS_WGS_J18`.

6. List all datasets associated with this project.

```bash
bs list datasets --terse -f csv --filter-field=Project.Name --filter-term=CVOS_WGS_J18
```

7. We can also list attributes associated with each dataset.

```bash
bs list datasets --terse -f csv --filter-field=Project.Name --filter-term=CVOS_WGS_J18 | \
while read id; do echo ${id}; bs list attributes dataset -i ${id}; done
```

8. Let's list the contents of each dataset as well.

```bash
bs list datasets --terse -f csv --filter-field=Project.Name --filter-term=CVOS_WGS_J18 | \
while read id; do echo ${id}; bs contents dataset -i ${id}; done
```

Here we can see paths for each file, but not sure how to use that.

9. We can download a dataset by doing the following. We will focus on dataset `ds.bcb4c9f5b0d34ccaa44859af9a5fb5e1` since it is relatively small.

```bash
# will take a second to start downloading
cd /path/to/store/data
bs download dataset -i ds.bcb4c9f5b0d34ccaa44859af9a5fb5e1 --extension=fastq.gz -o ./
```

10. Unfortunately, `bs download` does not download a MD5 checksum to verify the integrity of the data. (But maybe it checks it automatically as part of the download - something to ask Illumina). Fortunately another tool, `bs cp`, does download checksums.

```bash
cd /path/to/store/data
bs cp --write-md5 <basespace_location> ./
```

Just need to figure out how to determine the path on BaseSpace.
