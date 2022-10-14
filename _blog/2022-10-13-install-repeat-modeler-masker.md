---
layout: posts
title: "Installing RepeatModeler & RepeatMasker"
date: 2022-10-13
excerpt: "A guide for installing popular repeat annotation software."
---

## Introduction
[RepeatModeler](https://www.repeatmasker.org/RepeatModeler/) and [RepeatMasker](https://www.repeatmasker.org/RepeatMasker/) are sister software packages for *de novo* identification of repeats in genomes and annotation / masking of repeats based on an existing library, respectively. These programs have been widely used for over a decade and are vital tools for genome annotation, though newer repeat pipelines are also emerging. Despite the popularity of these programs, they can often be difficult to use primarily because of the installation overhead. Because installation is a bottleneck and I often receive request for guidance in installing RepeatModeler / RepeatMasker, I decided to make a tutorial on installing both pieces of software. I have been using both for about 10 years and in thattime, the installation process has only changed in minor ways, so I'm confident this guide will be helpful well into the future, though there may be some subtle differences with future versions of the software.

## Software Information
Part of what makes RepeatModeler and RepeatMasker so hard to install is the dependencies. Here is a full list of the software that will be installed, the version information, and software sources. Keep in mind that software changes over time and some aspects of this tutorial may shift with future updates.

I will be installing all prerequisites and both RepeatMasker and RepeatModeler. For this demo, I am choosing to install RepeatModeler2 with the new LTR structural search pipeline, which adds several additional prerequisites. I am also going to install these software with both the [Dfam](https://dfam.org/) and [Repbase](https://www.girinst.org/repbase/) repeat libraries.

### Prerequisites
1. Perl 5.8.0 or higher - RepeatModeler was tested using Perl 5.8.8 so I suggest having at least that version. I will assume the user has met this prerequisite as most systems will probably have this installed by default. But you may need to install or update so work with local system administrators as needed. Software environments can get quite complex and this can lead to issues, which I cannot support.
2. Python 3 and `h5py` library - Again, I will assume the user can install this without issue but the same caveats I noted about Perl also apply here. You can probably use `pip` or `conda` to manage this installation (as well as Perl). Again, due to the complexity of software environments, which vary system to system, I cannot provide support for this installation. There is lots of available information if you Google.
3. [Tandem Repeat Finder (TRF)](https://github.com/Benson-Genomics-Lab/TRF) - version 4.09.1
4. [RECON](http://eddylab.org/software/recon/) - version 1.08
5. RepeatScout (no software website) - version 1.0.6
6. [RMBlast](http://www.repeatmasker.org/RMBlast.html) - version 2.11.0
7. [GenomeTools](http://genometools.org/pub/) - version 1.6.2 (required for `ltrharvest` for running new RepeatModeler2 LTR structural search pipeline)
8. [LTR_retriever](https://github.com/oushujun/LTR_retriever) - version 2.9.0 (required for running new RepeatModeler2 LTR structural search pipeline)
9. [MAFFT](https://mafft.cbrc.jp/alignment/software/) - version 7.505 (required for running new RepeatModeler2 LTR structural search pipeline)
10. [CD-HIT](http://cd-hit.org/) - version 4.8.1 (required for running new RepeatModeler2 LTR structural search pipeline)
11. [Ninja](https://github.com/TravisWheelerLab/NINJA) - version 0.95-cluster_only (required for running new RepeatModeler2 LTR structural search pipeline)
12. Note: I typically only use the BLAST search engine but there are also options for [Cross_Match](http://www.phrap.org/phredphrapconsed.html), [HMMER](http://hmmer.org/), and [ABBlast/WUBlast](https://blast.advbiocomp.com/licensing/). I will not describe installing these engines but the process is generally the same as below once you install those software prerequisites.

### Core Software
1. [RepeatMasker](https://www.repeatmasker.org/RepeatMasker/) - version 4.1.3-p1
2. [RepeatModeler](https://www.repeatmasker.org/RepeatModeler/) - version 2.0.3

### Databases
1. [Dfam](https://dfam.org/home) - release 3.6
2. [Repbase](https://www.girinst.org/repbase/) - release 20181026. Note: This was the final release of Repbase before GIRI, unfortunately, put this resource behind a paywall. For those who already have this release, you should be able to install the software with Repbase. If you have a lucky person with individual or institutional subscriptions, you may have access to updated Repbase releases and can use those instead. I will describe installing release 20181026 but the process should be the same for other releases. You can always proceed without Repbase and only use Dfam if access is an issue.

## Installation

I will generally be installing the prerequisite software first before proceeding with installing RepeatMasker and RepeatModeler (in that order).

Let's perform this installation in our home directory. You can always perform the same process somewhere else.

```bash
cd $HOME
```

I will make a directory called `repeat-annotation` to hold the installation, but you could name it alternatively.

```bash
mkdir repeat-annotation
cd repeat-annotation
```

Now let's download the prerequisite software we need. Generally, I can do this via the command line but for some software/database files, you will need to download to your local computer and transfer to the computer you are installing on, if you are working with a remote system. I am installing using the [Harvard HPC Canon](https://www.rc.fas.harvard.edu/about/cluster-architecture/), which is a 64-bit Unix system. If you are using another operating system, some of the decisions you need to make may vary. I will go through the prerequisite installations one-by-one.

### Prerequisites

#### TRF

```bash
# download 64-bit Linux binaries of TRF
wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
# given this is a binary, no code needs to be compiled. just need to make file executable
chmod +x trf409.linux64
# create a symlink to 'trf'
ln -s trf409.linux64 trf
# test installation works correctly (should see help details)
./trf -h
```

#### RECON

```bash
# download source code of RECON
wget http://www.repeatmasker.org/RepeatModeler/RECON-1.08.tar.gz
# extract file contents
tar xvf RECON-1.08.tar.gz
# navigate to source code
cd RECON-1.08/src
# compile code
make
make install
# compiled software is now available in RECON-1.08/bin
# move back to parent directory
cd $HOME/repeat-annotation
```

#### RepeatScout

```bash
# download source code of RepeatScout
wget http://www.repeatmasker.org/RepeatScout-1.0.6.tar.gz
# extract file contents
tar xvf RepeatScout-1.0.6.tar.gz
# navigate to source code
cd RepeatScout-1.0.6
# compile code
make
# compiled software is now available in RepeatScout-1.0.6
# move back to parent directory
cd $HOME/repeat-annotation
```

#### RMBlast

```bash
# download 64-bit Linux binaries of RMBlast
wget http://www.repeatmasker.org/rmblast-2.11.0+-x64-linux.tar.gz
# extract file contents
tar xvf rmblast-2.11.0+-x64-linux.tar.gz
# given we downloaded binaries, software is available in rmblast-2.11.0/bin
# move back to parent directory
cd $HOME/repeat-annotation
```

#### GenomeTools (ltrharvest)

```bash
# download source code of genometools
wget http://genometools.org/pub/genometools-1.6.2.tar.gz
# extract file contents
tar xvf genometools-1.6.2.tar.gz
# navigate to source code
cd genometools-1.6.2
# compile with multithreading
make threads=yes
# compiled software is now available in ???
# move back to parent directory
cd $HOME/repeat-annotation
```

#### LTR_retriever

```bash
# download source code and binaries of LTR_retriever
wget -qO- https://github.com/oushujun/LTR_retriever/archive/refs/tags/v2.9.0.tar.gz > ltr_retriever_v2.9.0.tar.gz
# extract file contents
tar xvf ltr_retriever_v2.9.0.tar.gz
# compiled software is now available in LTR_retriever-2.9.0
# move back to parent directory
cd $HOME/repeat-annotation
```

#### MAFFT

```bash
# download source code of MAFFT
wget --no-check-certificate https://mafft.cbrc.jp/alignment/software/mafft-7.490-with-extensions-src.tgz
# extract file contents
tar xvf mafft-7.490-with-extensions-src.tgz
# navigate to source code
cd mafft-7.490-with-extensions
# compile core and extension code
cd core
# first you must modify 'Makefile' and change the first line from 'PREFIX = /usr/local' to 'PREFIX = /n/home13/dcard/repeat-annotation/mafft-7.490-with-extensions' (modify accordingly)
make clean
make
cd ..
cd extensions
# first you must modify 'Makefile' and change the first line from 'PREFIX = /usr/local' to 'PREFIX = /n/home13/dcard/repeat-annotation/mafft-7.490-with-extensions' (modify accordingly)
make clean
make
# compiled software is now available in mafft-7.490-with-extensions/bin
# move back to parent directory
cd $HOME/repeat-annotation
```

#### CD-HIT

```bash
# download source code of CD-HIT
wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
# extract file contents
tar xvf cd-hit-v4.8.1-2019-0228.tar.gz
# navigate to source code
cd cd-hit-v4.8.1-2019-0228
# compile code
make
# might as well compile the auxtools as well
cd cd-hit-auxtools
make
# compiled software is now available in cd-hit-v4.8.1-2019-0228
# move back to parent directory
cd $HOME/repeat-annotation
```

#### Ninja

```bash
# download source code of Ninja
wget -qO- https://github.com/TravisWheelerLab/NINJA/archive/refs/tags/0.95-cluster_only.tar.gz > Ninja_0.95-cluster_only.tar.gz
# extract file contents
tar xvf Ninja_0.95-cluster_only.tar.gz
# navigate to source code
cd NINJA-0.95-cluster_only/NINJA
# compile code
make
# compiled software is now available in cd-hit-v4.8.1-2019-0228
# move back to parent directory
cd $HOME/repeat-annotation
```

### RepeatMasker & RepeatModeler

#### RepeatMasker

```bash
# ensure we are in correct location
cd $HOME/repeat-annotation
# download source code of RepeatMasker
wget https://www.repeatmasker.org/RepeatMasker/RepeatMasker-4.1.3-p1.tar.gz
# extract file contents
tar xvf RepeatMasker-4.1.3-p1.tar.gz
```

Before we proceed with installation, we need to prepare the repeat annotation libraries from Dfam and Repbase. Release 3.6 of Dfam should be distributed with RepeatMasker v. 4.1.3-p1, but let's double check.

```bash
head -25 RepeatMasker/Libraries/Dfam.h5
# you should see the release number noted - ensure it is 3.6
```

If you are viewing this in the future, you can always download the newest version of Dfam and replace the `Dfam.h5` file in `RepeatMasker/Libraries/` with the newly downloaded release file.

I will not describe how to retrieve Repbase since it is behind a paywall and many may not have access. You should download the correct release fill formatted for RepeatMasker. For our release, that is `RepBaseRepeatMaskerEdition-20181026.tar.gz`, but the file name may vary if you are downloading other releases. We need to extract the Repbase release files to `RepeatMasker/Libraries/`.

```bash
# move RepBaseRepeatMaskerEdition-20181026.tar.gz to the extracted RepeatMasker directory
mv /path/to/RepBaseRepeatMaskerEdition-20181026.tar.gz RepeatMasker/
# navigate to source code
cd RepeatMasker
# extract file contents to overwrite existing release distributed with RepeatMasker
tar xvf RepBaseRepeatMaskerEdition-20181026.tar.gz
```

Now that all prerequisites have been installed and libraries have been properly staged, we can proceed with installation of RepeatMasker. This is an interactive process and the installer will prompt for the locations of installed software.

```bash
# should be located in $HOME/repeat-annotation/RepeatMasker/
perl ./configure
```

You will see some output like the following as RepeatMasker prepares the repeat libraries from Dfam and Repbase. Note that both are reflected in the output below, which is what we want. Keep in mind that `$HOME` for me is `/n/home13/dcard/` and adjust your paths accordingly below.

```
RepeatMasker Configuration Program


Checking for libraries...

Rebuilding RepeatMaskerLib.h5 master library
  - Read in 49011 sequences from /n/home13/dcard/repeat-annotation/RepeatMasker/Libraries/RMRBSeqs.embl
  - Read in 49011 annotations from /n/home13/dcard/repeat-annotation/RepeatMasker/Libraries/RMRBMeta.embl
  Merging Dfam + RepBase into RepeatMaskerLib.h5 library............................
```

Once that is complete, we will begin receiving prompts from the installer for the location of critical software.

```
The full path including the name for the TRF program.
TRF_PRGM:
```

All we must do is provide the path to TRF. However, we are stuck in the installer, which makes that difficult. You can always stop the job controlling the installer by executing `CTRL+X`, navigate to find the information you need, and then resume the installation using `fg` ('foreground' the job) to enter the path. You can always start the installation from the beginning again if you are having problems or make a mistake. If you are following my tutorial exactly, then the following paths should also work for you.

```
TRF_PRGM: /n/home13/dcard/repeat-annotation/trf
```

Then we receive options for the search engine that will be used. We are going to setup RMBlast.

```
Add a Search Engine:
   1. Crossmatch: [ Un-configured ]
   2. RMBlast: [ Un-configured ]
   3. HMMER3.1 & DFAM: [ Un-configured ]
   4. ABBlast: [ Un-configured ]

   5. Done


Enter Selection: 2
```

Now we are prompted to enter the path to the RMBlast software. A suggested default path may be displayed but you can enter a new one.

```
The path to the installation of the RMBLAST sequence alignment program.
RMBLAST_DIR [/n/home13/dcard/bin]: /n/home13/dcard/repeat-annotation/rmblast-2.11.0/bin
```

Then we are asked if we want to make this our default search engine, which we do.

```
Do you want RMBlast to be your default
search engine for Repeatmasker? (Y/N)  [ Y ]: Y
```

Since we are not installing any other search engines, we can end here. If you wanted to install others, this is where you would do so by pointing the installer to the path of the software like we did with RMBlast.

```
Add a Search Engine:
   1. Crossmatch: [ Un-configured ]
   2. RMBlast: [ Configured, Default ]
   3. HMMER3.1 & DFAM: [ Un-configured ]
   4. ABBlast: [ Un-configured ]

   5. Done


Enter Selection: 5
```

The installer will continue to work building the repeat libraries before finishing the installation

```
Building FASTA version of RepeatMasker.lib .................................
Building RMBlast frozen libraries..
The program is installed with a the following repeat libraries:
File: /n/home13/dcard/repeat-annotation/RepeatMasker/Libraries/RepeatMaskerLib.h5
Database: Dfam withRBRM
Version: 3.6
Date: 2022-04-12

Dfam - A database of transposable element (TE) sequence alignments and HMMs.
RBRM - RepBase RepeatMasker Edition - version 20181026

Total consensus sequences: 63852
Total HMMs: 18987


Further documentation on the program may be found here:
  /n/home13/dcard/repeat-annotation/RepeatMasker/repeatmasker.help
```

Now you should be able to execute `./RepeatMasker -h` from ` ` and see the help page, which verifies the installation is complete. Unfortunately for me, a required Perl module was not initially available and so I was getting the following error.

```
Can't locate EMBL.pm in @INC (you may need to install the EMBL module) (@INC contains: ...
```

This results from the fact that the `EMBL.pm` module distributed with RepeatMasker is not automatically available in the Perl module path. This underscores the difficulty of dealing with software environments for Perl/Python/etc. I overcome the issue by adding the path to the RepeatMasker files to the `PERL5LIB` path. The same fix may work if there are other Perl modules missing that are distributed with RepeatMasker (check for `<NAME>.pm` in directory).

```bash
export PERL5LIB=$PERL5LIB/:/n/home13/dcard/repeat-annotation/RepeatMasker
```

To have this done automatically when you login, you will need to add this line to your existing `$HOME/.bash_profile` file or similar, depending on your system. Consult Google or a knowledgeable colleague for assistance.

Great! Now that RepeatMasker is installed, we can try RepeatModeler!

#### RepeatModeler

```bash
# ensure we are in correct location
cd $HOME/repeat-annotation
# download source code of RepeatModeler
wget https://www.repeatmasker.org/RepeatModeler/RepeatModeler-2.0.3.tar.gz
# extract file contents
tar xvf RepeatModeler-2.0.3.tar.gz
# navigate to source code
cd RepeatModeler-2.0.3
```

All prerequisite software has been installed and the repeat libraries are properly configured through RepeatMasker, so we can proceed directly with installing RepeatModeler.

```bash
# should be located in $HOME/repeat-annotation/RepeatModeler-2.0.3/
perl ./configure
```

We are first met with a generic introductory prompt, where we can press `<Enter>` to proceed.

```
RepeatModeler Configuration Program

This program assists with the configuration of the
RepeatModeler program.  The next set of screens will ask
you to enter information pertaining to your system
configuration.  At the end of the program your RepeatModeler
installation will be ready to use.



<PRESS ENTER TO CONTINUE>
```

Next, the installer prompts for the installation location for Perl. In my experience, this does a pretty good job of automatically identifying whatever version of Perl is installed, so I almost always go with the suggested path by hitting `<Enter>`. Just ensure the Perl version is appropriate based on the prerequisites above.

```
**PERL INSTALLATION PATH**

  This is the full path to the Perl interpreter.
  ie. /usr/local/bin/perl

Enter path [ /n/sw/helmod/apps/centos7/Core/perl/5.26.1-fasrc01/bin/perl ]: <Enter>
```

Now we are prompted to direct the installer to the installation location for RepeatMasker, which we just installed. It may find this installation automatically like it did for me, but you can always enter the path (replacing `/n/home13/dcard` with your `$HOME`: `/n/home13/dcard/repeat-annotation/RepeatMasker`.

```
The path to the installation of RepeatMasker.
REPEATMASKER_DIR [/n/home13/dcard/repeat-annotation/RepeatMasker]:
```

Next, we must enter the path to our installation of RECON, as below.

```
The path to the installation of the RECON de-novo repeatfinding program.
RECON_DIR: /n/home13/dcard/repeat-annotation/RECON-1.08/bin
```

And the same with RepeatScout, as below.

```
The path to the installation of the RepeatScout ( 1.0.6 or higher ) de-novo repeatfinding program.
RSCOUT_DIR: /n/home13/dcard/repeat-annotation/RepeatScout-1.0.6
```

And with TRF. However, instead of the full path to the software like we provided to RepeatMasker, we only provide the path to where `trf` is installed.

```
The full path to TRF program.  TRF must be named \"trf\". ( 4.0.9 or higher )
TRF_DIR: /n/home13/dcard/repeat-annotation
```

And with CD-HIT.

```
The path to the installation of the CD-Hit sequence clustering package.
CDHIT_DIR: /n/home13/dcard/repeat-annotation/cd-hit-v4.8.1-2019-0228
```

Next the installer asks for the location of various 'TwoBit' tools from UCSC (sometimes referred to as the Kent Utilities). I do not know why these prerequesites were not mentioned anywhere in the documentation, so we will have to install these on the fly. Part of the joy of installing software! Fortunately, this is pretty easy to do since binaries for all UCSC tools are available to download directly from `https://hgdownload.soe.ucsc.edu/admin/exe/`. I am downloading those in the directory `linux.x86_64/` but you may need to select another option. To perform the installation, I would login or open a second terminal so you do not have to interupt the RepeatModeler installation in the existing terminal.

```bash
# navigate to correct location
cd $HOME/repeat-annotation
# download the tools with "twobit" in their name
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToTwoBit
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitDup
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitInfo
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitMask
wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
# make these tools executable
chmod +x faToTwoBit twoBitDup twoBitInfo twoBitMask twoBitToFa
# you can test these work as follows - you should see the help information for each tool
./faToTwoBit -h
```

Now, back in our other terminal with the RepeatModeler installer, we can provide the path to these tools.

```
The path to the installation directory of the UCSC TwoBit Tools (twoBitToFa, faToTwoBit, twoBitInfo etc).
UCSCTOOLS_DIR: /n/home13/dcard/repeat-annotation
```

Next, we can add our search engine, which will be RMBlast again.

```
Add a Search Engine:
   1. RMBlast - NCBI Blast with RepeatMasker extensions: [ Un-configured ]
   2. WUBlast/ABBlast: [ Un-configured ]

   3. Done


Enter Selection: 1
```

And enter the correct path.

```
The path to the installation of the RMBLAST sequence alignment program.
RMBLAST_DIR: /n/home13/dcard/repeat-annotation/rmblast-2.11.0/bin
```

Assuming that works properly, we can move on.

```
Add a Search Engine:
   1. RMBlast - NCBI Blast with RepeatMasker extensions: [ Configured ]
   2. WUBlast/ABBlast: [ Un-configured ]

   3. Done


Enter Selection: 3
```

Now we are prompted for the prerequisite software necessary to use the LTR Structure Identification Pipeline.

```
LTR Structural Identication Pipeline [optional]

In addition to RECON/RepeatScout this version of RepeatModeler
has the option of running an additional analysis to identify
structural features of LTR sequences.

Do you wish to configure RepeatModeler for this type
of analysis [y] or n?: y
```

We begin with the path to GenomeTools.

```
The path to the installation of the GenomeTools package.
GENOMETOOLS_DIR: /n/home13/dcard/repeat-annotation/genometools-1.6.2/bin
```

And then to LTR_retriever.

```
The path to the installation of the LTR_Retriever (v2.9.0 and higher) structural LTR analysis package.
LTR_RETRIEVER_DIR: /n/home13/dcard/repeat-annotation/LTR_retriever-2.9.0
```

And then to MAFFT.

```
The path to the installation of the MAFFT multiple alignment program.
MAFFT_DIR: /n/home13/dcard/repeat-annotation/mafft-7.490-with-extensions/bin
```

And, finally, Ninja.

```
The path to the installation of the Ninja phylogenetic analysis package.
NINJA_DIR: /n/home13/dcard/repeat-annotation/NINJA-0.95-cluster_only/NINJA
```

If all went well, the installation should finish with the following.

```
 -- Setting perl interpreter and version...

Congratulations!  RepeatModeler is now ready to use.
```

You can test that your installation works with the following command. You should see the RepeatModeler help information if all is well!

```bash
./RepeatModeler -h
```

And that's it! Now both RepeatMasker and RepeatModeler, along with necessary software dependencies and repeat libraries, should be installed properly and available to use. If you need any guidance for running RepeatModeler or RepeatMasker, you can visit some of my other [blog tutorials](https://darencard.net/blog/2022-07-09-genome-repeat-annotation/).

I also recommend adding the directories with software to your `$PATH` so that the software is accessible from anywhere. Otherwise, you must use the full path to each program when you want to run it. There is lots of guidance online for doing this. Both RepeatMasker and RepeatModeler also have a subdirectory `util` with useful support scripts/programs that you can also add to your `$PATH`.
