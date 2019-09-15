---
layout: posts
title: "Introductory RAD-seq Activity"
date: 2019-09-15
excerpt: "Experimental design activity to accompany an introductory lesson on RAD-seq"
---

## Introduction

Today, let's go through a brief experimental design activity to help demonstrate some of the features of RAD-seq and how they can be harnessed to explore numerous questions in ecology, evolutionary biology, and genetics. To perform this activity, we will use Harvard's [Odyssey compute cluster]().

## Setting up Our Analysis Environment

We are going to use the very useful environments feature of [Anaconda Python]() to manage the software dependencies we have for this activity. Anaconda is already installed on Odyssey and can be loaded as follows:

```bash
module load Anaconda3/5.0.1-fasrc02
```

Now we can create an analysis environment that will always be available for us to use for this activity, which will run using Python version 3.7.

```bash
conda create -n radseq python=3.7
```

This will take a few minutes to run and along the way, we will answer Yes (y) to allow installation. When it is done, we can activate our newly created environment.

```bash
conda activate radseq
```

We must also install some other important software in our environment, namely [Biopython](). For whatever reason, installing `biopython` through `conda` is now allowing later software to work for me, so we will install it using `pip` instead, which is built in with Python.

```bash
pip install --user biopython
```

Now let's create a directory for our activity called `radseq_activity`.

```bash
mkdir radseq_activity
```

We are also going to be using some software that I wrote for performing *in silico* RAD-seq called [RADis](). I have, unfortunately, not packaged this into easily-installable software, but we can still make use of it fairly easily. Let's download the GitHub repository for this purpose.

```bash
git clone https://github.com/darencard/RADis
```

One more thing we will need to do is add the software in this repository to our software `$PATH`, which is where the computer looks for software to run by default.

```bash
cd RADis
export PATH=$PATH:`pwd`
```
