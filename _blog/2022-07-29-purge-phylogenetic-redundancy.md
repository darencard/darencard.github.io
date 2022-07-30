---
layout: posts
title: "Purging Phylogenetic Redundancy"
date: 2022-07-29
excerpt: "Subsetting a phylogeny non-randomly based on genetic distance."
---

## Introduction

Large phylogenies with many tips can lead to downstream bottlenecks in computation, analysis, and visualization. I was recently working to extract data from a very large whole-genome alignment produced by the [Zoonomia Consortium](https://zoonomiaproject.org/). The original alignment data are described in [this paper](https://www.nature.com/articles/s41586-020-2876-6), were created using the [Progressive Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) method described in [this paper](https://www.nature.com/articles/s41586-020-2871-y), and are available online from [this webpage](https://cglgenomics.ucsc.edu/data/cactus/). Specifically, I was using the very large 605-way combined mammalian and avian alignment HAL file for some comparative genomics work I am doing. A really nice feature of the new Cactus method of aligning and the HAL format is the phylogenetic framework that the genomic data are incorporated into. 

While not as large as other phylogenies, a 605-tip tree is pretty large and extracting alignment data from the HAL file for all of these species has a very large computational burden. Specifically, I have been working on estimating a neutral-rates phylogeny for the species included in this tree so I could attempt to identify genomic loci with various patterns of evolutionary conservation in specific lineages of mammals or birds. However, performing this task on the full alignment of all species has not been possible for me, so I was interested in reducing my list of species so this process would be more computationally manageable. However, I do not want to simply remove tips randomly from the tree since this could result in a biased representation of mammalian/avian diversity and evolution that would make the resulting inference of neutral rates problematic, potentially leading to biased results.

Ideally, one subsets the tree in a non-random way instead, resulting in a smaller tree that still encompasses major lineages in the phylogeny that can be used to produce a non-biased inference of neutral rates of evolution across mammals and birds. The idea here is to prioritize long branches representing larger amounts of evolutionary distance, and thus "uniqueness", over shorter branches that simply fill in greater amounts of diversity within a lineage of organisms that are more closely related and probably more biologically similar. Hopefully, this makes intuitive sense. 

I could always manually produce a subsetted phylogeny based on my own knowledge and opinions of mammal and bird diversity by, say, taking a random representative species for each family. However, this could also lead to biases and it ends up being time consuming, so I was really interested in an automated way of reducing redundancy in my phylogeny in a non-biased way. This seemed like a simple task and I anticipated there would be some tools to perform this subsetting but, unfortunately, I had little luck finding anything. If anyone is aware of something I missed, please let me know, but the only tool I found that did what I envisioned is called [Treemer](https://github.com/fmenardo/Treemmer), which is published [here](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2164-8). What Treemer attempts to do is to subset a phylogeny by identifying the pairs of tips that have the smallest evolutionary distance and randomly removing one of those tips. Performing this process iteratively allows one to subset the phylogeny while maintaining the more diverged/unique taxa and, thus, a more realistic representation of all diversity in the dataset.

While Treemer seemed to be what I wanted, I quickly noticed that the documentation is minimal. I saw there were dependencies, but it would have been up to me to look into how to install these properly. I guess I'm lazy because this proved to be too much for me - a lesson in ensuring your software is well documented and easy to setup/use if you wish for others to use it! Given the subsetting task seemed pretty easy and intuitive, I decided to implement the approach myself. It turned out to be straight-forward, and I document how I did this below.

## Estimating Evolutionary Distances

The key to subsetting a phylogeny in this use case is estimating evolutionary distances from an existing phylogeny. My phylogeny is stored as a Newick file, which I exported from the 605-way HAL alignment file. I used [HAL software](https://github.com/ComparativeGenomicsToolkit/hal) `halStats` command to do this.

```bash
singularity exec --cleanenv /n/singularity_images/informatics/cactus/cactus_v2.0.1.sif halStats --tree 605-vertebrate-2020.hal > 605-vertebrate-2020_species_phylo.tree
```

I ran this software out of a Singularity container, but you could install the standalone software, in which case the above command would be:

```bash
halStats --tree 605-vertebrate-2020.hal > 605-vertebrate-2020_species_phylo.tree
```

There are many options for calculating evolutionary distances from a tree so there are many ways someone could accomplish this goal. Given I do most of my data processing in the Unix shell, I decided to use [gotree](https://github.com/evolbioinfo/gotree), which is newer software written in Go for summarizing and manipulating phylogenetic trees. I have used `gotree` a lot lately after discovering it and highly recommend it. There is another useful piece of software from the same developer called [goalign](https://github.com/evolbioinfo/goalign) for working with molecular sequence alignments that I also recommend. We can use `gotree` to quickly visualize our phylogeny directly in the terminal.

```bash
cat 04_605-vertebrate-2020_species_phylo.tree | gotree draw text -w 100
```

It is also possible to render the tree with the node labels displayed.

```bash
cat 605-vertebrate-2020_species_phylo.tree | gotree draw text --with-node-labels -w 100
```

gotree can do a lot with phylogenetic trees, but I will not demo all the capability. However, as a useful example, gotree can be used to extract a subphylogeny from this full mammal/bird phylogeny.

```bash
# extract the bird subphylogeny from the mammal/bird phylogeny
# do so by specifying the node label for the root of the bird phylogeny
# unfortunately, the bird root is poorly named 'root' in the 605-way alignment
# sed command was used to convert 'root' to 'birdroot'
# gotree subtree command was used to grab the bird phylogeny
# and then another sed command was used to convert 'birdroot' back to 'root'
singularity exec --cleanenv /n/singularity_images/informatics/cactus/cactus_v2.0.1.sif halStats --tree 605-vertebrate-2020.hal | sed 's/root\:/birdroot:/g' | gotree subtree -n birdroot | sed 's/birdroot\;/root;/g' > bird_phylo.tree
```

With some introduction out of the way, let's look at how we extract pairwise evolutionary distances (i.e., the branch lengths) between taxa/tips in the phylogeny.

```bash
cat 605-vertebrate-2020_species_phylo.tree | gotree matrix
```

This spits out a giant matrix of distances, so let's look at a part of that to visualize the format. Here, we will take the first 10 lines and the first 10 columns of output.

```bash
cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | head | cut -f 1-10
```

Here is the output.

```bash
605
Nycticebus_coucang	0.000000000000	0.092236200000	0.191195200000	0.212361990000	0.213949390000	0.215183790000	0.2288136900000.227023190000	0.211468000000
Otolemur_garnettii	0.092236200000	0.000000000000	0.194483800000	0.215650590000	0.217237990000	0.218472390000	0.2321022900000.230311790000	0.214756600000
Daubentonia_madagascariensis	0.191195200000	0.194483800000	0.000000000000	0.119177790000	0.120765190000	0.121999590000	0.135629490000	0.133838990000	0.118283800000
Propithecus_coquereli	0.212361990000	0.215650590000	0.119177790000	0.000000000000	0.047710400000	0.084886600000	0.0985165000000.096726000000	0.086588990000
Indri_indri	0.213949390000	0.217237990000	0.120765190000	0.047710400000	0.000000000000	0.086474000000	0.100103900000	0.098313400000	0.088176390000
Cheirogaleus_medius	0.215183790000	0.218472390000	0.121999590000	0.084886600000	0.086474000000	0.000000000000	0.0734099000000.071619400000	0.089410790000
Microcebus_murinus	0.228813690000	0.232102290000	0.135629490000	0.098516500000	0.100103900000	0.073409900000	0.0000000000000.044459300000	0.103040690000
Mirza_coquereli	0.227023190000	0.230311790000	0.133838990000	0.096726000000	0.098313400000	0.071619400000	0.044459300000	0.000000000000	0.101250190000
Lemur_catta	0.211468000000	0.214756600000	0.118283800000	0.086588990000	0.088176390000	0.089410790000	0.103040690000	0.101250190000	0.000000000000
```

As you can see, the first line prints the number of tips in the phylogeny. Then there are individual lines for each tip/species in the tree, with the first column displaying the taxon ID from the phylogeny. Then there are N additional columns that display the pairwise evolutionary distance (branch lengths between these two tips/species), where N is the number of tips in the phylogeny, which is 605 in this example. The ordering of the columns is exactly the same as the ordering of the rows, which is why we see a diagonal of `0.000000000000` values running through the matrix for when a given tip is compared to itself. This is the foundation for everything we need to do downstream, but the format of this matrix is not immediately useful to extract the information that we need.

## Gathering and Subsetting the Pair of Taxa with the Smallest Evolutionary/Phylogenetic Distance

We could extract the information we need in many ways, but here is how I did it. First, I trimmed off the top line with the number of taxa and then I extracted the first row of species/tip IDs.

```bash
cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | tail -n +2 | awk '{ print $1 }'
```

And for each of these rows, I want to extract both the index (i.e., the column number) and the pairwise distance estimate from the matrix. I do this using a `while` loop where for each line after the 1st line, I:

A. transpose the line so it runs down instead of right-to-left using `datamash transpose`;

B. remove the first line, which holds the species/tip ID, using `tail -n +2`; 

C. print each distance and a sequential index from 1 to N using `awk -v OFS="\t" '{ print NR, $0 }'`; 

D. sort the resulting data by genetic/phylogenetic distance from low to high using `sort -k2,2n`; and 

E. extract the 2nd row of the sorted data as the shortest evolutionary distance involving this species/tip, since the first row reflects the `0.000000000000` value of the given tip/species vs. itself. 

Here is a look at the full command.

```bash
cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | tail -n +2 | while read line;
do echo -e "${line}" | datamash transpose | tail -n +2 |
awk -v OFS="\t" '{ print NR, $0 }' | sort -k2,2n | sed -n '2,2p';
done
```

I can paste the two above commands together.

```bash
paste <(cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | tail -n +2 | awk '{ print $1 }') \
<(cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | tail -n +2 | while read line;
do echo -e "${line}" | datamash transpose | tail -n +2 |
awk -v OFS="\t" '{ print NR, $0 }' | sort -k2,2n | sed -n '2,2p';
done)
```

This gives us output in the following form, with the first column being the taxon name for each row, the 2nd column being the column index of the taxon with the lowest pairwise phylogenetic/evolutionary distance with the taxon in column 1, and the 3rd column being the measured lowest pairwise phylogenetic/evolutionary distance between these two taxa.

```
Nycticebus_coucang	2	0.092236200000
Otolemur_garnettii	1	0.092236200000
Daubentonia_madagascariensis	9	0.118283800000
Propithecus_coquereli	5	0.047710400000
Indri_indri	4	0.047710400000
Cheirogaleus_medius	8	0.071619400000
Microcebus_murinus	8	0.044459300000
Mirza_coquereli	7	0.044459300000
Lemur_catta	10	0.049048900000
Eulemur_fulvus	11	0.014128410000
```

The first two taxa, `Nycticebus_coucang` and `Otolemur_garnettii` are great to focus on as an example, as they are each other's closest relatives in the phylogeny. To see what I mean, take a look at the output of:


```bash
singularity exec --cleanenv /n/singularity_images/informatics/cactus/cactus_v2.0.1.sif halStats --tree ../605-vertebrate-2020.hal | gotree draw text -w 100 | head
```

This is why the column index for `Nycticebus_coucang` is `2`, which is `Otolemur_garnettii`, and vice versa. From here, we can use the column index to extract the species/tip ID for the closest species in the phylogeny. I do this by looping through the lines in the resulting data using another `while` loop and for each line I:

A. extract the index for the closest species/tip using `taxindex=$(echo -e "${distance}" | awk '{ print $2 }')`; 

B. calculate the pairwise species matrix and gather the first column of species/tip IDs before removing the first row (number of taxa) using `cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | cut -f 1 | tail -n +2`; 

C. excise a given row of the resulting list of species/tips using the row index to gather the species/tip ID corresponding with the closest species/tip using the variable `${taxindex}` and `sed -n "${taxindex},${taxindex}p"`; and 

D. print this species/tip ID alongside the original distance information using `awk -v OFS="\t" -v distance="${distance}" '{ print distance, $0 }'`.

```bash
paste <(cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | tail -n +2 | awk '{ print $1 }') \
<(cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | tail -n +2 | while read line;
do echo -e "${line}" | datamash transpose | tail -n +2 |
awk -v OFS="\t" '{ print NR, $0 }' | sort -k2,2n | sed -n '2,2p';
done) |
while read distance;
do taxindex=$(echo -e "${distance}" | awk '{ print $2 }');
cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | cut -f 1 | tail -n +2 | sed -n "${taxindex},${taxindex}p" |
awk -v OFS="\t" -v distance="${distance}" '{ print distance, $0 }';
done
```

This results in an output that looks like this, where the first 3 column are the output from the above command and the final column is the inferred closest species/tip ID based on the column matrix in column 2. For the examples of `Nycticebus_coucang` and `Otolemur_garnettii`, we now see the corresponding species/tip IDs for the closest species/tip ID of `Otolemur_garnettii` and `Nycticebus_coucang`, respectively.

```
Nycticebus_coucang	2	0.092236200000	Otolemur_garnettii
Otolemur_garnettii	1	0.092236200000	Nycticebus_coucang
Daubentonia_madagascariensis	9	0.118283800000	Lemur_catta
Propithecus_coquereli	5	0.047710400000	Indri_indri
Indri_indri	4	0.047710400000	Propithecus_coquereli
Cheirogaleus_medius	8	0.071619400000	Mirza_coquereli
Microcebus_murinus	8	0.044459300000	Mirza_coquereli
Mirza_coquereli	7	0.044459300000	Microcebus_murinus
Lemur_catta	10	0.049048900000	Eulemur_fulvus
Eulemur_fulvus	11	0.014128410000	Eulemur_flavifrons
```

Finally, let's rearrange the output a little and also add the index for each species in the first column, which is simply the row index since the species/tip IDs are written in order. We can use a little awk to do this: `awk -v OFS="\t" '{ print NR, $1, $2, $4, $3 }'`.

```bash
paste <(cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | tail -n +2 | awk '{ print $1 }') \
<(cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | tail -n +2 | while read line;
do echo -e "${line}" | datamash transpose | tail -n +2 |
awk -v OFS="\t" '{ print NR, $0 }' | sort -k2,2n | sed -n '2,2p';
done) |
while read distance;
do taxindex=$(echo -e "${distance}" | awk '{ print $2 }');
cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | cut -f 1 | tail -n +2 | sed -n "${taxindex},${taxindex}p" |
awk -v OFS="\t" -v distance="${distance}" '{ print distance, $0 }';
done | 
awk -v OFS="\t" '{ print NR, $1, $2, $4, $3 }'
```

Which results in the following output with the fields:

1. Taxon index for first species in pair,

2. taxon/species ID for the first species in the pair, 

3. taxon index for second species in pair, 

4. taxon/species ID for the second species in the pair, and 

5. the evolutionary/phylogenetic distance between these two taxa. 
 
It is good to note that we will see two lines for some of the pairs of species since distances were estimated reciprocally, which is why the first two lines are essentially duplicates of one another.

```
1	Nycticebus_coucang	2	Otolemur_garnettii	0.092236200000
2	Otolemur_garnettii	1	Nycticebus_coucang	0.092236200000
3	Daubentonia_madagascariensis	9	Lemur_catta	0.118283800000
4	Propithecus_coquereli	5	Indri_indri	0.047710400000
5	Indri_indri	4	Propithecus_coquereli	0.047710400000
6	Cheirogaleus_medius	8	Mirza_coquereli	0.071619400000
7	Microcebus_murinus	8	Mirza_coquereli	0.044459300000
8	Mirza_coquereli	7	Microcebus_murinus	0.044459300000
9	Lemur_catta	10	Eulemur_fulvus	0.049048900000
10	Eulemur_fulvus	11	Eulemur_flavifrons	0.014128410000
```

Finally, we get to another important part of this task, extracting the pair of species with the smallest evolutionary/phylogenetic distance. To do so, we will sort based on the 5th column of our new dataset and then extract the first row, which will be the pair of taxa with the smallest evolutionary/phylogenetic distance.

```bash
paste <(cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | tail -n +2 | awk '{ print $1 }') \
<(cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | tail -n +2 | while read line;
do echo -e "${line}" | datamash transpose | tail -n +2 |
awk -v OFS="\t" '{ print NR, $0 }' | sort -k2,2n | sed -n '2,2p';
done) |
while read distance;
do taxindex=$(echo -e "${distance}" | awk '{ print $2 }');
cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | cut -f 1 | tail -n +2 | sed -n "${taxindex},${taxindex}p" |
awk -v OFS="\t" -v distance="${distance}" '{ print distance, $0 }';
done | 
awk -v OFS="\t" '{ print NR, $1, $2, $4, $3 }' | 
sort -k5,5n | head -1
```

Which gives us the following final output of the taxon pair with the smallest evolutionary/phylogenetic distance.

```
58	Cavia_tschudii	59	Cavia_porcellus	0.001111819000
```

It is a nice sanity check to see that the species with the smallest evolutionary/phylogenetic distance are in the same genus. What I have chosen to do is to write the above output to a file.

```bash
paste <(cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | tail -n +2 | awk '{ print $1 }') \
<(cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | tail -n +2 | while read line;
do echo -e "${line}" | datamash transpose | tail -n +2 |
awk -v OFS="\t" '{ print NR, $0 }' | sort -k2,2n | sed -n '2,2p';
done) |
while read distance;
do taxindex=$(echo -e "${distance}" | awk '{ print $2 }');
cat 605-vertebrate-2020_species_phylo.tree | gotree matrix | cut -f 1 | tail -n +2 | sed -n "${taxindex},${taxindex}p" |
awk -v OFS="\t" -v distance="${distance}" '{ print distance, $0 }';
done | 
awk -v OFS="\t" '{ print NR, $1, $2, $4, $3 }' | 
sort -k5,5n | head -1 > smallest_pairwise_distance.txt
```

And now we can use this output to trim one of this tips away from the tree and reduce redundancy. First, we should randomly select one of the two taxa to purge from the phylogeny. With the following command, I print the two taxa on separate lines, randomly sort, and then grab the top line to randomly select a species/tip ID, which I set equal to the environmental variable `taxon2purge`.

```bash
taxon2purge=$(cat smallest_pairwise_distance | awk -v OFS="\n" '{ print $2, $4 }' | sort -R | head -1)
```

Finally, we can purge this taxon from the phylogeny using another handy feature of `gotree`.

```bash
cat 605-vertebrate-2020_species_phylo.tree | gotree prune ${taxon2purge} > 605-vertebrate-2020_species_phylo.phyloreduce.001.tree
```

## Protecting a Given Taxon from being Purged from the Phylogeny

The above code should work to purge one of the two taxa with the smallest evolutionary/phylogenetic distance. However, for a lot of comparative genomics research, results are often contextualized based on a given reference species or genome. For example, someone studying humans may wish to use my approach to purge redundancy from a large phylogeny that includes humans, but they probably *do not* want to purge human from the phylogeny since a lot of their downstream work will rely on that reference organism. Therefore, it is useful to have the ability to "protect" a given species/tip from being purged so that whenever that species has the smallest evolutionary/phylogenetic distance to another species, the other species ends up being purged from the phylogeny. To accomplish this, I modified the commands that create the `taxon2purge` variable as follows.

```bash
# set a "protected" taxon - Mus musculus in this case
protect="Mus_musculus"

# if one of the two taxa are protected, print the other taxon so it gets purged
# otherwise, print both and have the random sort extract one to purge (as above)
taxon2purge=$(cat ${3} | awk -v OFS="\n" -v protect="${protect}" '{ if ($2 == protect) print $4; else if ($4 == protect) print $2; else print $2, $4 }' | sort -R | head -1)
```

## Pulling Everything Together into a Usable Script

All the magic is now described in the above sections, but to run this code iteratively, it makes sense to build a more flexible script that can take arbitrary inputs and write appropriate output files. I am a fan of using [positional arguments](https://www.computerhope.com/jargon/p/positional-parameter.htm) in Unix to provide appropriate file names, which can be coded into the script as `${1}`, `${2}`, etc. I also modified the above code a little so that the distance matrix inferred by gotree is only inferred once. It can then be stored as a variable and provided to different portions of the code using a simple `echo -e` command. Otherwise, you will see the following code emmulates what we walked through in detail above. We can call this script `phyloreduce`.

```bash
# parameters
# ${1} = input phylogeny to reduce
# ${2} = name of taxon to "protect" (avoid purging from phylogeny)
#        useful if a certain reference taxon/genome is needed
# ${3} = output text file for information on shortest neighbor distance
#        the taxon in either column 2 or 4 will be randomly selected and pruned
# ${4} = output reduced phylogeny with selected taxon pruned

# gather pairwise taxon distance matrix from gotree and ${1}
# use index of minimum value per row to extract appropriate taxon
# format into taxon1_index, taxon1, taxon2_index, taxon2, distance
# sort numerically by distance and keep first record (lowest distance)
# write this first record to ${2} file
matrix=`cat ${1} | gotree matrix`;
paste <(echo -e "${matrix}" | tail -n +2 | awk '{ print $1 }') \
<(echo -e "${matrix}" | tail -n +2 |
while read line;
do echo -e "${line}" | datamash transpose | tail -n +2 |
awk -v OFS="\t" '{ print NR, $0 }' | sort -k2,2n | sed -n '2,2p';
done) |
while read distance;
do taxindex=$(echo -e "${distance}" | awk '{ print $2 }');
echo -e "${matrix}" | cut -f 1 | tail -n +2 | sed -n "${taxindex},${taxindex}p" |
awk -v OFS="\t" -v distance="${distance}" '{ print distance, $0 }';
done |
awk -v OFS="\t" '{ print NR, $1, $2, $4, $3 }' |
sort -k5,5n | head -1 > ${3}

# gather "protected" taxon
protect=$(echo ${2})

# based on output, gather the two taxa names for closest species and randomly select 1 to prune
taxon2purge=$(cat ${3} | awk -v OFS="\n" -v protect="${protect}" '{ if ($2 == protect) print $4; else if ($4 == protect) print $2; else print $2, $4 }' | sort -R | head -1)

# prune the randomly selected taxon from original tree and write to ${3}
cat ${1} | gotree prune ${taxon2purge} > ${4}
```

One can run this script as follows.

```bash
bash phyloreduce 605-vertebrate-2020_species_phylo.tree Mus_musculus 605-vertebrate-2020_species_phylo.phyloreduce.001.txt 605-vertebrate-2020_species_phylo.phyloreduce.001.tree
```

## Running `phyloreduce` Iteratively

Now that we have a flexible script, it is much easier to run `phyloreduce` iteratively, which is one of the major goals of this approach. It will seldom make sense to purge just a single tip from the tree. Rather, one probably wants to trim dozens or more species/tips from the phylogeny. With a 605-way alignment, we have a lot of room for purging redundancy. I was interested in getting down to about 25-50 species, so I decided to 575 interations with my script, each one purging a given tip of the tree based on whether it is one of the pairwise taxa with the smallest evolutionary/phylogenetic distance, with care to protect a given taxon from being purged.

Let's first standardize names a little bit so we can more readily keep track of them. We can create a symlink between our original phylogeny and a file with this standardized name.

```bash
ln -s 04_605-vertebrate-2020_species_phylo.tree 04_605-vertebrate-2020_species_phylo.phyloreduce.000.tree
```

Now we will have this common suffix `.phyloreduce.000` that keeps track of which iteration we are on. Given this is the original, full, input phylogeny, it made sense to use `000` here.

Now we can use a simple script to run `phyloreduce` iteratively.

```bash
# for iterations (i) 1 to 575
# gather input iteration number as ${iteration} - 1
# gather output iteration number as ${iteration}
# produce the command calling the script with the appropriate inputs and outputs
# protect "Mus_musculus" in each iteration
# echo and eval the command so it runs
for i in {1..575};
do input=$(printf "%03d\n" `echo "$((i-1))"`);
output=$(printf "%03d\n" ${i});
cmd="bash phyloreduce.jobscript 04_605-vertebrate-2020_species_phylo.phyloreduce.${input}.tree Mus_musculus 04_605-vertebrate-2020_species_phylo.phyloreduce.${output}.txt 04_605-vertebrate-2020_species_phylo.phyloreduce.${output}.tree";
echo $cmd;
eval $cmd;
done
```

This script took several hours to run for me. There are probably much more efficient ways of doing this and someone with some real skills could speed this up significantly, but it is good enough for me for the rare occassions where I might need to perform this task. After finishing, we now see 575 sets of output files from each iteration. The `.txt` file keeps track of the pairwise species/tips with the smallest evolutionary/phylogenetic distance in a given iteration and the `.tree` file is the resulting phylogeny from this iteration with one of the species/tips purged.

And that's it! Each iteration results in a phylogeny with one less tip than the iteration before where the tip that is purged is one of those that has the smallest pairwise evolutionary/phylogenetic distance. So by the 575th iteration, we are down to 30, which we can confirm.

```bash
cat 605-vertebrate-2020_species_phylo.phyloreduce.575.tree | gotree labels | wc -l
# 30
```

We can also confirm that `Mus_musculus` is retained in all of our 576 phylogenies (1 input + 575 iterations), as it should be based on the code I wrote.

```bash
grep "Mus_musculus" *phyloreduce.???.tree | wc -l
# 576
```

I am going to try to work on a visualization so one can see how redundancy is reduced over iterations, but I do not have this ready yet. In the original Treemer README there is a link to [another GitHub repo](https://github.com/thackl/treemmer-animate) where someone was able to create an animation using [`ggtree`](http://yulab-smu.top/treedata-book/) and [`gganimate`](https://gganimate.com/), but this does not appear to work with the new API for `gganimate`, so I did not have luck with this. Stop back and perhaps I will have an animation I can share in the near future!
