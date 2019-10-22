---
layout: posts
title: "A Brief Introduction to msprime"
date: 2019-09-27
excerpt: "Notes from a quick introduction session I attended on msprime"
---

Coalescent simulators are commonly used in evolutionary biology for a variety of reasons. They nicely complement forward-time simulators, which evaluate evolution forward in time, keep better track of individuals within populations, and allow users to implement more "biological" scenarios. Coalescent simulators, on the other hand, are more model-based, with simulations coming from draws of distributions of parameters, and they also simulate backward in time.

I sat through a brief introductory session on newer coalescent simulation software called `msprime`. This session was kindly provided by [Michael Miyagi](https://wakeleylab.oeb.harvard.edu/people/michael-miyagiin) the [Wakeley lab](https://wakeleylab.oeb.harvard.edu/) Harvard. `msprime` is published by [Kelleher et al. (2016)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004842) and has [detailed documentation](https://msprime.readthedocs.io/en/stable/). Some major advantages to `msprime` are that it is fast, accounts for numerous population genetic processes, and uses a new "tree sequence" data type to store geneologies, which more efficiently stores trees across genomes by taking advantage of the fact that the relationships between trees are often not random. While this leads to a more specialized data structure, it works well and common data formats (e.g., Newick) can be exported. `msprime` can also do forward-time simulations.

The tutorial I went through a couple relatively simple examples. `msprime` has a nice Python interface, which makes it quite easy to put together arbitrary simulation models that can be run effectively. The first scenario is a simple 4 population model with fixed population sizes, mutation rates, and recombination rates. The output that results are geneologies for each locus in Newick format and a binary genotype matrix.

```python
import msprime

tree_sequence=msprime.simulate(sample_size=4,Ne=1000000,mutation_rate=3e-7,recombination_rate=3e-7)
treelist=tree_sequence.trees()
for tree in treelist:
	print(tree.newick())
print(tree_sequence.genotype_matrix())
```

This script can be saved as the file `msprime_simple.py` and run simply using the command `python msprime_simple.py` (assuming `msprime` is properly installed).

We also ran through a second, more complex example with a 4-tip tree with a migration event between tips 1 and 3. Here is an image of the tree model we were considering. The image should be rotated clockwise 90 degrees but GitHub is not cooperating properly to do that, so please pardon that issue.

![msprime model](https://github.com/darencard/darencard.github.io/raw/master/assets/images/blog/msprime_model.jpg)

The script to create and run such a model, which includes more details on parameters, etc. is as follows. I have added some annotation comments to help make sense of different parts and of the output being produced.

```python
import msprime
import numpy as np

mutation_rate=3e-8
recombination_rate=3e-8
sample_per_pop=1
pop_size=1e6
length=1e3
intro_prob=0.5
mig_rate=0

def introgression_simulation(length, mutation_rate, recombination_rate, sample_size, pop_size, intro_prob, mig_rate):
	## initialize populations with equal initial sizes
	N_1=N_2=N_3=N_4=pop_size
	population_configurations = [
		msprime.PopulationConfiguration(sample_size=sample_size, initial_size = N_1),
       		msprime.PopulationConfiguration(sample_size=sample_size, initial_size = N_2),
        	msprime.PopulationConfiguration(sample_size=sample_size, initial_size = N_3),
		msprime.PopulationConfiguration(sample_size=sample_size, initial_size = N_4)]
	## number of demes (populations)
	d=4
	## migration rates
	m=mig_rate/(2*(d-1))
	## migration matrix
	migration_matrix = [
        	[0, m, m, m],
        	[m, 0, m, m],
        	[m, m, 0, m],
			[m, m, m, 0]]

	## time at which there are 3 populations
	T_3 = 6*pop_size
	## time at which there are 4 populations
	T_4 = 4*pop_size
	## time at which there is introgression
	T_intro = 2*pop_size
	## time at which there are 2 populations
	T_2 = 20*pop_size

	demographic_events = [
		## migration between 1 and 3 before populations split (backward in time)
		msprime.MassMigration(
			time = T_intro, source = 0, destination = 2, proportion = intro_prob),
		## migration from 1 to 2 all the time (represents the node between tips 1 and 2)
		msprime.MassMigration(
			time = T_4, source = 0, destination = 1, proportion = 1.0),
		## migration from 2 (inclusive) to 3 all the time (represents the node between tips [1,2] and 3)
		msprime.MassMigration(
			time = T_3, source = 1, destination = 2, proportion = 1.0),
		## migration from 3 (inclusive) to 4 all the time (represents the node between tips [][1,2],3] and 4)
		msprime.MassMigration(
			time = T_2, source = 2, destination = 3, proportion = 1.0)
	]
	dd = msprime.DemographyDebugger(
    		Ne=N_1,
	       	population_configurations=population_configurations,
        	migration_matrix=migration_matrix,
        	demographic_events=demographic_events)

	dd.print_history() #can comment out when we are happy demography is correct

	output = msprime.simulate(
		population_configurations = population_configurations,
		migration_matrix = migration_matrix,
		demographic_events = demographic_events,
		mutation_rate = mutation_rate,
		length = length,
		recombination_rate = recombination_rate)
	return output

def treePrinter(treeseq):
	output=''
	runningTot=0
	for index,tr in enumerate(treeseq.trees()):
		lengthQ=-int(np.round(tr.get_interval()[0]))+int(np.round(tr.get_interval()[1]))
		if lengthQ>0:
			runningTot+=lengthQ
			output=output+'['+str(lengthQ)+']'+str(tr.newick())+'\n'
	print(output)

tree_sequence=introgression_simulation(length, mutation_rate, recombination_rate, sample_per_pop, pop_size, intro_prob, mig_rate)

# for tree in tree_sequence.trees():
# 	print(tree.newick())
# print(tree_sequence.genotype_matrix())

##Looking at mutations
## output = position and tree location of the mutation - [2] = tip 2; [0, 1, 2] = ancestor to 0, 1, & 2
# tree=tree_sequence.first()
# for site in tree.sites():
# 	for mutation in site.mutations:
# 		print("Mutation at position {:.2f} over node {}".format(site.position,[n for n in tree.leaves(mutation.node)]))

## output = length of geneology (in bp) and newick of the geneology
## can be used to simulate sequences in seqgen
print(treePrinter(tree_sequence))
```

his script can be saved as the file `msprime_migration.py` and run simply using the command `python msprime_migration.py` (assuming `msprime` is properly installed).
