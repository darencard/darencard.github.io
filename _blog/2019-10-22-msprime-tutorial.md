---
layout: posts
title: "A Brief Introduction to `msprime`"
date: 2019-09-27
excerpt: "Notes from a quick introduction session I attended on `msprime`"
---

Coalescent simulators are commonly used in evolutionary biology for a variety of reasons. They nicely complement forward-time simulators, which evaluate evolution forward in time, keep better track of individuals within populations, and allow users to implement more "biological" scenarios. Coalescent simulators, on the other hand, are more model-based, with simulations coming from draws of distributions of parameters, and they also simulate backward in time.

I sat through a brief introductory session on newer coalescent simulation software called `msprime`. This session was kindly provided by [Michael Miyagi](https://wakeleylab.oeb.harvard.edu/people/michael-miyagiin) the [Wakeley lab](https://wakeleylab.oeb.harvard.edu/) Harvard. `msprime` is published by [Kelleher et al.](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004842) and documentation can be found [here](https://msprime.readthedocs.io/en/stable/). Some major advantages to `msprime` are that it is fast, accounts for numerous population genetic processes, and uses a new "tree sequence" data type to store geneologies, which more efficiently stores trees across genomes by taking advantage of the fact that the relationships between trees are often not random. While this leads to a more specialized data structure, it works well and common data formats (e.g., Newick) can be exported. `msprime` can also do forward-time simulations.

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
