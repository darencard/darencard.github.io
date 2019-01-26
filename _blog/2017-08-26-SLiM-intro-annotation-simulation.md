---
layout: posts
title: "Introductory walkthrough of SLiM"
date: 2017-08-26
excerpt: "Some notes I made as I began learning forward-time simulation using SLiM."
---

[SLiM](https://messerlab.org/slim/) is a newer, powerful piece of population genetic simulation software that is well documented, user-friendly, flexible, and has a pretty sweet GUI interface (plus command-line capability). The following script represents an initial dummy simulation situation I created as I got my feet wet with `SLiM`, and I added many commented notes to make it clear what each command was doing.

`//` in SLiM context are comments.

```
// set up a simple neutral simulation
initialize() {
	initializeMutationRate(1e-7);

	// m1 mutation type: neutral
	// arguments: 	arbitrary id
	//					dominance coefficient (0-1, no dominance)
	//					mutation distribution (f=fixed)
	//					mutation distribution argument (sel. coefficient)
	initializeMutationType("m1", 0.5, "f", 0.0);

	// g1 genomic element type: uses m1 for all mutations
	// arguments:	arbitrary id
	//					mutation type to use (from mutation type specified above)
	//					proportion mutations drawn from this mutation type
	initializeGenomicElementType("g1", m1, 1.0);

	// uniform chromosome of length 100 kb with uniform recombination
	// arguments:	genomic element type to use
	//					first base position
	//					last base position (=1kb length)
	initializeGenomicElement(g1, 0, 1000);

	// recombination rate across simulated chromosomes
	// arguments:	recombination rate per base
	initializeRecombinationRate(1e-8);
}

// perform action at generation 1
1 {
	// create a population of 500 individuals
	// arguments:	arbitrary id
	//					population size
	sim.addSubpop("p1", 500);
}

// run for burn-in number of generations to have population reach equilibrium
1000000 late() {
	p1.outputSample(100);
}

// output samples of 10 genomes periodically, all fixed mutations at end
// at generation N (= 1001, here)
// late() = at end of generation (after offspring)
1000001 late() {
	// output sample of 100 genomes (haploid) = 2N
	// from population 1 that was initialized above
	// p1.outputSample(100);
	sim.outputFull();
}

/// DESCRIPTION OF FULL OUTPUT (/// = description)
// #OUT: 2 A
/// #output prefix: generation type of output
//
// Populations:
/// Population-level output
// p1 500 H
/// subpopulation size reproduction_type (H=hermaphroditic)
//
// Mutations:
/// Mutation-level output: 1 mutation per line
// 13 4 m1 89111 0 0.5 p1 1 2
/// temp_id permanent_id mutation_type bp_position sel_coeff dominance_coeff
/// pop_id generation_mutation_arose prevalence_num_genomes_possess
//
// Individuals:
/// Individual-level ouptput: 1 individual per line, starting at 0
// p1:i0 H p1:0 p1:1
/// pop_id:ind_id reproduction_type pop_id:genome1_id pop_id:genome2_id
//
// Genomes:
/// Genome-level output: 1 genome per line, starting at 0
// p1:0 A
/// pop_id:genome_num autosome(A)/sex-linked <1+ mutations ordered sequentially)
```
