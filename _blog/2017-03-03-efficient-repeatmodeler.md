---
layout: posts
title: "Running RepeatModeler more efficiently on sequencing reads"
date: 2017-03-03
excerpt: "Instructions for modifying RepeatModeler to run more efficiently on sequencing reads"
---

RepeatModeler isn't very well suited for sample sequencing data, taking a long time and creating copious amounts of intermediate data files. It obviously wasn't designed for small fragments and reads, which are what we get with sample sequencing data, and here are the main difficulties.

1. The subsampling steps for each round take a long time (hours in later rounds) and are done using a single core, which is wasteful and inefficient. However, the full script depends on this subsetting to run properly, so there isn't really a way around this.
2. Parallelization occurs during the RECON analyses of rounds 2 to N, so overall, it makes little sense to parallelize heavily since a major bottleneck is the subsetting step (see #1).
3. Huge amounts of intermediate files are produced, which grow rapidly with each round. Most of these are the `batch-*` files that are used for parallelization during the RECON rounds. In later rounds (5+, the output size inflates to over 200GB, mostly as a result of these files. Therefore, if a user wants to run several RepeatModeler runs concurrently with a few cores each, they run into a huge issue with storage space to store all the files, which can't be deleted until a round ends.

I've come up with a few ways to try to circumvent these problems while not totally reimplementing RepeatModeler.

1. Sampling time is a bottleneck that we largely cannot avoid, but we can cut out the sampling from Round 1, in my experience, as RepeatScout never appears to find any repeat models. Therefore, RepeatModeler subsets for 3 hours, and RepeatScout runs a few minutes and returns nothing. Might as well cut out this step, but the later rounds rely on it completing in some way. Therefore, the shortcut is to simply minimize the amount of subsampling that RepeatModeler does in this round. The default is 40 Mbp and it is easy to just change this to e.g. 1000 bp by changing the `$rsSampleSize` variable on line 255, as follows:
```
my $rsSampleSize               = 1000; # 40000000;     # The size of the sequence given
                                               #   to RepeatScout. ( round #1 )
```
Note that this change may sometimes cause RepeatScout not to find any kmers in the sample it takes, leading to the following error that ends RepeatModeler:
```
OOPS no good lmers
build_lmer_table failed. Exit code 256
Command exited with non-zero status 1
```
If this happens, try running the same RepeatModeler command again a couple times and the new random sampling may product something with kmers, allowing the run to proceed. If this issue continues, increase the `$rsSampleSize` some arbitrary amount and it should go away with larger sample sizes.
2. Given that so much of RepeatModeler is being run with one core (i.e., the subsampling and many processing/refining steps), it makes little sense, in my opinion, to give RepeatModeler a ton of cores to run on. It speeds up some portions of the analysis but others are not changed at all, and are wasting valuable computational resources (especially on servers where this is tracked). In my limited experience playing around with this, I suggest using 4 cores so help a bit with speedup while not breaking the bank on overhead. The run time will vary in each run and per sample, but based on a trial run I round that 4 cores enables for a 7-8 day runtime (5 rounds).
3. The file size issue is quite a pain, especially if a user wanted to run several 4-core RepeatModeler runs concurrently. He or she would run out of disk space with just a couple full RepeatModeler runs, which prohibits parallelization. Fortunately, the file sizes really ramp up in rounds 5 and later and the first 4 rounds are under ~40 GB total size. Therefore, it is possible, on a 500 GB volume (for example), to run 10 RepeatModeler jobs concurrently if only the first 4 rounds are run. The user still has to hassle with rounds 5+, but at least they are running the first 4 efficiently. To implement this, I modified the RepeatModeler script in two different ways to produce `RepeatModelerStart` and `RepeatModelerFinish` flavors, which run rounds 1-4 and then 5+, respectively. Here are the changes to make for `RepeatModelerStart`:
```
# lines 261-262: set max sample size to 27M, which is the sampling amount for round 4, so the script will stop after this round
my $genomeSampleSizeMax        = 27000000; # 100000000;    # The max sample size for
                                               #   RECON analysis
```
```
# lines 483-484: do not reset the sample size to the maximum genome sample size when it is larger (this allows the script
# to stop at the desired point and not keep going)
  # $sampleSize = $genomeSampleSizeMax
      # if ( $sampleSize > $genomeSampleSizeMax );
```
```
# lines 512-514: constrain the while loop to only run when the sample size is below the maximum sample size by adding
# another term to statement
while ( $totSampledBP < $dbSize
        && ( $dbSize - $totSampledBP ) > $fragmentSize
	&& $sampleSize <= $genomeSampleSizeMax )
```
```
Lines 1292-1298: comment to prevent repeat classifier from running, as there is no point in running it now
# if ( $numModels > 0 )
# {
#   print "Classifying Repeats...\n";
#   system(   "$RepModelConfig::REPEATMODELER_DIR/RepeatClassifier "
#           . "-consensi $tmpDir/consensi.fa" );
#   print "Classification Time: " . elapsedTime( 1 ) . "\n";
# }
```
```
Line 1328: change command to find working directory as original doesn't work when restarting (for some reason?)
  my $origDir = abs_path(".");
```
With these changes, you can now run RepeatModeler like you normally would and it should stop after round 4.
```
RepeatModelerStart -pa 4 -engine ncbi -database <DB>
```
This gives you time to clean up files (see below) and run jobs concurrently without eating up the entire disk. When appropriate, `RepeatModelerFinish` can be run to restart the job and finish the later rounds. Using a fresh `RepeatModeler` script, the following change must be made:
```
# line 120: add the following line to avoid working directory error upon RepeatModeler
# note that this change has to be made anytime you want to restart a run, regardless of what is outlined in this Gist
use Cwd 'abs_path';
```
Now the RepeatModeler run can be restarted by first creating an empty `round-5` directory to trigger the program to continue and running RepeatModeler with the `-recoverDir` set to the existing directory of round 1-4 output. ((Note that the consensi.fa files for rounds 2-4 must also contain content. This should not be a problem with normal RepeatModeler runs, but if they happen to be empty you can manually fill them with any artificial text and this will trigger RepeatModeler to infer that a round is complete.))
```
# create empty round-5 directory (change according to output directory name, which will vary)
mkdir <output_dir>/round-5
RepeatModelerFinish -pa 4 -engine ncbi -recoverDir <output_dir> -database <DB>
```
This will run round-5 and potentially further rounds (usually only 6 rounds total, and honestly not sure why the 6th round is sometimes run???). If one wanted to ensure only 5 rounds were run, one could also make the changes outlined above for `RepeatModelerStart`, but set the `$genomeSampleSizeMax` variable to 81000000 (the number of Mbp sampled in round 5). This step will be less efficient, as less runs can be parallelized, but overall this scheme if far more efficient than just running RepeatModeler in series.

## RepeatModeler Intermediate File Cleanup

As outlined above, the huge amount of intermediate files created when RepeatModeler runs can be an issue, especially in the later couple of rounds. Turns out, most of the storage space appears to be used up by the `batch-*gilist*` files, which are used to parallelize the RECON searching. These really aren't interpetable and aren't needed, so I suggest deleting them from each of rounds 2-N. This alone will save massive amounts of space in each round and facilitate better longterm storage of RepeatModeler outputs that are much smaller. A user could also elect to delete more intermediate files as only two files are necessary to restart RepeatModeler, though I think other intermediate files are useful to have. Those two files must be in `round-(N-1)` directory (and contain content), where N-1 is the round prior to the round the user wants to run. For example, to restart and run round-5, the round-4 directory should include these two files. The two files are: `consensi.fa` and `sampleDB-N.fa` where N is the round number.

If wanting to delete all the `batch-*gilist*` files as they are finished being used, you can make an additional change to the RepeatModeler program. The following is based on a fresh RepeatModeler script (so adjust line numbers as necessary).
```
# line 1165 (after 'close OUT;'): add the following line to delete all batch-*gilist* files
system( "rm $roundTmpDir/batch-$batchNum-gilist*" );
```
This will save lots of space, but is throwing away a possibly valuable intermediate file that could be useful to keep. An alternative thing to do is to compress the data, which can help a lot. In this case, we can still delete `batch-*-gilist` files, as they are redundant with `batch-*-gilist.txt` files and could be regenerated. Then we'll zip these `batch-*-gilist.txt` files with our desired compression software (xz > bzip2 > gzip as far as storage savings).
```
# line 1165 (after 'close OUT;'): add the following lines to delete all batch-*-gilist files
# and compress batch-*-gilist.txt files.
# I'm using xz compression here, which saves the most space (output is about 1/9 size of input)
system( "rm $roundTmpDir/batch-$batchNum-gilist" );
system( "xz $roundTmpDir/batch-$batchNum-gilist.txt" );
```
