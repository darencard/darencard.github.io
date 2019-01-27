---
layout: posts
title: "gnuplot quickstart guide"
date: 2017-08-31
excerpt: "A quick-start guide for using `gnuplot` for in-terminal plotting."
---

Sometimes it is really nice to just take a quick look at some data. However, when working on remote computers, it is a bit of a burden to move data files to a local computer to create a plot in something like `R`. One solution is to use `gnuplot` and make a quick plot that is rendered in the terminal. It isn't very pretty by default, but it gets the job done quickly and easily. There are also advanced `gnuplot` capabilities that aren't covered here at all.

`gnuplot` has it's own internal syntax that can be fed in as a script, which I won't get into. Here is the very simplified `gnuplot` code we'll be using:

```gnuplot
set terminal dumb size 120, 30; set autoscale; plot '-' using 1:3 with lines notitle
```
Let's break this down:
- `set terminal dumb size 120, 30`: `gnuplot` has 'terminals', which is essentially the output format for the plot. Here we are using `dumb` which renders the plot in ASCII characters in the terminal. You can also specify size parameters for the plot, in this case we're using `size 120, 30` to make it a bit larger than default (you can play around with this).
- `set autoscale`: This just makes it so that the axes are automatically scaled, which is normally desireable.
- `plot '-' using 1:3 with lines notitle`: Performs the plotting magic. The `'-'` is the file from which to take the data and plot, which is being read in from STDIN (`'-'`), but you could specify a file name instead. Next, `using 1:3` tells `gnuplot` to use columns 1 and 3 from the data file for plotting (x = 1, y = 3). Change accordingly to plot any column combination you desire. Finally, `with lines notitle` just makes the plot a line plot with no title.

This should allow for many basic plots (but note the lack of axis labels!). This plotting script is called as follows:

```bash
gnuplot -p -e "set terminal dumb size 120, 30; set autoscale; plot '-' using 1:3 with lines notitle"
```

This basically feeds the script from above into the `gnuplot` command-line call. The `-p` flag just allows the plot to persist beyond the command call (otherwise it disappears) and `-e` tells `gnuplot` to expect the following script, which is surrounded by quotes.

This makes it pretty easy to call plotting in a loop, as follows:

```bash
for i in *.summary.txt; do echo $i; \
  cat $i | gnuplot -p -e "set terminal dumb size 120, 30; set autoscale; plot '-' using 1:3 with lines notitle"; done
```

Here, we are taking each file that ends in `.summary.txt`, and piping it's content into `gnuplot` to produce a plot. Pretty cool, right!?
