---
layout: posts
title: "Creating a scrolling DNA sequence visualization"
date: 2018-01-29
excerpt: "Colored DNA sequence scrolling animation."
---

Are you a scientist working in genomics? Have you ever given an interview to your university communications staff or local press? Did you ever find yourself wishing you could have a continuously scrolling nucleotide sequence running in the background during one of these interviews? Well look no further, because here is exactly what you need to create such an effect, which will really wow those who watch your interview.

We will create two scripts - one that creates a random string of nucleotides and one that colors each nucleotide a different color. Here is the first script, in `python`, named `random_seq.py`.

```python
#!/usr/bin/env python

import random
import sys

def DNA(length):
    return ''.join(random.choice('CGTA') for _ in xrange(length))

print DNA(int(sys.argv[1]))
```

This script just prints a long letter string of As, Cs, Gs, and Ts as STDOUT. There is a single argument that needs to be passed, which is the length of this string. To make the DNA scroll for a while you will need a pretty big number, but note that the larger the number the longer the script will take to get going and print anything.

Here is teh second script, in `perl`, named `color_characters.pl`, which I found somewhere online (can't remember where now).

```perl
#!/usr/bin/env perl

use Getopt::Std;
use strict;
use Term::ANSIColor;

my %opts;
getopts('hic:l:',\%opts);
    if ($opts{h}){
      print<<EoF;
Use -l to specify the pattern(s) to highlight. To specify more than one
pattern use commas.

-l : A Perl regular expression to be colored. Multiple expressions can be
     passed as comma separated values: -l foo,bar,baz
-i : makes the search case sensitive
-c : comma separated list of colors;

EoF
      exit(0);
    }

my $case_sensitive=$opts{i}||undef;
my @color=('bold red','bold blue', 'bold yellow', 'bold green',
           'bold magenta', 'bold cyan', 'yellow on_magenta',
           'bright_white on_red', 'bright_yellow on_red', 'white on_black');
if ($opts{c}) {
   @color=split(/,/,$opts{c});
}
my @patterns;
if($opts{l}){
     @patterns=split(/,/,$opts{l});
}
else{
    $patterns[0]='\*';
}


# Setting $| to non-zero forces a flush right away and after
# every write or print on the currently selected output channel.
$|=1;

while (my $line=<>)
{
    for (my $c=0; $c<=$#patterns; $c++){
    if($case_sensitive){
        if($line=~/$patterns[$c]/){
           $line=~s/($patterns[$c])/color("$color[$c]").$1.color("reset")/ge;
        }
    }
    else{
        if($line=~/$patterns[$c]/i){
          $line=~s/($patterns[$c])/color("$color[$c]").$1.color("reset")/ige;
        }
      }
    }
    print STDOUT $line;
}
```

This one has some options, but the `-l` one is what you are most likely to use (but you can control colors to some degree using `-c`). Now let's put these together to get the desired effect.

```bash
python random_seq.py 100000 | perl color_characters.pl -l A,C,G,T
```

That's it! It is pretty easy to run. Just vary your random DNA sequence length (it is set to 100000 here) to control how long you want the sequence to scroll. But note it will take longer to "warm up" and start generating the actual text.
