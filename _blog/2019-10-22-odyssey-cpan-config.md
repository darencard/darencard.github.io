---
layout: posts
title: "Configuring CPAN on Harvard Odyssey cluster"
date: 2019-10-26
excerpt: "Notes for how I configured CPAN on Harvard Odyssey cluster"
---

I have always found [Perl](https://www.perl.org/) to be a pain to manage on computers, but gradually I have improved at figuring out how to wrangle it to get what I want.

Today I was getting the error `Can't locate Bio/SeqIO.pm in @INC` when I was trying to run a script I downloaded somewhere. I tried installing [BioPerl](https://bioperl.org/) previously and it seemed to install properly. I did this using the [CPAN](https://www.cpan.org/) package manager built into Perl, but apparently my installation was not being found.

I was able to overcome this using a combination of the guidance from Harvard Research Computing - found [here](https://www.rc.fas.harvard.edu/resources/documentation/software-on-the-cluster/perl/) - and some nice advice on how to properly configure the CPAN installation directories found [here](https://linuxgazette.net/139/okopnik.html).

I basically following the directions in the 2nd source, but used the install directory provided by Harvard RC.

So my `.bashrc` file has the following lines:

```
if [ -z "$PERL5LIB" ]
then
	# If PERL5LIB wasn't previously defined, set it...
	PERL5LIB=~/apps/perl/lib
else
	# ...otherwise, extend it.
	PERL5LIB=$PERL5LIB:~/apps/perl/lib
fi

MANPATH=$MANPATH:~/apps/perl/man

export PERL5LIB MANPATH
```

And when I configured CPAN again from scratch, I used the following commands at each `cpan> ` prompt.

```
o conf makepl_arg "LIB=~/apps/perl/lib INSTALLSITEMAN1DIR=~/apps/perl/man/man1 INSTALLSITEMAN3DIR=~/apps/perl/man/man3"
o conf make_install_arg UNINST=0
o conf commit
```

Then I was able to properly install BioPerl's SeqIO module (`cpan Bio::SeqIO`), though with some difficulty. A dependency or two was required that was easy to install using `cpan` (as above). However, one dependency kept giving me trouble. I kept getting an error when trying to install `XML::DOM::XPath` that led the entire installation of the dependency, and thus SeqIO, to fail.

```
t/test_non_ascii.t .................... The encoding pragma is no longer supported at t/test_non_ascii.t line 10.
BEGIN failed--compilation aborted at t/test_non_ascii.t line 10.
# Looks like your test exited with 2 before it could output anything.
t/test_non_ascii.t .................... Dubious, test returned 2 (wstat 512, 0x200)
Failed 10/10 subtests
```

In looking around, I discovered that it is because the line 10 of the installation file `t/test_non_ascii.t` should be `use utf8;` instead of `use encoding 'utf8';` due to updates in Perl (hard to believe this error is persisting). Following the directions [here](https://stackoverflow.com/questions/47966512/error-installing-xmldomxpath), I found the path for my `t/test_non_ascii.t` file to be `~/.cpan/build/XML-DOM-XPath-0.14-0/t/test_non_ascii.t`. I modified line 10 as above and installed this module manually using the command `cpanm -force ~/.cpan/build/XML-DOM-XPath-0.14-0`. After doing this, I was then able to install `Bio::SeqIO` without issue! And the script I was trying to run is now working successfully! What a pain in the ass!
