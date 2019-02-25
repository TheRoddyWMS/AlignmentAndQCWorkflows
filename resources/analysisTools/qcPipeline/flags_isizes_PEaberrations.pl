#!/usr/bin/perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

# extended flagstats with number of secondary and supplementary alignments,
# paired end reads mapping to different chromosomes,
# bucketsort insert sizes, calculate average, median, standard deviation and "SDpercent" as standard deviation/median*100
# print out the insert size frequencies on stdout for R:
# insertsize_bin\tnumber_reads
# optional, print out median and "SDpercent" to a file (mimicking the text output of insertSizeDistribution.r)

use strict;
use warnings;
use Getopt::Std;

my %opts = (p=>1000,n=>1);
getopts("a:b:c:f:i:m:o:p:q:h", \%opts);

my $samfile = $opts{i};
my $maxisize = $opts{p};
my $isizebinout = $opts{b};
my $isizeoutfile = $opts{m};
my $matrixout = $opts{a};
my $peaberrstat = $opts{o};
my $chromfile = $opts{c};
my $flagsout = $opts{f};
my $minmapq = $opts{n};

if (defined $opts{h} || ! defined $samfile || ! defined $isizebinout  || ! defined $isizeoutfile || ! defined $matrixout || ! defined $chromfile || ! defined $flagsout)
{
	die "USAGE: $0 [options]
	-a FILE	file to write the matrix with different chromosomes out to (_DiffChroms.txt)
	-b FILE	file to write the insert size bins out to (_insertsizes.txt)
	-c FILE file with the chromosomes that are wanted (must fit w.r.t. pre- and suffix!)
	-f FILE	file to write the extended flagstats out to (_flagstats_extended.txt)
	-i FILE	SAM file (set to - for pipe from STDIN)
	-m FILE	file to write the insert size median and standard deviation/median*100 out to (_plot.png_qcValues.txt)
	-n INT minimal mapping quality for improper pair (default 1)
	-o FILE file to write the percentage of improper pairs out to (_DiffChroms.png_qcValues.txt)
	-p INT maximal insert size for a proper pair (default 1000)
	-h help\n";
}

open (IN, $samfile) or die "Cannot open $samfile: $!\n";
open (CF, $chromfile) or die "Cannot open $chromfile: $!\n";
open (IB, ">$isizebinout") or die "could not open $isizebinout for writing: $!\n";
open (IS, ">$isizeoutfile")  or die "could not open $isizeoutfile for writing: $!\n";
open (MA, ">$matrixout") or die "could not open $matrixout for writing: $!\n";
open (FL, ">$flagsout") or die "could not open $flagsout for writing: $!\n";
open (PE, ">$peaberrstat") or die "could not open $peaberrstat for writing: $!\n";

if ($maxisize !~ /^\d+$/)
{
	print STDERR "maximal insert size $maxisize is not a number, setting to default 1000\n";
	$maxisize = 1000;
}
if ($maxisize < 1000)
{
	print STDERR "maximal insert size $maxisize is unreasonably small, setting to default 1000\n";
	$maxisize = 1000;
}
if ($minmapq !~ /^\d+$/)
{
	print STDERR "minimal mapping quality $minmapq is not a number, setting to default 1\n";
	$minmapq = 1;
}

# first read in which chromosomes should be regarded for the statistics
my @help = ();
my %chroms = ();	# to look up later if the chromosome is wanted
my @chromarray = ();	# for the output in matrix form
my $chrnum = 0;
while (<CF>)
{
	if ($_ =~ /^#/)
	{
		next;
	}
	chomp;
	@help = split ("\t", $_);
	$chroms{$help[0]} = 1;
	$chromarray[$chrnum] = $help[0];
	$chrnum++;
}
close CF;
print STDERR "$chrnum chromosomes to keep: @chromarray\n";

# extended flagstats
# all reads
my $all = 0;
# only the non-dup, non-secondary, non-supplementary reads
my $uniq = 0;
# of these, with min. mapping quality
my $minmapuniq = 0;
# of these, those on wanted chroms
my $onchr = 0;
# of these, both on wanted chroms - these will be the denominator
my $both = 0;
# the PE aberrations among them
my $aberrant = 0;
my %chrompairs = ();	# hash of hashes: chromosomes of read => chr. of mate => number
my $flag = 0;

# insert size stats
my @isizes = ();
my $entry = 0;	# absolute value of insert size
my $min = 10000000000000000000;
my $max = 0;
my $ctr = 0;	# number of read 1 of proper pairs for insert sizes
my $sum = 0;	# sum of insert sizes
my $pp_strange = 0;
my $average = 0;

while (<IN>)
{
	if ($_ =~ /^\@/)	# there might be a SAM header
	{next;}
	$all++;
	@help = split ("\t", $_);
	$flag = $help[1];
	# not unmapped, no duplicate and no secondary/supplementary alignment, and mapqual >= X
	if (!($flag & 4) && !($flag & 1024) && !($flag & 256) && !($flag & 2048))
	{
		$uniq++;
		if ($help[4] >= $minmapq)
		{
			# of these, with mapqual >= X
			$minmapuniq++;
		}
		else
		{
			next;
		}
		# is the read itself on a wanted chromosome?
		if (defined $chroms{$help[2]})
		{
			$onchr++;
			# and the mate on a wanted chromosome
			if (defined $chroms{$help[6]} || $help[6] eq "=")	# same chrom is usually indicated by "=" instead of repeating the name
			{
				$both++;
				# paired end aberration:  mate also has to be mapped, on a different chrom
				if (!($flag & 8) && $help[6] ne "=" && ($help[2] ne $help[6]))
				# keep matrix symmetrical to see whether there is a bias, e.g. more 1->10 than 10->1
				{
					$aberrant++;
					# only use read1 info, since read2 might have mapq 0, and the info of having bias w.r.t. which read
					# is more interesting
					if ($flag & 64)
					{
						$chrompairs{$help[2]}{$help[6]}++;
					}
				}
				# for insert sizes, take first read of a proper pair (-f 67 = 64 (first in pair) + 2 (proper pair) + 1 (paired));
				# discarding duplicates (-F 1024) is already done further up
				if ($flag & 64 && $flag & 2 && $flag & 1)
				{
					$entry = abs($help[8]);	# insert size
					if ($entry < $min)
					{
						$min = $entry;
					}
					if ($entry > $maxisize)	# bwa can give "proper pair" to ones with insert size over 200 million bp, makes no sense at all!
					{
						$pp_strange++;
						next;
					}
					if ($entry > $max)
					{
						$max = $entry;
					}
					$isizes[$entry]++;
					$ctr++;
					$sum+=$entry;
				}
			}
		}
	}
}
close IN;
# extended flagstats
print FL "total alignments\t$all\nnon-duplicate, non-secondary, non-supplementary reads\t$uniq\nsuch with mapping quality >=$minmapq\t$minmapuniq\nsuch on regarded chromosomes\t$onchr\nsuch with both reads on regarded chromosomes\t$both\nsuch mapping on different chromosomes\t$aberrant\nproper pairs read 1\t$ctr\n";
close FL;


# in case we have single end reads, there will be no counted ones for isizes ($ctr < 1)
# and PE aberrations ($aberrant < 1) => fill the files with placeholder "NA"

if ($aberrant < 1)
{
	print STDERR "no aberrant paired reads found, single end reads?\n";
	print PE "NA\n";
}
else
{
	$average = $sum/$ctr;
	print STDERR "insert sizes: $pp_strange proper pairs with insert size > $maxisize; $ctr values in range. minimum: $min; maximum: $max; average: $average\n";
	print STDERR "from $all reads, $aberrant are paired end aberrations\n";
	my $percentage = sprintf ("%2.2f\n", ($aberrant/$both*100));
	print PE $percentage;
}
close PE;

# print out matrix
if ($aberrant < 1)
{
	print MA "NA\n";
}
else
{
	my $chromstring = join ("\t", @chromarray);
	# $chrnum is the number of chromosomes wanted = length of the @chromarray
	print MA "\t$chromstring\n";
	my $chrom = "";
	my $chrmate = "";
	for (my $c1 = 0; $c1 < $chrnum; $c1++)
	{
		$chrom = $chromarray[$c1];	# the actual chromosome name of the read
		print MA "$chrom";
		for (my $c2 = 0; $c2 < $chrnum; $c2++)
		{
			$chrmate = $chromarray[$c2];	# ... of the mate
			if (exists $chrompairs{$chrom} && exists $chrompairs{$chrom}{$chrmate})
			{
				print MA "\t", $chrompairs{$chrom}{$chrmate};
			}
			else
			{
				print MA "\t0";
			}
		}
		print MA "\n";
	}
}
close MA;

# insert size statistics
# The median is in a sorted list the value in the middle. The list is actually sorted but contains the same values several times.
# I know how many entries are there. So with some smart calculation, I can "jump" to it. In case of a even number, the median is
# the average of the 2 middle values. But for isizes that will be the same value as middle-1 anyways -
# here for simplicity, just take the value when read number is > 1/2 of all
# Also get the standard deviation by the sqrt of the deviations from the average by going through the list efficiently

if ($ctr < 1)
{
	print STDERR "no properly paired reads found, single end reads?\n";
	print IB "NA\n";
	print IS "NA\nNA\nNA\n";
}
else
{
	my $isize_mindiff = $ctr/50000;	#10000;	# for 30x genome, but less reads means less difference!
	print STDERR "for insert size slope detection, set min difference to $isize_mindiff\n";
	my $median = 0;
	my $middle = int($ctr/2);
	my $howmany = 0;	# number of entries for each isize bin
	my $ctr2 = 0;
	my $docollect = 1;	# a flag for the median determination
	my $devsqsum = 0;		# the sum of the squared deviations from the average
	# the modes = peaks; should be only 1 (ideally, the median), but might be 2 in case of adapter dimers, then it's bad!
	# a real peak must have several 1000 reads to really stand out
	my $mode1 = 0;
	my $mode = 0;
	my $modnr = 0;
	my $oldnr = 0;	# number of reads in previous bin
	my $newmode = 1;	# a flag for setting the first mode

	# print for plotting with R:
	# insertsize_bin\tnumber_reads
	for (my $i = $min; $i <= $max; $i++)
	{
		$howmany = $isizes[$i] || 0;
		print IB "$i\t$howmany\n";
		# find highest peak
		if ($howmany > 1000 && $howmany > $modnr)
		{
			$modnr = $howmany;
			$mode1 = $i;
		}
		# drop after a peak => keep this peak as the first mode
		# this is not optimal because after a small dip, there could be the real first peak
		# if ($howmany > 1000 && $howmany < $oldnr && $newmode > 0)
		# try to avoid detecting noise: the slope after the peak must be really steep
		if ($howmany > 1000 && $howmany < $oldnr && $newmode > 0 && (abs($howmany - $isizes[$i+1]) > $isize_mindiff))
		{
			$mode = $mode1;
			$mode1 = 0;
			$modnr = 0;
			$newmode = 0;
		}
		$ctr2+=$howmany;
		$oldnr = $howmany;
		if ($ctr2 >= $middle && $docollect > 0)	# we have hit the middle, fix the median!
		{
			$median = $i;
			print STDERR "hit middle ($middle) at $ctr2 entries: median is $i\n";
			$docollect = 0;
		}
		# squared deviations from the average
		$devsqsum += $howmany * (($i - $average)**2);
	}

	my $stddev;
	my $sdrd;
	my $sdperc;
	my $rounded;
	if ($ctr - 1 == 0) { # Cannot calculate SD and derived values due to lack of data.
		$stddev = "NA";
		$sdrd = "NA";
		$sdperc = "NA";
		$rounded = "NA";
	} else {
		$stddev = sqrt($devsqsum/($ctr-1));
		$sdrd = int($stddev+0.5);
		$sdperc = $stddev/$median*100;
		$rounded = int($sdperc+0.5);
	}
	print STDERR "stddev: $stddev; rounded: $sdrd\n";
	print STDERR "median: $median; SDpercent: $sdperc; rounded: $rounded\n";
	print IS "$median\n$rounded\n$sdrd\n";

	# mode(s) often there is an intitial "peak" at 0 that does not count
	if ($mode1 > 0 && $mode > 0)
	{
		print STDERR "2 peaks: $mode and $mode1";
		# if the difference between the 2 peaks is small, it was just a random fluctuation,
		# i.e. a small dip before the real peak
		if (($mode1 - $mode) > 20)
		{
			print STDERR " - caution, bimodal distribution, might be adapter dimers!\n";
			print IS "INSERT SIZE DISTRIBUTION POSSIBLY BIMODAL: (first) 2 peaks at $mode and $mode1\n";
		}
		else
		{
			print STDERR " - small difference does not hint at critical binomial distribution.\n";
		}
	}
	else
	{
		print STDERR "1 peak: $mode\n";
	}
}
close IS;
close IB;
exit;
