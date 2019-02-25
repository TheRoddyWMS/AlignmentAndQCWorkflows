#!/usr/bin/perl
#
# Copyright (c) 2018 German Cancer Research Center (DKFZ).
#
# Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
#

use strict;
use warnings;

if (@ARGV < 2)
{
        die "$0 needs two arguments: 1.sam file (if piped in from samtools view, use -), 2. library name\n";
}

my $samfile = shift;
my $libname = shift;    # sample_pid

open (IN, $samfile) or die "Cannot open $samfile: $!\n";

# duplicates are UNPAIRED_READ_DUPLICATES + 2*READ_PAIR_DUPLICATES
# the denominator is mapped reads: UNPAIRED_READS_EXAMINED + 2*READ_PAIRS_EXAMINED

# obviously secondary/supplementary alignments are not counted
# because the read sum is lower than the number of alignments in flagstats
# does unpaired really mean unpaired? According to flagstats, all are "paired in sequencing"
# probably that's the singletons with mate unmapped (number of alignments comes close to this)

# optical duplicates cannot be determined and biobambam usually has 0 => 0
# library size estimation is anyways off => NA

my $all = 0;
my $unpaired = 0;       # singletons
my $pairs = 0;
my $unmapped = 0;
my $singledups = 0;
my $pairdups = 0;

while (<IN>)
{
        if ($_ =~ /^\@/)        # there might be a SAM header
        {next;}
        $all++;
        my @help = split ("\t", $_);
        my $flag = $help[1];
        if (($flag & 256 || ($flag & 2048))             ) # secondary or supplementary alignment
        {next;}
        if ($flag & 4) # unmapped
        {
                $unmapped++;
        }
        elsif ($flag & 8) # mapped, but mate unmapped = singleton
        {
                $unpaired++;    # here I count many more when I do if instead of elsif
                if ($flag & 1024)       #duplicate
                {
                        $singledups++;
                }
        }
        else    # both mapped
        {
                if ($flag & 64)  # first in pair bc. we count pairs, not reads => ignore second read
                {
                        $pairs++;       # here I cound some more when I do if instead of elsif above
                        if ($flag & 1024)       #duplicate
                        {
                                $pairdups++;
                        }
                }
        }
}
close IN;
# duplicates are UNPAIRED_READ_DUPLICATES + 2*READ_PAIR_DUPLICATES
# the denominator is mapped reads: UNPAIRED_READS_EXAMINED + 2*READ_PAIRS_EXAMINED

my $alldups=$singledups+2*$pairdups;
my $mappedreads=$unpaired+2*$pairs;
my $duprate;
if ($mappedreads != 0) {
        $duprate = sprintf ("%.3f", ($alldups/$mappedreads));
} else {   # No mapped reads, e.g. with test-data.
        $duprate = "NA";
}
print "##METRICS\nLIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE\n";
print "$libname\t$unpaired\t$pairs\t$unmapped\t$singledups\t$pairdups\t0\t$duprate\tNA\n";

exit;
