#!/usr/bin/env perl
# partition_genome.pl
# January 24, 2012
# Kamil Slowikowski
#
# Given a genome file (described below), print a BED format to stdout
# with each chromosome partitioned into bins of equal size.
#
# A genome file has 2 tab-delimited columns: chrom, size

use strict;
use warnings;

use List::Util qw(min);

if (scalar @ARGV ne 2) {
    die "Usage: ./bin/partition_genome.pl genome 100000\n";
}

my ($genome, $binsize) = @ARGV;

open my $fh, '<', $genome or die $!;

while (<$fh>) {
    chomp;
    my ($chrom, $size) = split;
    
    for (my $i = 1; $i < $size; $i += $binsize) {
        print join("\t",
            $chrom,
            $i,
            min($size, $i + $binsize - 1)
        ), "\n";
    }
}
