#!/usr/bin/perl
# blasttable2gff3.pl
# January 24, 2012
# Kamil Slowikowski
#
# Convert a blasttable (-outfmt 7) to GFF3.
#
# See: http://gmod.org/wiki/GFF3
#
# Usage: ./bin/blasttable2gff3.pl < out/file.blasttable > out/file.gff3

use strict;
use warnings;

# hash for labelling every HSP with an id number
my %count;

while (<STDIN>) {
    chomp;
    
    # skip comments
    next if /^#/;
    
    # remove all spaces, because bit_score sometimes starts with a space
    s/ //g;
    
    # split on tab
    my ($query_id, $subject_id, $percent_identity, $alignment_length, 
    $mismatches, $gap_opens, $query_start, $query_end, $subject_start, 
    $subject_end, $evalue, $bit_score) = split /\t/;
    
    # sanity check
    if ($query_start > $query_end) {
        die "$.: query_start > query_end: $query_start > $query_end\n";
    }
    
    # the coordinates define which strand was hit
    my $strand;
    if ($subject_end - $subject_start >= 0) {
        $strand = '+';
    } else {
        $strand = '-';
        # swap start and end
        my $temp       = $subject_start;
        $subject_start = $subject_end;
        $subject_end   = $temp;
    }
    
    # create key-value pairs
    my $attributes = join(';', (
            "query_id=$query_id",
            "percent_identity=$percent_identity",
            "alignment_length=$alignment_length",
            "mismatches=$mismatches",
            "gap_opens=$gap_opens",
            "query_start=$query_start",
            "query_end=$query_end",
            "bit_score=$bit_score",
        )
    );
    
    # increment the count for this query_id, so each HSP gets a unique name
    $count{$query_id}++;
    
    my $name = $query_id . '-' . $subject_id . '-' . $count{$query_id};

    # see http://gmod.org/wiki/GFF3 for details about the GFF3 format
    # seqid, source, type, start, end, score, strand, phase, attributes
    print
    join("\t", (
            $subject_id,
            'BLASTN',
            $name,            # this is useful to see tandems with mergeBed
            $subject_start,
            $subject_end,
            $evalue,
            $strand,
            '.',
            $attributes,
        )
    ), "\n";
}
