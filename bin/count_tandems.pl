#!/usr/bin/perl
# count_tandems.pl
# January 24, 2012
# Kamil Slowikowski
#
# Read output from mergBed, find the most frequent word in the tandem
# (contiguous) monomer repeats, and output a tab-delimited report.
#
# Usage: mergeBed -i <(bzcat out/minisatellites.gff3.bz2) -d 5 -nms | \
#        cut -f4 | ./bin/count_tandems.pl > out/tandem_counts.tab

use strict;
use warnings;

# Count the frequency of each tandem: { ABA => 12, ABBA => 25, ... }
my %tandems;

# Read the output from mergeBed that lists the names of merged entries.
while (<STDIN>) {
    chomp;
    # Extract the single letter for the minisatellite.
    my @monomers = map { /Gm_ms_(\w)/; uc $1 } split /;/;
    # Count this tandem.
    $tandems{ join('', @monomers) }++;
}

# These are the monomers in each tandem.
my @mons = qw/A B C D E/;

# Print the output table, sorted by descending tandem length.
for my $tandem (sort { length($b) <=> length($a) } keys %tandems) {
    # Skip empty tandems.
    next unless $tandem;
    my $tandemlen = length $tandem;
    
    # Get the most frequent 3-7 character word in the tandem.
    my ($word, $count) = get_most_frequent_word($tandem);
    
    print join("\t",
        $tandems{$tandem},
        $word,
        $count,
        length($tandem),
        $tandem
    ), "\n";
}

# Input a string of letters (monomers).
# Output an array with two elements:
# * the most frequent word (lengths 3-7)
# * the number of times it occurs in the string
sub get_most_frequent_word {
    my $tandem = shift;
    my $tandemlen = length $tandem;
    
    my @noresult = ('-', 0);
    
    if ($tandemlen <= 3) {
        return @noresult;
    }
    
    my %words;
    
    # Find all possible words in this string with lengths 3-7.
    for my $wordlen (3 .. 7) {
        my $lasti = $tandemlen - $wordlen;
        
        for (my $starti = 0; $starti <= $lasti; $starti++) {
            $words{ substr($tandem, $starti, $wordlen) } ||= 0;
        }
    }
    
    # Count the occurrences of all our words.
    for my $word (keys %words) {
        while ($tandem =~ /$word/g) {
            $words{$word}++;
        }
    }
    
    my $best_word = (sort {
        # word with most coverage comes first
        $words{$b} * length($b) <=> $words{$a} * length($a) ||
        
        # more frequent word comes first
        # $words{$b} <=> $words{$a} ||
        
        # longer word comes first
        # length($b) <=> length($a) ||
        
        # word that sorts first alphabetically comes first
        $a cmp $b
    } keys %words)[0];
    
    # When reporting the word, sometimes it's nice to report the 
    # rotated version that would come first when all the rotations are 
    # sorted alphabetically.
    # This makes output consistent and easier to interpret.
    # my $rotated_word = (sort map {
        # substr($best_word, $_) . substr($best_word, 0, $_)
    # } 0 .. length($best_word) - 1)[0];
    
    if ($words{$best_word} > 1) {
        return ($best_word, $words{$best_word});
    }
    
    return @noresult;
}
