#!/usr/bin/perl
# soytedb2gff3.pl
# January 24, 2012
# Kamil Slowikowski
#
# Convert the soyTEdb tab-delimited file to GFF3. Omit transposable 
# elements anchored to scaffolds, only the elements anchored to 
# chromosomes will be output.
#
# See: http://gmod.org/wiki/GFF3
#
# Usage: ./bin/soytedb2gff3.pl < <(bzcat data/soytedb.tab.bz2) | \
#        bzip2 > data/soytedb.gff3.bz2

use strict;
use warnings;

my $header = join("\t",
    'Element Name',
    'Reference',
    'Class',
    'Sub_Class',
    'Order',
    'Super_Family',
    'Family',
    'Description',
    'Chromosome',
    'Unanchored_Scaffold',
    'Start Position',
    'End Position'
);

while (<STDIN>) {
    chomp;
    
    # Future-proof this script by dying if the format has changed.
    if ($. == 1) {
        if ($_ ne $header) {
            die "Header does not match! SoyTEdb must have changed the file format!\n";
        }
    }
    
    my %F;
    
    # Split on tab.
    @F{qw/ name reference class
        subclass order super_family
        family description chromosome
        unanchored_scaffold start_pos end_pos /} = split /\t/;
    
    if (scalar values %F != 12) {
        die "Number of values in row $. not equal to 12!\n";
    }
    
    # It should be on a chromosome, not on a scaffold.
    next unless $F{chromosome} =~ /(\d+)/;
    
    my ($id) = $F{name} =~ /(\d+)$/;
    my $attributes = join(';', (
            "desc=$F{description}",
            "order=$F{order}",
            "sfamily=$F{super_family}",
            "family=$F{family}",
            "id=$id",
        )
    );
    
    # see http://gmod.org/wiki/GFF3 for details about the GFF3 format
    # seqid, source, type, start, end, score, strand, phase, attributes
    print join("\t", (
            $F{chromosome},
            'soyTEdb',
            $F{name},
            $F{start_pos},
            $F{end_pos},
            '.',
            '.',
            '.',
            $attributes,
        )
    ), "\n";
}
