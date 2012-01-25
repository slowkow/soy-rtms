#!/usr/bin/perl
# analyze_putative_element_blasttable.pl
# January 24, 2012
# Kamil Slowikowski
#
# Read the blasttable for each putative element and determine if it 
# belongs to a known transposable element family. Alignments against 
# the minisatellite cluster are counted, but they do not influence the 
# decision about the family to which the putative element is assigned. 
# Alignments that overlap at least one flank and meet the criteria of 
# $MIN_LEN and $MIN_PERC determine the family to which the putative 
# element is assigned.
#
# Usage: cat file.blasttable | ./bin/analyze_putative_element_blasttable.pl

use strict;
use warnings;

use Getopt::Long::Descriptive;

my ($opt, $usage) = describe_options(
    '%c [options...] < file.gff3.pair',
    [
        'flank-length|f:i',
        'length in nt of putative element up and downstream flanks',
        { default => 500 }
    ],
    [
        'min-length|l:i',
        'minimum hit length to count a soytedb family',
        { default => 400 }
    ],
    [
        'min-percent|p:i',
        'minimum hit percent identity to count a soytedb family',
        { default => 80 }
    ],
    [
        'gff3|g',
        'output GFF3 format instead of tab-delimited family counts'
    ],
    [ 'help|h', 'print this usage message and exit' ],
);

$usage->die if $opt->help;

# length of upstream and downstream flanks in each putative element
my $FLANK_LEN = $opt->flank_length;

# minimum length of BLAST alignment that overlaps a flank
# to be considered a good match
my $MIN_LEN   = $opt->min_length;

# minimum percent identity of BLAST alignment to be considered a good match
my $MIN_PERC  = $opt->min_percent;

# these are the headers of a blasttable tab-delimited file (-outfmt 7)
my @header = qw/ query_id subject_id percent_identity alignment_length 
    mismatches gap_opens q_start q_end s_start s_end evalue bitscore /;

# hash with all information for each putative element, looks like this:
# %pseudos = {
#    "query_id01" => {
#        ms_hits => 10,
#        ms_matches => 1,
#        flank_hits => 25,
#        flank_matches => 2,
#        length => 5891,
#        families => { NA => 0, "family01" => 3},
#    },
#    "query_id02" => {
#        ms_hits => 20,
#        ms_matches => 5,
#        flank_hits => 20,
#        flank_matches => 5,
#        length => 1000,
#        families => { NA => 1, "family04" => 2, "family07" => 15},
#    },
# }
my %pseudos;

# grab the query id from the comments
my $query_id;

# read a BLAST table (-outfmt 7)
while (<STDIN>) {
    chomp;
    
    if (/^#\s*Query:/) {
        ($query_id) = /Query:\s+(\S+)/;
        
        if ($pseudos{$query_id}) {
            die "Query $query_id found more than once in blasttable!\n";
        }
        
        $pseudos{$query_id} = {
            # a hit is counted as a match when the identity > $MIN_PERC
            
            # hits and matches against the minisatellite cluster
            ms_hits    => 0, ms_matches    => 0,
            
            # hits and matches against the up and downstream flanks
            flank_hits => 0, flank_matches => 0,
            
            # get the pseudoelement length from its id
            length     => do { $query_id =~ /:(\d+)-(\d+)$/; 1 + abs $1 - $2 },
            
            # give it a default family 'NA' to avoid undef later
            families   => { NA => 0 },
        };
    }
    
    # skip all other comments and blank lines
    next if (/^#/ || /^\s*$/);
    
    # split on tab
    my %H;
    @H{@header} = split /\t/;
    
    if ($query_id ne $H{query_id}) {
        die "$query_id should equal $H{query_id}\n";
    }
    
    # convenient reference to this member of the hash
    my $pseudo = $pseudos{ $H{query_id} };
    
    # hit aligned against just the minisatellite cluster, not against
    # the upstream or downstream flanks
    if ($FLANK_LEN < $H{q_start}
    && $H{q_end} < $pseudo->{length} - $FLANK_LEN) {
        
        $pseudo->{ms_hits}++;
        
        # the minisatellite cluster may be short,
        # so don't enforce MIN_LEN
        if ($H{percent_identity} > $MIN_PERC) {
            $pseudo->{ms_matches}++;
        }
    } 
    # hit aligned against the flanks of the pseudoelement
    else {
        $pseudo->{flank_hits}++;
        
        if ($H{alignment_length} > $MIN_LEN
        && $H{percent_identity} > $MIN_PERC) {
            
            $pseudo->{flank_matches}++;
            
            # grab the family name from the subject_id of the alignment
            my ($subject_family) = $H{subject_id} =~ /^.+_(.+)_.+$/;
            $pseudo->{families}->{$subject_family}++;
        }
    }
}

# Now we'll assign a family for each putative element.

# This is the GFF3 format...
# Gm01	Perl	pe_Gmr9_Gm01-1	6539725	6542107	.	.	.	.
# Gm01	Perl	pe_Gmr9_Gm01-2	6229053	6230507	.	.	.	.
# Gm01	Perl	pe_Gmr9_Gm01-3	6225800	6226826	.	.	.	.
# Gm01	Perl	pe_Gmr9_Gm01-4	6187853	6188902	.	.	.	.
# Gm01	Perl	pe_NA_Gm01-1	2860010	2861046	.	.	.	.
if ($opt->gff3) {
    # number of putative elements assigned to this family
    my %familypes;
    
    for my $query_id (keys %pseudos) {
        
        # get a convenient reference to this putative element's hash
        my %pseudo = %{ $pseudos{$query_id} };
        
        # sort the element's families by flank_matches count, high to low
        my @families =
            sort { $pseudo{families}{$b} <=> $pseudo{families}{$a} }
                keys %{ $pseudo{families} };
        
        my $max_family = shift @families;
        
        my ($chrom, $start, $end) = split(/[:-]/, $query_id);
        
        print join("\t",
            $chrom,
            'Perl',
            "pe_${max_family}_$chrom-".++$familypes{$max_family},
            $start,
            $end,
            ('.')x4,
        ), "\n";
    }
}
# This is the family tab-delimited format...
# query_id                ms_hits  >80%  fl_hits  >80% && >400nt  families  max_family  fl_matches
# Gm01:13898374-13900167  2970     1.00  563      0.46            2         Gmr9        260         Gmr6   1    NA     0
# Gm01:10369013-10370090  0        0.00  460      0.93            4         Gmr9        427         Gmr37  1    Gmr19  1   Gmr35   1  NA    0
# Gm01:11353495-11354520  0        0.00  140      0.45            5         Gmr139      59          Gmr4   1    Gmr3   1   Gmr431  1  Gmr7  1  NA     0
# Gm01:14146368-14148640  2006     0.97  687      0.72            2         Gmr9        493         Gmr21  2    NA     0
# Gm01:10184067-10185118  0        0.00  401      0.68            9         Gmr9        148         Gmr37  106  Gmr5   10  Gmr4    3  uuu   1  Gmr35  1  Gmr243  1  Gmr22  1  Gmr25  1  NA  0
# Gm01:14012504-14013531  0        0.00  292      0.03            4         Gmr9        4           Gmr4   2    Gmr1   1   Gmr5    1  NA    0
# Gm01:11777945-11779050  2        1.00  440      0.87            4         Gmr9        378         Gmr37  3    Gmr19  2   Gmr5    1  NA    0
# Gm01:10659549-10660836  2        1.00  82       0.05            2         Gmr338      2           Gmr9   2    NA     0
# Gm01:12491791-12494252  5190     0.96  1467     0.47            2         Gmr9        684         Gmr6   4    NA     0
else {
    # print header
    print join("\t",
        'query_id',
        'ms_hits',
        ">$MIN_PERC%",
        'fl_hits',
        ">$MIN_PERC% && >${MIN_LEN}nt",
        'families',
        'max_family',
        'fl_matches',
    ), "\n";
    for my $query_id (keys %pseudos) {
        
        # get a convenient reference to this putative element's hash
        my %pseudo = %{ $pseudos{$query_id} };
        
        # sort the element's families by flank_matches count, high to low
        my @families =
            sort { $pseudo{families}{$b} <=> $pseudo{families}{$a} }
                keys %{ $pseudo{families} };
        
        my $num_families = scalar @families - 1;
        
        # print the most frequent families for each putative element
        my $max_family = shift @families;
        my $previous_max = $pseudo{families}{$max_family};
        
        my @other_families;
        for my $family (@families) {
            my $next_max = $pseudo{families}{$family};
            
            push @other_families, "$family\t$next_max";
            $previous_max = $next_max;
        }
        # calculate the percent of hits that are matches
        my $ms_perc = sprintf('%.2f',
            $pseudo{ms_hits} > 0
            ? $pseudo{ms_matches} / $pseudo{ms_hits}
            : 0);
        my $flank_perc = sprintf('%.2f',
            $pseudo{flank_hits} > 0
            ? $pseudo{flank_matches} / $pseudo{flank_hits}
            : 0);
        
        print join("\t",
            (
                $query_id,
                $pseudo{ms_hits},
                $ms_perc,
                $pseudo{flank_hits},
                $flank_perc,
                $num_families,
                $max_family,
                $pseudo{families}{$max_family},
                join("\t", @other_families),
            )
        ), "\n";
    }
}
