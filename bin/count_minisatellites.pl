#!/usr/bin/perl
# count_minisatellites.pl
# January 24, 2012
# Kamil Slowikowski
#
# Read output from intersectBed. Count minisatellites in each 
# transposable element family.

use strict;
use warnings;

use List::Util qw/sum/;
use Getopt::Long::Descriptive;
use Sort::Versions qw/versioncmp/;

my ($opt, $usage) = describe_options(
    'intersectBed -a ms.gff3 -b te.gff3 -wo | %c [options...]',
    [
        'percent-identity|p:i',
        'minimum percent identity required to count a minisatellite HSP',
        { default => 90 }
    ],
    [ 'help|h', 'print this usage message and exit' ],
);

$usage->die if $opt->help;

# Hash of all minisatellites, looks like:
# { 'Gm_ms_a-Gm20-1' =>
#       { elem => 'RLC_Gmr7_Gm1-1', desc => 'SOLO', len => 455 }
#   'Gm_ms_a-Gm20-2' =>
#       { ... }
# }
my %minis;

while (<STDIN>) {
    chomp;
    
    my @F = split /\t/;
    
    my %ms_attributes = split /[=;]/, $F[8];
    
    # Get attributes for transposable elements or empty array 
    # if attributes are missing.
    my %te_attributes = $F[17] =~ /[=;]/ ? split /[=;]/, $F[17] : ();
    
    # We are interested in counting only the BLAST hits
    # that have a high percenty identity.
    next if $ms_attributes{percent_identity} < $opt->percent_identity;
    
    my ($mini, $elem, $beg, $end) = @F[2,11,12,13];
    
    my $len = abs($end - $beg) + 1;
    
    # Check if the element is SOLO or INTACT,
    # or else assume INTACT if the description is undefined.
    my $desc = $te_attributes{desc} ? $te_attributes{desc} : 'INTACT';
    
    # SoyTEdb sometimes uses LTR_INTACT instead of INTACT, but
    # we'll treat both the same way.
    $desc = 'SOLO' if $desc =~ /solo/i;
    $desc = 'INTACT' if $desc =~ /intact/i;
    
    # The element's id says if it is a putative element.
    $desc = 'PUTATIVE' if $elem =~ /^pe/i;
    
    # Make sure we have a valid element, length, and description.
    die unless
        $elem && $len && grep { /^$desc$/ } qw/PUTATIVE SOLO INTACT/;
    
    # Want innermost nested element that contains this minisatellite.
    if (!$minis{$mini} || $len < $minis{$mini}->{len}) {
        $minis{$mini} = {
            elem => $elem,
            desc => $desc,
            len  => $len
        };
    }
}

# Count of minisatellites in each family, key is family.
my %fams;

for my $mini (keys %minis) {
    # Get one of: qw/ a b c d e /
    my ($m)   = $mini =~ /Gm_ms_(\w)/i;
    
    # Get the family name: Gmr1 Gmr2 Gmr3 ...
    my $elem  = $minis{$mini}->{elem};
    my ($fam) = $elem =~ /_(\w+)_/;
    
    # Initialize the family's hash to avoid undefined values later.
    $fams{$fam} ||= {
        a => { minis => 0, PUTATIVE => {}, SOLO => {}, INTACT => {} },
        b => { minis => 0, PUTATIVE => {}, SOLO => {}, INTACT => {} },
        c => { minis => 0, PUTATIVE => {}, SOLO => {}, INTACT => {} },
        d => { minis => 0, PUTATIVE => {}, SOLO => {}, INTACT => {} },
        e => { minis => 0, PUTATIVE => {}, SOLO => {}, INTACT => {} },
    };
    
    # Increment count of this minisatellite in this family.
    $fams{$fam}->{$m}->{minis}++;
    
    # Put this MS's parent element in PUTATIVE, SOLO, or INTACT hash.
    $fams{$fam}->{$m}->{ $minis{$mini}->{desc} }->{$elem}++;
}

# Build each line of output incrementally.
my @lines;

for my $fam (keys %fams) {
    # Build one line of output per family.
    my @line;
    push @line, $fam;
    
    # Loop over each MS in alphabetical order: a, b, c, d, e
    for my $m (qw/a b c d e/) {
        # Make a hashref for this family and this minisatellite.
        my $href = $fams{$fam}->{$m};
        
        # The number of PUTATIVE, SOLO, and INTACT elements
        # that contain this minisatellite.
        my $pelems  = scalar keys %{ $href->{PUTATIVE} };
        my $solos   = scalar keys %{ $href->{SOLO} };
        my $intacts = scalar keys %{ $href->{INTACT} };
        
        # The total number of non-putative elements that
        # contain this minisatellite.
        my $elems = $solos + $intacts;
        
        my $perc_solo = sprintf "%.2f", $elems > 0
            ? $solos / $elems 
            : 0;
        
        #             a          b          c             d              e
        # @line index (1,2,3,4), (5,6,7,8), (9,10,11,12), (13,14,15,16), (17,18,19,20)
        push @line, ($href->{minis}, $pelems, $elems, $perc_solo);
    }
    
    # build the table line by line
    push @lines, \@line;
}

# Print the header.
print join("\t",
    qw/ family
    A Ape Ae Aeps
    B Bpe Be Beps
    C Cpe Ce Ceps
    D Dpe De Deps
    E Epe Ee Eeps /
), "\n";

# Print the table, one line per family.
map { print join("\t", @$_), "\n" }
# Sort lines by descending total MS count,
# then alphabetically by family name.
sort {
    sum(@$b[1,5,9,13,17]) <=> sum(@$a[1,5,9,13,17]) ||
    versioncmp(@$a[0], @$b[0])
} @lines;
