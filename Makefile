SHELL = /bin/bash

# # # #

# Input data files
# (data/minisatellites.fa is already included)
te_tab = data/soytedb.tab.bz2
te_gff3 = data/soytedb.gff3.bz2
te_fa = data/soytedb.fa.bz2
# chroms = $(shell seq -f"data/Gm%02g.fa.bz2" 20)
soy_fa = data/soy.fa

# BLAST databases
soy_blastdb = blast/soy.nhr blast/soy.nin blast/soy.nsq
te_blastdb = blast/soytedb.nhr blast/soytedb.nin blast/soytedb.nsq

# minisatellite mappings
ms_blast = out/minisatellites.blasttable.bz2
ms_gff3 = out/minisatellites.gff3.bz2

# the soybean genome partitioned into 100k bins
soy_bins = out/soy.hist.bed.bz2

# histograms of MS and TE frequency
soy_ms_hist = out/soy.minisatellites.hist.bz2
soy_te_hist = out/soy.soytedb.hist.bz2

# pngs
pngs = png/minisatellites_length_histogram.png \
png/minisatellites_identity_histogram.png \
$(shell seq -f"png/Gm%02g_histogram.png" 20)

# putative element annotations and sequences
pe_bed = out/putative_elements.bed.bz2
pe_fa = out/putative_elements.fa.bz2

# BLAST putative elements against soytedb, assign TE families
pe_te_blast = out/putative_elements.soytedb.blasttable.bz2
pe_te_tab = out/putative_elements.soytedb.blasttable.tab.bz2
pe_te_gff3 = out/putative_elements.soytedb.blasttable.gff3.bz2

# BLAST putative elements that are not assigned a family against soy
pe_soy_blast = out/putative_elements.soy.blasttable.bz2
pe_gss_blast = out/putative_elements.gss.blasttable.bz2

# final table with MS counts for each TE family
ms_count = out/minisatellite_counts.tab

# final table with monomers and tandem counts
tandem_count = out/tandem_counts.tab

# # # #

.PHONY : all
all : \
blast out png \
$(te_tab) $(te_fa) $(soy_fa) \
$(soy_blastdb) $(te_blastdb) \
$(ms_blast) $(ms_gff3) \
$(soy_bins) \
$(soy_ms_hist) $(soy_te_hist) \
$(pngs) \
$(pe_bed) $(pe_fa) $(pe_te_blast) $(pe_te_tab) $(pe_te_gff3) \
$(pe_soy_blast) $(pe_gss_blast) \
$(ms_count) $(tandem_count)

.PHONY : clean
clean :
	@echo "clean is not implemented yet!"

# # # #

# Directories

blast:
	mkdir blast

out:
	mkdir out

png:
	mkdir png

# # # #

# Input data files

$(soy_fa) :
	wget -O - "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&rettype=fasta&id=283570573,283570572,283570571,283570570,283570569,283570568,283570567,283570566,283570565,283570564,283570563,283570562,283570561,283570560,283570559,283570558,283570557,283570556,283570555,283570554" | \
	perl -pe 'BEGIN{ $$i="01" } s/^>.*/>Gm$$i/ && $$i++' > $@

## transposable element annotations, 423K
$(te_tab) :
	wget -O - "http://www.soybase.org/soytedb/Dump.php?type=summary" | \
	bzip2 > $@

## convert to GFF3 so BEDTools can read it
$(te_gff3) : $(te_tab)
	./bin/soytedb2gff3.pl < <(bzcat $(te_tab)) | bzip2 > $@

## transposable element sequences, 42M
$(te_fa) :
	wget -O - "http://www.soybase.org/soytedb/Dump.php?type=fasta" | \
	bzip2 > $@


# BLAST databases

## soy blast database
$(soy_blastdb) : $(soy_fa)
	cat $(soy_fa) | \
	makeblastdb -dbtype nucl -title soy -out blast/soy

## soytedb blast database
$(te_blastdb) : $(te_fa)
	bzcat $(te_fa) | \
	makeblastdb -dbtype nucl -title soytedb -out blast/soytedb

## map the minisatellite sequences to the soybean genome
$(ms_blast) : $(soy_blastdb)
	blastn \
	-query data/minisatellites.fa \
	-db blast/soy \
	-word_size 4 \
	-num_threads 2 \
	-outfmt 7 | \
	 bzip2 > $@

## convert the mappings to GFF3
$(ms_gff3) : $(ms_blast)
	./bin/blasttable2gff3.pl < <(bzcat $(ms_blast)) | bzip2 > $@


# Histograms

## partition the genome into a BED file with 100k nt bins
$(soy_bins) : data/genome
	./bin/partition_genome.pl data/genome 100000 | bzip2 > $@

## histogram data file for MS frequency across chromosomes
$(soy_ms_hist) : $(ms_gff3) $(soy_bins)
	coverageBed \
	-a <(bzcat $(ms_gff3)) \
	-b <(bzcat $(soy_bins)) | \
	bzip2 > $@

## histogram data file for TE frequency across chromosomes
$(soy_te_hist) : $(te_gff3) $(soy_bins)
	coverageBed \
	-a <(bzcat $(te_gff3)) \
	-b <(bzcat $(soy_bins)) | \
	bzip2 > $@

## png images of the histograms
$(pngs) : $(ms_blast) $(soy_ms_hist) $(soy_te_hist)
	R --slave --no-save < bin/create_histograms.r

# putative element annotations in BED format
$(pe_bed) : $(ms_gff3) $(te_gff3)
	intersectBed -v \
	-a <(bzcat $(ms_gff3)) \
	-b <(bzcat $(te_gff3)) | \
	mergeBed -d 1000 | \
	slopBed -g data/genome -b 500 | \
	bzip2 > $@

# putative element sequences
$(pe_fa) : $(pe_bed)
	bzcat $(pe_bed) | \
	fastaFromBed -bed stdin -fi data/soy.fa -fo stdout | \
	bzip2 > $@

# blasttable of putative elements against soytedb
$(pe_te_blast) : $(pe_fa) $(soy_blastdb)
	bzcat $(pe_fa) | \
	parallel --pipe --recstart '>' -N1 \
	blastn -db blast/soytedb -outfmt 7 | \
	bzip2 > $@

$(pe_te_tab) : $(pe_te_blast)
	bzcat $(pe_te_blast) | \
	./bin/analyze_putative_element_blasttable.pl | \
	bzip2 > $@

$(pe_te_gff3) : $(pe_te_blast)
	bzcat $(pe_te_blast) | \
	./bin/analyze_putative_element_blasttable.pl --gff3 | \
	bzip2 > $@

$(pe_soy_blast) : $(pe_te_gff3) $(soy_blastdb) $(soy_fa)
	bzcat $(pe_te_gff3) | \
	grep _NA_ | \
	fastaFromBed -bed stdin -fi $(soy_fa) -fo stdout | \
	blastn -num_threads 2 -db blast/soy -outfmt 7 | \
	bzip2 > $@

$(pe_gss_blast) : $(pe_soy_blast) $(soy_fa)
	bzcat $(pe_soy_blast) | \
	grep -v '^#' | \
	perl -wanE '$$H{$$F[0]}++; END{ do{ say join "\t", split /[:-]/ if $$H{$$_} < 10 } for keys %H }' | \
	fastaFromBed -bed stdin -fi $(soy_fa) -fo stdout | \
	blastn -remote -db gss -entrez_query '"Glycine max"[porgn:__txid3847]' -outfmt 7 | \
	bzip2 > $@

$(ms_count) : $(te_gff3) $(pe_te_gff3) $(ms_gff3)
	bzcat $(te_gff3) $(pe_te_gff3) | \
	intersectBed -a <(bzcat $(ms_gff3)) -b stdin -wo | \
	./bin/count_minisatellites.pl > $@

$(tandem_count) : $(ms_gff3)
	mergeBed -i <(bzcat $(ms_gff3)) -d 5 -nms | \
	cut -f4 | ./bin/count_tandems.pl > $@
