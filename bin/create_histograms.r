# create_histograms.r
# January 24, 2012
# Kamil Slowikowski
#
# Create the following histograms:
#       - count of minisatellites by length in nt
#       - count of minisatellites by percent identity
#       - count of minisatellites and soyTEdb elements on each chr
#
# Usage: R --no-save < bin/create_histograms.r

library(ggplot2)

blasttable = read.delim(
    file          = "out/minisatellites.blasttable.bz2",
    header        = F,
    col.names     = c("query_id",
                    "subject_id",
                    "percent_identity",
                    "alignment_length",
                    "mismatches",
                    "gap_opens",
                    "query_start",
                    "query_end",
                    "subject_start",
                    "subject_end",
                    "evalue",
                    "bit_score"),
    comment.char = "#"
)

minisatellites = read.delim(
    file      = "out/soy.minisatellites.hist.bz2",
    header    = F,
    col.names = c("chr", "start", "end", "count", "nt", "length", "perc")
)
minisatellites = transform(minisatellites, type = "minisatellite")

soytedb = read.delim(
    file      = "out/soy.soytedb.hist.bz2",
    header    = F,
    col.names = c("chr", "start", "end", "count", "nt", "length", "perc")
)
soytedb = transform(soytedb, type = "soytedb")

ms_and_te = rbind(minisatellites, soytedb)

# tab-delimited file with length of each chromosome
genome = read.delim(
    file      = "data/genome",
    header    = F,
    col.names = c("chr", "length")
)

create_histograms = function() {
    width   = 800 # px
    height  = 600 # px
    res     = 100 # dpi
    
    # a histogram for each minisatellite showing count of HSPs with each length
    file    = "png/minisatellites_length_histogram.png"
    title   = "Count of each minisatellite by length"
    
    png(file = file, width = width, height = height, res = res)
    print(
        qplot(
            data   = blasttable,
            x      = alignment_length,
            geom   = "histogram",
            facets = . ~ query_id,
            xlab   = "length (nt)",
            main   = title
        )
    )
    dev.off()
    
    # a histogram for each minisatellite showing count of HSPs with each perc
    # identity
    file = "png/minisatellites_identity_histogram.png"
    title = "Count of each minisatellite by percent identity"
    
    png(file = file, width = width, height = height, res = res)
    print(
        qplot(
            data = blasttable,
            x = percent_identity,
            geom = "histogram",
            facets = . ~ query_id,
            xlab = "percent identity",
            main = title
        )
    )
    dev.off()
    
    # a histogram of each chromosome showing the count of soyTEdb elements
    # and minisatellites per 1e5 nucleotides
    for (i in c(1:20)) {
        chrom = sprintf("Gm%02i",i)
        length = genome[i,]$length
        data = subset(ms_and_te, chr == chrom)
        file = paste(
            sep = "",
            "png/", chrom, "_histogram.png"
        )
        title = paste(
            chrom,
            "Count of minisatellites and transposable elements per 1e5 nt"
        )
        
        png(file = file, width = width, height = height, res = res)
        print(
            qplot(
                data = data,
                x = end,
                weight = count,
                binwidth = 100000,
                facets = type ~ .,
                main = title,
                xlab = "position (nt)"
            ) +
            coord_cartesian(xlim = c(1, length)) +
            scale_x_continuous(breaks = seq(0, length, 5000000))
        )
        dev.off()
    }
}
create_histograms()
