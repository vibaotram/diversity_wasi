library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(Rsamtools)
library(karyoploteR)

snp_info <- read.delim("./raw_genotype/snp_261_info.tsv")
snp_gr <- makeGRangesFromDataFrame(snp_info, 
                                   seqnames.field = "CHROM",
                                   start.field = "POS",
                                   end.field = "POS")

annot <- readGFFAsGRanges("/shared/projects/most_kmer/afterkiss/gea/annotation/CC1.8_v2_named_annotation_0.5_BTI_unique.gff3")

## number of SNPs on genes
length(snp_gr[snp_gr %within% annot])

##number of SNPs not on genes
length(snp_gr[!snp_gr %within% annot])

## mean distance
distanceToNearest(snp_gr)
mean(distance(snp_gr[-length(snp_gr)], snp_gr[-1]), na.rm = T)


genome_fa <- "/shared/projects/vietcaf/data/reference/CC1.8_v2_pseudomolecule_cat.fa"
genome <- scanFa(genome_fa)
genome <- genome[grepl("Chr", names(genome))] # remove contigs
names(genome) <- gsub(" .+", "", names(genome))
genome_GR <- GRanges(seqnames = names(genome), ranges = IRanges(start = 1, end = width(genome)))

# png(file = "./genotyping/lgc_results/snp_distributions.png", width = 1200, height = 700)
gkp <- plotKaryotype(genome = genome_GR, plot.type = 1) 
kpDataBackground(gkp, data.panel ="ideogram", color = "gray80")
kpAddBaseNumbers(gkp, tick.dist = 2e7, add.units = T)
kpPlotRegions(gkp, data = snp_gr, data.panel = "ideogram", col = "black", avoid.overlapping = T)
kpAddMainTitle(gkp, "Distribution of 261 genotyped SNPs in the genome")
