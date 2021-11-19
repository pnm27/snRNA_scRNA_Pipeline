#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(ggsci)
library(movAPA)
library(umap)
library(knitr)
library(scAPAtrap)
library(argparse)
library(movAPA)


parser <- ArgumentParser(description='Rscript to run scAPAtrap on bam files produced by STARsolo(default) on 10x data(default)')
parser$add_argument('samtools_path', help='Path to the samtools executable (linux)')
parser$add_argument('fC_path', help='Path to the featureCounts executable (linux)')
parser$add_argument('bam_file_path', help='Path to the bam file')
parser$add_argument('barcodes_file_path', help='Path to the barcodes file')

parser$add_argument('-t', '--threads', type=int, default=8, help='Specify number of threads to use (Default: 8)'
parser$add_argument('-x', '--seqmethod', action='store_true', help='Is the data from 10X squencing? (Default: True)')
parser$add_argument('-o', '--outdir', default=getwd(),help='Output dir (Default: Current Working Directory)')
parser$add_argument('-c', '--chrpref', action='store_false', help='Are the chromosomes named like "chr1", "chr2", etc? (Default: False)')
parser$add_argument('-l', '--readlen', type=int, default=90, help='Read Length (Default: 90)')
parser$add_argument('-m', '--maxwidth', type=int, default=1000, help='Max width beyond which th peak will be iteratively split')
parser$add_argument('-c', '--mincells', type=int, default=10, help='Peak is expressed by at least these number of cells (Default: 10)')
parser$add_argument('-e', '--minpeaks', type=int, default=0, help='Minimum peak expression (Default: 0)')
parser$add_argument('-p', '--pacs', help='PACs file')


# Parse command arguments
args <-parser$parse_args()

# positional arguments
sam_p <- args$samtools_path
fC_p <- args$fC_path
bam_f <-args$bam_file_path
bc_f <- args$barcodes_file_path

# optional arguments
outputdir <- args$outdir

nextinput <- findUniqueMap(sam_p, input = bam_f, thread = args$threads, index = T, sort = T)
#nextinput <- dedupByPos(umitools.path, nextinput, TenX = args$seqmethod)
nextinput <- separateBamBystrand(sam_p, nextinput, args$threads)

chrs <- c(as.character(1:22),'X','Y')
if (args$c == 'True'){
chrs <- paste('chr', chrs, sep='')
}
maxwidth <- args$maxwidth
readlength <- args$readlen
outputdir <- args$outdir
fullcovF <- loadBpCoverages(nextinput[1],chrs)
fullcovR <- loadBpCoverages(nextinput[2],chrs)

forwardPeaks <-findPeaks(fullcovF, '+', readlength, maxwidth)
reversePeaks <-findPeaks(fullcovR, '-', readlength, maxwidth)

head(forwardPeaks)
head(reversePeaks)

peaksfile <- generateSAF(forwardPeaks, reversePeaks, outputdir)

final_bam <- generateFinalBam(fC_p, sam_p, bam_f, peaksfile, args$threads)
counts_tsv <- countPeaks('umi_tools', final_bam, outputdir, TenX=args$seqmethod)

tails <- findTails(bamfile=bam_f)

barcode <- read.delim2(bc_f,header = F)
barcode <- gsub('-[0-9]','',barcode$V1)
expma<- generatescExpMa(counts_tsv, peaksfile, barcode, tails, min.cells = 2,min.count = 0)

coldata <- data.frame(group = colnames(expma)[7:ncol(expma)], row.names = colnames(expma)[7:ncol(expma)])
scPACds <- movAPA::readPACds(pacfile = expma, coldata = coldata)

saveRDS(scPACds, file=args$pacs)
