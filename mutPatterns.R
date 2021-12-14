library(ggplot2)
library(reshape2)
library(MutationalPatterns)
library(BSgenome)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
library(ref_genome, character.only = TRUE)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)
snvFile <- args[1]
sampleName <- args[2]
outputPdf <- args[3]

setwd("/mnt/codon/DATA/DavidOConnor_tALL/snv_benchmarking/sage_100/")
snvFile <- "P003_diag_sage100.vcf.gz"
sampleName <- "P003_diag"
outputPdf <- "test.pdf"

vcf_files <- list(snvFile)
sample_names <- list(sampleName)

grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
type_occurrences <- mut_type_occurrences(grl, ref_genome)
plot.muts <- plot_spectrum(type_occurrences, CT = TRUE, legend = TRUE, error_bars = 'none')

mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)

signatures = get_known_signatures()
fit_res <- fit_to_signatures(mut_mat, signatures)
plot.cosmic <- plot_contribution(fit_res$contribution,
                  coord_flip = FALSE,
                  mode = "absolute")
plot.cosine <- plot_original_vs_reconstructed(mut_mat, fit_res$reconstructed, 
                               y_intercept = 0.95)

plot.profile <- plot_96_profile(mut_mat)
chromosomes <- seqnames(get(ref_genome))
plot.rainfall <- plot_rainfall(grl[[1]], chromosomes=chromosomes,  cex = 1.5, 
                               ylim = 1e+09, type = "snv")

pdf(outputPdf, width=12, height=10)
plot_grid(plot_grid(plot.muts, plot.cosmic, nrow=1), plot.profile, plot.rainfall, ncol=1)
dev.off()
