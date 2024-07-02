log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(dplyr)
        library(readr)
        library(tibble)
})

# ---------- Snakemake parsing ---------- #


merged_ctsTSV <- snakemake@input[[1]]

tpm_ctsTSV <- snakemake@output[[1]]

# ---------- Functions ---------- #

tpm_funct <- function(counts,len) {
                  x <- counts/len
                  return(t(t(x)*1e6/colSums(x)))
              }

# ---------- Main code ---------- #

cat("Reading input data...", sep="\n")
merged_cts <- read.table(merged_ctsTSV, header = TRUE,  sep = '\t',  stringsAsFactors = FALSE, check.names = FALSE) # initialize merged_cts with gene info
cat("\n")

cat("Calculating TPM...", sep="\n")
genes_metadata <- merged_cts %>% select(1:8) # select only cols with gene info
cts <- merged_cts %>% select(9:ncol(merged_cts)) # select only cols with counts
tpm_cts <- tpm_funct(cts, merged_cts$gene_length)
head(tpm_cts)
tpm_cts <- cbind(genes_metadata, tpm_cts)
head(tpm_cts)
cat("\n")

cat(paste0("Saving outputs:\n\t- ", tpm_ctsTSV), sep="\n")
write.table(tpm_cts, file=tpm_ctsTSV, sep="\t", quote=FALSE, row.names=FALSE)
cat("DONE!", sep="\n")