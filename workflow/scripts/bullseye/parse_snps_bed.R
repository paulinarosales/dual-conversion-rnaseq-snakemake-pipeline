log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(dplyr)
        library(readr)
        library(tibble)
        library(tidyr)
})

# ---------- Snakemake parsing ---------- #

snpTXT <- snakemake@input[[1]]

snpBED <- snakemake@output[[1]]

# ---------- Main code ---------- #

cat("Reading input data...", sep="\n")
snp_tab <- read.table(snpTXT, header = FALSE,  sep = '\t', comment.char = "")
cat("\n")

cat("Parsing input table...", sep="\n")
parse_tab <- separate(snp_tab, 
                    col=V1, into=c("ref", "mut", "chr", "start"), sep='\\:')

head(parse_tab)
parse_tab$conv <- paste0(parse_tab$ref, parse_tab$mut)
head(parse_tab)  

cat("Extracting C-to-T SNPs...", sep="\n")
parse_tab <- parse_tab %>% filter(conv == "CT")
parse_tab$start <- as.numeric(parse_tab$start)
head(parse_tab)
parse_tab <- parse_tab %>% select(c("chr", "start"))
parse_tab$end <- parse_tab$start + 1
head(parse_tab)
cat("\n")

cat(paste0("Saving outputs:\n\t- ", snpBED), sep="\n")
write.table(parse_tab, file=snpBED, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
cat("DONE!", sep="\n")