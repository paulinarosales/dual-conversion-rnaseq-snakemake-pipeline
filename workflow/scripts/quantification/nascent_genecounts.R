log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(dplyr)
        library(readr)
        library(tibble)
})

# ---------- Snakemake parsing ---------- #
bakr_metaTSV <- snakemake@input[["bakr_metaTSV"]]
feature_ctsTSV <- snakemake@input[["feature_ctsTSV"]]
genesetTSV <- snakemake@input[["genesetTSV"]]

nascent_gene_ctsTSV <- snakemake@output[["nascent_gene_ctsTSV"]]

min_reads_th <- snakemake@params[["min_reads"]]
min_conv_th <- snakemake@params[["min_conv"]]

# ---------- Main code ---------- #

cat("Reading input data...", sep="\n")
# Gene annotation
geneset <- read.table(genesetTSV, header = TRUE,  sep = '\t',  stringsAsFactors = FALSE)
geneset <- geneset %>% select(c("gene_id", "gene_name"))

# BakR mut table per read
conv_tab <- read.table(bakr_metaTSV, stringsAsFactors=FALSE, header = TRUE, sep="\t")
conv_tab <- conv_tab %>% select(c("GF", "T_C"))
names(conv_tab) <- c("gene_id", "T_C")

# Gene readcounts
feature_cts <- read.table(feature_ctsTSV, stringsAsFactors=FALSE, header = TRUE, sep="\t", comment.char = "#")
feature_cts <- feature_cts %>% select(c(1, 7))
head(feature_cts)
names(feature_cts) <- c("gene_id", "total_counts")

# initialize output table
# cts_tab <- data.frame(gene_id = unique(geneset$gene_id),
#                       total_counts = numeric(length(unique(geneset$gene_id))),
#                       nascent_counts = numeric(length(unique(geneset$gene_id))))
# cat("\n")

cat("Assigning total reads from featureCounts...", sep="\n")
cts_tab <- feature_cts
cts_tab$nascent_counts <- 0
# cts_tab$total_counts <- feature_cts[match(cts_tab$gene_id, feature_cts$gene_id, nomatch = 0), "counts"]
cat("\n")



# low conversion distribution (for printing only)
tc_counts_dist <- data.frame("T_to_C_conv" = c("0", "1", "2", ">2"),
                             "n_reads" = numeric(4))

tc_counts_dist$n_reads <- c(sum(conv_tab$T_C == 0),
                            sum(conv_tab$T_C == 1),
                            sum(conv_tab$T_C == 2),
                            sum(conv_tab$T_C > 2))

cat(paste0("Filtering out reads with <", min_conv_th, " T-to-C converions nascent reads.\nNascent reads distribution:"), sep="\n")
print(tc_counts_dist)
# count filtering according to min conversions per read threshold
nascent_cts <- conv_tab %>% filter(T_C >= min_conv_th) # remove reads with < min_conv_th T-to-C conversions
cat("\n")

cat("Counting nascent reads per gene...", sep="\n")
nascent_cts <- as.data.frame(table(nascent_cts$gene_id)) # Count gene_id instances
names(nascent_cts) <- c("gene_id", "nascent_counts")
head(nascent_cts)
nrow(nascent_cts)

cat(paste("Mergening nascent counts with total counts from featureCounts..."), sep="\n") 
cts_tab$nascent_counts <- nascent_cts[match(cts_tab$gene_id, nascent_cts$gene_id), "nascent_counts"] # ensure right order
cts_tab[is.na(cts_tab)] <- 0  # if gene_id not found assign 0
cat("\n")

cat(paste("Filtering out gene entries with 0 read counts."), sep="\n")
removed <- nrow(cts_tab)
cts_tab <- cts_tab %>% filter(total_counts > 0)
removed <- removed - nrow(cts_tab)
cat(paste("Removed", removed, "entries"), sep="\n")
cat("\n")

cat("Calculating nascent fraction per gene...", sep="\n")
cts_tab$nascent_fraction <- round(cts_tab$nascent_counts/cts_tab$total_counts,2)
# cts_tab <- cts_tab[order(cts_tab$nascent_fraction, decreasing = TRUE),]
cat("\n")

cat(paste("Adding gene annotation..."), sep="\n") 
cts_tab <- inner_join(geneset, cts_tab, by="gene_id")
cat("\n")

cat(paste0("Saving outputs:\n\t- ", nascent_gene_ctsTSV), sep="\n")
write.table(cts_tab, file=nascent_gene_ctsTSV, sep="\t", quote=FALSE, row.names=FALSE) # Gene level
cat("DONE!", sep="\n")