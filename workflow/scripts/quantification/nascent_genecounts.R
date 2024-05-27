log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(dplyr)
        library(readr)
        library(tibble)
})

# ---------- Snakemake parsing ---------- #
bakr_mergedTSV <- snakemake@input[["bakr_mergedTSV"]]
tx_infoTSV <- snakemake@input[["tx_infoTSV"]]
genesetTSV <- snakemake@input[["genesetTSV"]]

nascent_ctsTSV <- snakemake@output[["nascent_ctsTSV"]]

min_reads_th <- snakemake@params[["min_reads"]]
min_conv_th <- snakemake@params[["min_conv"]]

# ---------- Main code ---------- #

cat("Reading input data...", sep="\n")
conv_tab <- read.table(bakr_mergedTSV, stringsAsFactors=FALSE, header = TRUE, sep="\t")
conv_tab <- conv_tab %>% select(c("GF", "T_C"))

tx_info <- read.table(tx_infoTSV, header = TRUE,  sep = '\t',  stringsAsFactors = FALSE)
tx_info <- tx_info %>% select(c("transcript_id", "gene_id"))

geneset <- read.table(genesetTSV, header = TRUE,  sep = '\t',  stringsAsFactors = FALSE)

# initialize output table
cts_tab <- data.frame(transcript_id = unique(conv_tab$GF),
                      total_counts = numeric(length(unique(conv_tab$GF))),
                      nascent_counts = numeric(length(unique(conv_tab$GF)))
                      )
cat("\n")

cat("Counting trasncripts total reads...", sep="\n")
total_cts <- as.data.frame(table(conv_tab$GF))
total_cts <- total_cts %>% column_to_rownames("Var1")
cts_tab$total_counts <- total_cts[cts_tab$transcript_id, "Freq"]
cat("\n")

cat("Counting transcripts nascent reads...", sep="\n")
# low conversion distribution (for printing only)
tc_counts_dist <- data.frame("T_to_C_conv" = c("0", "1", "2", ">2"),
                             "n_reads" = numeric(4))

tc_counts_dist$n_reads <- c(sum(conv_tab$T_C == 0),
                            sum(conv_tab$T_C == 1),
                            sum(conv_tab$T_C == 2),
                            sum(conv_tab$T_C > 2))

# count filtering according to min conversions per read threshold
nascent_cts <- conv_tab %>% filter(T_C >= min_conv_th)
nascent_cts <- as.data.frame(table(nascent_cts$GF))
nascent_cts <- nascent_cts %>% column_to_rownames("Var1")
cts_tab$nascent_counts <- nascent_cts[cts_tab$transcript_id, "Freq"]

cat(paste("Ignoring reads with <", min_conv_th, "T-to-C converions nascent reads.\nNascent reads distribution:"), sep="\n")
print(tc_counts_dist)
cat("\n")

cat(paste("Collapsing counts according to gene ID..."), sep="\n") 
cts_tab <- inner_join(tx_info, cts_tab, by="transcript_id")
cts_tab$transcript_id <- NULL

cts_tab <- aggregate(.~gene_id, data = cts_tab, FUN = sum)

cat(paste("Filtering out gene entries with <", min_reads_th, "total reads."), sep="\n")
removed <- nrow(cts_tab)
cts_tab <- cts_tab %>% filter(total_counts >= min_reads_th)
removed <- removed - nrow(cts_tab)
cat(paste("Removed", removed, "gene entries"), sep="\n")

cts_tab <- inner_join(geneset, cts_tab, by="gene_id")
cat("\n")

cat("Computing total/nascent ratios...", sep="\n")
cts_tab$nascent_to_total_ratio <- round(cts_tab$nascent_counts/cts_tab$total_counts, 2)
cts_tab <- cts_tab[order(cts_tab$nascent_to_total_ratio, decreasing = TRUE),]
cat("\n")

cat(paste0("Saving outputs:\n\t- ", nascent_ctsTSV), sep="\n")
write.table(cts_tab, file=nascent_ctsTSV, sep="\t", quote=FALSE, row.names=FALSE)
cat("DONE!", sep="\n")