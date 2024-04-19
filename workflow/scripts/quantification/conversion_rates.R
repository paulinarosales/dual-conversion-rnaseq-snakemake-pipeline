
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(data.table)
        library(readr)
        library(dplyr)
        })

# ---------- Snakemake parsing ---------- #
collapsedTSV <- snakemake@input[["collapsedTSV"]]
transcriptBED <- snakemake@input[["transcriptBED"]]
outputTSV <- snakemake@output[[1]]

# ---------- Functions ---------- #
computeRates <- function(all_rates_tab, ref_base) {
        cat(paste("Calculating conversion rates for mutations of", ref_base, "as reference base..."), sep="\n")
        base_content <- paste0("n", ref_base)
        regex <- paste0("^", ref_base, "_")
        n_rates <- collapsed_conv %>% select(matches(regex)) # select all combinations of mutations with same original base mutated 
                                                                                # (e.g.  T_C, T_A, T_G for T as original)
        # n_rates$Ncontent <- rowSums(n_rates)
        print(head(n_rates))
        n_rates <- (n_rates/collapsed_conv[[base_content]]) * 100 # number of mutations of a given base / total number of the original base detected (e.g. T_C/nT)
        print(head(n_rates))
        # n_rates$Ncontent <- NULL
        n_rates <-  round(n_rates, 2)
        return(n_rates)
        cat("\n")
}

# ---------- Main code ---------- #
cat("Reading input data...", sep="\n")
collapsed_conv <- read.table(collapsedTSV, header = TRUE, sep = "\t")
# collapsed_conv <- collapsed_conv %>% select(-ends_with("_N"))

tx_info <- read.table(transcriptBED, header = TRUE, sep = "\t")
rownames(tx_info) <- tx_info$transcript_id

all_rates_tab <- tx_info[collapsed_conv$GF,]
cat("\n")


cat("Merging genecounts assignment with conversion counts...", sep="\n")
ref_code <- c("A", "T", "G", "C") # ignore Ns
for(ref_base in ref_code){
        n_rates_tab <- computeRates(all_rates_tab, ref_base)
        all_rates_tab <- cbind(all_rates_tab, n_rates_tab)
}
cat("\n")

cat(paste0("Saving outputs:\n\t- ", outputTSV), sep="\n")
write.table(all_rates_tab, file = outputTSV, row.names= FALSE, quote = FALSE, sep = "\t")
cat("\n")

cat("DONE!", sep="\n")