log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(dplyr)
        library(readr)
        library(tibble)
})

# ---------- Snakemake parsing ---------- #


sample_manifestTSV <- snakemake@input[["sample_manifestTSV"]]
fraction_ctsTSV <- snakemake@input[["fraction_ctsTSV"]]

decayTSV <- snakemake@output[[1]]

# ---------- Main code ---------- #

cat("Reading input data...", sep="\n")

sample_man <- read.table(sample_manifestTSV, header=TRUE, sep="\t", check.names=FALSE)
sample_man <- sample_man %>% filter(Batch == 2023)
sample_man$Sample <- paste(sample_man$Sample_type, sample_man$Treatment, "Chase-time", sample_man$Chase_time_h, "Bio-rep", sample_man$Bio_rep, sep="_")
sample_man$Group <- paste(sample_man$Sample_type, sample_man$Treatment, sep="_")

fract_t <- read.table(fraction_ctsTSV, header=TRUE, sep="\t", check.names=FALSE)
colnames(fract_t)
cat("\n")

decay_d <- fract_t %>% select(c(1:8))
decay_d[ ,unique(sample_man$Group)] <- 0
head(decay_d)

cat("Calculating decay fold change...", sep="\n")
for(group in unique(sample_man$Group)){
        sample_0h <- paste(group, "Chase-time_0h_Bio-rep_1", sep="_")
        sample_3.5h <- paste(group, "Chase-time_3-5h_Bio-rep_1", sep="_")
        print(sample_0h)
        print(sample_3.5h)
        decay_d[[group]] <- fract_t[[sample_0h]] - fract_t[[sample_3.5h]]
        print(summary(decay_d[[group]]))
}

head(decay_d)
cat("\n")

cat(paste0("Saving outputs:\n\t- ", decayTSV), sep="\n")
write.table(decay_d, file=decayTSV, sep="\t", quote=FALSE, row.names=FALSE)
cat("DONE!", sep="\n")