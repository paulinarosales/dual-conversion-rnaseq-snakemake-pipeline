log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(dplyr)
        library(ggplot2)
        library(RColorBrewer)
        library(tidyr)
})

# ---------- Snakemake parsing ---------- #
tpmTSV <- snakemake@input[["tpmTSV"]]
sample_manifestTSV <- snakemake@input[["sample_manifestTSV"]]
goisTSV <- snakemake@input[["goisTSV"]]
# treatment_order <- c(snakemake@params[["treatment_order"]])

outputPDF <- snakemake@output[[1]]

filter_batch <- snakemake@params[["filter_batch"]]

# ---------- Color palette ---------- #
pawlette <- c("#e18727", # orange
		          "#0072b5", # dark blue
		          "#20854e", # dark green
		          "#bc3c29", # red
		          "#ffdb24", # yellow
		          "#7876b1", # purple
		          "#cd6090", # pink
		          "#6ca6cd", # light blue
		          "#698b22", # light green
		          "#8b5a2b") # brown


pawlette_light <- c("#F0B77A", # orange
		                "#77B4DB", # dark blue
		                "#6AB38D", # dark green
		                "#BF7C71", # red
		                "#FDE988", # yellow
		                "#A09FBC", # purple
		                "#CC9AB0", # pink
		                "#A5C2D4", # light blue
		                "#ADC183", # light green
		                "#AE8F73") # brown

pawlette_paired <- c(rbind(pawlette, pawlette_light))

# ---------- Snakemake parsing ---------- #
cat("Reading input data...", sep="\n")
gois <- read.table(goisTSV, header = TRUE, sep = "\t")
gois$gene_name <- NULL

sample_man <- read.table(sample_manifestTSV, header=TRUE, sep="\t", check.names=FALSE)
sample_man <- sample_man %>% filter(Batch == filter_batch)
sample_man$Sample <- paste(sample_man$Sample_type, sample_man$Treatment, "Chase-time", sample_man$Chase_time_h, "Bio-rep", sample_man$Bio_rep, sep="_")
sample_man$Group <- paste(sample_man$Sample_type, sample_man$Treatment, sep="_")
print(sample_man)


tpm_d <- read.table(tpmTSV, header=TRUE, sep="\t", check.names=FALSE)
tpm_d <- inner_join(gois, tpm_d, by="gene_id")
tpm_d <- tpm_d %>% select(c("gene_id", "gene_name", sample_man$Sample))
head(tpm_d)


if(filter_batch == 2023){
	cols <- pawlette_paired[1:length(unique(sample_man$Group))]
}else{
	cols <- pawlette[1:length(unique(sample_man$Group))]
}

cat("\n")

# ---------- Build dataframe ---------- #

pivot_tpm <- tpm_d %>% pivot_longer(cols = 3:ncol(tpm_d), names_to = "Sample", values_to = "tpm")
# gathered_genes <- tpm_d %>% gather(colnames(subset(tpm_d, select=-c(gene_id, gene_name))), key = "Sample", value = "tpm")
print(head(pivot_tpm))
pivot_tpm <- inner_join(sample_man, pivot_tpm, by= "Sample")
# gathered_genes$Treatment <- factor(gathered_genes$Chase_time_h, levels= c("0h", "3-5h"))
print(head(pivot_tpm))


p_timegrid <- ggplot(pivot_tpm, aes(x = Chase_time_h, 
                                         y = tpm, 
                                         group = interaction(gene_id, Sample_type, Treatment), 
                                         shape = Treatment, 
                                         color = Group)) +
                      geom_line(aes(linetype = Treatment), alpha=0.8) + 
                      geom_point() +
                      scale_colour_manual(values = cols) +
                      scale_linetype_manual(values = c(1:2)) +
                	  labs(x = "Chase time (h)", y = "Nascent fraction per gene") +
                      theme_bw()  + 
                      facet_wrap(~ gene_name, scales="free", ncol = 5) 

cat(paste0("Saving output:\n\t- ", outputPDF), sep="\n")

pdf(file=outputPDF, width=8, height=6)
                print(p_timegrid)
        invisible(dev.off())
cat("\n")
cat("DONE!", sep="\n")
