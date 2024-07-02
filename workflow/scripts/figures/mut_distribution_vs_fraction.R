log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(readr)
        library(dplyr)
        library(tidyr)
	library(gridExtra)
        library(ggplot2)
        library(ggrepel)
})

# ---------- Snakemake parsing ---------- #
gene_mutTSV <- snakemake@input[["gene_mutTSV"]]
fraction_ctsTSV <- snakemake@input[["fraction_ctsTSV"]]
sample_manifestTSV <- snakemake@input[["sample_manifestTSV"]]
genesetTSV <- snakemake@input[["genesetTSV"]]

outPDF <- snakemake@output[[1]]

target_mut <- snakemake@params[["target_mut"]]

# ---------- Color palette ---------- #
pawlette <- c("#e18727", # orange
		"#0072b5", # dark blue
		"#20854e", # dark green
		"#bc3c29", # red
		"#7876b1", # purple
		"#cd6090", # pink
		"#ffdb24", # yellow
		"#6ca6cd", # light blue
		"#cd6090", # pink
		"#698b22", # light green
		"#8b5a2b") # brown


# ---------- Main code ---------- #
cat("Reading input data...", sep="\n")
geneset <- read.table(genesetTSV, header = TRUE, sep = "\t")
geneset <- geneset %>% select(c("gene_id", "gene_name"))

sample_man <- read.table(sample_manifestTSV, header = TRUE, sep = "\t")
rownames(sample_man) <- paste(sample_man$Sample_type, sample_man$Treatment, "Chase-time", sample_man$Chase_time_h, "Bio-rep", sample_man$Bio_rep, sep="_")

pawlette <- pawlette[1:nrow(sample_man)]
names(pawlette) <- rownames(sample_man)
print(pawlette)

# Read nascent transxript fraction (per gene)
fract_t <- read.table(fraction_ctsTSV, header = TRUE, sep = "\t", check.names = FALSE)
colnames(fract_t)
fract_t <- fract_t %>% select(c("gene_id", basename(dirname(gene_mutTSV))))
colnames(fract_t) <- c("gene_id", "fraction")
head(fract_t)

# Read mutation count table per gene
mut_t <- read.table(gene_mutTSV, header = TRUE, sep = "\t")
mut_t <- mut_t %>% select(c("GF", target_mut))
colnames(mut_t) <- c("gene_id", "mut")
mut_t <- inner_join(geneset, mut_t, by="gene_id")
mut_t$gene_id <- NULL
head(mut_t)

gene_tab <- inner_join(fract_t, mut_t, by="gene_id") 
head(gene_tab)
tail(gene_tab)
cat("\n")

gene_tab <- gene_tab[order(gene_tab$mut, decreasing = T),] 

top_mut <- gene_tab$gene_id[1:15]

cat("Plotting...", sep="\n")
dist_p <- ggplot(gene_tab, aes(x=mut, y=fraction)) +
                labs(x = "Mutation count", y = "Nascent fraction per gene") +
                geom_point(color=dplyr::case_when(gene_tab$gene_id %in% top_mut ~ "#bc3c29"
                                                !(gene_tab$gene_id %in% top_mut) ~ "gray"), 
                    alpha=dplyr::case_when(gene_tab$gene_id %in% top_mut)  ~ 0.8, 
                                            !(gene_tab$gene_id %in% top_muts) ~ 0.2),
                    size=dplyr::case_when(gene_tab$gene_id %in% top_mut  ~ 1.5, 
                                            !(gene_tab$gene_id %in% top_mut ~ 1)) +
            geom_text_repel(aes(label=dplyr::case_when(gene_tab$gene_id %in% top_mut ~ gene_tab$gene_name)), size=3) +
                theme_bw()


cat(paste0("Saving outputs:\n\t- ", outPDF), sep="\n")

pdf(file=outPDF, width=8, height=6)
                print(dist_p)
        invisible(dev.off())

cat("\n")
cat("DONE!", sep="\n")