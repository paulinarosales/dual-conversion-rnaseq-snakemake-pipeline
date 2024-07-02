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
fraction_ctsTSV <- snakemake@input[["fraction_ctsTSV"]]
ref_halflifeCSV <- snakemake@input[["ref_halflifeCSV"]]
goisTSV <- snakemake@input[["goisTSV"]]

outPDF <- snakemake@output[[1]]

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
ref_hf <- read.table(ref_halflifeCSV, header = TRUE, sep = ",", check.names = FALSE)
ref_hf <- ref_hf %>% select(c("Gene", "Pluripotent_hf_hour"))
colnames(ref_hf) <-  c("gene_name", "hf_hour")
head(ref_hf)

gois_t <- read.table(goisTSV, header = TRUE, sep = "\t")
head(gois_t)

# Read nascent fraction (per gene)
decay_t <- read.table(fraction_ctsTSV, header = TRUE, sep = "\t", check.names = FALSE)
colnames(decay_t)
decay_t <- decay_t %>% select(c("gene_name", basename(dirname(outPDF))))
colnames(decay_t) <- c("gene_name", "decay")
head(decay_t)

hf_tab <- inner_join(decay_t, ref_hf, by="gene_name") 
tail(hf_tab)
cat("\n")


cat("Plotting...", sep="\n")
dist_p <- ggplot(hf_tab, aes(x=hf_hour, y=decay)) +
                geom_point(color=dplyr::case_when(hf_tab$gene_name %in% gois_t$gene_name ~ "#bc3c29",
                                                !(hf_tab$gene_name %in% gois_t$gene_name) ~ "gray"), 
                    alpha=dplyr::case_when(hf_tab$gene_name %in% gois_t$gene_name  ~ 1, 
                                            !(hf_tab$gene_name %in% gois_t$gene_name) ~ 0.3),
                    size=dplyr::case_when(hf_tab$gene_name %in% gois_t$gene_name  ~ 1.5, 
                                            !(hf_tab$gene_name %in% gois_t$gene_name) ~ 1)) +
            geom_text_repel(aes(label=dplyr::case_when(hf_tab$gene_name %in% gois_t$gene_name ~ hf_tab$gene_name)), size=3) +
                labs(title = basename(dirname(outPDF)), x = "Half-life (h)", y = "Decay") +
                theme_bw()


cat(paste0("Saving outputs:\n\t- ", outPDF), sep="\n")

pdf(file=outPDF, width=8, height=6)
                print(dist_p)
        invisible(dev.off())

cat("\n")
cat("DONE!", sep="\n")