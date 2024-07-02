log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

suppressPackageStartupMessages({
        library(readr)
        library(dplyr)
        library(tidyr)
	library(gridExtra)
        library(ggplot2)
})

# ---------- Snakemake parsing ---------- #
bakR_mutTSV <- snakemake@input[["bakR_mutTSV"]]
sample_manifestTSV <- snakemake@input[["sample_manifestTSV"]]

raw_distPDF <- snakemake@output[["raw_distPDF"]]
filter_distPDF <- snakemake@output[["filter_distPDF"]]

target_mut <- snakemake@params[["target_mut"]]

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


# ---------- Main code ---------- #
cat("Reading input data...", sep="\n")
sample_man <- read.table(sample_manifestTSV, header = TRUE, sep = "\t")
rownames(sample_man) <- paste(sample_man$Sample_type, sample_man$Treatment, "Chase-time", sample_man$Chase_time_h, "Bio-rep", sample_man$Bio_rep, sep="_")

pawlette_paired <- pawlette_paired[1:length(unique(sample_man$Group))]
names(pawlette_paired) <- rownames(unique(sample_man$Group))
sample_man$col <- sample_man$Group %in% pawlette_paired
print(sample_man$col)
cat("\n")

cat(paste("Extracting ", target_mut, "mutation counts per read..."), sep="\n")
mut_cts <- read.table(bakR_mutTSV, header = TRUE, sep = "\t")
mut_cts <- mut_cts %>% select(c("qname", target_mut))
print(head(mut_cts))
colnames(mut_cts) <- c("qname", "mut")
mut_cts$mut <- as.numeric(mut_cts$mut)
cat("\n")

print(max(mut_cts[,"mut"]))
print(mut_cts[mut_cts[,"mut"] == max(mut_cts[,"mut"]),])
mean_ct <- mean(mut_cts[,"mut"])
head(mut_cts)

cat("Plotting...", sep="\n")

raw_dist_p <- ggplot(mut_cts, aes(x=mut)) + 
                geom_histogram(aes(y=(after_stat(count)/sum(after_stat(count))*100))) +
                scale_fill_manual(values = pawlette[basename(dirname(bakR_mutTSV))]) +
                geom_vline(xintercept=mean_ct, linetype="dashed", color = "red") +
                # labs(x = "Mutation", y = "Global conversion rate [%]") +
                ylim(c(0,100)) +
                theme_bw()

mean_ct <- mean_ct %>% filter(mut>0)
mean_ct <- mean(mut_cts[,"mut"])

filter_dist_p <- ggplot(mut_cts, aes(x=mut)) + 
                geom_histogram(aes(y=(after_stat(count)/sum(after_stat(count))*100))) +
                scale_fill_manual(values = pawlette[basename(dirname(bakR_mutTSV))]) +
                geom_vline(xintercept=mean_ct, linetype="dashed", color = "red") +
                # labs(x = "Mutation", y = "Global conversion rate [%]") +
                ylim(c(0,100)) +
                theme_bw()

cat(paste0("Saving outputs:\n\t- ", raw_distPDF, "\n\t- ", filter_distPDF), sep="\n")

pdf(file=outPDF, width=8, height=6)
                print(dist_p)
        invisible(dev.off())

cat("\n")
cat("DONE!", sep="\n")