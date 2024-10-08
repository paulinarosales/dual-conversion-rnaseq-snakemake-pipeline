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
convFiles <- snakemake@input[["convFiles"]]
countFiles <- snakemake@input[["countFiles"]]
sample_manifestTSV <- snakemake@input[["sample_manifestTSV"]]

all_globalRatesTSV <- snakemake@output[["all_globalRatesTSV"]]
barplot_newPDF <- snakemake@output[["barplot_newPDF"]]
barplot_oldPDF <- snakemake@output[["barplot_oldPDF"]]
# barplot_norm_newPDF <- snakemake@output[["barplot_norm_newPDF"]]
# barplot_norm_oldPDF <- snakemake@output[["barplot_norm_oldPDF"]]

y_lim <- 2.2
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
sample_man <- read.table(sample_manifestTSV, header = TRUE, sep = "\t")
rownames(sample_man) <- paste(sample_man$Sample_type, sample_man$Treatment, "Chase_time", sample_man$Chase_time_h, "Bio-rep", sample_man$Bio_rep, sep="_")
sample_man$convFiles <- convFiles
sample_man$countFiles <- countFiles
cat("\n")

cat("Merging conversion rates tables for all samples...", sep="\n")
rates_t <- data.frame("T_A" = numeric(nrow(sample_man)), "C_A" = numeric(nrow(sample_man)), "G_A" = numeric(nrow(sample_man)),                                                 
                      "A_T" = numeric(nrow(sample_man)), "C_T" = numeric(nrow(sample_man)), "G_T" = numeric(nrow(sample_man)), 
                      "A_C" = numeric(nrow(sample_man)), "T_C" = numeric(nrow(sample_man)), "G_C" = numeric(nrow(sample_man)), 
                      "A_G" = numeric(nrow(sample_man)), "T_G" = numeric(nrow(sample_man)), "C_G" = numeric(nrow(sample_man)),
                      "Lib_size" = numeric(nrow(sample_man)),
                      "Sample_type" = sample_man$Sample_type, 
                      "Treatment" = sample_man$Treatment,
                      "Chase_time_h" = sample_man$Chase_time_h,
                      "Batch" = sample_man$Batch,
                      "Sample" = rownames(sample_man))



for(i in 1:nrow(sample_man)){
        conv_rates <- read.table(sample_man$convFiles[i], header = TRUE, sep = "\t")
        cts <- read.table(sample_man$countFiles[i], header = TRUE, sep = "\t", comment.char = "#")

        rates_t[i, colnames(conv_rates)] <- conv_rates
        rates_t[i, "Lib_size"] <- sum(cts[,ncol(cts)])
}

rates_t_new <- rates_t %>% filter(Batch == 2023)
head(rates_t_new)
rates_t_old <- rates_t %>% filter(Batch == 2022)
head(rates_t_old)

# pivot_rates_t$mut <- sub("_", "-to-", pivot_rates_t$mut) # chage mutation label from N_N --> N-to-N for clarity

rates_t <- rates_t %>% select(c("Sample", 1:12)) # select only rate to save plot

cat("\n")


cat("Plotting for new batch...", sep="\n")
pivot_rates_new <- rates_t_new %>% pivot_longer(cols = 1:12, names_to = "mut", values_to = "rate")

# pivot_rates_new$norm_rate <- (pivot_rates_new$rate/pivot_rates_new$Lib_size)*1000000
# pivot_rates_new$norm_rate <- round(pivot_rates_new$norm_rate, 2)

pivot_rates_new$rate <- round(pivot_rates_new$rate, 2)

rates_p_new <- ggplot(pivot_rates_new, aes(x=mut, y=rate, fill=Sample_type)) + 
                geom_bar(position="dodge", stat="identity") +
                geom_text(aes(label = rate, group=Sample_type), 
                angle = 90, position = position_dodge(width = 1), 
                hjust = -0.2, vjust = 0.5, size = 2) +
                scale_fill_manual(values = pawlette[1:length(unique(pivot_rates_new$Sample_type))]) +
                labs(x = "Mutation", y = "Global conversion rate [%]") +
                facet_grid(Treatment~Chase_time_h) +
                ylim(c(0,y_lim)) +
                theme_bw(base_size = 10)

pdf(file=barplot_newPDF, width=8, height=6)
                print(rates_p_new)
        invisible(dev.off())

# rates_p_new_norm <- ggplot(pivot_rates_new, aes(x=mut, y=norm_rate, fill=Sample_type)) + 
#                 geom_bar(position="dodge", stat="identity") +
#                 geom_text(aes(label = norm_rate, group=Sample_type), 
#                 angle = 90, position = position_dodge(width = 1), 
#                 hjust = -0.2, vjust = 0.5, size = 2) +
#                 scale_fill_manual(values = pawlette[1:length(unique(pivot_rates_new$Sample_type))]) +
#                 labs(x = "Mutation", y = "Global conversion rate [%]") +
#                 facet_grid(Chase_time_h~Treatment) +
#                 ylim(c(0,1)) +
#                 theme_bw()

# pdf(file=barplot_norm_newPDF, width=8, height=6)
#                 print(rates_p_new_norm)
#         invisible(dev.off())

cat("\n")

cat("Plotting for old batch...", sep="\n")
pivot_rates_old <- rates_t_old %>% pivot_longer(cols = 1:12, names_to = "mut", values_to = "rate")

# pivot_rates_old$norm_rate <- (pivot_rates_old$rate/pivot_rates_old$Lib_size)*1000000
# pivot_rates_old$norm_rate <- round(pivot_rates_old$norm_rate, 2)

pivot_rates_old$rate <- round(pivot_rates_old$rate, 2)

rates_p_old <- ggplot(pivot_rates_old, aes(x=mut, y=rate, fill=Sample_type)) + 
                geom_bar(position="dodge", stat="identity") +
                geom_text(aes(label = rate, group=Sample_type), 
                angle = 90, position = position_dodge(width = 1), 
                hjust = -0.2, vjust = 0.5, size = 3) +
                scale_fill_manual(values = pawlette[1:length(unique(pivot_rates_old$Sample_type))]) +
                labs(x = "Mutation", y = "Global conversion rate [%]") +
                # facet_grid(Chase_time_h~Treatment) +
                ylim(c(0,y_lim)) +
                theme_bw() 

pdf(file=barplot_oldPDF, width=8, height=6)
                print(rates_p_old)
        invisible(dev.off())


# rates_p_old_norm <- ggplot(pivot_rates_old, aes(x=mut, y=norm_rate, fill=Sample_type)) + 
#                 geom_bar(position="dodge", stat="identity") +
#                 geom_text(aes(label = norm_rate, group=Sample_type), 
#                 angle = 90, position = position_dodge(width = 1), 
#                 hjust = -0.2, vjust = 0.5, size = 2) +
#                 scale_fill_manual(values = pawlette[1:length(unique(pivot_rates_old$Sample_type))]) +
#                 labs(x = "Mutation", y = "Global conversion rate [%]") +
#                 # facet_grid(Chase_time_h~Treatment) +
#                 # ylim(c(0,2)) +
#                 theme_bw()

# pdf(file=barplot_norm_oldPDF, width=8, height=6)
#                 print(rates_p_old_norm)
#         invisible(dev.off())

cat("\n")

cat(paste0("Saving outputs:\n\t- ", all_globalRatesTSV, "\n\t-", barplot_newPDF, "\n\t-", barplot_oldPDF), sep="\n")

write.table(rates_t, file = all_globalRatesTSV, row.names= FALSE, quote = FALSE, sep = "\t")


cat("\n")
cat("DONE!", sep="\n")