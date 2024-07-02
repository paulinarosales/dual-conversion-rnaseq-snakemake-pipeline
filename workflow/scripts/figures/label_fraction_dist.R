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
sample_manifestTSV <- snakemake@input[["sample_manifestTSV"]]
nascentFractTSV <- snakemake@input[["nascentFractTSV"]]
outputPDF <- snakemake@output[[1]]

# ---------- Main code ---------- #
cat("Reading input data...", sep="\n")
sample_man <- read.table(sample_manifestTSV, header = TRUE, sep = "\t") #%>% filter(c("Sample_type", "Treatment", "Chase_time_h", "Batch"))
sample_man$sample <- paste(sample_man$Sample_type, sample_man$Treatment, "Chase-time", sample_man$Chase_time_h, "Bio-rep", sample_man$Bio_rep, sep="_")

fract_t <- read.table(nascentFractTSV, header = TRUE, sep = "\t", check.names = FALSE)

fract_t <- fract_t %>% pivot_longer(cols = 9:ncol(fract_t), names_to = "sample", values_to = "fraction")

fract_t <- inner_join(fract_t, sample_man, by="sample")
fract_t$group <- paste(fract_t$Sample_type, fract_t$Chase_time_h, sep="_")
cat("\n")

cat("Ploting...", sep="\n")

rates_p <- ggplot(fract_t, aes(x=fraction, fill=group)) + 
                geom_histogram(aes(y=after_stat(density)), position = "dodge") +
                scale_fill_manual(values = pawlette_paired[1:length(unique(fract_t$group))]) +
				# stat_summary(aes(x = fraction, y=after_stat(density)), fun.x=mean, colour="red", geom="line", linetype = "dashed") +
                labs(y = "Density", x = "Nascent fraction per gene") +
                facet_grid(Treatment~Sample_type, scales="free", space="free") + 
                # ylim(c(0,1)) +
                theme_bw()


# rates_p <- ggplot(fract_t, aes(x=sample,y=fraction, fill=group)) + 
#                 geom_boxplot(trim=FALSE) + 
# 				stat_summary(fun.y=mean, geom="point", size=2, color="black") +
# 				stat_summary(fun.y=median, geom="point", size=2, color="red") +
# 				# stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red") +
#                 # scale_color_manual(values = pawlette[1:length(unique(fract_t$Sample_type))]) +
#                 scale_fill_manual(values = pawlette_paired[1:length(unique(fract_t$group))]) +
#                 # scale_alpha_manual(values = c(1, 0.1, 0.7)) +
#                 # facet_wrap(~Treatment, ncol=3, scales = "free") +
#                 facet_grid(~Treatment, scales="free", space="free") + 
#                 labs(x = "Sample", y = "Nascent fraction") +
# 				guides(x = "none") +
#                 # coord_cartesian(ylim=c(0, ymax))
#                 # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
#                 theme_bw()
                # theme(axis.ticks.x = element_blank(), ) 
cat("\n")
# cat("ANOVA test:", sep="\n")
# anovaTest <- aov(values ~ class, data = rates_t)
# print(TukeyHSD(x = anovaTest, 'class', conf.level = 0.95)$class)
# cat("\n")
cat(paste0("Saving output:\n\t- ", outputPDF), sep="\n")

pdf(file=outputPDF, width=8, height=6)
                print(rates_p)
        invisible(dev.off())
cat("\n")
cat("DONE!", sep="\n")



	# curPlot = ggplot(plotTab, aes(x=class,y=values,fill=highlight,col=highlight)) + 
        # stat_boxplot(geom ='errorbar') + 
        # geom_boxplot(outlier.shape = NA,lwd=0.8,fatten=2) + 
        # facet_grid(~group, scales="free", space="free") + 
        # xlab("") + ylab("Mutation rate per UTR base [%]") +
	# scale_fill_manual(values=c("white","white")) + 
        # scale_color_manual(values=c("black", "red")) + 
        # theme(axis.ticks.x = element_blank(), legend.position = "none") + 
        # coord_cartesian(ylim=c(0, ymax))

	# plotList[[length(plotList)+1]] <- curPlot + ggtitle(rates$sample[i])


