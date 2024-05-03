#!/usr/bin/env Rscript

# The script is used to visualize the flu segment depth per base file
# The depth file is produced by bedtool genomecov with -d option
# ARGUMENTS: 1 -> genome depth file in bed format, which have three columns
# 2 -> output file prefix

#import libraries
library(ggplot2)
#library(showtext)
# Set up variable to control command line arguments
args <- commandArgs(TRUE)
depth_file <- args[1]
prefix <- args[2]

C20 <- c(
    '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'
)
test <- read.table(file=depth_file)

#otherwise the text is not showing in the png file
#showtext::showtext_auto()
#showtext::showtext_opts(dpi = 300)

#test <- read.table(file='/nfs/Genomics_DEV/projects/xdong/deve/nf-fluAB/test/illumina/results/231222_S_I_390-S3/mapping/231222_S_I_390-S3-bwa.perbase.bed', sep = '\t', header = FALSE, fill = TRUE)
gg <-  ggplot(test, aes(x = V2, y = V3, color=V1, group=V1)) 
gg <- gg + ggtitle("Flu segment position vs depth plot") 
gg <- gg + theme(plot.title = element_text(hjust=0.5))
gg <- gg + xlab("Segment Position") 
gg <- gg + ylab("Sequence depth") 
gg <- gg + guides(color = guide_legend(title = "Segments"))
gg <- gg + scale_x_continuous(limits=c(0, 2500), breaks=seq(0,2500,200))
#gg <- gg + scale_y_continuous(limits=c(0, max(test$V3)), breaks=seq(0,max(test$V3),500))


gg <- gg + geom_line(linetype = "solid", color="grey1", linewidth=0.2)
gg <- gg + geom_point(size=1)
gg <- gg + scale_color_manual(values=C20)
#ggsave(paste(prefix, "_plot.jpeg", sep = ""), plot = gg,  width = 30, height = 12, units = "cm",dpi = 300)
ggsave(paste(prefix, "_plot.pdf", sep = ""), plot = gg,  dpi = 300, width=11, height=8.5)
