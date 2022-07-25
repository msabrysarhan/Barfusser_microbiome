#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

data_tab = read.table(args, sep = ";", header = FALSE)

#data_tab = read.csv("C:/Users/Sabry/OneDrive - Scientific Network South Tyrol/Desktop/Body_water.seq_stats.txt")

pdf(file = paste(args, "pdf", sep = ".", collapse = ""), width = 7, height = 7)

plot.new()
for (i in 1:nrow(data_tab)){
  title(sub=data_tab[i, 1], line = -(nrow(data_tab)-i), adj = 0)
}

dev.off()
