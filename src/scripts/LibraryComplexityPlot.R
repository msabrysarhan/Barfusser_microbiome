#!/usr/bin/env Rscript
#Author: Mohamed Sabry Sarhan
#Email: mohamed.sarhan@eurac.edu, m.sabrysarhan@gmail.com 
#Eurac Research, Institute for Mummy Studies, Bolzano, Italy

args = commandArgs(trailingOnly=TRUE)

data_tab = read.table(args[1])

x = data_tab$V1

y = data_tab$V2

pdf(file = args[2], width = 5, height = 5)

plot(x, y, lwd=2, pch=16, xlab="", ylab="")

fit3 <- lm(y~poly(x,3,raw=TRUE), data=data_tab)
lines(x, predict(fit3, data.frame(x=data_tab)), col="red", lwd=2)
title(main="Library complexity", xlab="Number of reads", ylab="Unique reads (%)") 

dev.off()
