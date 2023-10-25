library(tidyverse)
library(pheatmap)
library(ade4)
library(gridExtra)


DATA<-read.table("/media/inter/mkapun/projects/Montenegrina_2023/data/distance.txt",
header=T,
na.string="NA")

Geo<-as.dist(xtabs(DATA[, 3] ~ DATA[, 2] + DATA[, 1]))
DNA<-as.dist(xtabs(DATA[, 4] ~ DATA[, 2] + DATA[, 1]))
Morph<-as.dist(xtabs(DATA[, 5] ~ DATA[, 2] + DATA[, 1]))

sink("/media/inter/mkapun/projects/Montenegrina_2023/analyses/distances/distance.stat")
print("Geo vs. DNA\n")
mantel.rtest(Geo,DNA,nrepet=10000)
print("DNA vs. Morph\n")
mantel.rtest(DNA,Morph,nrepet=10000)
print("Geo vs. Morph\n")
mantel.rtest(Geo,Morph,nrepet=10000)
sink()

pdf("/media/inter/mkapun/projects/Montenegrina_2023/analyses/distances/distance_DNA.pdf")
pheatmap(as.matrix(DNA))
dev.off()
pdf("/media/inter/mkapun/projects/Montenegrina_2023/analyses/distances/distance_Geo.pdf")
pheatmap(as.matrix(Geo))
dev.off()
pdf("/media/inter/mkapun/projects/Montenegrina_2023/analyses/distances/distance_Morph.pdf")
pheatmap(as.matrix(Morph))
dev.off()