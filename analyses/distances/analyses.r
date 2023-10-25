library(tidyverse)
library(pheatmap)
library(ecodist, help, pos = 2, lib.loc = NULL)

### read Data
DATA<-read.table("/media/inter/mkapun/projects/Montenegrina_2023/data/distance.tsv",
header=T,
na.string="NA")

## 
Geo<-as.dist(xtabs(DATA[, 3] ~ DATA[, 2] + DATA[, 1]))
DNA<-as.dist(xtabs(DATA[, 4] ~ DATA[, 2] + DATA[, 1]))
Morph<-as.dist(xtabs(DATA[, 5] ~ DATA[, 2] + DATA[, 1]))

sink("/media/inter/mkapun/projects/Montenegrina_2023/analyses/distances/distance.stat")
print("DNA vs. Geo")
mantel(DNA ~ Geo,nperm=1000000)
print("Morph vs. DNA")
mantel(Morph ~ DNA,nperm=1000000)
print("Morph vs. Geo")
mantel(Morph ~Geo,nperm=1000000)
#mantel(Morph ~Geo+DNA,nperm=10000)
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