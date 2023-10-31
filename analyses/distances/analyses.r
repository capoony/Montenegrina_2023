library(tidyverse)
library(pheatmap)
library(ecodist)

### read Data
DATA<-read.table("/media/inter/mkapun/projects/Montenegrina_2023/data/distance.txt",
header=T,
na.string="NA")

## convert to distance matrices
Geo<-as.dist(xtabs(DATA[, 3] ~ DATA[, 2] + DATA[, 1]))
DNA<-as.dist(xtabs(DATA[, 4] ~ DATA[, 2] + DATA[, 1]))
Morph<-as.dist(xtabs(DATA[, 5] ~ DATA[, 2] + DATA[, 1]))


## do Mantel tests in all combinations 
sink("/media/inter/mkapun/projects/Montenegrina_2023/analyses/distances/distance.stat")
print("DNA vs. Geo")
mantel(DNA ~ Geo,nperm=1000000)
print("Morph vs. DNA")
mantel(Morph ~ DNA,nperm=1000000)
print("Morph vs. Geo")
mantel(Morph ~Geo,nperm=1000000)
sink()

## plot heatmaps
pdf("/media/inter/mkapun/projects/Montenegrina_2023/analyses/distances/distance_DNA.pdf")
pheatmap(as.matrix(DNA))
dev.off()
pdf("/media/inter/mkapun/projects/Montenegrina_2023/analyses/distances/distance_Geo.pdf")
pheatmap(as.matrix(Geo))
dev.off()
pdf("/media/inter/mkapun/projects/Montenegrina_2023/analyses/distances/distance_Morph.pdf")
pheatmap(as.matrix(Morph))
dev.off()


## plot pairwise regression plots
reg <- function(x, y, col) abline(lm(y~x), col=col) 

panel.lm =  function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
    cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...)  {
    points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) reg(x[ok], y[ok], col.smooth)
}

panel.cor <- function(x, y, digits = 2, prefix = "R²: ", cex.cor, ...)
{
 usr <- par("usr"); on.exit(par(usr))
 par(usr = c(0, 1, 0, 1))
 r <- 100*abs(cor(x, y))**2
 txt <- format(c(r, 0.123456789), digits = digits)[1]
 txt <- paste0(prefix, txt,"%")
 text(0.5, 0.5, txt, cex = 1.1, font = 4)
 }

pdf("/media/inter/mkapun/projects/Montenegrina_2023/analyses/distances/PairwiseRegression.pdf")

pairs(na.omit(DATA[,c(5,4,3)]), panel = panel.lm,
    cex = 1.5, pch = 19, col = adjustcolor(4, .4), cex.labels = 2, 
    font.labels = 2, lower.panel = panel.cor)

dev.off()

png("/media/inter/mkapun/projects/Montenegrina_2023/analyses/distances/PairwiseRegression.png",
    width=2000,
    height=2000,res=300)

pairs(na.omit(DATA[,c(5,4,3)]), panel = panel.lm,
    cex = 1.5, pch = 19, col = adjustcolor(4, .4), cex.labels = 2, 
    font.labels = 2, lower.panel = panel.cor)

dev.off()