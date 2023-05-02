library(tidyverse) ## for ggplot
library(cultevo) ## for Hamming distance
library(ape) ## for Phylogenetics in general
library(phangorn) ## for NJ tree
library(treeio) ## as.treedata
library(phytools) ## for midpoint rooting
library(usedist) ## to add labels
library(ggtree) ## for plotting "nice" trees
library(readxl) ## read EXCEL input

## set working directory
setwd("D:/GitHub/Montenegrina_2023/")

################## Functions ###################

## 100X Bootstrapping with the ape package
BOOT <- function(X, Y, B) {
    f <- function(x) midpoint.root(nj(x))
    tr <- f(X)

    ## do bootstrapping
    bp <- boot.phylo(tr, as.matrix(X), f, B = B) / B
    bp[is.na(bp)] <- paste0("<", 1 / B)
    bp2 <- tibble(node = 1:Nnode(tr) + Ntip(tr), bootstrap = bp)
    tr <- full_join(tr, bp2, by = "node")

    ## rename labels and convert shape IDs to factors
    x <- full_join(as_tibble(tr), Y, by = c("label" = "ID2"))
    x$Species <- as.factor(x$Species)
    tr <- as.treedata(x, boot = bootstrap)

    return(tr)
}

## create codes for plotting
LISTs <- function(X) {
    ## Clade = Color
    C <- list()
    for (i in seq(1, nrow(X), 1)) {
        C[[X$Clade[i]]] <- X$Clade[i]
    }

    ## Clade = Labels
    L <- list()
    for (i in seq(1, nrow(X), 1)) {
        L[[X$Clade[i]]] <- X$Subclade[i]
    }

    ## Clade = Labels
    S <- list()
    for (i in seq(1, nrow(X), 1)) {
        S[[as.character(X$Species[i])]] <- X$Species[i]
    }
    return(list(C, L, S))
}

## plot tree
TREE <- function(X, Y, Xlim, Xlab, width, height) {
    ## get codes for plotting
    LIST <- LISTs(Y)
    Cf <- LIST[[1]]
    Lf <- LIST[[2]]
    Sf <- LIST[[3]]

    ## plot tree
    PLOT <- ggtree(X, layout = "roundrect") +
        ## Set scales manually
        scale_color_manual(values = Cf, labels = Lf) +
        scale_fill_manual(values = Cf, labels = Lf) +
        scale_shape_manual(values = Sf) +

        ## add bootstrap
        geom_text(aes(label = bootstrap),
            hjust = 1.25,
            vjust = -0.75,
            size = 2,
            col = "blue"
        ) +

        ## add tip labels
        geom_tiplab(
            aes(
                label = label,
                angle = 0
            ),
            hjust = -0.1,
            size = 2,
            show.legend = FALSE
        ) +

        ## add colored symbols as tip-points according to clades
        geom_tippoint(aes(color = Clade, fill = Clade, shape = Species),
            size = 2,
            stroke = 0.75,
            show.legend = FALSE
        ) +
        theme_tree2() +
        theme_bw() +
        xlim(0, Xlim) +
        xlab(Xlab) +
        ylab("Number of Individuals")

    ## save trees as PNG and PDF
    PDF <- paste0("analyses/PLOT.", deparse(substitute(X)), ".pdf")
    ggsave(PDF,
        PLOT,
        width = width,
        height = height
    )
    PNG <- paste0("analyses/PLOT.", deparse(substitute(X)), ".png")
    ggsave(PNG,
        PLOT,
        width = width,
        height = height
    )
}


####################### Code ###################

## read dataset
DATA <- read_excel("data/Montenegrina_rawdata.xlsx")

## rename
DATA <- DATA %>%
    rename(
        Clade = CladeColor,
        ID2 = ID_short,
        Species = SpeciesShape
    )

## restrict to Measurment columns
DATA.matrix <- as.data.frame(DATA[, 6:ncol(DATA)])
rownames(DATA.matrix) <- DATA$ID2

## Qual/Quant Data
Quant.full <- DATA.matrix[, 1:24]
Qual.full <- DATA.matrix[, 25:ncol(DATA.matrix)]

## subset D1 D2
DATA.d1d2 <- DATA %>%
    filter(Subclade %in% c("D1", "D2"))

DATA.matrix.d1d2 <- as.data.frame(DATA.d1d2[, 6:ncol(DATA)])
rownames(DATA.matrix.d1d2) <- DATA.d1d2$ID2

Quant.d1d2 <- DATA.matrix.d1d2[, 1:24]
Qual.d1d2 <- DATA.matrix.d1d2[, 25:ncol(DATA.matrix)]

################ Analysis of Full Dataset ###################

## Make Neighborjoining tree based on euclidean distance and root @ midpoint
Quant.full.dist <- dist(Quant.full)

## make NJ tree with bootstrapping
Quant.full.tree <- BOOT(Quant.full.dist, DATA, 100)

## plot tree
TREE(
    X = Quant.full.tree,
    Y = DATA,
    Xlim = 80,
    Xlab = "Euclidean Distance",
    width = 10,
    height = 25
)

## calculate Hamming Distances
Qual.full.dist <- hammingdists(Qual.full)
Qual.full.dist <- dist_setNames(
    Qual.full.dist, DATA$ID2
)

## make NJ tree with bootstrapping
Qual.full.tree <- BOOT(Qual.full.dist, DATA, 100)

## plot tree
TREE(
    X = Qual.full.tree,
    Y = DATA,
    Xlim = 8.5,
    Xlab = "Hamming Distance",
    width = 10,
    height = 25
)

############## D1/D2 only

## Make Neighborjoining tree based on euclidean distance and root @ midpoint
Quant.d1d2.dist <- dist(Quant.d1d2)
Quant.d1d2.tree <- BOOT(Quant.d1d2.dist, DATA.d1d2, 100)

TREE(
    X = Quant.d1d2.tree,
    Y = DATA.d1d2,
    Xlim = 45,
    Xlab = "Euclidean Distance",
    width = 5,
    height = 7
)

## calculate Hamming Distances and make trees
Qual.d1d2.dist <- hammingdists(Qual.d1d2)
Qual.d1d2.dist <- dist_setNames(
    Qual.d1d2.dist, DATA.d1d2$ID2
)
Qual.d1d2.tree <- BOOT(Qual.d1d2.dist, DATA.d1d2, 100)

TREE(
    X = Qual.d1d2.tree,
    Y = DATA.d1d2,
    Xlim = 8.5,
    Xlab = "Hamming Distance",
    width = 5,
    height = 7
)
