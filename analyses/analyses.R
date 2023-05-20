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
BOOT <- function(X, Y, B, NAME) {
    f <- function(x) midpoint.root(nj(x))
    tr <- f(X)

    ## do bootstrapping
    bp <- boot.phylo(tr, as.matrix(X), f, B = B) / B
    bp[is.na(bp)] <- paste0("<", 1 / B)
    bp2 <- tibble(node = 1:Nnode(tr) + Ntip(tr), bootstrap = bp)
    tr <- full_join(tr, bp2, by = "node")

    ## rename labels and convert shape IDs to factors
    x <- full_join(as_tibble(tr), Y, by = c("label" = NAME))
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
            size = 2.5,
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
DATA.matrix <- as.data.frame(DATA[, 7:ncol(DATA)])

## Qual/Quant Data
Quant.full <- DATA.matrix[, 1:24]
rownames(Quant.full) <- DATA$ID2

Qual.full <- DATA.matrix[, 25:ncol(DATA.matrix)]
Qual.full$ID <- DATA$LocationID

Qual.full <- aggregate(. ~ ID, data = Qual.full, mean)
rownames(Qual.full) <- Qual.full$ID

DATA.loc <- distinct(DATA[, c(1, 3, 4, 5, 6)])

Qual.full <- Qual.full %>%
    select(-ID)
## subset D1 D2
DATA.d1d2 <- DATA %>%
    filter(Subclade %in% c("D1", "D2"))

DATA.matrix.d1d2 <- as.data.frame(DATA.d1d2[, 7:ncol(DATA)])
rownames(DATA.matrix.d1d2) <- DATA.d1d2$ID2

Quant.d1d2 <- DATA.matrix.d1d2[, 1:24]

Qual.d1d2 <- DATA.matrix.d1d2[, 25:ncol(DATA.matrix.d1d2)]
Qual.d1d2$ID <- DATA.d1d2$LocationID

Qual.d1d2 <- aggregate(. ~ ID, data = Qual.d1d2, mean)
rownames(Qual.d1d2) <- Qual.d1d2$ID

Qual.d1d2 <- Qual.d1d2 %>%
    select(-ID)

DATA.loc.d1d2 <- distinct(DATA[, 3:6]) %>%
    filter(Subclade %in% c("D1", "D2"))


################ Analysis of Full Dataset ###################

## Make Neighborjoining tree based on euclidean distance and root @ midpoint
Quant.full.dist <- dist(Quant.full)

## make NJ tree with bootstrapping
Quant.full.tree <- BOOT(Quant.full.dist, DATA, 100, "ID2")

## plot tree
TREE(
    X = Quant.full.tree,
    Y = DATA,
    Xlim = 65,
    Xlab = "Euclidean Distance",
    width = 10,
    height = 25
)

## calculate Hamming Distances
Qual.full.dist <- hammingdists(Qual.full)
Qual.full.dist <- dist_setNames(
    Qual.full.dist, rownames(Qual.full)
)

## make NJ tree with bootstrapping
Qual.full.tree <- BOOT(Qual.full.dist, DATA.loc, 100, "LocationID")

## plot tree
TREE(
    X = Qual.full.tree,
    Y = DATA.loc,
    Xlim = 10,
    Xlab = "Hamming Distance",
    width = 5,
    height = 12
)

############## D1/D2 only

## Make Neighborjoining tree based on euclidean distance and root @ midpoint
Quant.d1d2.dist <- dist(Quant.d1d2)
Quant.d1d2.tree <- BOOT(Quant.d1d2.dist, DATA.d1d2, 100, "ID2")

TREE(
    X = Quant.d1d2.tree,
    Y = DATA.d1d2,
    Xlim = 55,
    Xlab = "Euclidean Distance",
    width = 5,
    height = 7
)

## calculate Hamming Distances and make trees
Qual.d1d2.dist <- hammingdists(Qual.d1d2)
Qual.d1d2.dist <- dist_setNames(
    Qual.d1d2.dist, rownames(Qual.d1d2)
)
Qual.d1d2.tree <- BOOT(Qual.d1d2.dist, DATA.loc.d1d2, 100, "LocationID")

TREE(
    X = Qual.d1d2.tree,
    Y = DATA.loc.d1d2,
    Xlim = 10,
    Xlab = "Hamming Distance",
    width = 5,
    height = 3
)

############### Finally, plot Gerhards MP tree

SynTable <- read_excel("data/Montenegrina_rawdata.xlsx",
    sheet = "Tabelle2"
)
Qual.mp.tree <- read.nexus("data/MP_tree.tre")
Qual.mp.tree <- as_tibble(Qual.mp.tree)

Qual.mp.tree <- full_join(Qual.mp.tree, SynTable, by = c("label" = "MP"))
Qual.mp.tree$label <- Qual.mp.tree$LocationID
Qual.mp.tree$Species <- as.factor(Qual.mp.tree$Species)
DATA.loc2 <- as.data.frame(na.omit(Qual.mp.tree[, 5:8]))
Qual.mp.tree <- as.treedata(Qual.mp.tree)

## plot tree
TREE.mp <- function(X, Y, Xlim, Xlab, width, height) {
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

        ## add tip labels
        geom_tiplab(
            aes(
                label = label,
                angle = 0
            ),
            hjust = -0.1,
            size = 2.5,
            show.legend = FALSE
        ) +

        ## add colored symbols as tip-points according to clades
        geom_tippoint(aes(color = Clade, fill = Clade, shape = Species),
            size = 1,
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

TREE.mp(
    X = Qual.mp.tree,
    Y = DATA.loc2,
    Xlim = 15,
    Xlab = "Maximum Parsimony",
    width = 5,
    height = 10
)
