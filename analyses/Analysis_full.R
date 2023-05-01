library(tidyverse) ## for ggplot
library(cultevo) ## for Hamming distance
library(ape) ## for Phylogenetics in general
library(phangorn) ## for NJ tree
library(treeio) ## as.treedata
library(phytools) ## for midpoint rooting
library(usedist) ## to add labels
library(ggtree) ## for plotting "nice" trees
library(poppr)

setwd("D:/GitHub/Montenegrina_2023/")

### 100X Bootstrapping with the ape package
BOOT <- function(X) {
    f <- function(x) midpoint.root(nj(x))
    tr <- f(X)

    bp <- boot.phylo(tr, as.matrix(X), f) / 100

    bp2 <- tibble(node = 1:Nnode(tr) + Ntip(tr), bootstrap = bp)
    tr <- full_join(tr, bp2, by = "node")
    return(tr)
}

DATA <- read.table("data/mainmatrix copy with colours and symbols.dat",
    header = TRUE,
    sep = "\t",
    comment.char = ""
)

DATA$Num <- seq(1, nrow(DATA), 1)
DATA$ID2 <- paste(DATA$ID, DATA$Num, sep = "_")


Colorcode <- read.table("data/ColorCode.txt",
    header = TRUE,
    sep = "\t",
    comment.char = ""
)

ShapeData <- read.table("data/Montenegrina qualitative_data copy.dat",
    header = TRUE,
    sep = "\t",
    comment.char = ""
)

ShapeData$Num <- seq(1, nrow(ShapeData), 1)
ShapeData$ID2 <- paste(ShapeData$ID, DATA$Num, sep = "_")

Shapes <- ShapeData %>%
    select(ID2, Species)

DATA$ID2 <- Shapes$ID2

# build grouping by combination of variables
DATA <- DATA %>%
    select(-Color, -Shape) %>%
    left_join(Colorcode, by = c("Subclade" = "subclade")) %>%
    left_join(Shapes, by = c("ID2" = "ID2")) %>%
    mutate(ID2 = str_replace(ID2, "Montenegrina", "M."))


## subset D1 D2
DATA.d1d2 <- DATA %>%
    filter(Subclade %in% c("D1", "D2"))

DATA.matrix <- DATA[, 9:ncol(DATA) - 5]
rownames(DATA.matrix) <- DATA$ID2

DATA.matrix.d1d2 <- DATA.d1d2[, 9:ncol(DATA) - 5]
rownames(DATA.matrix.d1d2) <- DATA.d1d2$ID2

Quant.full <- DATA.matrix[, 1:24]
Quant.d1d2 <- DATA.matrix.d1d2[, 1:24]
Qual.full <- DATA.matrix[, 25:ncol(DATA.matrix)]
Qual.d1d2 <- DATA.matrix.d1d2[, 25:ncol(DATA.matrix)]

## Make Neighborjoining tree based on euclidean distance and root @ midpoint
Quant.full.dist <- dist(Quant.full)
Quant.full.tree <- BOOT(Quant.full.dist)

## rename labels and convert shape IDs to factors
x <- full_join(as_tibble(Quant.full.tree), DATA, by = c("label" = "ID2"))
x$Species <- as.factor(x$Species)
Quant.full.tree <- as.treedata(x, boot = bootstrap)


## calculate Hamming Distances
Qual.full.dist <- hammingdists(Qual.full)
Qual.full.dist <- dist_setNames(
    Qual.full.dist, DATA$ID2
)
Qual.full.tree <- BOOT(Qual.full.dist)

## rename labels and convert shape IDs to factors
x <- full_join(as_tibble(Qual.full.tree), DATA, by = c("label" = "ID2"))
x$Species <- as.factor(x$Species)
Qual.full.tree <- as.treedata(x, boot = bootstrap)

##

## Clade = Color
C.full <- list()
for (i in seq(1, nrow(DATA), 1)) {
    C.full[[DATA$Clade[i]]] <- DATA$Clade[i]
}

## Clade = Labels
L.full <- list()
for (i in seq(1, nrow(DATA), 1)) {
    L.full[[DATA$Clade[i]]] <- DATA$Subclade[i]
}

## Clade = Labels
S.full <- list()
for (i in seq(1, nrow(DATA), 1)) {
    S.full[[as.character(DATA$Species[i])]] <- DATA$Species[i]
}

## plot tree
PLOT.Qual.full.tree <- ggtree(Qual.full.tree, layout = "roundrect") +
    ## Set scales manually
    scale_color_manual(values = C.full, labels = L.full) +
    scale_fill_manual(values = C.full, labels = L.full) +
    scale_shape_manual(values = S.full) +

    ## add bootstrap
    geom_text(aes(label = bootstrap), hjust = 1.25, vjust = -0.75, size = 2, col = "blue") +


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
        size = 1,
        stroke = 1,
        show.legend = FALSE
    ) +
    theme_tree2() +
    theme_bw() +
    xlim(0, 11) +
    xlab("Hamming Distance") +
    ylab("Number of Individuals")


PLOT.Qual.full.tree

## save trees as PNG and PDF
ggsave("analyses/PLOT.Qual.full.tree.pdf",
    PLOT.Qual.full.tree,
    width = 10,
    height = 25
)

ggsave("analyses/PLOT.Qual.full.tree.png",
    PLOT.Qual.full.tree,
    width = 25,
    height = 25
)


## plot tree
PLOT.Quant.full.tree <- ggtree(Quant.full.tree, layout = "roundrect") +
    ## Set scales manually
    scale_color_manual(values = C.full, labels = L.full) +
    scale_fill_manual(values = C.full, labels = L.full) +
    scale_shape_manual(values = S.full) +

    ## add bootstrap
    geom_text(aes(label = bootstrap), hjust = 1.25, vjust = -0.75, size = 2, col = "blue") +

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
        size = 1,
        stroke = 1,
        show.legend = FALSE
    ) +
    theme_tree2() +
    theme_bw() +
    xlim(0, 80) +
    xlab("Euclidean Distance") +
    ylab("Number of Individuals")


## save trees as PNG and PDF
ggsave("analyses/PLOT.Quant.full.tree.pdf",
    PLOT.Quant.full.tree,
    width = 10,
    height = 25
)

ggsave("analyses/PLOT.Quant.full.tree.png",
    PLOT.Quant.full.tree,
    width = 10,
    height = 25
)


## Make Neighborjoining tree based on euclidean distance and root @ midpoint
Quant.d1d2.dist <- dist(Quant.d1d2)
Quant.d1d2.tree <- BOOT(Quant.d1d2.dist)

## rename labels and convert shape IDs to factors
x <- full_join(as_tibble(Quant.d1d2.tree), DATA.d1d2, by = c("label" = "ID2"))
x$Species <- as.factor(x$Species)
Quant.d1d2.tree <- as.treedata(x, boot = bootstrap)


## calculate Hamming Distances
Qual.d1d2.dist <- hammingdists(Qual.d1d2)
Qual.d1d2.dist <- dist_setNames(
    Qual.d1d2.dist, DATA.d1d2$ID2
)
Qual.d1d2.tree <- BOOT(Qual.d1d2.dist)

## rename labels and convert shape IDs to factors
x <- full_join(as_tibble(Qual.d1d2.tree), DATA.d1d2, by = c("label" = "ID2"))
x$Species <- as.factor(x$Species)
Qual.d1d2.tree <- as.treedata(x, boot = bootstrap)

##
## Clade = Color
C.d1d2 <- list()
for (i in seq(1, nrow(DATA.d1d2), 1)) {
    C.d1d2[[DATA.d1d2$Clade[i]]] <- DATA.d1d2$Clade[i]
}

## Clade = Labels
L.d1d2 <- list()
for (i in seq(1, nrow(DATA.d1d2), 1)) {
    L.d1d2[[DATA.d1d2$Clade[i]]] <- DATA.d1d2$Subclade[i]
}

## Clade = Labels
S.d1d2 <- list()
for (i in seq(1, nrow(DATA.d1d2), 1)) {
    S.d1d2[[as.character(DATA.d1d2$Species[i])]] <- DATA.d1d2$Species[i]
}

## plot tree
PLOT.Qual.d1d2.tree <- ggtree(Qual.d1d2.tree, layout = "roundrect") +
    ## Set scales manually
    scale_color_manual(values = C.d1d2, labels = L.d1d2) +
    scale_fill_manual(values = C.d1d2, labels = L.d1d2) +
    scale_shape_manual(values = S.d1d2) +

    ## add bootstrap
    geom_text(aes(label = bootstrap), hjust = 1.25, vjust = -0.75, size = 2, col = "blue") +

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
        size = 1,
        stroke = 1,
        show.legend = FALSE
    ) +
    theme_tree2() +
    theme_bw() +
    xlim(0, 11) +
    xlab("Hamming Distance") +
    ylab("Number of Individuals")


PLOT.Qual.d1d2.tree

## save trees as PNG and PDF
ggsave("analyses/PLOT.Qual.d1d2.tree.pdf",
    PLOT.Qual.d1d2.tree,
    width = 5,
    height = 7
)

ggsave("analyses/PLOT.Qual.d1d2.tree.png",
    PLOT.Qual.d1d2.tree,
    width = 5,
    height = 7
)


## plot tree
PLOT.Quant.d1d2.tree <- ggtree(Quant.d1d2.tree, layout = "roundrect") +
    ## Set scales manually
    scale_color_manual(values = C.d1d2, labels = L.d1d2) +
    scale_fill_manual(values = C.d1d2, labels = L.d1d2) +
    scale_shape_manual(values = S.d1d2) +

    ## add bootstrap
    geom_text(aes(label = bootstrap), hjust = 1.25, vjust = -0.75, size = 2, col = "blue") +

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
        size = 1,
        stroke = 1,
        show.legend = FALSE
    ) +
    theme_tree2() +
    theme_bw() +
    xlim(0, 80) +
    xlab("Euclidean Distance") +
    ylab("Number of Individuals")


## save trees as PNG and PDF
ggsave("analyses/PLOT.Quant.d1d2.tree.pdf",
    PLOT.Quant.d1d2.tree,
    width = 5,
    height = 7
)

ggsave("analyses/PLOT.Quant.d1d2.tree.png",
    PLOT.Quant.d1d2.tree,
    width = 5,
    height = 7
)

write.table(DATA,
    file = "data/RawData.txt",
    quote = F,
    row.names = F
)

aboot(as.matrix(Qual.d1d2), dist = "diss.dist")


ggtree(tr) +
    geom_text(aes(label = bootstrap), hjust = -.25, size = 3) +
    geom_tiplab(aes(label = label))
