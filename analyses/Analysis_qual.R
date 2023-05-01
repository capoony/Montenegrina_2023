library(tidyverse) ## for ggplot
library(cultevo) ## for Hamming distance
library(ape) ## for Phylogenetics in general
library(phangorn) ## for NJ tree
library(treeio) ## as.treedata
library(phytools) ## for midpoint rooting
library(usedist) ## to add labels
library(ggtree) ## for plotting "nice" trees

setwd("D:/GitHub/Montenegrina_2023/")

DATA <- read.table("data/D1_D2_FD_speciesnames copy.dat",
    header = TRUE,
    sep = "\t",
    comment.char = ""
)

# build grouping by combination of variables
DATA <- DATA %>%
    select(-X) %>%
    group_by(Species) %>%
    # add row number which works per group due to prior grouping
    dplyr::mutate(duplicateID = dplyr::row_number()) %>%
    # ungroup to prevent unexpected behaviour down stream
    dplyr::ungroup() %>%
    mutate(Sample = paste0(Species, "_", duplicateID))
## replace Color Column based on Colorcode from Lisi

DATA.matrix <- as.data.frame(DATA)[, 6:ncol(DATA) - 2]
rownames(DATA.matrix) <- DATA$Sample

DATA2 <- as.data.frame(DATA)
rownames(DATA2) <- DATA$Sample

DATA2$Clade <- DATA2$Color

DATA.dist <- dist(DATA.matrix)
tree <- nj(DATA.dist)

## Make Neighborjoining tree and root @ midpoint
my_nj <- midpoint.root(tree)

Tree.tibble <- as_tibble(my_nj)

## rename labels and convert shape IDs to factors
x <- full_join(as_tibble(my_nj), DATA2, by = c("label" = "Sample"))
x$Shape <- as.factor(x$Shape)
tree2 <- as.treedata(x)

##

## Clade = Color
C <- list()
for (i in seq(1, nrow(DATA2), 1)) {
    C[[DATA2$Color[i]]] <- DATA2$Color[i]
}

## Clade = Labels
L <- list()
for (i in seq(1, nrow(DATA2), 1)) {
    L[[DATA2$Color[i]]] <- DATA2$Species[i]
}

## Clade = Labels
S <- list()
for (i in seq(1, nrow(DATA2), 1)) {
    S[[as.character(DATA2$Shape[i])]] <- DATA2$Shape[i]
}



# Species <- c("rugilabris rugilabris", "pallida", "lambdaformis welterschultesi", "rugilabris irmengardis", "klemmi", "edmundi", "fuchsi", "janinensis crassilabris", "golikutensis", "gregoi")

## plot tree
PLOT.tree <- ggtree(tree2, layout = "roundrect") +
    ## Set scales manually
    scale_color_manual(values = C, labels = c("D1", "D2")) +
    scale_fill_manual(values = C, labels = c("D1", "D2")) +
    scale_shape_manual(values = S, labels = L) +


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
    geom_tippoint(aes(color = Clade, fill = Clade, shape = Shape),
        size = 1,
        stroke = 1,
        show.legend = FALSE
    ) +
    theme_tree2() +
    theme_bw() +
    xlim(0, 27) +
    xlab("Euclidean Distance") +
    ylab("Number of Individuals")


PLOT.tree

## save trees as PNG and PDF
ggsave("analyses/tree_D1D2.pdf",
    PLOT.tree,
    width = 5,
    height = 7
)

ggsave("analyses/tree_D1D2.png",
    PLOT.tree,
    width = 5,
    height = 7
)
