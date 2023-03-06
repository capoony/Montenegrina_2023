library(tidyverse) ## for ggplot
library(cultevo) ## for Hamming distance
library(ape) ## for Phylogenetics in general
library(phangorn) ## for NJ tree
library(treeio) ## as.treedata
library(phytools) ## for midpoint rooting
library(usedist) ## to add labels
library(ggtree) ## for plotting "nice" trees


## read Data
DATA <- read.table("H:/Work/Projects/HaringL/Montenegrina_2023/data/Montenegrina qualitative_data copy.dat",
    header = TRUE,
    sep = "\t",
    comment.char = ""
)

Colorcode <- read.table("H:/Work/Projects/HaringL/Montenegrina_2023/data/ColorCode.txt",
    header = TRUE,
    sep = "\t",
    comment.char = ""
)

Colorcode$color[Colorcode$subclade %in% DATA$Subclade]

## replace Color Column based on Colorcode from Lisi
DATA <- DATA %>%
    select(-Color) %>%
    mutate(ID = str_replace(ID, "Montenegrina", "M.")) %>%
    left_join(Colorcode, by = c("Subclade" = "subclade"))

## focus on columns with Qualitative data and convert to factors
DATA.matrix <- as.data.frame(DATA[, 7:ncol(DATA) - 1], stringsAsFactors = TRUE)
DATA$Num <- seq(1, nrow(DATA), 1)
DATA$ID2 <- paste(DATA$ID, DATA$Num, sep = "_")

DATA$Shape <- DATA$PCH
DATA$Clade <- DATA$Color
## add rownames with unique IDs
rownames(DATA.matrix) <- paste(DATA$ID, DATA$Num, sep = "_")

## calculate Hamming Distances
Hamming <- hammingdists(DATA.matrix)
Hamming <- dist_setNames(
    Hamming, DATA$ID2
)

## Make Neighborjoining tree and root @ midpoint
my_nj <- midpoint.root(nj(Hamming))

## rename labels and convert shape IDs to factors
x <- full_join(as_tibble(my_nj), DATA, by = c("label" = "ID2"))
x$Shape <- as.factor(x$Shape)
tree2 <- as.treedata(x)

## make lists with aesthetics for plotting colors and shapes

## Clade = Color
G <- list()
for (i in seq(1, nrow(DATA), 1)) {
    G[[DATA$Clade[i]]] <- DATA$Clade[i]
}

## Clade = Labels
L <- list()
for (i in seq(1, nrow(DATA), 1)) {
    L[[DATA$Clade[i]]] <- DATA$Subclade[i]
}

## Shape symbols
S <- list()
for (i in seq(1, nrow(DATA), 1)) {
    S[[as.character(DATA$Shape[i])]] <- DATA$Shape[i]
}

## plot tree
PLOT.tree <- ggtree(tree2, layout = "roundrect") +
    ## Set scales manually
    scale_color_manual(values = G, labels = L) +
    scale_fill_manual(values = G, labels = L) +
    scale_shape_manual(values = S) +

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
        show.legend = TRUE
    ) +
    theme_tree2() +
    theme_bw() +
    xlim(0, 10) +
    guides(
        shape = guide_legend(override.aes = list(fill = "black"))
    )
# PLOT.tree

## save trees as PNG and PDF
ggsave("H:/Work/Projects/HaringL/Montenegrina_2023/analyses/tree.pdf",
    PLOT.tree,
    width = 10,
    height = 20
)

ggsave("H:/Work/Projects/HaringL/Montenegrina_2023/analyses/tree.png",
    PLOT.tree,
    width = 10,
    height = 20
)

## plot tree
PLOT.tree <- ggtree(tree2, layout = "circular") +
    ## Set scales manually
    scale_color_manual(values = G, labels = L) +
    scale_fill_manual(values = G, labels = L) +
    scale_shape_manual(values = S) +

    ## add tip labels
    geom_tiplab2() +

    ## add colored symbols as tip-points according to clades
    geom_tippoint(aes(color = Clade, fill = Clade, shape = Shape),
        size = 2,
        stroke = 1,
        show.legend = TRUE
    ) +
    theme_tree2() +
    theme_bw() +
    guides(
        shape = guide_legend(override.aes = list(fill = "black"))
    )
PLOT.tree

## save trees as PNG and PDF
ggsave("H:/Work/Projects/HaringL/Montenegrina_2023/analyses/tree_circular.pdf",
    PLOT.tree,
    width = 25,
    height = 25
)

ggsave("H:/Work/Projects/HaringL/Montenegrina_2023/analyses/tree_circular.png",
    PLOT.tree,
    width = 25,
    height = 25
)
