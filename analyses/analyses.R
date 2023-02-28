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

## focus on columns with Qualitative data and convert to factors
DATA.matrix <- as.data.frame(DATA[, 7:ncol(DATA) - 1], stringsAsFactors = TRUE)
DATA$Num <- seq(1, nrow(DATA), 1)
DATA$ID2 <- paste(DATA$ID, DATA$Num, sep = "_")

## add rownames with unique IDs
rownames(DATA.matrix) <- paste(DATA$ID, DATA$Num, sep = "_")

## calculate Hamming Distances
Hamming <- hammingdists(DATA.matrix)
Hamming <- dist_setNames(
    Hamming, DATA$ID2
)
my_nj <- midpoint.root(nj(Hamming))

rownames(DATA) <- DATA$ID2
x <- full_join(as_tibble(my_nj), DATA, by = c("label" = "ID2"))
tree2 <- as.treedata(x)

G <- list()
for (i in seq(1, nrow(DATA), 1)) {
    G[[DATA$Color[i]]] <- DATA$Color[i]
}

PLOT.tree <- ggtree(tree2, layout = "roundrect") +
    scale_color_manual(values = G) +
    # geom_tiplab(aes(color = Color,
    #     	label=label,
    #         angle=0),
    #     show.legend = FALSE
    # ) +
    geom_tiplab(
        aes(
            label = label,
            angle = 0
        ),
        hjust = -0.1,
        show.legend = FALSE
    ) +
    geom_tippoint(
        pch = DATA$PCH,
        size = 2,
        stroke = 1,
        fill = "black"
    ) +
    theme_tree2() +
    theme_bw() +
    xlim(0, 20)
PLOT.tree

ggsave("H:/Work/Projects/HaringL/Montenegrina_2023/analyses/black_tree.pdf",
    PLOT.tree,
    width = 10,
    height = 40
)

ggsave("H:/Work/Projects/HaringL/Montenegrina_2023/analyses/black_tree.png",
    PLOT.tree,
    width = 10,
    height = 40
)
