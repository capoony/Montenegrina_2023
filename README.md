### Montenegrina NJ Tree based on qualitative morphological characterization

Here, I constructed a distance matrix based on hamming distances across all 12 characters using the _cultevo_ package in _R_. After that, I used the Neighbor-Joining Method as implemented in the package _phangorn_ to reconstruct a tree and applied midpoint rooting using the package _phytools_ prior to plotting the tree with _ggtree_. 
The raw input data can be found [here](data\Montenegrina\ qualitative_data\ copy.dat) and the R code for the analysis in [here](analyses/analyses.R). I recoded the shape descriptions to PCH parameters as used in R to encode different shapes.

The resulting tree can be found here: 

![tree](analyses/black_tree.png)
