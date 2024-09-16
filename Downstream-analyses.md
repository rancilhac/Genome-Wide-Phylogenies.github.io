---
title: Trees exploration and downstream analyses
---

Now that we have calculated the trees, we can explore them and use them for downstream phylogenetic analyses. Here, we will see how to filter the trees based on the windows metadata, visualize the trees and run downstream phylogenetic analyses.

## Filtering the trees

First, let's import the trees and windows metadata in R. 

```R:
library(ape)

# the metadata are read as a table, with one genomic window per line
metadata <- read.table("test_windows_stats.tsv", header=T)

# the trees are imported as a multiPhylo object of the ape package, which functions like a list
trees <- read.tree("test_NJ_trees.trees")
```
The two files should have the same number of lines, corresponding to the number of windows. However, the lines containing `NA` in the trees file (i.e., windows where the phylogenetic analysis failed) are not read by the `read.tree()` function, hence these windows must be removed from metadata table as well.

```R:
metadata <- na.omit(metadata)
nrows(metadata)
length(trees)
```
If after removing the `NA` the number of windows and trees still don't match, then something went wrong in the previous step. If they match, we can continue.

We can now filter the windows according to criteria of our choice. For exemple here we will keep only windows that span less than 50 kbp, extract the corresponding trees and write them to a new file.

```R:
windows.less.50kbp <- which(metadata$CHR.SIZE < 50000)
trees.less.50kbp <- trees[windows.less.50kbp]
write.tree(trees.less.50kbp, "trees_less_50kbp.trees")
```

A similar approach can be used to filter the windows based on the different columns of the `metadata` table (proportion of parsimony informative sites and missing genotypes, number of tips in the tree).

## Manipulating and visualizing the trees

Visualy inspecting the trees can be useful to detect problems or get a first idea of the signal in the data. It can also be interesting to edit the trees, for exemple to drop some samples before further analyses. `ape` provides some useful tools for that.

- Rooting the trees: 
```R:
outgroup <- c("01", "02")
trees.rooted <- root.multiPhylo(trees)
```

- Plotting the trees:
```R:
# plot a single tree:
plot(trees.rooted[[1]])

# plot all trees one after the other (press "enter" to see the next and "escape" to exit)
for(i in 1:length(trees.rooted)){
  plot(trees.rooted[[i]])
  readline(prompt="Press [enter] to continue")
}
```
Further plotting options are available in the `plot.phylo()` function, and in the [ggphylo R package](https://github.com/gjuggler/ggphylo).

- Removing tips from the trees:
```R:
tips.to.remove <- c("01", "02")
trees.pruned <- drop.tip.multiphylo(trees.rooted, tips.to.remove)
```

