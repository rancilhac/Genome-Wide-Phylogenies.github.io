---
title: Trees manipulation and exploration
---

Now that we have a set of trees, we can manipulate, filter and vizualize them before using them for downstream analyses. The ape package offers many functionalities for that.

## Filtering the trees

First, let's import the windows metadata in R: 

```R:
# the metadata are read as a table, with one genomic window per line
metadata <- read.table("TW_tutorial_windows_stats.tsv", header=T)

dim(metadata)
> [1] 2808    9

head(metadata)
>    CHR CHR.START CHR.END CHR.SIZE NSITES PROP.MISS PROP.PIS TREE NTIPS
1 chr20      2972  164561   161589    500    0.1147    0.148  YES    58
2 chr20    164563  169043     4480    500    0.0312    0.098  YES    58
3 chr20    169062  173762     4700    500    0.0292    0.082  YES    58
4 chr20    173834  179341     5507    500    0.0523    0.096  YES    58
5 chr20    179345  184187     4842    500    0.0356    0.084  YES    58
6 chr20    184191  189102     4911    500    0.0420    0.092  YES    58
```
Next, we can import the trees with the `read.tree()` function. This creates a `multiPhylo` object, which is essentially a list of phylogenetic trees and can be manipulated as a regular list in R:
```R:
library(ape)

trees <- read.tree("TW_tutorial.trees")
trees
>2808 phylogenetic trees
```
There should be as many lines in the `metadata` table as items in the `trees` object. However, phylogenetic inference may fail in some windows, in which case `NA` is written in both the `TREE` column of the `metadata` file and in the corresponding line of the `trees` file. These lines are not recognised by `read.tree()`, meaning that `trees` will have less items than rows in `metadata`. This can be corrected using the following:

```R:
metadata <- na.omit(metadata)
nrow(metadata) == length(trees)
>[1] TRUE
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
for(i in 1:length(trees.rooted))
  plot(trees.rooted[[i]])
  readline(prompt="Press [enter] to continue")
```
Further plotting options are available in the `plot.phylo()` function, and in the [ggphylo R package](https://github.com/gjuggler/ggphylo).

- Removing tips from the trees:
```R:
tips.to.remove <- c("01", "02")
trees.pruned <- drop.tip.multiphylo(trees.rooted, tips.to.remove)
```

