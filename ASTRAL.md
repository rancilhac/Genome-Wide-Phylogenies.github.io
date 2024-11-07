---
title: Species tree inference with ASTRAL
---

In this page, we will see how to use the trees calculated at the previous step to infer a _species tree_. In this context, _species tree_ refers to a phylogeny representing the history of splits across individuals in the dataset, as opposed to a _gene tree_, which represents the evolutionary relationships of the DNA sequences carried by these individuals at a given locus of the genome. In that sense, the sliding windows trees that we inferred earlier are _gene trees_, with loci defined by the windows. The distinction between _species tree_ and _gene tree_ is important, because the relationship at any given locus may be discordant from the actual evolutionary history of the individuals, because of variations in the order of coalescence events that can be caused for exemple by introgression or incomplete lineage sorting. Despite the name, a _species tree_ does not necessarily describe relationships between species, and the tips can represent any taxonomic level. Assuming an absence of introgression, the _species tree_ is a function of the levels of variation across _gene trees_, and several methods have been developped to estimate the former from the latter (in the presence of introgression things get trickier, but we will talk about it in another page). Here we will use a software called ASTRAL. More details about species trees and related softwares can be found in the litterature.

## Preparing the input trees

ASTRAL takes as input a set of gene trees in newick format, concatenated in a single file (1 tree per line). The trick is that the loci used to calculate gene trees must satisfy several conditions to avoid violating the statistical model underlying ASTRAL:
* ASTRAL assumes unlinked loci, i.e. the topology of a given gene tree is independent of the previous and next ones. In other words, the windows used in this analysis should be under linkage equilibrium. To take that into account, we can thin our tree set so that the windows considered are separated by enough physical distance to ensure free recombination in most cases. This distance can be informed by Linkage Disequilibrium decay analyses, and additionnal knowledge on the recombination landscape.
* ASTRAL assumes an absence of within-locus recombination, i.e. that each window corresponds to a single haplotype. This is difficult to take into account, especially when recombination rate varies strongly. What we will do is to filter windows based on their physical size, as chances of recombination events increase with size. Again, LD decay analyses may inform on the best size to use.
* One may want to filter based on the proportion of missing data and parsimony informative sites to avoid errors in gene trees. Whether this is important or not depends on the filters applied to SNPs before calculating the trees.

We can filter the trees as shown in the previous page:

```R:

metadata <- read.table("TW_tutorial_windows_stats.tsv", header=T)
trees <- read.tree("TW_tutorial.trees")

# this is a custom function to thin the windows based on the desired physical size
thin.pos <- function(pos, start, interval){
  thinned.pos <- c(start)
  p <- pos[start]
  i <- 2
  while(i <= length(pos)){
    if(pos[i]-p >= interval){ thinned.pos <- c(thinned.pos, i)
    p <- pos[i]
    i <- i+1 }
    else if(pos[i]-p < interval){ i <- i+1 }
  }
  return(thinned.pos)
}

# first we select windows with less than 40% missing data and a physical size of 15 kbp
metadata.ASTRAL <- metadata[which(metadata$PROP.MISS < 0.4 & metadata$CHR.SIZE < 15000), ]
# next, we thin the windows so that they are separated by at least 25kbp
metadata.ASTRAL <- metadata.ASTRAL[thin.ranges(metadata.ASTRAL[,2:3], 25000), ]

trees.ASTRAL <- trees[as.numeric(riow.names(metadata.ASTRAL))]
write.tree(trees.ASTRAL, "TW_tutorial_ASTRAL.trees")
```

