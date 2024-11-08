---
title:Topology weighting with TWISST
---

An interesting question in phylogenomics is how phylogenetic signal is distributed along the genome. Knowing the amount and location of genomic regions supporting a given topology may inform on evolutionary processes such as introgression. In principle, this is a straightforward thing to do: one simply need to calculate gene trees in sliding windows and match them to each possible species tree topology. However things get tricky in the presence of allele sharing across species (i.e., species that are not reciprocally-monophyletic at a gene tree), which is often the case when considering small genomic regions and recently diverged lineages. Topology Weighting overcomes this issue by supsampling tips of the gene tree to quantify its support to all possible species tree topologies. For more detailed explanations on this method, I recommend reading the [TWISST paper by Martin & Van Belleghem (2017)](https://academic.oup.com/genetics/article/206/1/429/6064218), as well as [this paper by Stankowski et al. (2024)](https://www.science.org/doi/full/10.1126/science.adi2982) (especially Figure 2).

## Preparing the input trees

Similarly to ASTRAL, TWISST takes as input a list of phylogenetic trees in newick format, in a single file. Because TWISST is a purely descriptive method, it does not need filters to retain loci meeting certain underlying assumptions. Especially, we do not want to thin the windows here, because we want to have the densest coverage of the genome possible to see fine scale changes in topological support. The only requirement is that all samples are present in all input gene trees. Because Topowindows sometimes remove sequences if they hve too many missing genotypes, this may not always be the case. Thus, we first need to remove gene trees with missing tips.

```R:

metadata <- read.table("TW_tutorial_windows_stats.tsv", header=T)
trees <- read.tree("TW_tutorial.trees")

# The number of tips is indicated in the NTIPS column of the metadata table. Here we have a maximum of 58 tips (29 phased samples)
trees.TWISST <- trees[which(metadata$NTIPS == 58)]
write.tree(trees.TWISST, "TW_tutorial_TWISST.trees")
```

## Running TWISST

TWISST is distributed as a python script and just requires a few python packages to be installed. Download and installation instructions can be found [here](https://github.com/simonhmartin/twisst).
To run the analysis, we need the input trees and to assign samples to groups. The groups will be used to define the possible species tree topologies to be tested. For exemple, if I define three group A, B and C, the following rooted species tree are possible: (A,(B,C)); (B,(A,C)); and (C(A,B));. A set of outgroup samples must be included to root the gene trees, and TWISST will quantify the support of each gene tree to all possible species trees.
