---
title: Topology weighting with TWISST
---

An interesting question in phylogenomics is how phylogenetic signal is distributed along the genome. Knowing the amount and location of genomic regions supporting a given topology may inform on evolutionary processes such as introgression. In principle, this is a straightforward thing to do: one simply need to calculate gene trees in sliding windows and match them to each possible species tree topology. However things get tricky in the presence of allele sharing across species (i.e., species that are not reciprocally-monophyletic at a gene tree), which is often the case when considering small genomic regions and recently diverged lineages. Topology Weighting overcomes this issue by supsampling tips of the gene tree to quantify its support to all possible species tree topologies. For more detailed explanations on this method, I recommend reading the [TWISST paper by Martin & Van Belleghem (2017)](https://academic.oup.com/genetics/article/206/1/429/6064218), as well as [this paper by Stankowski et al. (2024)](https://www.science.org/doi/full/10.1126/science.adi2982) (especially Figure 2).

## Preparing the input trees

Similarly to ASTRAL, TWISST takes as input a list of phylogenetic trees in newick format, in a single file. Because TWISST is a purely descriptive method, it does not need filters to retain loci meeting certain underlying assumptions. Especially, we do not want to thin the windows here, because we want to have the densest coverage of the genome possible to see fine scale changes in topological support. The only requirement is that all samples are present in all input gene trees. Because Topowindows sometimes remove sequences if they have too many missing genotypes, this may not always be the case. Thus, we first need to remove gene trees with missing tips.

```R:

metadata <- read.table("TW_tutorial_windows_stats.tsv", header=T)
trees <- read.tree("TW_tutorial.trees")

# The number of tips is indicated in the NTIPS column of the metadata table. Here we have a maximum of 58 tips (29 phased samples)
trees.TWISST <- trees[which(metadata$NTIPS == 58)]
write.tree(trees.TWISST, "TW_tutorial_TWISST.trees")
```

## Running TWISST

TWISST is distributed as a python script and just requires a few python packages to be installed. Download and installation instructions can be found [here](https://github.com/simonhmartin/twisst).
To run the analysis, we need the input trees and to assign samples to groups. The groups will be used to define the possible species tree topologies to be tested. For exemple, if I define three group A, B and C, the following rooted species tree are possible: (A,(B,C)); (B,(A,C)); and (C(A,B));. A set of outgroup samples must be included to root the gene trees, and TWISST will quantify the support of each gene tree to all possible species trees. Samples can be mapped to groups using a tab separated table (here the samples are included twice with _0 adn _1 suffixes because the data is phased):

```bash:
head samples_map.txt
P12604_106_S6_0 Mot_agu
P12604_106_S6_1 Mot_agu
P12604_107_S7_0 Mot_agu
P12604_107_S7_1 Mot_agu
P12604_108_S8_0 Mot_agu
P12604_108_S8_1 Mot_agu
P12604_109_S9_0 Mot_agu
P12604_109_S9_1 Mot_agu
P12604_110_S10_0        Mot_agu
P12604_110_S10_1        Mot_agu
```
The complete sample map for the exemple data set can be downloaded [here](https://drive.google.com/file/d/1XUbxjjQfhvEqIsqSvVDafs0FhK2ylKVE/view?usp=sharing).

Now let's run TWISST:
```bash:
python twisst.py -t TW_tutorial_TWISST.trees -w TWISST.weights --outputTopos TWISST.topologies --outgroup Mot_cin_outgroup --method complete -g Mot_agu -g Mot_alb -g Mot_gra -g Mot_sam -g Mot_mad -g Mot_cin_outgroup --groupsFile samples_map.txt
```
## Vizualizing the results: summaries of all gene trees

TWISST outputs two files: `TWISST.topologies` contains a list of all possible species tree topologies, in parenthetical format; `TWISST.weights` is a table containing the weights (as well as the list of topologies at the beggining of the file). In the weights table each column is a topology and each line a gene tree. The value in a cell is the number of subtree of a given gene tree that match a given species tree. Here is a simple example with three gene trees and three species trees:
```bash:
T1  T2  T3
100 0   0
0   50  50
100 0   0
```
This tells us that the first and third gene tree match species tree T1, while the second supports both T2 and T3, indicating allele sharing across species.
Let's now have a look at the real results. As usual, we will import this data into R to explore the data and do some vizualization:

```R:
weights <- read.table("TWISST.weights", header=T)

# TWISST gives raw weights, but proportions are easier to interpret (i.e., the proportion of subtrees supporting a topology, rather than the raw number of subtrees). Let's first convert counts into proportions:

prop.weight <- function(raw){raw/apply(raw, 1, sum)}
weights <- prop.weight(weights)

# Now we can calculate the mean weight for each topology

mean.weight <- apply(weights, 2, mean)
barplot(mean.weight[order(mean.weight)], ylab="Mean weight")
```
![barplot](/Genome-Wide-Phylogenies.github.io/assets/Mean_weight_barplot.png)

As you can see, the results are fairly difficult to interpret. This is because the number of possible species tree topologies increases very fast with the number of species. Here, with 5 ingroup species there are already 105 possible rooted topologies. To keep things manageable, it is best to consider trees with less species. Three ingroup species is the best possible scenario (=3 species tree topologies), but four species is also manageable (=15 species tree topologies). Another observation is that the support is low even for the most supported topology. This is because most topologies are largely redundant. However, if there are high levels of allele sharing (as is the case with these species), one can expect every topology to receive a bit of support. With a large number of species tree topologies like here, this means that even the support for the commonest topology will remain low.
For the sake of the example, and simpler visual representations, we will now consider the three commonest topologies:

```R:
library(ggplot2)

# Names of the three commonest topologies
topo.names <- names(mean.weight[order(mean.weight, decreasing=T)][1:3])

# Let's have a look at the three most supported topologies
topologies <- read.tree("TWISST.topologies")

par(mfrow=c(2,2))
plot(topologies[[order(mean.weight, decreasing=T)[1]]], main=topo.names[1])
plot(topologies[[order(mean.weight, decreasing=T)[2]]], main=topo.names[2])   
plot(topologies[[order(mean.weight, decreasing=T)[3]]], main=topo.names[3])

# We can also redo the barplot with only the three most common topologies
barplot(mean.weight[order(mean.weight, decreasing=T)[1:3]], ylab="Mean weight")
```
![topos](/Genome-Wide-Phylogenies.github.io/assets/Plot_3_commonest_topos.png)

A more informative alternative to barplots is to use violin plots, which will not only give us information about the mean weights of a given species tree topology but also the weights distribution across gene trees.

```R:
library(ggplot2)

ggplot(weights) + geom_violin(aes(x = "T13", y = topo13)) + geom_boxplot(aes(x = "T13", y = topo13), width=0.05) +
  geom_violin(aes(x = "T101", y = topo101)) + geom_boxplot(aes(x = "T101", y = topo101), width=0.05) +
  geom_violin(aes(x = "T48", y = topo48)) + geom_boxplot(aes(x = "T48", y = topo48), width=0.05) +
  theme_classic() + ylab("Weight") + xlab("Topology") +
  scale_x_discrete(limits=c("T13", "T101", "T48"))
```
![violin](/Genome-Wide-Phylogenies.github.io/assets/Violin_plots_TW.png)

As expected with high levels of allele sharing and large number of species trees, the weight distribution is skewed towards small values. There are few intermediates (~0.5), but also some gene trees that fully support these topologies (weight = 1).

One final strategy to visualize all gene trees jointly is to use a ternary plot. This vizualization works only when there are three species tree topologies, and was introduced by [Stankowski et al. (2024)](https://www.science.org/doi/full/10.1126/science.adi2982). It essentially encompass the same information as violin plots but in a more compelling and esay to get way. I also allows to test for assymetries which indicate introgression (see Stankowski et al. for descriptions of this test). The supplements of this paper also includes lots of simulations showing expected ternary distributions under a wide range of demographic scenarios.

```R:
library(ggtern)

ggtern(data=weights, aes(x=topo13, y=topo101, z=topo48)) + geom_hex_tern(bins=150) + scale_fill_gradientn(colours = heat.colors(5))
```

