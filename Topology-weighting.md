---
title:Topology weighting with TWISST
---

An interesting question in phylogenomics is how phylogenetic signal is distributed along the genome. Knowing the amount and location of genomic regions supporting a given topology may inform on evolutionary processes such as introgression. In principle, this is a straightforward thing to do: one simply need to calculate gene trees in sliding windows and match them to each possible species tree topology. However things get tricky in the presence of allele sharing across species (i.e., species that are not reciprocally-monophyletic at a gene tree), which is often the case when considering small genomic regions and recently diverged lineages. Topology Weighting overcomes this issue by supsampling tips of the gene tree to quantify its support to all possible species tree topologies. For more detailed explanations on this method, I recommend reading the [TWISST paper by Martin & Van Belleghem (2017)](https://academic.oup.com/genetics/article/206/1/429/6064218), as well as [this paper by Stankowski et al. (2024)](https://www.science.org/doi/full/10.1126/science.adi2982) (especially Figure 2).


