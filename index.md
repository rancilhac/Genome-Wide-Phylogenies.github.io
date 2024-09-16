---
title: Calculating local phylogenies from whole-genome data - A tutorial
---

In this page, you will find guidelines to calculate local phylogenetic trees from Whole-genome data based on a vcf file and perform downstream analyses and vizualizations. Here, local phylogenies refer to phylogenetic trees inferred from specific genomic regions, either arbitrary windows or specific genomic features (genes, exons, ect...). 

## Why calculating local phylogenies?

Local phylogenetic signal can be used to investigate many evolutionary questions based on genomic data. Most prominently, this approach is useful to estimate phylogenetic trees and networks using methods that take gene trees as input (i.e., based on the Multispecies Coalescent), or to describe variations in phylogenetic signal along a genome using a Topology Weighting approach. Local phylogenies may be used to investigate many processes such as selection and introgression, and to study the evolutionary history of specific genomic regions. An important limitation of this approach is that recombination breakpoints need to be known to accurately study variation in phylogenetic relationships along a genome, which is practically impossible in most cases. Statistical approaches based on the Ancestral Recombination Graph are able to reconstruct local relationships while accounting for recombination, but these are computationaly very intensive (and potentially sensitive to demographic parameters), making them difficult to implement. Thankfully, local phylogenetic reconstruction, even though blind to recombination, is often sufficient to uncover evolutionary patterns.

The approach presented here is not appropriate for reduced-representation data (e.g., RADseq), which do not yield variants continuously distributed across chromosomes.

## Contact

If you have any questions or suggestions, you can contact me at loisrancilhac@gmail.com.
