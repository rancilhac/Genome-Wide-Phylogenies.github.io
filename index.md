---
title: Calculating local phylogenies from whole-genome data: A tutorial
---

In this page, you will find guidelines to calculate local phylogenetic trees based on a vcf from Whole-genome data and perform downstream analyses and vizualizations. Here, local phylogenies refer to phylogenetic trees inferred from specific genomic regions, either arbitrary windows or specific genomic features (genes, exons, ect...). 

## Why calculating local phylogenies?

Local phylogenetic signal can be used to investigate many evolutionary questions based on genomic data. Most prominently, this approach is useful to estimate phylogenetic trees and networks using methods that take gene trees as input (i.e., based on the Multispecies Coalescent), or to describe variations in phylogenetic signal along a genome using a Topology Weighting approach. Local phylogenies may be used to investigate many processes such as selection and introgression, and to study the evolutionary history of specific genomic regions. An important limitation of this approach is that recombination breakpoints need to be known to accurately study variation in phylogenetic relationships along a genome, which is practically impossible in most cases. Statistical approaches based on the Ancestral Recombination Graph are able to reconstruct local relationships while accounting for recombination, but these are computationaly very intensive (and potentially sensitive to demographic parameters), making them difficult to implement. Thankfully, local phylogenetic reconstruction, even though blind to recombination, is often sufficient to uncover evolutionary patterns.

## Prerequisites

In this tutorial I use a custom R script to calculate phylogenies using a vcf file as input: https://github.com/rancilhac/TopoWindows. The input vcf should contain data from a single chromosome/contig, and the analysis can be run in parallel on several chromosomes. In case of a very fragmented reference assembly, it is probably better to filter out the smaller contigs beforehand, to make the analysis more manageable. 
