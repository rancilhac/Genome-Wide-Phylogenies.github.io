# Calulating phylogenetic trees in sliding windows

In the first part of this tutorial, we will calculate phylogenetic trees in sliding windows along a chromosome.

## Getting started

This analysis can be run on any system with R installed. Download and install TopoWindows as detailed here: https://github.com/rancilhac/TopoWindows. Create a working directory containing the TopoWindows scripts and your input vcf (or vcfs). The output files will saved there as well.

## Inferring the trees

TopoWindows provides two functions for sliding windows `topo.windows.sites` and `topo.windows.coord`. The difference is that in the latter the size of the windows is defined in terms of coordinates on the chromosome, while in the former it is defined based on a fixed number of sites. Besides this difference both work essentially in the same way. Whether to use fixed coordinates or number of sites will depend on the aim of the analysis and the dataset (i.e., variations in SNP density). For this exemple, we will use `topo.windows.sites`.

To infer neighbor-joining phylogenies in non-overlapping 500 SNPs windows, run the following:
```R:
setwd("TopoWindows_tutorial")
source("Topo_windows_v03.R")
topo.windows.sites(vcf = "YW_ASTRAL_PHASED_57i_CHR23.recode.vcf", size = 500, incr = 0, phased = T, prefix = "test", 
                   write.seq = T, nj = T, dna.dist = "JC69")
```
Here are what the options do:
`vcf`: name of the vcf file (or full path if it is not in the working directory)
`size`: size of the window (number of sites). Note that invariant sites are removed when the vcf is read in (see the documentation of vcfR for more details)
`incr`: an increment defining overlap across windows (number of sites). 0 means no overlap
`phased`: whether the data is phased or not. If the data is phased, two sequences are used for each (diploid) samples. If the data is not phased, a consensus sequence is generated with IUPAC ambiguity codes used at heterozygous positions
`prefix`: a prefix for the output files
`write.seq`: wether to write multispecies alignments in fasta format for each window
`dna.dist`: the model to use to calculate distances across sequences. For a complete list of models, see the documentation of the function dist.dna() in the ape package

Most of these parameters are relatively straighforward to set. The size of the window depends a lot on the aim of the analysis, divergence times across species, the proportion of allele sharing and LD decay. I recommend to try different values and vizualy explore the resulting trees before deciding on a value. Whether to use overlapping windows depend on the question addressed. Typically, phylogenetic inference based on gene trees requires independent markers, and thus overlapping windows won't do. On the other hand, genome scans aimed at detecting fine-scale changes in topological support may benefit from using overlapping windows.

## Exploring the results



