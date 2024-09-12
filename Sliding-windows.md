# Calulating phylogenetic trees in sliding windows

In the first part of this tutorial, we will calculate phylogenetic trees in sliding windows along a chromosome.

## Getting started

This analysis can be run on any system with R installed. Download and install TopoWindows as detailed here: https://github.com/rancilhac/TopoWindows. Create a working directory containing the TopoWindows scripts and your input vcf (or vcfs). The output files will saved there as well.

## Inferring the trees

TopoWindows provides two functions for sliding windows `topo.windows.sites` and `topo.windows.coord`. The difference is that in the latter the size of the windows is defined in terms of coordinates on the chromosome, while in the former it is defined based on a fixed number of sites. Besides this difference both work essentially in the same way. Whether to use fixed coordinates or number of sites will depend on the aim of the analysis and the dataset (i.e., variations in SNP densities). For this exemple, we will use `topo.windows.sites`.

To infer phylogenies in non-overlapping 500 SNPs windows, run the following:
```R:
setwd("TopoWindows_tutorial")
source("Topo_windows_v03.R")
topo.windows.sites("YW_ASTRAL_PHASED_57i_CHR23.recode.vcf", size = 1000, incr = 0, phased = T, prefix = "test", 
                   write.seq = T, nj = T, dna.dist = "JC69")
```

