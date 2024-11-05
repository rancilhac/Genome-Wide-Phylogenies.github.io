---
title: Calulating phylogenetic trees in sliding windows
---

## Getting started

This tutorial shows how to use the [TopoWindows R functions](https://github.com/rancilhac/TopoWindows) to calculate phylogenetic trees in sliding windows and user-defined genomic regions. These functions can be run on any system with R installed, and rely on the vcfR package to parse the data and the ape and phangorn packages for phylogenetic analyses and tree output. Further information can be found in the TopoWindows github repository.

The only input file needed is a vcf file containing single nucleotide variants from a single chromosome/contig. Variants can be phased or not, and the vcf can also contain invariant sites and indels, but those will be removed when the vcf is imported into R. The vcf file can be gziped. In the following exemples I will use a vcf ([download link](https://drive.google.com/file/d/1JGOdDQLBWibT8xaSoz4d31fKFJBDQUhF/view?usp=sharing)) corresponding to chromosome 20 in [Rancilhac et al. (2024)](https://academic.oup.com/sysbio/article/73/1/12/7294611), which includes 1,403,718 phased SNPs from 6 species of wagtails (songbirds).

To get started, place the vcf file and Topowindows scripts in a working directory and open an R terminal.

## Inferring the trees

In the first part of this tutorial, we will calculate phylogenetic trees in sliding windows along a chromosome.

TopoWindows provides two functions for sliding windows: `topo.windows.sites` and `topo.windows.coord`. The difference is that in the latter the windows are defined in terms of coordinates on the chromosome, while in the former they are defined based on a fixed number of variants. Besides this difference both work essentially in the same way. Whether to use fixed coordinates or number of variants will depend on the aim of the analysis and the dataset (i.e., variations in SNP density). For this exemple, we will use `topo.windows.sites`.

To infer neighbor-joining phylogenies in non-overlapping 500 SNPs windows, run the following:
```R:

setwd("TopoWindows_tutorial")

source("Topo_windows_v03.R")

topo.windows.sites(vcf = "BWW_chr20.recode.vcf.gz", size = 500, incr = 0, phased = T, prefix = "TW_tutorial", 
                   write.seq = T, tree = "NJ", dna.model = "JC", missing.thresh=0.7, force=F)
```
With the exemple vcf, it takes a couple of minutes. Here are what the options do: \
`vcf`: name of the vcf file (or full path if it is not in the working directory), can be gziped \
`size`: size of the window (number of sites). Note that invariant sites and indels are removed when the vcf is read in (see the documentation of vcfR for more details) \
`incr`: an increment defining overlap across windows (number of sites). 0 means no overlap \
`phased`: whether the data is phased or not. If the data is phased, two sequences are used for each (diploid) samples. If the data is not phased, a consensus sequence is generated with IUPAC ambiguity codes used at heterozygous positions \
`prefix`: a prefix for the output files \
`write.seq`: whether to write multispecies alignments in fasta format for each window. For phased data, the two haplotypes are named `sample_name_0` and `sample_name_1`. \
`tree`: whether to calulate trees. Can be "NJ" (neighbor-joining trees), "ML" (Maximum-likelihood trees) or "N" (no trees). \
`dna.model`: the nucleotide substitution model to use. When calculating NJ trees, this correspond to the models available to the dist.dna() function in ape. \
`missing.thresh`: At a given window, sequences will be remove if their proportion of missing SNPs is above the specified threshold. \
`force`: whether to overwrite pre-existing output files with the same prefix.

To use Maximum-likelihood instead of Neighbor-joining, simply change `tree = "NJ"` to `tree = "ML"`. The ML approach implemented is very simple (no bootstraping, no optimization of substitution model), but more elaborate inference can also be performed if needed, as discussed at the end of this page.

Most of these parameters are relatively straighforward to set. The size of the window depends a lot on the aim of the analysis, divergence times across species, the proportion of allele sharing and LD decay. I recommend to try different values and vizualy explore the resulting trees before deciding on a value. Whether to use overlapping windows depend on the question addressed. Typically, phylogenetic inference based on gene trees requires independent markers, and thus overlapping windows won't do. On the other hand, genome scans aimed at detecting fine-scale changes in topological support may benefit from using overlapping windows.

*Note*: the vcf importation into R can take a lot of memory. For very large chromosomes, it may be helpful to split the vcfs.

### Outputs

Two output files were created:
`TW_tutorial.trees`: this file contains the neighbor-joining trees in Newick format, one per line. In case the calculation of the tree fails, NA is printed instead of a tree. It can be opened in a regular phylogenetic tree viewer (e.g., FigTree)
`TW_tutorial_windows_stats.tsv`: this file is a tab-separated table containing metadata for each window on a separate line

```
CHR    CHR.START  CHR.END  CHR.SIZE  NSITES  PROP.MISS  PROP.PIS  TREE  NTIPS
chr22  67         20215    20148    500      0.0772     1         YES   114
chr22  20275      36961    16686    500      0.0742     1         YES   114
chr22  37118      61175    24057    500      0.0871     1         YES   114
```
`CHR`: chromosome name \
`CHR.START`: start coordinate \
`CHR.END`: end coordinate \
`CHR.SIZE`: actual size on the chromosome (i.e., CHR.END - CHR.START) \
`NSITES`: number of SNPs in the window \
`PROP.MISS`: proportion of missing genotypes \
`PROP.PIS`: proportion of parcimony-informative SNPs \
`TREE`: whether a NJ tree could be calculated (if not, `NA` in this column and in the `.trees` file) \
`NTIPS`: the number of tips in the tree (here twice the number of samples since the vcf contains phased diploid individuals, for unphased data it should be the number of samples). Note that this number can be smaller than the number of sequences, because some sequences are droped when they contain too much missing data. \

In the next pages of this tutorial, we will see how to manipulate these two files to subset the trees and use them in downstream analyses.

### Inferring trees using more elaborate approaches

The phylogenetic methods implemented in `TopoWindows` are fairly simple and most appropriate for time efficient inference of a large number of trees in closely related species. In some cases, more elaborated inference may be needed. This can be done based on the fasta alignments written for each windows with the `write.seq = T` option. Alignements are written in the `TW_tutorial_sequences` directory with names like `chr22-67-20215.fasta`(CHR-CHR.START-CHR.END.fasta), and can be fed into any phylogenetic inference softwares that accept fasta as input format. Here is an exemple with IQTREE 2 (ML inference) with two threads and 1000 bootstrap replicates:

```bash:
cd test_sequences

for i in *fasta
do
iqtree2 -s ${i} -m TEST+ASC -T 2 -B 1000
if [ -f ${i}.treefile ]
then
cat ${i}.treefile >> ML_trees_IQTREE.trees
else
echo "NA" >> ML_trees_IQTREE.trees
fi
done
```

The `.trees` file can now be used together with the `.tsv` file produced in the previous step.

### Working in a High-Performance Cluster environment (or parallelizing in a UNIX command line environment)

Working directly in an R environment is not very handy because analyses of large contigs quickly become memory- and time-consuming. Here we will see how to parallelize the runs on a chromosome/contig basis in a High-Performance Cluster with the Slurm job manager. A similar approach can probably be achieved on a UNIX working station with the `parallel` command, although I have never tried.

First, we need to prepare single VCFs for each chromosome/contig and use a consistent naming convention such as `chr1.vcf.gz`, `chr2.vcf.gz`, `chr3.vcf.gz` ect... This can easily be achieved with `bcftools (v1.18)`:
```bash:
# generate a list of all chromosome from original vcf
bcftools query -f '%CHROM\n' file.vcf | sort | uniq > chromosome_list.txt

# 
while read CHR
do
bcftools view -r ${CHR} -O z -o ${CHR}.vcf.gz file.vcf
done < chromosome_list.txt
```

Next, we create a text file listing all the input files (one per line):
```bash:
ls *.vcf.gz > all_vcfs.txt
```

We will use Slurm job arrays to run `topo.windows.sites` on all vcfs simultaneously. To do that, create a slurm job script with the specifications corresponding to your cluster. At the end of the header, add the following line:
```bash:
#SBATCH -a 0-2
```

This tells slurm to run an array of three jobs with indices 0, 1 and 2. You need to specify as many jobs as you have chromosomes, counting on a 0-basis (here I have three chromosomes, hence). For each job, the index is stored in a variable named ${SLURM_ARRAY_TASK_ID}. We can now use the command line wrapper `Topo_windows_v03_cl_wrapper.R` to run the R script from a UNIX command line environment. Here is what jour job file should look like:
```bash:
#SBATCH [regular Slurm specifications]
#SBATCH -a 0-2

#Read the list of single-chromosome vcf files
VCFLIST=$(<all_vcfs.txt)
#Identify the current chromosome based on the array ID (starting with 0 because indices are 0-based in bash contrary to e.g. R)
VCF=${VCFLIST}[${SLURM_ARRAY_TASK_ID}]

#Define a prefix for the output files (adapt to your naming convention, here it takes all the vcf file name before .vcf.gz)
PREF=$(echo ${CHR} | cut -d'.' -f1)

Rscript Topo_windows_v03_cl_wrapper.R --prefix ${PREF} --vcf ${VCF} --type s --size 500 --incr 0 --phased T --tree NJ --ali T --dna-model NJ --force F --missingness 0.7
```
The options are essentially the same as earlier. `type` defines whether to use `topo.windows.sites` (`--type s`) or `topo.windows.coord` (`--type c`).
This will result in the same output files as previously. 

*Note*: Slurm arrays are fairly flexible and can be adapted to different cases. Here I use them as indices to pick file names in a bash array, but if your files are named with numbers you can refer to them directly with the ${SLURM_ARRAY_TASK_ID} variable. For exemple, here my files are named `chr1.vcf.gz`, `chr2.vcf.gz` and `chr3.vcf.gz` so I could use the following script:
```bash:
#SBATCH [regular Slurm specifications]
#SBATCH -a 1-3

#Define a prefix for the output files (adapt to your naming convention, here it takes all the vcf file name before .vcf.gz)
PREF=$(echo chr${SLURM_ARRAY_TASK_ID}.vcf.gz | cut -d'.' -f1)

Rscript Topo_windows_v03_cl_wrapper.R --prefix ${PREF} --vcf chr${SLURM_ARRAY_TASK_ID}.vcf.gz --type s --size 500 --incr 0 --phased T --tree NJ --ali T --dna-model NJ --force F --missingness 0.7
```

*Note 2*: Memory usage will greatly vary depending on chromosome size. To optimize the runs, I recommend to split the small and large chromosome in different runs (because memory specification will be the same for all jobs in the array), and maybe splitting very large chromosomes if run time are too long (it should happen only for small windows).

*Note 3*: If you have a mix of phased and unphased chromosomes (e.g., when sex chromosomes are unphased), you will need to run them separately.

Finally, a similar prallelization strategy may be used to calculate ML trees from the fasta sequences. First, move to the `TW_tutorial_sequences` directory and create a list of the fasta files:
```bash:
ls *.fasta > all_fastas.txt
```

Now, we can use a similar job script as previously (again, assuming three fasta files):
```bash:
#SBATCH [regular Slurm specifications]
#SBATCH -a 0-2

#Read the list of fastas into a variable
FASTALIST=$(<all_fastas.txt)
#Identify the current chromosome based on the array ID (starting with 0 because indices are 0-based in bash contrary to e.g. R)
FASTA=${FASTALIST}[${SLURM_ARRAY_TASK_ID}]

iqtree2 -s ${FASTA} -m TEST+ASC -T 1
```
Once the job is finished, the trees can be concatenated in a single file:
```bash:
cat *.treefile > test_ML_trees.trees
```
*Note*: Slurm arrays are limited to 999, so this will require a bit of adaptation for larger numbers of windows.


