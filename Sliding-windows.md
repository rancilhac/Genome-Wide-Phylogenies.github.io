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
For a vcf with 87000 SNPs, it takes a couple of minutes. Here are what the options do: \
`vcf`: name of the vcf file (or full path if it is not in the working directory), can be gziped \
`size`: size of the window (number of sites). Note that invariant sites and indels are removed when the vcf is read in (see the documentation of vcfR for more details) \
`incr`: an increment defining overlap across windows (number of sites). 0 means no overlap \
`phased`: whether the data is phased or not. If the data is phased, two sequences are used for each (diploid) samples. If the data is not phased, a consensus sequence is generated with IUPAC ambiguity codes used at heterozygous positions \
`prefix`: a prefix for the output files \
`write.seq`: whether to write multispecies alignments in fasta format for each window. For phased data, the two haplotypes are named `sample_name_0` and `sample_name_1` \
`dna.dist`: the model to use to calculate distances across sequences. For a complete list of models, see the documentation of the function dist.dna() in the ape package \

Most of these parameters are relatively straighforward to set. The size of the window depends a lot on the aim of the analysis, divergence times across species, the proportion of allele sharing and LD decay. I recommend to try different values and vizualy explore the resulting trees before deciding on a value. Whether to use overlapping windows depend on the question addressed. Typically, phylogenetic inference based on gene trees requires independent markers, and thus overlapping windows won't do. On the other hand, genome scans aimed at detecting fine-scale changes in topological support may benefit from using overlapping windows.

Practical note: the vcf importation into R can take a lot of memory. For very large chromosomes, I sometime had to split the vcfs in several parts.

## Outputs

Two output files were created:
`test_NJ_trees.trees`: this file contains the neighbor-joining trees in Newick format, one per line. In case the calculation of the tree fails, NA is printed instead of a tree. It can be opened in a regular phylogenetic tree viewer (e.g., FigTree)
`test_windows_stats.tsv`: this file is a tab-separated table containing metadata for each window on a separate line

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

Next we will see how to use these two files to filter and subset the trees.

## Inferring trees using likelihood approaches

The `TopoWindows` R functions can only calculate NJ trees internaly. However, when `write.seq = T`, fasta alignments are written for each windows that can be analysed with an external software. The alignements are written in `test_sequences` with names like `chr22-67-20215.fasta`(CHR-CHR.START-CHR.END.fasta).
I tend to use IQTREE for maximum-likelihood phylogenetic inference, but other softwares can also be used (e.g., RAxML) depending on preference, including Bayesian inference softwares. Here is how to calculate the trees with IQTREE in a UNIX comand line environment:
```bash:
cd test_sequences

for i in *fasta
do
iqtree2 -s ${i} -m TEST+ASC -T 1
done

cat *.treefile > test_ML_trees.trees
```

The `.trees` file can now be used together with the `.tsv` file produced in the previous step. If you are not interested in getting nj trees, only ML trees, you can run `topo.windows.sites` with `write.seq = T` and `nj = F`. 

## Working in a High-Performance Cluster environment (or parallelizing in a UNIX command line environment)

Working directly in an R environment is not very handy because analyses of large contigs quickly become memory- and time-consuming. Here we will see how to parallelize the runs on a chromosome/contig basis in a High-Performance Cluster with the Slurm job manager. A similar approach can probably be achieved on a UNIX working station with the `parallel` command, although I have never tried.

First, we need to prepare single VCFs for each chromosome/contig and use a consistent naming convention such as `chr1.vcf.gz`, `chr2.vcf.gz`, `chr3.vcf.gz` ect...

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

#Read the list of chromosomes/contigs into a variable
CHRLIST=$(<all_vcfs.txt)
#Identify the current chromosome based on the array ID (starting with 0 because indices are 0-based in bash contrary to e.g. R)
CHR=${CHRLIST}[${SLURM_ARRAY_TASK_ID}]

#Define a prefix for the output files (adapt to your naming convention, here it takes all the vcf file name before .vcf.gz)
PREF=$(echo ${CHR} | cut -d'.' -f1)

Rscript Topo_windows_v03_cl_wrapper.R --prefix ${PREF} --vcf ${CHR} --type s --size 500 --incr 0 --phased T --nj T --ali T --dist NJ69
```
The options are essentially the same as earlier. `type` defines whether to use `topo.windows.sites` (`--type s`) or `topo.windows.coord` (`--type c`).
This will result in the same output files as previously. 

*Note*: Memory usage will greatly vary depending on chromosome size. To optimize the runs, I recommend to split the small and large chromosome in different runs (because memory specification will be the same for all jobs in the array), and maybe splitting very large chromosomes if run time are too long (only for small windows).

*Note 2*: If you have a mix of phased and unphased chromosomes (e.g., when sex chromosomes are unphased), you will need to run them separately.

Finally, a similar prallelization strategy may be used to calculate ML trees from the fasta sequences. First, move to the `test_sequences` directory and create a list of the fasta files:
```bash:
ls *.fasta > all_fastas.txt
```

Now, we can use a similar job script as previously (again, assuming three fasta files):
```bash:
#SBATCH [regular Slurm specifications]
#SBATCH -a 0-2

#Read the list of fastas into a variable
FASTALIST=$(<all_vcfs.txt)
#Identify the current chromosome based on the array ID (starting with 0 because indices are 0-based in bash contrary to e.g. R)
FASTA=${FASTALIST}[${SLURM_ARRAY_TASK_ID}]

iqtree2 -s ${FASTA} -m TEST+ASC -T 1
```
Once the job is finished, the trees can be concatenated in a single file:
```bash:
cat *.treefile > test_ML_trees.trees
```


