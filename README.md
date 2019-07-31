# SNPs HumanOrigins Recipes
## Commands and R scripts to analyse SNP data

This project contains different folders (sections) for some of the most used analysis methods to be performed on SNP chip databases. It is particularly tailored to the Human Origins Affymetrix array, and for analysis of human population history.

Most of the commands are part of the analysis pipeline for the paper **The current genomic landscape of western South America: Andes, Amazonia and Pacific Coast** https://10.1093/molbev/msz174

The commands use publicity available softwares and packages such as:
 
### PLINK https://www.cog-genomics.org/plink2 

for data manipulation and several useful analysis tools. Uses the plink file format(s). Extensive documentation available online. The general syntax on the terminal is:
```
plink --file nameofyourfile --flagAction1 --flagAction2 --etc.
```
PLINK  parses each command line as a collection of flags (each of which starts with two dashes), plus parameters (which immediately follow a flag), such as nameofyourfile in the example above.
File format include usually two files for the same dataset: one is the tables of samples and one is the variant calls. In my commands I work with the regular text format, including the two files nameofyourfile.ped and nameofyourfile.map , The .ped includes ID, pedigree(optional) + genotype table, the .map is basically the list of SNPs with chromosome position and alleles. The two file have to be called with the same name. 
Because PLINK was developed for GWAS medical studies, many basic informations will be not used in our analysis, such as pedigree or phenotype. 


### EIGENSOFT package - EIGENSTRAT programs https://github.com/DReichLab/EIG

It includes CONVERTF to manipulate files and POPGEN for running PCA, between other functions available. The general syntax on the terminal is:
```
convertf -p parfile
```
Where convertf is the command chosen and parfile is a file that includes the path and names of the input files, path and names of the output files to be generated, and additional files to call during the command.
The input files in the Eigenstrat format are nameofyourfile.geno, nameofyourfile.snp, nameofyourfile.ind. the .geno file contains genotype data for each individual (1 line per SNP); the .snp file contains information about each SNP; the .ind file contains information about each individual. Other formats are also supported. 

Other commands from software packages such as bcftools https://samtools.github.io/bcftools/ are sporadicly used. Bash scripts and some awk simple commands come in help.


## CONTENT
Start with the General Commands for how to first manipulate your dataset, and then check inside each folder for the relative command.md file, which will contain instructions to perform the analysis proposed. In most cases the commands refer to a separate R script to further elaborate and plot data. In some cases also examples of files to prepare separately for the analysis are included.

The General Commands and its related R script visualize file include format conversion, data screening for relatives and outliers, Runs of Homozygosity (ROH), merging datasets, visualizing data on a map (Figure S1), correlation between percentage of native ancestry and consanguinity (Figure S6), plotting ROH in bins (Figure 3, Figure S5).

The uniparental commands in R include some script to elaborate the haplogrup assignation for Y chromosome and mtDNA and plot frequency pie charts on a map (Figure S8). it also includes one script to plot percentage of native ancestry in autosomal (from supervised ADMIXTURE), y ch and mtdna (Figure S10).

The folder BEAGLEandIBD includes the terminal commands to run Beagle and refindedIBD. its relative R script includes commands to plot IBD sharing between individuals, between populations (Figure 4a, Figure S12), within populations, and between populations over a map (Figure 4B).

The folder BEAGLE includes old scripts and should not be considered.

The folder PCA includes commands to run smartPCA and R scripts to visualize PCA (Figure 2) and calculate Neighbor Joining trees, MDS and heatplots from Fst Distances (Figure S3), and plot it.  

The folder ADMIXTURE includes commands to run ADMIXTURE and R scripts to elaborate the results of ADMIXTURE and plot them (Figure S9, and similar to Figure S2, which in the paper is elaborated with PONG https://github.com/ramachandran-lab/pong). It also include an R script to prepare for f3 run (from the Eigensoft package) and elaborate the results (Figure S11).

The folder ALDER includes commands to run MALDER and the R script to plot the results (Figure 5).





### Example data 
Some examples of elaborated data and plots are based on the datasets published in [Patterson et al. 2012](https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/2012_Patterson_AncientAdmixture_Genetics.pdf) 
and [Lazaridis et al. 2014](https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/2014_Nature_Lazaridis_EuropeThreeAncestries.pdf)


_____________________________________________
## DISCLAIMER and biography notes:
I am not a bioinformatician, and I am aware that my R scripting is sometimes not very "elegant". I learned bioinformatic pipelines and methods for SNP analysis when R alone proved to be not sufficient anymore for my genetic data. These commands were put together with the help of colleagues and friends and can be used by anyone new to this SNP platform and anyone who wants to find inspiration for simple data visualization in R. I am glad to improve this material with suggestions and corrections.

Contact me by [email](mailto:barbieri.chiara@gmail.com)



_________________________________________
Developed at University of Zurich - Department of Evolutionary Biology and Environmental Studies. 


![alt text](https://upload.wikimedia.org/wikipedia/de/thumb/8/89/Universit%C3%A4t_Z%C3%BCrich_logo.svg/200px-Universit%C3%A4t_Z%C3%BCrich_logo.svg.png)

_____________________________________



### Credits:
External help from Irina Pugach (Plink, Admixture, Beagle, Alder), Fabrizio Mafessoni (Bash, awk) and colleagues from the Max Planck Institute for Evolutionary Anthropology, Leipzig; Hiba Babiker (convertf, smartpca) and colleagues from the Max Planck Institute for the Science of Human History, Jena; Luca Pagani (plink tips) from University of Padova. 

