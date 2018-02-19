# SNPs HumanOrigins Recipes
## Commands and R scripts to analyse SNP data.

The folder contains different sections for most used analysis methods to perform on SNP chip databases. It is particularly tailored to the use of the Human Origins Affymetrix array, and for analysis of human population history.

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
Start with the General Commands for data overview , and then check each folder for instruction on how perform the analysis proposed. Each folder contains a command .md file and in most case a dedicated R script for visualization purposes.

### Example data 
Some example of elaborated data and plots are based on the dataset published in [Patterson et al. 2012](https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/2012_Patterson_AncientAdmixture_Genetics.pdf
and [Lazaridis et al. 2014](https://reich.hms.harvard.edu/sites/reich.hms.harvard.edu/files/inline-files/2014_Nature_Lazaridis_EuropeThreeAncestries.pdf)


_____________________________________________
## DISCLAIMER and biography notes:
I am not a bioinformatician, I am a molecular anthropologist. I recently started to explore and learn these methods when R alone proved to be not sufficient anymore for my genetic analysis. These commands were put together with the help of colleagues and friends and can be used by colleagues who are new to this SNP platform and who want to find inspiration for data visualization in R. I am looking forward to improve this material with suggestions and corrections.

Contact me with [email](mailto:barbieri.chiara@gmail.com)



_________________________________________
Developed at University of Zurich - Department of Evolutionary Biology and Environmental Studies. 


![alt text](https://upload.wikimedia.org/wikipedia/de/thumb/8/89/Universit%C3%A4t_Z%C3%BCrich_logo.svg/200px-Universit%C3%A4t_Z%C3%BCrich_logo.svg.png)



