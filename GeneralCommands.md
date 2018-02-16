
# Prepare files for subsequent analysis

## Format conversion
My new genotype file came in plink format. I have two files: *yourfile.ped* and *yourfile.map* I will convert them to eigenstrat format to merge it with other public dataset. Some examples here https://reich.hms.harvard.edu/datasets

```
convertf -p par.EIGENSTRATmakefiles
```
Note: the parfile used for the Eigensoft package are stored in the [parfiles](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/tree/master/parfiles) folder


Explore your old and new files with a quick look with bash command

```
less yourpath/Eigenstrat/yourfileEigenformat.geno
```

## Explore your dataset for errors and samples to be removed

A good idea is to perform some exploratory analysis before you start with the real analysis.

### Count the number of missing call per sample

```
cat yourfile.ped |  awk '{co=0;for (i=6;i<=NF;i++){if ($i=="0"){co=co+1}};print co/2}' > missing_data.txt
```

### Check for relatives
You might have duplicate samples or close relatives that escaped your sample filtering in the lab. Close relatives impact some of your future population analysis. A function of plink checks for this: we use the flag --genome with a threshold for minimum allele frequency of 0.05.

```
plink --file yourfile --genome --allow-no-sex --maf 0.05
```
--genome invokes an IBS/IBD computation, and then writes a report to *plink.genome*.
We are interested in the value *PI_HAT* :   Proportion IBD, i.e. P(IBD=2) + 0.5*P(IBD=1)

Now we run another plink command to explore F, the degree of consanguineity, and eventually delete outliers with a very high F.

```
plink --file yourfile --het
```
--het computes observed and expected autosomal homozygous genotype counts for each sample, and reports method-of-moments F coefficient estimates. It writes a report to *plink.het*


go to R to visualize the results of the two plink analysis with this [script](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/visualizeGeneralCommand.r)

### Check for outliers with a PCA

We search for obvious outliers with a quick PCA analysis --> see PCA folder and [instructions](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/PCA/CommandsPCA.md)

______________________________

## Merge with other dataset

Copy your Eigen files in a new directory MERGE
```
Mkdir MERGE
cp Eigenstrat/yourfileEigenformat* MERGE/
```

Check that the snp files of your dataset and of the dataset you want to merge have the same format.
For example, I had to change my .snp file to adapt to the format of the file I wanted to merge with.
A simple awk command to merge and reorder the columns of my .snp file

```
awk '{print $2 "_" $4 "    " $2 "   " $3"    " $4"   " $5 "   "$6}' yourfileEigenformat.snp > yourfileEigenformat.simple.snp
```

I obtain this kind of .snp file, where column $1 is the name of the snp, $2 is the chromosome, $3 is the position of the SNP in centimorgans on the chromosome, $3 is the position in bases, $4 and $5 are the reference and variant alleles.

* 1_565596   1        0.006200        565596        G        A
* 1_567137        1        0.006200        567137        C        T

Manually change the .ind file to include the populations. With awk, again. I create a separate list of populations per each individual (in one column). I call the file "addStringPopsToIndFile.txt". i want to replace the column $3 of the .ind file with the column $1 (the only column) of my "addStringPopsToIndFile.txt".

```
awk 'FNR==NR{a[NR]=$1;next}{$3=a[FNR]}1' addStringPopsToIndFile.txt yourfileEigenformat.ind > yourfileEigenformat_Familypop_EIGEN.ind
```
And now merge, with the eigensoft command *mergeit* and appropriate [parfile](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/parfiles/par.MERGE). Your dataset is the yourfileEigenformat, to be merged with a set called downloadedDatasetToMergeEigenformat, let's place them all in the folder /MERGE

```
mergeit -p par.MERGE
```
______________________________


## Subset a list of populations from your merged dataset

After creating the large merged dataset *MyDataset.merged*, it might be convenient to run some of the next analysis only with a subset of populations, for the purpose of the analysis or simply to speed up the calculations. I create a file *poplistSelection.txt* with only one column including a list of populations I want to extract from the original dataset. I make the [par file](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/parfiles/par.convertf_subsetsamples) to keep the Eigenstrat format (i use this format to run smartPCA) but i could also convert to PED plink format (I use this format to run Admixture)

```
convertf -p par.convertf_subsetsamples
```

Note: after running convertf, always check the output files to see if the format, number of individuals and number of SNPs have been preserved during the conversion. Make your tricks in bash or R. some examples

```
wc - l yourfile  #count the rows of the file
cmp --silent yourfile.simple.snp yourfileAfterConvertf.simple.snp || echo "files are different"
```
______________________________



