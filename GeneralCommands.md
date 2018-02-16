
# Prepare files for subsequent analysis

## format conversion
My new genotype file came in plink format. I will convert them to eigenstrat format to merge it with other public dataset. Some examples here https://reich.hms.harvard.edu/datasets

```
convertf -p par.EIGENSTRATmakefiles
```
Note: the parfile used for the Eigensoft package are stored in the [parfiles](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/tree/master/parfiles) folder


Explore your old and new files with a quick look with bash command

```
less yourpath/Eigenstrat/yourfileEigenformat.geno
```

## explore your dataset for errors and samples to be removed

A good idea is to perform some exploratory analysis before you start with the real analysis.

We search for obvious outliers with a quick PCA analysis --> see PCA folder.


# Merge with other dataset
 
Copy your Eigen files in a new directory MERGE
```
Mkdir MERGE
cp Eigenstrat/yourfileEigenformat* MERGE/
```

Check that the snp files of your dataset and of the dataset you want to merge have the same format.
For example, I had to change my .snp file to adapt to the format of the file I wanted to merge with.
A simple awk command to merge and reorder the columns of my .snp file

```
   awk '{print $2 "_" $4 "        " $2 "        " $3"        " $4"        " $5 "        "$6}' yourfileEigenformat.snp > yourfileEigenformat.simple.snp
```

I obtain this kind of .snp file, where column $1 is the name of the snp, $2 is the chromosome, $3 is the position of the SNP in centimorgans on the chromosome, $3 is the position in bases, $4 and $5 are the reference and variant alleles.

* 1_565596   1        0.006200        565596        G        A
* 1_567137        1        0.006200        567137        C        T
    
Manually change the .ind file to include the populations. With awk, again. I create a separate list of populations per each individual (in one column). I call the file "addStringPopsToIndFile.txt"

```
awk 'FNR==NR{a[NR]=$1;next}{$3=a[FNR]}1' addStringPopsToIndFile.txt yourfileEigenformat.ind > yourfileEigenformat_Familypop_EIGEN.ind
```
And now merge, with appropriate [parfile](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/parfiles/par.MERGE). Your dataset is the yourfileEigenformat, to be merged with a set called downloadedDatasetToMergeEigenformat, let's place them all in the folder /MERGE

```
mergeit -p par.MERGE
```
______________________________

