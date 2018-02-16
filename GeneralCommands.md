
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

* 1_565596        1        0.006200        565596        G        A
* 1_567137        1        0.006200        567137        C        T
    
Manually change the .ind file to include the populations. With awk, again. I create a separate list of populations per each individual (in one column). I call the file "addStringPopsToIndFile.txt"

```
awk 'FNR==NR{a[NR]=$1;next}{$3=a[FNR]}1' addStringPopsToIndFile.txt yourfileEigenformat.ind > yourfileEigenformat_Familypop_EIGEN.ind
```
And now merge, with appropriate [parfile](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/parfiles/par.MERGE). Your dataset is the yourfileEigenformat, to be merged with a set called downloadedDatasetToMergeEigenformat, let's place them all in the folder /MERGE

```
mergeit -p par.MERGE
```
   
 -----------------------------#### "par.convertf_subsetsamples2" 
 ## paramfile for convertf
 
genotypename: /projects1/users/barbieri/MERGE/SouthAme_PlusALLReich.merged.geno
snpname: /projects1/users/barbieri/MERGE/SouthAme_PlusALLReich.merged.simple.snp
indivname: /projects1/users/barbieri/MERGE/SouthAme_PlusALLReich.merged.ind

outputformat:   EIGENSTRAT

genotypeoutname: /projects1/users/barbieri/PCA/CompAndMyAme/set2Comp.geno
snpoutname: /projects1/users/barbieri/PCA/CompAndMyAme/set2Comp.snp
indivoutname: /projects1/users/barbieri/PCA/CompAndMyAme/set2Comp.ind

poplistname: /projects1/users/barbieri/PCA/CompAndMyAme/listpopRefAndMyAme.txt


---------------------

597569 SNPs left in the set1 with ancient

--- now SMARTPCA    with lsqproject: YES, and selected list of pops set1 (ancient)

        scp par.SMARTPCA_global cdag1.cdag.shh.mpg.de:/projects1/users/barbieri/


 smartpca -p par.SMARTPCA_global
 
 --------------------------### par.SMARTPCA_global
 ## paramfile for smartpca
genotypename: /projects1/users/barbieri/MERGE/SouthAme_PlusALLReich.merged.geno
snpname: /projects1/users/barbieri/MERGE/SouthAme_PlusALLReich.merged.simple.snp
indivname: /projects1/users/barbieri/MERGE/SouthAme_PlusALLReich.merged.ind


evecoutname:     /projects1/users/barbieri/PCA/EIGENVECTOROUT_SouthAme_PlusALLReich_set1.evec
evaloutname:     /projects1/users/barbieri/PCA/EIGENVECTOROUT_SouthAme_PlusALLReich_set1.eval
altnormstyle:    NO
familynames:     NO
grmoutname:      grmjunk
lsqproject: YES
#outliermode: 2
numoutlieriter: 0
numoutevec: 3
poplistname: /projects1/users/barbieri/PCA/AmeModernAncient/poplistAmeAndAncient.txt
_____________________________________________________________________________________


and smartpca local with set 1

 smartpca -p par.SMARTPCA_local1


//////////////////////////////////////////
relatives check
plink --file SouthAmericaATLAS_180ExcludeRelatives_plink --genome
plink --file SouthAmericaATLAS_180ExcludeRelatives --genome --allow-no-sex --maf 0.05

///////////////////////////////////////////////
and go local
scp -r cdag1.cdag.shh.mpg.de:/projects1/users/barbieri/PCA/* 









make ped files so i can run Admixture

convertf -p par.PEDmakefiles


------------------------------ ## paramfile for convertf into ped

genotypename: /projects1/users/barbieri/PCA/CompAndMyAme/set2Comp.geno
snpname: /projects1/users/barbieri/PCA/CompAndMyAme/set2Comp.snp
indivname: /projects1/users/barbieri/PCA/CompAndMyAme/set2Comp.ind

outputformat:   PED

genotypeoutname: /projects1/users/barbieri/PED.EIGENSTRAT/set2CompPED.ped
snpoutname: /projects1/users/barbieri/PED.EIGENSTRAT/set2CompPED.pedsnp
indivoutname: /projects1/users/barbieri/PED.EIGENSTRAT/set2CompPED.pedind
familynames: NO
___________________________________________________________________________




=========================== random commands ==============================================


grep -f -v /projects1/users/barbieri/PCA/AmeModernAncient/poplistAmeAndAncient.txt <file>
grep -o /projects1/users/barbieri/PCA/AmeModernAncient/poplistAmeAndAncient.txt indmerge.txt | wc 

    scp poplistAmeAndAncient.txt cdag1.cdag.shh.mpg.de:/projects1/users/barbieri/PCA/AmeModernAncient
    
scp listpopRefAndMyAme.txt cdag1.cdag.shh.mpg.de:/projects1/users/barbieri/PCA/CompAndMyAme

awk '{ print $3 }' /projects1/users/barbieri/MERGE/SouthAme_PlusALLReich.merged.ind > indmerge.txt

cat indmerge.txt | grep -Fvf /projects1/users/barbieri/PCA/AmeModernAncient/poplistAmeAndAncient.txt | grep '^Q'

cant open file /projects1/users/barbieri/PCA/AmeModernAncient/AmeModernAncientEIGEN.pedsnp of type r
