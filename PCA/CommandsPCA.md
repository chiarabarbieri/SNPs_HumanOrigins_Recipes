
# Running PCA with smartPCA


Running Principal Component Analysis with smartPCA of the Eigensoft package is very easy.
You will use the command smartpca and the relative [parfile](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/parfiles/par.SMARTPCA_global).

```
smartpca -p par.SMARTPCA_global
```

In the parfile i specify the location of a poplist. The PCA will run over the populations selected, but all the individuals of the initial files will be projected onto them. For visualization in R, you can decide which populations to plot (only the population of the poplist file, or all the populations of the initial database).

The [documentation available](https://github.com/DReichLab/EIG/tree/master/POPGEN) is exhaustive, especially for all the functions one can specify in the parfile.

Here some extracts:

smartpca runs Principal Components Analysis on input genotype data and outputs principal components (eigenvectors) and eigenvalues.  Eigenvalue_k/(Sum of eigenvalues) is the proportion of variance explained by eigenvector_k.  The method assumes that samples are unrelated - however, a small number of cryptically related individuals is usually not a problem in practice as they will typically be discarded as outliers.

* *poplistname*:   If wishing to infer eigenvectors using only individuals from a subset of populations, and then project individuals from all populations onto those eigenvectors, this input file contains a list of population names, one population name per line, which will be used to infer eigenvectors. It is assumed that the population of each individual is specified in the indiv file.
* *fstonly*: If set to YES, then skip PCA and just calculate FST values.
* *lsqproject*:  If set to YES, PCA projections is carried out by solving least squares equations rather than an orthogonal projection step.  This is approriate if PCs are calculated using samples with little missing data but it is desired to project samples with much missing data onto the top PCs (for example, ancient DNA samples).
* *numoutevec*:     number of eigenvectors to output. I chose 3.


The output file will be one .evec file (eigenvectors) and one .eval file (eigenvalues).

Go to R and plot your results with the [scripts](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/PCA/PlottingPCA.r)

________________________________________

Some examples:

![alt text](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/PCA/FIGURE1.png)

more beautiful plot, but less control over colors and symbols

![alt text](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/PCA/FIGURE2.png)


