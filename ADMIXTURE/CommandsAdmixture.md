# The "classic" Admixture
Identifying ancestry components and their proportion  for each individual in a set of populations.

*note*: commands adapted from an original script by Irina Pugach.
see [Pugach et al. 2018](https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msx333/4782510) for more inputs on Admixture analysis, and other analysis with SNP chip datasets.

ADMIXTURE is a software that works as Structure, but with faster computation. Download from https://www.genetics.ucla.edu/software/admixture/. It takes plink format input files.


___________________________

We run ADMIXTURE on the autosomal chromosomes (from 1 to 22) to avoid sex biased effects.

First, prune the dataset for LD with Plink (settings which define window size, step and the r2 threshold described in the publication above).
```
plink --file yourfile --indep-pairwise 200 25 0.4 --out x.tmp
plink --file yourfile --extract x.tmp.prune.in --recode12 --out yourfile.pruned
```

Now the proper ADMIXTURE run. The following commands will run ADMIXTURE with 10 runs for each *K* (number of ancestry blocks) desired. The initial exploratory tests can be done with less runs per more *K*, to explore the diversity of many Ks, and then reduce them to a meaningful number for the analysis. There will be one *K* more supported by the analysis: this is the setup with the best representation of the structure in the data.

```
typeset -i run=0
while (( run < 10 )); do  ##  with 10 runs for each K
run=$(( run + 1 ));
for K in 4 5 6 7 8 9 10; do  # select a meaningful series of K - the more Ks, the longer the run obviously
admixture -s time --cv yourfile.pruned.ped $K -j6 | tee log.K${K}.RUN$run.out;
mv yourfile.pruned.$K.P K$K.Run$run.P;
mv yourfile.pruned.$K.Q K$K.Run$run.Q;
done;
done
```
For each run there are three output: .out, .P, and .Q

Now we will elaborate the outputs with a mix of bash commands and *R*.

Cross-validation error to identify the best value of *K*

```
grep -h CV log*out > CV.txt
```
plot the distribution of values associated to each K in R.

![alt text](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/ADMIXTURE/amixture1.png)

For each K, determine which of the 10 runs has the highest likelihood.

```
for Kvalue in {4..10}; do
for Try in {1..10}; do grep -h "^Loglikelihood:" log.K$Kvalue.RUN$Try.out; done  > LL.K$Kvalue.txt;
done

```


Elaborate this likelihood file in R, and visualize the ADMIXTURE results: follow this [script](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/ADMIXTURE/plotting_Admixture.r)

Example of admixture for a set of population with K=7 (seven ancestry components)

![alt text](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/ADMIXTURE/admixture2.png)

