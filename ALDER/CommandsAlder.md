
# Running Alder to test for Admixture times

Alder and Malder are software designed to detect admixture events from parental populations into a target population, and date those events. They use Linkage Disequilibrium decay.
Malder allows for muldiple admixture test.
ALDER stands for Admixture Linkage Disequilibrium for Evolutionary Relationships

[Original Publication](http://www.genetics.org/content/193/4/1233)



Download from https://github.com/joepickrell/malder/tree/master/MALDER

syntax like
```
malder -p par.MALDER | tee NameofTargetPop.logfile
```

With this [parfile](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/parfiles/par.MALDER) i will ask Alder to take one population from MyDataset.merged, and test for admixture with three parental populations chosen. In this case, Spanish;Yoruba;Xavante are chosen as close enough to the putative parent populations, if and when the admixture have occurred with a source from Europe, Africa or South America.

Now imagine that I want to run Malder on a list of target populations, with the same setting of three potential parents . I have to create different parfiles, and run Malder on each of them. Here's a bash script to automatize part of the process - with a different [parfile](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/parfiles/par.MALDER_Xavante):

```
for target in  Pima  Mixe Aymara Quechua Piapoco Guarani  Karitiana   ; do
sed s/TARGET/${target}/g par.MALDER_Xavante >  par.MALDER_Xavante_${target}
malder -p par.MALDER_Xavante_${target} | tee M_Xavante_${target}.logfile
grep 'success' M_Xavante_${target}.logfile > ${target}_Xavante_success_out.txt
done

```
I chose a list of target population to make a concrete example: we are testing admixture in American populations between the three sources described above. Some of them might not be admixed with the parents chosen: the run won't be marked as "success". The loop creates a file with the successfull values for each target population.

```
cat *Xavante_success_out.txt > outputSuccessMalder.cat.txt
```

I put together all the relevant information for each Target population that passed the admixture test. This includes the set of target and two parents, the p-value, the Z score, the time of the admixture (in generations ago) with confidence interval.

## LD decay values for each target population

Take all the LD decay output files (they are called *Targetpopulation_LDout*) and add one extra column with the name of the target population (so you don't confuse the values after you merge them all in a single file)

```
for target in  Pima  Mixe Aymara Quechua Piapoco Guarani  Karitiana   ; do
awk '$5+=${target}' ${target}_LDout > ${target}_newLDout
done
```

Then make a single file with all the LD decays
```
cat *_newLDout > LD_decaycurveALL.txt
```

Exclude the rows that begin with #
```
grep -v "#" LD_decaycurveALL.txt > LD_decaycurveALLfilter.txt
```

________________________
Now to R to plot the results of the file outputSuccessMalder.cat.txt. and LD_decaycurveALLfilter.txt with this [script](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/ALDER/PlottingALDER.R)

Examples of figures:


Admixture Times from Africa and Europe

![alt text](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/ALDER/Alder1.png)

Linkage Disequilibrium Decay in three populations, sources from Africa and Europe

![alt text](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/ALDER/Alder2.png)




