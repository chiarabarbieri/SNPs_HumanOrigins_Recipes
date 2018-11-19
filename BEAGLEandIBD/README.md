

# BEAGLE: phasing the data

Beagle is a software that performs several functions, including a good and fast algorithm for phasing chromosomes.
wget "https://faculty.washington.edu/browning/beagle/beagle.03Jul18.40b.jar"

Some commands are needed to prepare the data from the original plink format.
Beagle uses Variant Call Format (VCF) 4.2 for input and output file. The file must be gzipped (.gz format). I use *plink* to convert to .vcf format.

```
plink --bfile  yourinput  --allow-no-sex --recode vcf-iid --alleleACGT --out  yourinputVCF
```


I use a quick *bcftools* command to remove invariant sites. then check how many SNPs are left.

```
bcftools view -v snps -o yourinputVCF_noinvariant yourinputVCF.vcf
 wc -l yourinputVCF_noinvariant
```

I gzip the file and i use *tabix* (http://www.htslib.org/doc/tabix.html) to index it.

```
bgzip yourinputVCF_noinvariant
tabix -p vcf yourinputVCF_noinvariant.gz
```
To speed up the beagle run, I will split the whole vcf file into single chromosome files. I create a loop to perform the operation for each chromosome.
```

for chr in `bcftools view -h yourinputVCF_noinvariant.gz  | perl -ne '
 if (/^##contig=<ID=([^,]+)/) { if (length($1)<=2) {
   print "$1\n"
 } }'`; do
  bcftools view -Oz -r $chr yourinputVCF_noinvariant.gz  > splitted_$chr.vcf.gz &
done
```

Create a directory to run Beagle, and check how the files look like.

```
mkdir BEAGLE_files
mv *.vcf.gz BEAGLE_files
cd BEAGLE_files

zcat splitted_4.vcf.gz | less
```

## Provide a map file for accurate positioning of the distance between SNPs.
Almost ready. But before running beagle, I must setup a map file with the physical positions (in centimorgans) for all the SNPs I have in my vcf file. Usually this information is stored in the initial plink file before vcf conversion, specifically in the .map file. I also need it splitted into chromosomes, like the vcf input files i am using. I am using the .map [files indicated by beagle](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/) as a template.

I do the first steps in R:

-------------------------
### R

mapBeagle<-yourfileplinkformat.map
mapBeagle[,2]<-"."   # little tweak to add a dot between values.
mapBeagle<-mapBeagle[,1:4]

for (i in 1:24){
temp<-mapBeagle[which(mapBeagle[,1]==i),]
write.table(temp, paste(c("referenceFromMyData_Chr",i,".map"), collapse=""), sep="\t",  col.names = F, row.names = F)
}

----------------------------

now back to the shell for a final checkup.
I notice i have the map files in morgans, but beagle wants them in centimorgan. I multiply the column with the physical distance on the chromosome for 100.

```
for chromosome in {1..22}; do
awk '{$3 = $3 * 100; print}' referenceFromMyData_Chr${chromosome}.map > referenceFromMyDataCentimorgan_Chr${chromosome}.map
done
```

## Ready to run beagle

And now i create another  loop to run beagle on each chromosome. Note: i do not consider sex chromosomes and the mitochondrial DNA. I will run from chromosome 1 to 22.
This is repeated for 3 runs to later check for consistency.

```
wget "https://faculty.washington.edu/browning/beagle/beagle.03Jul18.40b.jar"

for run in {1..3}; do
for chromosome in {1..22}; do
	seed=$RANDOM
    java  -Xss5m  -Xmx4g -jar beagle.03Jul18.40b.jar gt=splitted_${chromosome}.vcf.gz  map=referenceFromMyDataCentimorgan_Chr${chromosome}.map window=20 seed=$seed out=RUN${run}_Beagle5Phased${chromosome}
done
done
```
*beagle.03Jul18.40b.jar* is the latesed version of Beagle i downloaded. I call java to run the .jar file of Beagle from my folder.
Xmx4g: 4 is the maximum permitted size of the memory pool in gigabytes (that I chose).
At the end of the run, I will have for each chromosome a logfile with the specifications of the run and a phased vcf file.

The resulting phased .vcf file can be used for several analysis which require phased data - like many analysis that work on admixture and ancestry block detection.




# REFINED IBD

This software works on the phased dataset to infer IBD (Identity by Descent) block sharing. This can be associated to recent layers of contact, a concept similarly to a haplotype sharing.
Three runs are performed on the previous three phasing runs from Beagle.

```
wget "http://faculty.washington.edu/browning/refined-ibd/refined-ibd.12Jul18.a0b.jar"

for run in {1..3}; do
for chromosome in {1..22}; do
java -Xss5m  -Xmx4g -jar refined-ibd.12Jul18.a0b.jar gt=RUN${run}_Beagle5Phased${chromosome}.vcf.gz map=/home/chiara_barbieri/Chad_personalAncestry/MAP/referenceFromMyData_Chr${chromosome}.map window=20 trim=0.3 out=RefinedIBD_RUN${run}.$chromosome
done
done
```

Merge and remove gaps with their utility *merge-ibd-segments.12Jul18.a0b.jar*. 
After unzipping the ibd results of the three runs, we are going to put together all the fragments from the three runs.
Then for each chromosome gaps between contiguous blocks are removed (resulting in one larger block) and only the blocks that are present in all three runs are kept.

In the final step, the results from each chromosome are merged in a single file.

```
wget "https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.12Jul18.a0b.jar"


gunzip *.ibd.gz

for chromosome in {1..22}; do
cat RefinedIBD_RUN1.${chromosome}.ibd RefinedIBD_RUN2.${chromosome}.ibd RefinedIBD_RUN3.${chromosome}.ibd > Chr${chromosome}.123.ibd
cat Chr${chromosome}.123.ibd | java -jar merge-ibd-segments.12Jul18.a0b.jar RUN2_Beagle5Phased${chromosome}.vcf.gz  /home/chiara_barbieri/GenotypeAtlasPlink_Andes_2plates/Beagle_newisBetter/MapFromMyData/referenceFromMyDataCentimorgan_Chr${chromosome}.map 2 1 > Chr${chromosome}.IBD.Merged
done


cat *.IBD.Merged > all.refinedIBD.Merged
```




## Visualization

The IBD block sharing data can be directly analyzed to reveal recent layers of contact, a concept similarly to a haplotype sharing.
[This script](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/BEAGLE/plotting_IBD_fromBeagle.r) in R provides some examples for data visualization: sharing within population, between populations, on a matrix and on a map.

**Section 1** will show you how to plot a matrix of IBD block sharing between populations and a length and LOD score for within population sharing.

examples:

![alt text](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/BEAGLE/1.1.png)

![alt text](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/BEAGLE/1.2.png)


Section 2 will show you how to plot barplots with different sharing patterns for each populations: average length of blocks, number of populations each pop shares with, frequency of sharing per pop, total cumulative lenght of the IBD blocks shared.

examples:

![alt text](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/BEAGLE/2.png)


**Section 3** visualizes the distribution (boxplot) of sharing for each population against two targets: in this case, one population to represent the sharing with Africa, and one population to represent the sharing with Europe.

![alt text](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/BEAGLE/3.png)


**Section 4** will create a distance matrix (dissimilarity matrix) from the number of sharing events between populations, and will visualize population distances with a Multi Dimensional Scalind (MDS) plot.

**Section 5** will display two different ways of visualizing data sharing on a map: with ggplot, each population will have its own plot, with dots over the populations who share with it, and size proportional to the amount of sharing. With library Maps, the network of sharing between populations will be plotted simultaneously, with thicker lines connecting populations who share more IBD blocks with each other.

A simple map
![alt text](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/BEAGLE/4.1.png)


![alt text](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/BEAGLE/4.png)

![alt text](https://media.springernature.com/lw900/springer-static/image/art%3A10.1038%2Fs41598-017-17728-w/MediaObjects/41598_2017_17728_Fig5_HTML.jpg)
from Barbieri et al. 2017 Scientific Reports https://www.nature.com/articles/s41598-017-17728-w
