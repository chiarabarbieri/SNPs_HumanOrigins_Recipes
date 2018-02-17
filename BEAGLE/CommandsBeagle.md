# Phasing and estimation of IBD blocks with Beagle

Beagle is a software that performs several functions, including a good and fast algorithm for phasing chromosome and a method to reconstruct blocks shared by identity by descent (IBD).
Download from https://faculty.washington.edu/browning/beagle/beagle.html.  Beagle runs with java (install java if needed). The [documentation](https://faculty.washington.edu/browning/beagle/beagle_4.1_21Jan17.pdf) is basic.

Some commands are needed to prepare the data from the original plink format.
The phasing can be performed with the whole merged dataset, or with a subset of representative populations, to save computational time. To select a subset of individuals from target populations follow the instructions with convertf in the [General commands](https://github.com/chiarabarbieri/SNPs_HumanOrigins_Recipes/blob/master/GeneralCommands.md)

Beagle uses Variant Call Format (VCF) 4.2 for input and output file. The file must be gzipped (.gz format). I use *plink* to convert to .vcf format.

```
plink --file  yourfileplinkformat  --allow-no-sex --recode vcf-iid --alleleACGT --out  yourfileBeagleVCF
```

I use a quick *bcftools* command to remove invariant sites.

```
bcftools view -v snps -o yourfileBeagleVCF_noinvariant yourfileBeagleVCF.vcf
```

I gzip the file and i use *tabix* (http://www.htslib.org/doc/tabix.html) to index it.
```
bgzip yourfileBeagleVCF_noinvariant
tabix -p vcf yourfileBeagleVCF_noinvariant.gz
```
To speed up the beagle run, I will split the whole vcf file into single chromosome files. I create a loop to perform the operation for each chromosome.

```
for chr in `bcftools view -h yourfileBeagleVCF_noinvariant.gz  | perl -ne '
if (/^##contig=<ID=([^,]+)/) { if (length($1)<=2) {
print "$1\n"
} }'`; do
bcftools view -Oz -r $chr yourfileBeagleVCF_noinvariant.gz  > splitted_$chr.vcf.gz &
done
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
for chromosome in {2..22}; do
awk '{$3 = $3 * 100; print}' referenceFromMyData_Chr${chromosome}.map > referenceFromMyDataCentimorgan_Chr${chromosome}.map
done
```

## Ready to run beagle

And now i create another  loop to run beagle on each chromosome. Note: i do not consider sex chromosomes and the mitochondrial DNA. I will run from chromosome 1 to 22.

```
for chromosome in {1..22}; do
java  -Xss5m  Xmx4g -jar beagle.08Jun17.d8b.jar gt=splitted_${chromosome}.vcf.gz  ibd=true impute=false map=MapFromMyData/referenceFromMyData_Chr${chromosome}.map out=BeaglePhased${chromosome}
done

```


*beagle.08Jun17.d8b.jar* is the latesed version of Beagle i downloaded. I call java to run the .jar file of Beagle from my folder.
Xmx4g: 4 is the maximum permitted size of the memory pool in gigabytes (that I chose).
With *ibd=true* i ask beagle to produce the IBD and HBD blocks as outputs for each chromosome.
At the end of the run, I will have for each chromosome a logfile with the specifications of the run, a phased vcf file, a file containing all the IBD blocks with the reference between the two samples which share each block, start and end point, and a measure of robustness, and a HBD file similar to the IBD file.

The beagle run might take between one hour and a whole day.

The resulting phased .vcf file can be used for several analysis which require phased data - like many analysis that work on admixture and ancestry block detection.

The IBD block sharing data can be directly analyzed to reveal recent layers of contact, a concept similarly to a haplotype sharing.
[This script] in R provides some examples for data visualization: sharing within population, between populations, on a matrix and on a map.




