### Generate haplogroup affiliation for the samples

# mtDNA: extract a .vcf file with the mitochondrial variants and run it in Haplogrep 
https://haplogrep.uibk.ac.at/
  
# Y chromosome: extract a .vcf file with the Y chromosome variants and run it with Yhaplo 
  https://github.com/23andMe/yhaplo

### uniparental visualization: pie charts on a map (Figure S8)

# info file with meta information for each population
info<-read.table("/Users/chiarabarbieri/Desktop/SNP_humOriginsAmericas/listpopInfo.txt", header=T, as.is=T,comment.char = "")

# put the haplogroup results in a table, each line corresponding to one individual, columns with sampleID, population, haplogroup affiliation mtDNA_hg Ych_hg
haplogroups<-read.table("uniparental.txt",header=T)

infoRED<-info[which(info$population%in%haplogroups$Population),]
mtfreq<-table(haplogroups$Population,haplogroups$mtDNA_hg)
write.table(mtfreq,"mtfreq.txt", sep="\t")
mtfreq<-read.table("mtfreq.txt", sep="\t", header=T)

mtfreq$size<-apply(mtfreq,1,sum)
mtfreq$population<-rownames(mtfreq)
haploY<-haplogroups[-which(haplogroups$sex=="F"),] # exclude females
Yfreq<-table(haploY$Population,haploY$Ych_hg)
Yfreq<-Yfreq[,-2]
write.table(Yfreq,"Yfreq", sep="\t")
Yfreq<-read.table("Yfreq", sep="\t", header=T)

Yfreq$sizeY<-apply(Yfreq,1,sum)
Yfreq$population<-rownames(Yfreq)

infomerg<-merge(infoRED,mtfreq)
infomerg<-merge(infomerg, Yfreq) # this object has the frequency of mt and Ych haplogroup per population, and the population sample size (different sample size for Ych and mt, Y ch has only males)
library(ggplot2)
library(scatterpie)
library(ggrepel)


infomerg$lat<-as.numeric(infomerg$lat)
infomerg$lon<-as.numeric(infomerg$lon)

columnshg<-c("A", "B", "C", "D", "others") #mtDNA

world <- map_data('world')
p <- ggplot(world, aes(long, lat)) +
  geom_map(data=world, map=world, aes(map_id=region), fill="white", colour="black", size=0.15)+
  coord_quickmap(ylim=c(-18,4), xlim=c(-90,-65))

p + geom_scatterpie(data=infomerg, aes(x=lon, y=lat, group=population, r=size/10),
                     cols=columnshg, color="gray20", alpha=.8) +
  geom_scatterpie_legend(infomerg$size/10, x=-87, y=-15, n=3)+
  geom_label_repel(data=infomerg, aes(x=lon, y=lat,label=shortPop), size=2, point.padding = 0.5)

ggsave("MapPieChart_mtDNA_all2.pdf", useDingbats=FALSE)

# repeat the above plot after specifying the Y chromosome columns
columnshg<-c("Q","C_", "othersY") #Y chromosome

ggsave("MapPieChart_Ych_all.pdf", useDingbats=FALSE)


#---------------------------------------
### proportion of non-Native ancestry in autosomal, mt and Ych (Figure S10)
#---------------------------------------
admixprop<-read.table("admixtureProportionsInd.txt", header=T, as.is=T, comment.char = "") # this file comes from the supervised ADMIXTURE run elaborated, it looks like this

----------
sample_ID	pop1 Africa	Europe	Native
HGDP00702	Piapoco	0.00001	0.00001	0.99998
HGDP00704	Piapoco	0.00001	0.00001	0.99998
HGDP00706	Piapoco	0.00001	0.00001	0.99998
HGDP00920	Yoruba	0.99998	0.00001	0.00001
HGDP00924	Yoruba	0.99998	0.00001	0.00001
---------

# mtDNA proportion per population
femaleNonNative<-c()
for (i in 1:nrow(mtfreq)){
  femaleNonNative[i]<-(mtfreq$others[i]/sum(mtfreq[i,]))
}

# Ych proportion per population

maleNonNative<-c()
for (i in 1:nrow(Yfreq)){
  maleNonNative[i]<-(Yfreq$othersY[i]/sum(Yfreq[i,]))
}

# Autosome mean proportion per population

autosomeNonNative<-c()
for (i in 1:nrow(Yfreq)){
  temp<-autsomeAdmix[which(autsomeAdmix$population==rownames(mtfreq)[i]),]
  autosomeNonNative[i]<-mean(temp$Africa+temp$Europe)
}


tabNonNativeProportion<-cbind(femaleNonNative,maleNonNative,autosomeNonNative)
rownames(tabNonNativeProportion)<-rownames(mtfreq)
library(reshape)
meltNativeProportion<-melt(tabNonNativeProportion) # melt dataset for ggplot
colnames(meltNativeProportion)<-c("population","geneticMarker","proportionNonNativeAncestry")

meltNativeProportion$alphaexternal<-c(rep(0.7,25),rep(0.5,25),rep(0.9,25))

info<-read.table("listpopInfo.txt", header=T, as.is=T,comment.char = "") #extra info, one line per population. includes the grouping scheme.
infoRED<-info[which(info$population%in%rownames(tabNonNativeProportion)),]
poporder<-as.character(infoRED$population)


coso<-merge(meltNativeProportion,infoRED, by = "population",  all.x = TRUE)


#trick to order the color patterns according to group populations belong to
nomibelli<-as.character(unique(coso$group))
riordinacolori<-as.character(unique(coso$group))
listpopInfo2<-infoRED
listpopInfo2<-listpopInfo2[!duplicated(listpopInfo2$group),]
rownames(listpopInfo2)<-listpopInfo2$group
listpopInfo2<-listpopInfo2[riordinacolori,]
colorelli<-as.character(listpopInfo2$color)


ggplot(coso, aes(geneticMarker, proportionNonNativeAncestry, fill=group,color=group,alpha=alphaexternal)) +   
    geom_bar( stat="identity") +    
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  theme(text = element_text(size=10))+
  scale_fill_manual(values = colorelli,
                    breaks=nomibelli,
                    name="Source Geographic Group")+
  scale_color_manual(values = colorelli,
                     breaks=nomibelli,
                     name="Source Geographic Group")+
  facet_wrap(~population)
ggsave("NONNativeAncestryProportion_3markers.pdf", useDingbats=FALSE)
