### Generate haplogroup affiliation for the samples

# mtDNA: extract a .vcf file with the mitochondrial variants and run it in Haplogrep 
https://haplogrep.uibk.ac.at/
  
# Y chromosome: extract a .vcf file with the Y chromosome variants and run it with Yhaplo 
  https://github.com/23andMe/yhaplo

### uniparental visualization

# info file with meta information for each population
info<-read.table("/Users/chiarabarbieri/Desktop/SNP_humOriginsAmericas/listpopInfo.txt", header=T, as.is=T,comment.char = "")

# put the haplogroup results in a table, each line corresponding to one individual.
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
infomerg<-merge(infomerg, Yfreq) # this object has the frequency of mt and Ych haplogroup per population, and the population sample size (different for Ych and mt)
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
