
####### in R

# #__________________________________________
# Visualize population locations on a map (Figure S1)
#__________________________________________

infoRed<-read.table("infoRed.toPlot.txt", header=T, as.is = T, comment.char = "") # the info file has information for each population, different groups, lat and lon coordinates
library(rworldmap)
library(mapproj)
library(ggmap)
library(ggplot2)

# trick for plotting the pop names with different colors according to the group they belong to

riordinacolori<-as.character(levels(as.factor(infoRed$group)))
nomibelli<-as.character(unique(infoRed$group))
infoRed2<-infoRed
infoRed2<-infoRed2[!duplicated(infoRed2$group),]
rownames(infoRed2)<-infoRed2$group
infoRed2<-infoRed2[riordinacolori,]
colorelli<-as.character(infoRed2$color)
#####

map.world <- map_data(map="world")

library(ggrepel)

gg <- ggplot()
gg <- gg + theme()
gg <- gg + geom_map(data=map.world, map=map.world, aes(map_id=region), fill="white", colour="black", size=0.15)

gg<- gg+coord_quickmap(ylim=c(-25,30), xlim=c(-111,-50)) # fix a window that includes the region of your interest

 gg + 
  geom_label_repel(data=infoRed, aes(x=lon, y=lat,label=population, color=group), size=2, label.padding=0.1)+
   geom_point(data=infoRed, aes(x=lon, y=lat, color=group), size=0.5 ) +
    scale_color_manual(values = colorelli,
                     breaks=nomibelli,
                     name="Source Geographic Group")+    
   theme(text = element_text(size=5))

 gg
ggsave("MapDotscolorelli.pdf", useDingbats=FALSE)



#__________________________________________
### search for relatives with plink --genome   (data screening)
#__________________________________________

filone<-read.table("plink.genome", header=T)
pdf("distributionPIHAT.pdf")
hist(filone$PI_HAT, breaks=100)
dev.off()

filone[which(filone$PI_HAT>0.25),]->Relat  # set the threshold to spot relatives. Second degree relative correspond to a cousin.

write.table(Relat, "couplesCousins.txt", sep="\t") # open this table to manually screen potential relative pairs

library(ggplot2)

# plot the distribution of pi hat for family groups
# remember when you replaced the .ind file with the names of the populations, before running plink? 
# FID1:	Family ID for first sample should correspond to the population

p <- ggplot(filone, aes(FID1,PI_HAT))
p + geom_boxplot()
dev2bitmap("boxplotDistribuPiHAT.pdf", type="pdfwrite")


## now something simila for the F value, a measure of heterozygozity and consanguineity 
het<-read.table("newHetCalculationsBeagleSet.het", header=T)
hetinfo<-merge(het,infoID[,c(1,14,15,16)],by.x = "FID",by.y = "population" ) # add info on population and group for each individual


# the next lines are for controlling the colors for populations and groups
riordinacolori<-as.character(levels(as.factor(hetinfo$group)))
listpopInfo2<-infoRED
listpopInfo2<-listpopInfo2[!duplicated(listpopInfo2$group),]
rownames(listpopInfo2)<-listpopInfo2$group
listpopInfo2<-listpopInfo2[riordinacolori,]
colorelli<-as.character(listpopInfo2$color)

# plot the distribution of F (consanguineity)  for family groups (Fig S6A)

p <- ggplot(hetinfo, aes(x=FID,y=F))+
  geom_boxplot(aes(y =  F, fill = group, alpha=0.7,color=group)) +
  labs(x = "", y = "F distribution",
       title = "distribution of F (consanguineity) within population")+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) +
  scale_x_discrete(limits=poporder)+
  scale_fill_manual(values = colorelli,
                    breaks=nomibelli,
                    name="Source Geographic Group")+
  scale_color_manual(values = colorelli,
                     breaks=nomibelli,
                     name="Source Geographic Group")+    theme(text = element_text(size=10))

p
ggsave("distributionF_boxplot.pdf", useDingbats=FALSE)


#__________________________________________
##### correlation between percentage of native ancestry and consanguineity (Figure S6B)
#__________________________________________

het<-read.table("newHetCalculationsBeagleSet.het", header=T)
hetinfo<-merge(het,infoID[,c(1,2,14,15,16)],by.x = "IID",by.y = "sample_ID" )
admixprop<-read.table("admixtureProportionsInd.txt", header=T, as.is=T, comment.char = "") # this file comes from the supervised ADMIXTURE run elaborated, it looks like this

----------
sample_ID	pop1 Africa	Europe	Native
HGDP00702	Piapoco	0.00001	0.00001	0.99998
HGDP00704	Piapoco	0.00001	0.00001	0.99998
HGDP00706	Piapoco	0.00001	0.00001	0.99998
HGDP00920	Yoruba	0.99998	0.00001	0.00001
HGDP00924	Yoruba	0.99998	0.00001	0.00001
---------

hetinfo<-merge(hetinfo,admixprop[,c(1,3:7)],by.x = "IID",by.y = "sample_ID" )
hetinfoAme<-hetinfo[-which(hetinfo$FID%in%c("Spanish","Yoruba","Italian_North")),]

infoRED<-info[which(info$population%in%hetinfoAme$FID),]

pdf(file="correlateAdmixtureConsanguineity.pdf", pointsize=12,width=12, useDingbats=FALSE)
plot(x=hetinfoAme$F,y=hetinfoAme$Native,col=as.character(hetinfoAme$color),pch=as.numeric(hetinfoAme$pch),xlab="consanguineity (F)", ylab="proportion of Native American ancestry")
abline(lm(hetinfoAme$Native ~ hetinfoAme$F))
legend("bottomright", legend=infoRED$population,  pch=as.numeric(infoRED$pch), col=as.character(infoRED$color), text.col="gray", ncol=5,bty="n", cex=0.5)
dev.off()



#__________________________________________
### plot Runs of Homozygosity (ROH)
#__________________________________________



# use an external file with information on population, group, plotting color etc for each individual
roh<-read.table("yourfile_ROH.hom.indiv", header=T, as.is=T)
myinfo<-read.table("infoData.txt",header=T)
roh<-roh[which(roh$IID%in%myinfo$sample_ID),]
rownames(myinfo)<-myinfo$sample_ID
myinfoSort<-myinfo[roh$IID,]

# an external file with characteristic same as above but for each population (is for the legend)
popinfo<-read.table("listpopInfoPlay.txt", header=T)

pdf(file="ROH_plot_numberofROH_length", pointsize=12,width=12)
plot(x=roh[,5]/1000,y=roh[,4],col=as.character(myinfoSort$color),pch=as.numeric(myinfoSort$pch),xlab="total length in ROH (Mb)", ylab="number of ROH")
legend("bottomright", legend=popinfo$population,  pch=as.numeric(popinfo$pch), col=as.character(popinfo$color), text.col="gray", ncol=5,bty="n", cex=0.5)
dev.off()


# play with ROH blocks, use the bin categories from Pemberton et al. 2013 and the six bin categories used in the Taino's paper

roh2<-read.table("newHetCalculationsBeagleSet.hom", header=T, as.is=T)
colnames(roh2)[2]<-"sample_ID"
roh2<-roh2[which(roh2$sample_ID%in%myinfoSort$sample_ID),]
rohinfo<-merge(roh2, myinfoSort)
rohinfo$KB<-rohinfo$KB/1000
colnames(rohinfo)[9]<-"MB"
categories<-matrix(NA,6,2)
categories[,1]<-c(0.5, 1,2,4,8,16)
categories[,2]<-c( 1,2,4,8,16,max(rohinfo$MB))
sketchrownames<-apply(categories,1,paste, collapse="-")
indROH<-matrix(NA,nrow(myinfoSort),length(sketchrownames))
colnames(indROH)<-sketchrownames
rownames(indROH)<-myinfoSort$sample_ID
for (i in 1:nrow(indROH)){
  targetID<-row.names(indROH)[i]
  tempID<-rohinfo[which(rohinfo$sample_ID==targetID),]
  for (k in 1:length(sketchrownames)){
    tempIDk<-tempID[which(tempID$MB>=categories[k,1]&tempID$MB<=categories[k,2]),]
    indROH[i,k]<- sum(tempIDk$MB)
  }
}

indROHinfo<-cbind(indROH,myinfoSort)
indROHinfo<-indROHinfo[,-which(colnames(indROHinfo)%in%c("order","pch","Beagle"))]
library(reshape)
library(tidybayes)
groupp<-table(indROHinfo$group)
rohgroup<-matrix(NA,length(groupp),length(sketchrownames)*3+2)
rownames(rohgroup)<-names(groupp)
colnames(rohgroup)<-c(sketchrownames,1:6,101:106,"color","set")
for (i in 1:length(groupp)){
  targettemp<-names(groupp)[i]
  blocktemp<-indROHinfo[which(indROHinfo$group==targettemp),]
  rohgroup[i,19]<-blocktemp$color[1]
  rohgroup[i,20]<-blocktemp$set[1]
  
  for (k in 1: length(sketchrownames)){
    rohgroup[i,k]<-mean(blocktemp[,k])  #/(groupp)[i]
    rohgroup[i,k+6]<-min(blocktemp[,k]) 
    rohgroup[i,k+12]<-max(blocktemp[,k])
  }
}
rohgroup<-as.table(rohgroup)
rohgroup[8:10,20]<-rep("Published",3)

indROHinfoMELT<-melt(rohgroup[,1:6],varnames=c("group","bin"))
indROHinfoMELT$minn<-melt(rohgroup[,7:12])[,3]
indROHinfoMELT$maxx<-melt(rohgroup[,13:18])[,3]
indROHinfoMELT$sett<-rep(rohgroup[,20],6)
write.table(indROHinfoMELT,"meltROHbinsXgroups.txt",sep="\t")
indROHinfoMELT<-read.table("meltROHbinsXgroups.txt",header=T)



riordinacolori<-as.character(levels(indROHinfoMELT$group))
listpopInfo2<-rohgroup
listpopInfo2<-listpopInfo2[riordinacolori,]
colorelli<-as.character(listpopInfo2[,19])


  p <- ggplot(indROHinfoMELT, aes(x=bin,y=value,group=group))+
    geom_ribbon(data=indROHinfoMELT,aes(ymin = minn , ymax = maxx,  fill = group), alpha=0.2)+
    geom_line(aes(color=group) )+
    scale_x_discrete(limits=sketchrownames)+
    scale_fill_manual(values = colorelli,
                      breaks=riordinacolori,
                      name="Source Geographic Group")+
    scale_color_manual(values = colorelli,
                      breaks=riordinacolori,
                      name="Source Geographic Group")+
    theme_bw()   +
    labs(x = "ROH length category (Mb)",y="Total length of ROH (Mb) per individual")
 p +facet_wrap(~sett, ncol = 1)
  
 ggsave("ROHbinsXgroupTEST.pdf", useDingbats=FALSE)  # (Figure 3B)


### now consider each population, to be plotted separately over the ribbons of the groups
 indROHinfoRED<-indROHinfo[which(indROHinfo$set=="presentstudy"),]
 popp<-table(indROHinfoRED$population)
 rohpops<-matrix(NA,length(popp),length(sketchrownames)*3+2)
 rownames(rohpops)<-names(popp)
 colnames(rohpops)<-c(sketchrownames,1:6,101:106,"color","set")
 for (i in 1:length(popp)){
   targettemp<-names(popp)[i]
   blocktemp<-indROHinfoRED[which(indROHinfoRED$population==targettemp),]
   rohpops[i,19]<-blocktemp$color[1]
   rohpops[i,20]<-blocktemp$set[1]
   
   for (k in 1: length(sketchrownames)){
     rohpops[i,k]<-mean(blocktemp[,k])  
     rohpops[i,k+6]<-min(blocktemp[,k]) 
     rohpops[i,k+12]<-max(blocktemp[,k])
   }
 }
 
 
# create a ROH melt file per population
 indROHinfoMELTpop<-melt(rohpops[,1:6],varnames=c("group","bin"))
 indROHinfoMELTpop$minn<-melt(rohpops[,7:12])[,3]
 indROHinfoMELTpop$maxx<-melt(rohpops[,13:18])[,3]
 write.table(indROHinfoMELTpop,"meltROHbinsXpops.txt",sep="\t")
 indROHinfoMELTpop<-read.table("meltROHbinsXpops.txt",header=T)

 # and now a loop to plot each pop (Figure S5)
 indROHinfoRED<-indROHinfo[which(indROHinfo$set=="presentstudy"),] #only the populations from the new dataset
popp<-table(indROHinfoRED$population)

 riordinacolori<-as.character(levels(indROHinfoMELTpop$group))
 listpopInfo3<-rohpops
 listpopInfo3<-listpopInfo3[riordinacolori,]
 colorelli<-as.character(listpopInfo3[,19])


infoREDD<-popinfo[which(popinfo$population%in%labels(popp)[[1]]),] # info per population
 for (i in 1:length(popp)){
   targettemp<-infoREDD$population[i]
   plottatemp<-indROHinfoMELTpop[which(indROHinfoMELTpop$group==targettemp),]
   ppp <- ggplot(indROHinfoMELTpresent, aes(x=bin,y=value,group=group))+
   geom_ribbon(aes(ymin = minn , ymax = maxx,  fill = group), alpha=0.2)+
   scale_x_discrete(limits=sketchrownames)+
   scale_fill_manual(values = colorelli,
                     breaks=riordinacolori,
                     name="")+
   theme_bw()   +
   theme(legend.position="none")+
 #    theme(legend.position="bottom")+
   labs(x = "ROH length category (Mb)",y="Total length of ROH (Mb) per individual",title=targettemp)
   
   
   print( ppp+  
   geom_ribbon(data=plottatemp,aes(ymin = minn , ymax = maxx),fill=as.character(infoREDD$color[i]), alpha=0.6)+
   geom_line(data=plottatemp,color=as.character(infoREDD$color[i]) )   )
 ggsave(paste(targettemp,"_ROHbinsXgroupTESTonlypresentstudy.pdf",collapse = "", sep=""), useDingbats=FALSE)
 }

 
### two ROH bins category: small and large ROH (less or more than 1.6 Mb), one against the other (Figure 3A)
 indROH2bins<-matrix(NA,nrow(myinfoSort),2)
 colnames(indROH2bins)<-c("small","large")
 rownames(indROH2bins)<-myinfoSort$sample_ID
 for (i in 1:nrow(indROH2bins)){
   targetID<-row.names(indROH2bins)[i]
   tempID<-rohinfo[which(rohinfo$sample_ID==targetID),]
   indROH2bins[i,1]<-sum(tempID[which(tempID$MB<1.6),]$MB)
   indROH2bins[i,2]<-sum(tempID[which(tempID$MB>1.6),]$MB)
    }
 indROHinfo2bins<-cbind(indROH2bins,myinfoSort)
 
 popinfo<-info[which(info$population%in%indROHinfo2bins$population),]
 
 pdf(file="ROH_largeAndSmall_individuals.pdf", useDingbats=FALSE)
 plot(x=indROHinfo2bins$small,y=indROHinfo2bins$large,col=as.character(indROHinfo2bins$color),pch=as.numeric(indROHinfo2bins$pch),cex=0.7,xlab="total length of small ROH (<1.6 Mb)", ylab="total length of large ROH (>1.6 Mb)")
 legend("topleft", legend=popinfo$population,  pch=as.numeric(popinfo$pch), col=as.character(popinfo$color), text.col="gray", ncol=3,bty="n", cex=0.5)
 dev.off()
 
