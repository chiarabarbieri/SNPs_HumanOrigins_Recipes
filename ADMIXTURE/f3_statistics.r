
## prepare a the population list to run f3
## with a combination of Target, Native American unadmixed, Spanish 
##  and with a combination of Target, Native American unadmixed, Yoruba 


pure<-as.character(read.table("poplistPureAme.txt")[,1]) # list of Native American unadmixed
all<-as.character(read.table("poplistAllAme.txt")[,1]) # all target populations

combinations<-matrix(NA,length(pure),length(all))
colnames(combinations)<-all
rownames(combinations)<-pure

library(reshape)
coso<-melt(combinations)
coso<-coso[-which(as.character(coso[,1])==as.character(coso[,2])),]

yor<-rep("Yoruba", nrow(coso))
spa<-rep("Spanish", nrow(coso))
source1<-c(yor,spa)
coso2<-rbind(coso,coso)

finalcomb<-cbind( source1, coso2)
finalcomb<-finalcomb[,-4]

write.table(finalcomb,"f3_poplist.txt", row.names=F, col.names=F, quote=FALSE)


# ------------------------------------------
# now run f3 with the command 
# qp3Pop  -p par.F3 > results_f3
# ------------------------------------------

# the par.F3 file will look as follows:

genotypename: 	yourfile.geno
indivname: 		yourfile.ind
snpname: 		yourfile.snp
popfilename:  f3_poplist.txt
# ------------------------------------------



### now elaborate the result file in R to be plottable 

coso<-read.table("results_f3.txt",header=T)

popoloni<-table(coso$Target)
poponames<-names(popoloni)

screma<-matrix(NA,length(poponames),6)
rownames(screma)<-poponames
for (i in 1:nrow(screma)){
  poptemp<-poponames[i]
  blocco<-coso[which(coso$Target==poptemp),]
  bloccoA<-blocco[which(blocco$Source1=="Spanish"),]
  bloccoB<-blocco[which(blocco$Source1=="Yoruba"),]
  screma[i,c(1:2)]<- as.numeric(bloccoA[order(bloccoA$f_3),][1,5:6])
  screma[i,3]<- as.character(bloccoA[order(bloccoA$f_3),][1,3])
  screma[i,c(4:5)]<- as.numeric(bloccoB[order(bloccoB$f_3),][1,5:6])
  screma[i,6]<- as.character(bloccoB[order(bloccoB$f_3),][1,3])
}

colnames(screma)<-c("f3Spanish","f3Spanishinterval","f3SpanishSource","f3African","f3Africaninterval","f3AfricanSource")
screma<-cbind(screma,poponames)
write.table(screma,"F3resultsclean.txt",sep="\t")
screma<-read.table("F3resultsclean.txt",sep="\t",header = T, as.is = T)

# time to plot, European admixture f3

scremaE<-screma[which(screma$f3Spanish<0),]

coso2<-table(scremaE$f3SpanishSource)
colori<-rainbow(length(coso2))
trick<-cbind(unlist(labels(coso2)),colori)
colnames(trick)<-c("f3SpanishSource","colors1")
trick[,1]<-unlist(labels(coso2))
scremaEurope<-merge(scremaE,trick)

  
library(Hmisc)
scremaOrdered = scremaEurope[order(scremaEurope$f3Spanish),]

pdf("F3_europe_Targetpops.pdf")
plot(scremaOrdered$f3Spanish)
errbar(1:nrow(scremaOrdered), scremaOrdered$f3Spanish,
       (scremaOrdered$f3Spanish+scremaOrdered$f3Spanishinterval),
       (scremaOrdered$f3Spanish-scremaOrdered$f3Spanishinterval), pch=20, las=2, cex.axis=0.4, xaxt='n', col=as.character(scremaOrdered$colors1),
       xlab="", ylab="F3")
axis(1,at=1:nrow(scremaOrdered), labels=scremaOrdered$poponames, las=2, cex.axis=0.6)
legend("topleft", legend=trick[,1],  pch=20, col=as.character(trick[,2]), text.col="gray", ncol=5,bty="n", cex=0.6)
title("F3 with European Source")
dev.off()

# time to plot, African admixture f3

scremaA<-screma[which(screma$f3African<0),]

coso2<-table(scremaA$f3AfricanSource)
colori<-rainbow(length(coso2))
trick<-cbind(unlist(labels(coso2)),colori)
colnames(trick)<-c("f3AfricanSource","colors1")
trick[,1]<-unlist(labels(coso2))
scremaAfrica<-merge(scremaA,trick)
scremaAfrica<-scremaAfrica[which(scremaAfrica$poponames%in%MYINFORED$populationAdmixture),]



library(Hmisc)
scremaOrdered = scremaAfrica[order(scremaAfrica$f3African),]

pdf("F3_africa_targetpops.pdf")
plot(scremaOrdered$f3African)
errbar(1:nrow(scremaOrdered), scremaOrdered$f3African,
       (scremaOrdered$f3African+scremaOrdered$f3Africaninterval),
       (scremaOrdered$f3African-scremaOrdered$f3Africaninterval), pch=20, las=2, cex.axis=0.4, xaxt='n', col=as.character(scremaOrdered$colors1),
       xlab="", ylab="F3")
axis(1,at=1:nrow(scremaOrdered), labels=scremaOrdered$poponames, las=2, cex.axis=0.6)
legend("topleft", legend=trick[,1],  pch=20, col=as.character(trick[,2]), text.col="gray", ncol=5,bty="n", cex=0.6)
title("F3 with African Source")
dev.off()
