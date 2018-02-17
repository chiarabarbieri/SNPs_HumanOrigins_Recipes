

### prepares the outfile to manipulate in R

maldersuccess<-read.table("outputSuccessMalder.cat.txt")
colnames(maldersuccess)<-c("test", "status",     "p-value", "test_pop","ref_A","ref_B", "Zscore",  "divergenceTime", "null", "interval")

library(ggplot2)
ThisYear<-2018  # to have the admixture times in calendar years
generationyeras<-30  # chose your generation time

maldersuccess$AgeCalendarYear<-ThisYear-(maldersuccess$divergenceTime*generationyeras)
maldersuccess$Min<-ThisYear-((maldersuccess$divergenceTime+maldersuccess$interval)*generationyeras)
maldersuccess$Max<-ThisYear-((maldersuccess$divergenceTime-maldersuccess$interval)*generationyeras)
maldersuccess$couple<-paste(maldersuccess$ref_A,maldersuccess$ref_B)

admixtureTarget<-c("Karitiana", "Surui", "Puertorico", "Lima", "Pima", "Aymara", "Mixe")
poporder<-admixtureTarget


pdf("divergenceAge_EuropeAfrica_Xavante_Malder.pdf")
gg<-ggplot(maldersuccess,aes(x=AgeCalendarYear,y=test_pop, group=couple,color= couple))+
            geom_errorbarh(aes(xmin=Min,xmax=Max), alpha=0.2) +
            geom_point(na.rm=TRUE)+
            ggtitle("Admixture Times from Africa and Europe") +
            theme_bw() +
            scale_fill_brewer(palette = "Accent")+
              scale_y_discrete(limits=poporder)     # this will plot all the Target admixed in the right order on the Y axis       
gg
dev.off()

# i want to exclude the admixture not supported by p value <0.01
maldersignificant<-maldersuccess[which(maldersuccess[,3]<0.001),]

pdf("divergenceAge_EuropeAfrica_Xavante_MalderSignificant0001.pdf")
gg<-ggplot(maldersignificant,aes(x=AgeCalendarYear,y=test_pop,color= couple))+
            geom_errorbarh(aes(xmin=Min,xmax=Max), alpha=0.2) +
            geom_point(na.rm=TRUE)+
            ggtitle("Admixture Times from Africa and Europe") +
            theme_bw() +
            scale_fill_brewer(palette = "Accent")+
              scale_y_discrete(limits=poporder)            
gg
dev.off()

# now color for the Z-score 
pdf("divergenceAge_EuropeAfrica_Xavante_MalderSignificant0001_zscore.pdf")
gg<-ggplot(maldersignificant,aes(x=AgeCalendarYear,y=test_pop,color= Zscore))+
            geom_errorbarh(aes(xmin=Min,xmax=Max), alpha=0.2) +
            geom_point(na.rm=TRUE)+
            scale_colour_gradientn(colours = c("darkred","red","orange"))+
             ggtitle("Admixture Times from Africa and Europe") +
            theme_bw() +
              scale_y_discrete(limits=poporder)            
gg
dev.off()




# _________ plot the admixture LD decay

LD<-read.table("LD_decaycurveALLfilter.txt",header=T)
LD<-LD[-which(LD$d=="Inf"),] # exclude values of d that are marked as Infinite

subset<-LD[which(LD$Target=="Karitiana"),]

library(ggplot2)
poporder1<-as.character(poporder[which(poporder%in%maldersignificant$test_pop)])

# i will plot two ways: with all the three possible combinations of parental populations, and with the two combinations that are more relevant.
totalLD_3admix<-matrix(NA,1,5)
colnames(totalLD_3admix)<-c(  "d", "Italian_North.Xavante", "Italian_North.Yoruba", "Yoruba.Xavante", "target")
totalLD_2admix<-matrix(NA,1,4)
colnames(totalLD_2admix)<-c(  "d", "Italian_North.Xavante", "Yoruba.Xavante", "target")


for (i in 1:length(poporder1)){
poptemp<-poporder1[i]
tempLD<-read.table(paste(poptemp,"_Xavante_LDout",sep="",collapse=""),header=T,as.is=T)
tempLD<-tempLD[-which(tempLD$d=="Inf"),]
tempLD$target<-poptemp
if(ncol(tempLD)==5){
totalLD_3admix<-rbind(totalLD_3admix,tempLD)
} else {
totalLD_2admix<-rbind(totalLD_2admix,tempLD)
}
}

totalLD_2admix<-totalLD_2admix[-1,]
totalLD_3admix<-totalLD_3admix[-1,]

gg<-ggplot(totalLD_2admix,aes(x=d, y=value, color=variable))+
geom_point(aes(x=d,y = Italian_North.Xavante, col="Italian_North.Xavante"), shape=3, size=0.02)+
geom_point(aes(x=d,y = Yoruba.Xavante, col="Yoruba.Xavante"), shape=3, size=0.02)+
 ggtitle("Admixture LD decay from Africa and Europe") +
  theme_bw() +
facet_wrap(~target, ncol=2)

gg

dev2bitmap("LDdecay_Xavante_2sources.pdf",type="pdfwrite")


gg<-ggplot(totalLD_3admix,aes(x=d, y=value, color=variable))+
geom_point(aes(x=d,y = Italian_North.Xavante, col="Italian_North.Xavante"), shape=3, size=0.02)+
geom_point(aes(x=d,y = Yoruba.Xavante, col="Yoruba.Xavante"), shape=3, size=0.02)+
geom_point(aes(x=d,y = Italian_North.Yoruba, col="Italian_North.Yoruba"), shape=3, size=0.02)+
 ggtitle("Admixture LD decay from Africa and Europe") +
  theme_bw() +
facet_wrap(~target, ncol=3)
gg

dev2bitmap("LDdecay_Xavante_3sources.pdf",type="pdfwrite")


