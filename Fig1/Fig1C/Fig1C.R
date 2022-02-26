##############################################################################################################################
#Description: generate Fig1C
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date:   Aug-2021
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
###############################################################################################################################
#--check if tumor purity differs between early and late passages in each Cancer
rm(list=ls())
library(ggpubr)

MINSAMPLES = 20 		#mimimal number of samples for a cancer to be included
MINPDX = 10				#mimimal number of PDX for a cancer to be included
MAX_PASSAGE= 10 		#max passage for a samples to be inclued 
MIN_SAMPLE2PLOT=30		#for individual cancer, minimal number of samples to be plotted
MIN_PASSAGE2PLOT=4		#for individual cancer, minimal number of passage difference to be plotted
MAX_SAMPLE_PER_PDX = 2

nod<-read.csv('nod.csv')

#--remove Cancer with few samples--#
a=as.data.frame(table(nod$Cancer)<MINSAMPLES)
a=cbind(a,a)
colnames(a) = c('t1','t2')
toRemove=rownames(a[a$t1==TRUE,])
nod=nod[!(nod$Cancer %in% toRemove),]

#--remove Cancer with few unique PDX models--#
Cancers = unique( nod$Cancer)
CancerKeep = c()
for(i in 1:length(Cancers)){
	nod2 = nod[nod$Cancer== Cancers[i],]
	if(length(unique(nod2$DuplicateModel))>MINPDX){
		CancerKeep = c(CancerKeep, Cancers[i])
	}
}
length(CancerKeep)

nod=nod[(nod$Cancer %in% CancerKeep),]

#--for each PDX, randomly keep at 2samples--#
Cancers = unique( nod$Cancer)
nod_new <-c()
for(i in 1:length(Cancers)){
	nod2 = nod[nod$Cancer== Cancers[i],]
	u<-unique(nod2$DuplicateModel)

	for(j in 1:length(u)){
		nod3 = nod2[nod2$DuplicateModel== u[j],]
		if(dim(nod3)[1]>=MAX_SAMPLE_PER_PDX){
			nod4 = nod3[sample(nrow(nod3), MAX_SAMPLE_PER_PDX), ]		
			nod_new<-rbind(nod_new,nod4)
		}#else{
		#	nod_new<-rbind(nod_new,nod3)
		#}		
	}
}

nod<-nod_new

#--for each Cancer, randomly pick two samples from a PDX and calculate correlation--#
TP1_vector_all<-c()
TP2_vector_all<-c()
#pdf(file='PurityCorrelation_by_Cancer.pdf')
Cancers = unique( nod$Cancer)
for(i in 1:length(Cancers)){
	nod2 = nod[nod$Cancer== Cancers[i],]
	u<-unique(nod2$DuplicateModel)

	TP1_vector<-c()
	TP2_vector<-c()
	
	for(j in 1:length(u)){
		nod3=nod2[nod2$TestModel2 == u[j],]
		
		tm<- sapply(strsplit(nod3$TestModel, "-"), `[`, 2)
		tm<-gsub("R[0-9]+","",tm)		
		tm<-gsub("^P","",tm)
		tm<-as.numeric(tm)
		if(length(tm)==1){next;}
		if(is.na(tm[1])==TRUE){next;}
		if(is.na(tm[2])==TRUE){next;}
		if(tm[1]<tm[2]){
			TP1_vector<-c(TP1_vector,nod3$TumorPurity[1])
			TP2_vector<-c(TP2_vector,nod3$TumorPurity[2])
		}else{		
			TP1_vector<-c(TP1_vector,nod3$TumorPurity[2])
			TP2_vector<-c(TP2_vector,nod3$TumorPurity[1])
		}
	}
	
	TP<-cbind(TP1_vector,TP2_vector)
	colnames(TP) = c('EP','LP')
	TP<-as.data.frame(TP)
	if(dim(TP)[1]<30){next;}
	
	gg<-ggscatter(TP, x="EP", y="LP", 
			  xlab="Early passage tumor purity (%)",
			  ylab="Late passage tumor purity (%)",
			  main=Cancers[i],
			  add="reg.line",
			  conf.int=TRUE,
			   ) +
			  stat_cor(size=6,label.x.npc = "left", label.y.npc = "top")+
			  theme_pubr(base_size = 20)
	
	print(gg)
	
	TP1_vector_all <- c(TP1_vector_all,TP1_vector)
	TP2_vector_all <- c(TP2_vector_all,TP2_vector)
	
}
#dev.off()

TP<-cbind(TP1_vector_all,TP2_vector_all)
colnames(TP) = c('EP','LP')
TP<-as.data.frame(TP)
if(dim(TP)[1]<30){next;}
	


pdf(file='Fig1C.pdf')
df<-TP
x <- densCols(df$EP,df$LP, colramp=colorRampPalette(c("black", "white")))
df$dens <- col2rgb(x)[1,] + 1L
cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                            "#FCFF00", "#FF9400", "#FF3100"))(256)
df$col <- cols[df$dens]

plot(LP~EP, data=df[order(df$dens),], pch=20, col=col, cex=1,
		cex.lab=1.667, 
		cex.axis=1.667,		
		bty = 'l',
		xlab="Early passage tumor purity (%)", ylab="Late passage tumor purity (%)")
text(x=20, y=100, paste("R=0.60"),cex=1.667)
	
	
abline(lm(LP ~ EP+0,df),lwd=3)
dev.off()

