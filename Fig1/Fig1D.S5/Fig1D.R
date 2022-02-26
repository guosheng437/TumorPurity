#################################################################################################
#Description: generate Fig1D and Fig S1 for PDX models
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date:   Aug-2021
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
#################################################################################################
#--get density plot of within and between PDX purity for each cancer and as a whole
library(ggpubr)

MINSAMPLES = 30
MINPDX = 10
MINPASSAGE=1
MAXPASSAGE=10
MAX_SAMPLE_PER_PDX=50
MIN_SAMPLE4PLOT=30				#for individual cancer, minimal number of samples to be plotted

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


#--keep only PDX <= MAXPASSAGE passage--#
tm<- sapply(strsplit(nod$TestModel, "-"), `[`, 2)
tm<-gsub("R[0-9]+","",tm)		
rn<-startsWith(tm, 'P')
tm<-gsub("^P","",tm)
tm<-as.numeric(tm)
rn<-rn & !is.na(tm) & tm<=MAXPASSAGE
table(rn)		
nod=nod[rn==TRUE,]

#--keep only PDX >= MINPASSAGE passage--#
tm<- sapply(strsplit(nod$TestModel, "-"), `[`, 2)
tm<-gsub("R[0-9]+","",tm)		
rn<-startsWith(tm, 'P')
tm<-gsub("^P","",tm)
tm<-as.numeric(tm)
rn<-rn & !is.na(tm) & tm>=MINPASSAGE
table(rn)		
nod=nod[rn==TRUE,]

#--keep only the first sample for a mouse with multiple samples: each mouse has unique 5-digit mouse number--#
Cancers = unique( nod$Cancer)
nod_new <-c()
for(i in 1:length(Cancers)){
	nod2 = nod[nod$Cancer== Cancers[i],]
	u<-unique(nod2$DuplicateModel)

	nod3=nod2
	
	#--get mouse number--#
	TestModel<-gsub('-(\\d{8})', '', nod3$TestModel) #--remove 8-digit yyyymmdd
	#nod3$TestModel
	p3<-sapply(strsplit(TestModel, "-"), `[`, 3)
	
	t1<-grepl("^[0-9]{5}$", p3, perl = T) #TRUE is 5-digit numbers		
	t2<-duplicated(p3) #duplicated mouse number will be marked TRUE
		
	nod3=nod3[t1 &(!t2)==TRUE,]	
	
	if(dim(nod3)[1]>1){
		nod_new<-rbind(nod_new,nod3)
	}	
}
nod<-nod_new

#--for each PDX, randomly keep at most 5 samples--#
Cancers = unique( nod$Cancer)
nod_new <-c()
for(i in 1:length(Cancers)){
	nod2 = nod[nod$Cancer== Cancers[i],]
	u<-unique(nod2$DuplicateModel)

	for(j in 1:length(u)){
		nod3 = nod2[nod2$DuplicateModel== u[j],]
		if(dim(nod3)[1]>MAX_SAMPLE_PER_PDX){
			nod4 = nod3[sample(nrow(nod3), MAX_SAMPLE_PER_PDX), ]		
			nod_new<-rbind(nod_new,nod4)
		}else{
			nod_new<-rbind(nod_new,nod3)
		}		
	}
}


#--for each cancer type, compute within and between PDX purity difference--#
wdist_all<-c()
bdist_all<-c()
Cancers = sort(unique( nod$Cancer))
pdf(file='FigS5.pdf')
for(i in 1:length(Cancers)){
	nod2 = nod[nod$Cancer== Cancers[i],]
	u<-unique(nod2$DuplicateModel)
	 nod2[order(nod2$TestModel),]


	#--compute within-PDX purity difference
	wdist<-c()
	for(j in 1:length(u)){
		nod3 = nod2[nod2$DuplicateModel== u[j],]
		np = dim(nod3)[1]
		if(np==1){next}
		
		wdist<-c(wdist, as.numeric(dist(nod3[,7])))
	}
	
	if(length(wdist)<MIN_SAMPLE4PLOT){next;}

	
	#--compute between-PDX purity difference
	bdist<-c()
	for(j in 1:(length(u)-1)){
		nod3 = nod2[nod2$DuplicateModel== u[j],]
	
		for(k in (j+1):length(u)){
			nod4 = nod2[nod2$DuplicateModel== u[k],]
			bdist<-c(bdist,abs(outer(nod3[,7],nod4[,7],"-")))
		}		
	}
	
	wdata<-data.frame( 	TumorPurity=c(wdist,bdist), TPtype=factor(c(rep('within-PDX',length(wdist)),rep('between-PDX',length(bdist)) )) ) 
	
	pval=wilcox.test(bdist, wdist)$p.value
	pval=formatC(pval, format = "e", digits = 2)

	
	gg<-ggdensity(wdata, x = "TumorPurity",
	xlab='Tumor purity difference (%)',
	ylab='Density',
	main=Cancers[i],
	add = "median", rug = TRUE,
	color = "TPtype", fill = "TPtype",
	legend = "top",                                 
    legend.title = paste('Mann–Whitney p-value=',pval,sep=''),
	palette = c("#00AFBB", "#E7B800"))
	
	gg<- gg+ geom_text(data = aggregate(TumorPurity~TPtype, data = wdata, FUN = median),
            aes(x = TumorPurity, y = Inf, color = TPtype, label = round(TumorPurity,2)), 
            vjust = "inward",  hjust = "inward")
 
	print(gg)
	
	wdist_all<-c(wdist_all, wdist)
	bdist_all<-c(bdist_all, bdist)
}
dev.off()

wdata<-data.frame( 	TumorPurity=c(wdist_all,bdist_all), TPtype=factor(c(rep('within-PDX',length(wdist_all)),rep('between-PDX',length(bdist_all)) )) ) 
pdf(file='Fig1D.pdf')
	
	gg<-ggdensity(wdata, x = "TumorPurity",
	xlab='Tumor purity difference (%)',
	ylab='Density',
	add = "median", 
	rug = FALSE,
	color = "TPtype", fill = "TPtype",
	#legend = "top",                                 
    #legend.title = paste('Mann–Whitney p-value=',pval,sep=''),
	palette = c("#00AFBB", "#E7B800"))
	
	gg<- gg+ geom_text(data = aggregate(TumorPurity~TPtype, data = wdata, FUN = median),
            aes(x = TumorPurity, y = Inf, color = TPtype, label = round(TumorPurity,2)),			
            vjust = "inward",  hjust = "inward")
			#vjust =1)
	gg<- gg+ theme_pubr(base_size = 20)

 
	print(gg)
dev.off()
