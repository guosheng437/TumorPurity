#################################################################################################
#Description: generate Fig2A-B
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date: Aug-2021
#Note: only use NOD and BALB/c mice with passage information
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
#################################################################################################
rm(list=ls())
library(ggpubr)

MINSAMPLES = 25
MAX_SAMPLE_PER_PDX = 5

a<-read.csv('nod_balbc.csv')

#--analysis done in NOD/SCID
nod=a[grep("NOD",a$Vendor),]
balb=a[grep("Balb",a$Vendor),]


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

nod<-nod_new

Cancers = unique( balb$Cancer)
balb_new <-c()
for(i in 1:length(Cancers)){
	balb2 = balb[balb$Cancer== Cancers[i],]
	u<-unique(balb2$DuplicateModel)

	for(j in 1:length(u)){
		balb3 = balb2[balb2$DuplicateModel== u[j],]
		if(dim(balb3)[1]>MAX_SAMPLE_PER_PDX){
			balb4 = balb3[sample(nrow(balb3), MAX_SAMPLE_PER_PDX), ]		
			balb_new<-rbind(balb_new,balb4)
		}else{
			balb_new<-rbind(balb_new,balb3)
		}		
	}
}

balb<-balb_new



#--select cancers with at least MINSAMPLES samples in both strains
table(balb$Cancer)>MINSAMPLES 
table(nod$Cancer)>MINSAMPLES
b=table(balb$Cancer)>MINSAMPLES
b=b[b==TRUE]
c=table(nod$Cancer)>MINSAMPLES
c=c[c==TRUE]
common=intersect(names(b), names(c))

nod_b<-c()
balb_b<-c()
for(i in 1:length(common)){
	
	nod2 = nod[nod$Cancer==common[i],]
	nod3 = balb[balb$Cancer==common[i],]
	
	n = min(dim(nod2)[1], dim(nod3)[1])
	#n = MINSAMPLES

	nod2 = nod2[sample(nrow(nod2), n), ]		
	nod3 = nod3[sample(nrow(nod3), n), ]		
	
	#nod_b<-c(nod_b,nod2$TumorPurity,nod2$Cancer)
	#balb_b<-c(balb_b,nod3$TumorPurity,nod2$Cancer)
	nod_b = rbind(nod_b, data.frame(nod2$TumorPurity,nod2$Cancer))
	balb_b = rbind(balb_b, data.frame(nod3$TumorPurity,nod3$Cancer))

}
colnames(nod_b) = c('TumorPurity','Cancer')
colnames(balb_b) = c('TumorPurity','Cancer')

nod_b = data.frame(nod_b$Cancer, rep('NOD/SCID',length(nod_b$Cancer)), nod_b$TumorPurity)
balb_b = data.frame(balb_b$Cancer, rep('BALB/c nude',length(balb_b$Cancer)), balb_b$TumorPurity)
colnames(nod_b)=c('Cancer', 'MouseStrain','TumorPurity')
colnames(balb_b)=c('Cancer', 'MouseStrain','TumorPurity')

nb = rbind(nod_b,balb_b)

#--a single boxplot between two mouse strains--#
pdf(file='Fig2A.pdf')
p <- ggboxplot(nb, x = "MouseStrain", y = "TumorPurity",
				color = "MouseStrain", palette =c("#E7B800","#00AFBB"),
                add = "jitter",
				ylab='Tumor purity (%)', 
				xlab='',
				notch=FALSE,
				order=c('BALB/c nude','NOD/SCID'),
				add.params = list(size = 0.5, jitter = 0.2)
				)
p<- p + stat_compare_means(method = "wilcox.test",  label.x.npc = "center", label.y.npc = "bottom", size=6, label = "p")
p<- p + theme_pubr(base_size = 20)				
ggpar(p, legend = "none")
dev.off()


#--a grouped boxplot for 9 cancers--#
pdf(file='Fig2B.pdf', width=8, height=4)
p<-ggboxplot(nb, "Cancer", "TumorPurity", color = "MouseStrain",
		palette = c("#E7B800","#00AFBB"),
		add = "jitter",
		ylab='Tumor purity (%)', 
		xlab='Cancer',
		size=0.25,
		notch=FALSE,
		#order=c('BALB/c nude','NOD/SCID'),
		add.params = list(size = 0.2, jitter = 0.2)
	)
p<- p + stat_compare_means(aes(group = MouseStrain), method = "wilcox.test",  label.x.npc = "center", label.y.npc = "bottom", size=3, label = "p")
p
dev.off()

