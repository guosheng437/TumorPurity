##############################################################################################################################
#Description: Generate Figure 4G-H
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date:   Feb-2022
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
###############################################################################################################################
library(ggpubr)
  
MINSAMPLES=10

a<-read.table('PDX_hybridESTIMATEscores.txt', sep='\t',header=TRUE, row.names=1)

#--keep only models in service: Cancer type is correct.
a=a[a$Is_in_service=='yes' | a$Cancer=='AL',]

#--remove MC and otherLU--#
a=a[a$Cancer!='MC',]
a=a[a$Cancer!='otherLU',]

#--merge LUAD and LUSC into NSCLC--#
a$Cancer[a$Cancer=='LUAD']='NSCLC'
a$Cancer[a$Cancer=='LUSC']='NSCLC'

#--remove Cancer with few samples--#
b=as.data.frame(table(a$Cancer)<MINSAMPLES)
b=cbind(b,b)
colnames(b) = c('t1','t2')
toRemove=rownames(b[b$t1==TRUE,])
a=a[!(a$Cancer %in% toRemove),]

medianStromalScore = aggregate(a$StromalScore, list(a$Cancer), median)
medianStromalScore = medianStromalScore[order(medianStromalScore$x),]
colnames(medianStromalScore) = c('Cancer', 'StromalScore')
pdf(file='Fig4G.pdf', width=14, height=4)
g<-ggboxplot(a, x = "Cancer", y = "StromalScore", color = "Cancer", 
		  xlab='', ylab="Stromal score",#notch =TRUE,
          add = "jitter", legend = "none" , size=0.25, order = c(medianStromalScore$Cancer),
		  add.params = list(size = 0.2, alpha = 0.25)
		  ) +
		rotate_x_text(angle = 45) 
ggpar(g, palette = "p1_aaas")  
dev.off()

medianImmuneScore = aggregate(a$ImmuneScore, list(a$Cancer), median)
medianImmuneScore = medianImmuneScore[order(medianImmuneScore$x),]
colnames(medianImmuneScore) = c('Cancer', 'ImmuneScore')
pdf(file='Fig4H.pdf', width=14, height=4)
g<-ggboxplot(a, x = "Cancer", y = "ImmuneScore", color = "Cancer", 
		  xlab='', ylab="Immune score",#notch =TRUE,
          add = "jitter", legend = "none" , size=0.25, order = c(medianImmuneScore$Cancer),
		  add.params = list(size = 0.2, alpha = 0.25)
		  ) +
		rotate_x_text(angle = 45) 
ggpar(g, palette = "p1_aaas")		  
dev.off()



