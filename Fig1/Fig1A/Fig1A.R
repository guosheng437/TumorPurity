##############################################################################################################################
#Description: generate Fig1A
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date:   Aug-2021
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
###############################################################################################################################
library(ggpubr)

a<-read.csv('PDXsomatic_mutation_freq_minread30_results.csv')
a=a[a$Passage<11,]
median(a$MedianMutationFrequency)

a$Passage = paste('P', a$Passage,sep='')
pdf(file='Fig1A.pdf', width=8, height=4)
g<-ggboxplot(a, x = "Passage", y = "MedianMutationFrequency", color = "Passage", 
		  xlab='Passage', ylab="Median variant allelic fraction (VAF) (%)",notch =FALSE,
          add = "jitter", legend = "none" , size=0.25,
		  add.params = list(size = 0.2, alpha = 0.25)
		  ) +		  
		rotate_x_text(angle = 45) 
g<-g+stat_compare_means(method = "anova")#, label.y = 40)
ggpar(g, palette = "p1_aaas") 
		  
dev.off()



