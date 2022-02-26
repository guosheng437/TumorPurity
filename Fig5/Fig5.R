#################################################################################################
#Description: generate Fig5A-D
#Author: Sheng Guo PhD  guosheng@crownbio.com
#Date:   Aug-2021
#Apache License
#Version 2.0, January 2004
#https://opensource.org/licenses/Apache-2.0
########################################################################################################################
library(ggpubr)

a<-read.table('SyngeneicModelTumorPurity.txt', head=TRUE)
a$TumorPurity=a$TumorPurity*100

medianPurity = aggregate(a[, 3], list(a$Model), median)
medianPurity = medianPurity[order(medianPurity$x),]
colnames(medianPurity) = c('Model', 'TumorPurity')
#write.csv(file='medianTumorPurityByModel.csv', medianPurity)

medianValidSNP = aggregate(a[, 4], list(a$Model), median)
#write.csv(file='medianValidSNPByModel.csv', medianValidSNP)

pdf(file='Fig5A.pdf', width=6, height=4)
g<-ggstripchart(a,"Model", "TumorPurity", 
			color = "Model",
			shape=17,
			xlab='',
			ylab="Tumor purity (%)",
			add='none',
			jitter=0.0,
			legend = "none",
			order = c(medianPurity$Model),
			add.params = list(size = 0.5, alpha = 0.25)) +
	rotate_x_text(angle = 45) +
	ylim(25,100)
ggpar(g, palette = "p1_aaas") # nature
dev.off()



medianStromalScore = aggregate(a[, 6], list(a$Model), median)
medianStromalScore = medianStromalScore[order(medianStromalScore$x),]
colnames(medianStromalScore) = c('Model', 'TumorStromalScore')

pdf(file='Fig5B.pdf', width=6, height=4)
g<-ggstripchart(a,"Model", "StromalScore", 
			color = "Model",
			shape=17,
			xlab='',
			ylab="Stromal score",
			add='none',
			jitter=0.0,
			legend = "none",
			order = c(medianStromalScore$Model),
			add.params = list(size = 0.25, alpha = 0.25)) +
	ylim(-2000,1000)+
	rotate_x_text(angle = 45) 
	
ggpar(g, palette = "p1_aaas") # nature
dev.off()


medianImmuneScore = aggregate(a[, 7], list(a$Model), median)
medianImmuneScore = medianImmuneScore[order(medianImmuneScore$x),]
colnames(medianImmuneScore) = c('Model', 'TumorImmuneScore')

pdf(file='Fig5C.pdf', width=6, height=4)
#pdf(file='TumorImmuneScoreByModel.pdf')
g<-ggstripchart(a,"Model", "ImmuneScore", 
			color = "Model",
			shape=17,
			xlab='',
			ylab="Immune score",
			add='none',
			jitter=0.0,
			legend = "none",
			order = c(medianImmuneScore$Model),
			add.params = list(size = 0.25, alpha = 0.25)) +
	ylim(-1200,2500)+
	rotate_x_text(angle = 45) 
	
ggpar(g, palette = "p1_aaas") # nature
dev.off()

pdf(file='Fig5D.pdf',width=6,height=4)
ggscatter(a, x="ESTIMATEScore", y="TumorPurity", 
		  ylab="Tumor purity (%)",
		  xlab="ESTIMATE score",		  
		  xlim=c(-3000,3000),
		  ylim=c(0,100),
		  size=1.5,
		  col='blue',
		  #add="reg.line",
		  #add.params = list(size = 1.0, alpha = 0.25),
		  conf.int=FALSE,
		  add.params = list(size = 0.5, alpha = 0.25)) +
		  #stat_cor(size=6)+
		  annotate("text", x = -2000, y = 5, label = "italic(R) == -0.74", parse = TRUE)
		  #theme_pubr(base_size = 20)
dev.off()



