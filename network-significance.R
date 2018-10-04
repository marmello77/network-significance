############################################################
#                                                          # 
#    SCRIPT FOR ESTIMATING P-VALUES FOR NETWORK METRICS    #
#                                                          # 
############################################################

##### Ecological Synthesis Lab (SintECO)
##### https://marcomellolab.wordpress.com
##### Authors: Renata Muylaert & Pavel Dodonov
##### E-mail: renatamuy@gmail.com
##### How to cite: Muylaert RL & Dodonov P. 2016. How to estimate
##### P-values of network metrics and compare pairs of networks
##### using Monte Carlo procedures. Ecological Synthesis Lab
##### at the University of São Paulo, Brazil.

##### Published on April 25th, 2017 (English version).
##### Run in R 3.3.3 (2017-03-06) -- "Another Canoe"

##### Disclaimer: You may use this script freely for non-comercial
##### purposes at your own risk. We assume no responsibility or
##### liability for the use of this software, convey no license
##### or title under any patent, copyright, or mask work right
##### to the product. We reserve the right to make changes in
##### the software without notification. 
##### We also make no representation or warranty that such
##### application will be suitable for the specified use without
##### further testing or modification. If this script helps you
##### produce any academic work (paper, book, chapter, 
##### dissertation etc.), please acknowledge the authors and
##### cite the source.


#############################################################


####SUMMARY##################################################


#1. Preparation
#2. P-value for one network (Dormann)
#3. P-value for one network (Muylaert & Dododonov)
#4. P-value for QuanBiMo modularity
#5. Comparing one metric between two networks
#6. Comparing QuanBiMo between two networks
#7. Suggested readings


#############################################################


####1. PREPARATION##########################################


#Load the packages.
library(bipartite)
library(reshape2)
library(sna)
library(igraph)
library(network)
library(tcltk)
library(vegan)
library(network)

#Set the working directory
setwd("path to the folder")


##############################################################


####2. P-value for one network (Dormann)####

#Delete all previous objects
rm(list= ls())

#Create the object to be analyzed
data<-read.table("rede1.txt", head=TRUE)

#Calculate the desired metric for the original network (observed)
obs <- unlist(networklevel(data, index="nestedness"))
obs

#Create randomized networks using a null model. 
#(It may take long, depending on your computer's processing power)
#The paramater "method" sets the null model and the parameter "N" sets the number of permutations
nulls <- nullmodel(data, N=10, method=5, autotransform="equiprobable")

#Calculate the same metric for all randomized networks
null <- unlist(sapply(nulls, networklevel, index="nestedness")) 
null

#Plot the observed value against the distribution of randomized values
plot(density(null), xlim=c(min(obs, min(null)), max(obs, max(null))), 
     main="Observed vs. randomized", xlab = "The chosen metric")
abline(v=obs, col="red", lwd=2)    

#Estimate the P-value
mean(null)
sd(null)
obs
praw <- sum(null>obs) / length(null)
ifelse(praw > 0.5, 1-praw, praw)    # P-value


#########################################################################


####3. P-value for one network (Muylaert & Dododonov)####


#Patefield (1981) null model

#Delete all previous objects
rm(list= ls())

#Create the object to be analyzed
data<-read.table("rede1.txt", head=TRUE)
data_neo=as.matrix(data)
data_neo

#Choose the metric
metrics=c("H2") 

#Save the observed value
orig=networklevel(data_neo, index=metrics)
orig
real<-data.frame(orig)
real

#Create the randomized networks. It is advised to set to at least 999.
#(It may take long, depending on your computer's processing power)
Nperm = 9
randomized.patef=matrix(nrow=length(orig),ncol=Nperm+1)
row.names(randomized.patef)=names(orig)
randomized.patef[,1]=orig 

#Repeat the procedure in a loop
i<-1
while(i<=Nperm){ 

  data_aleat=permatfull(data_neo,fixedmar="both",mtype="count",times=1)
  
  data_aleat=data_aleat$perm[[1]]
  
  linha<-networklevel(data_aleat, index=metrics)
  
  randomized.patef[,i+1]=linha
  print(i)
	i=i+1
}
 randomized.patef

#Plot the observed value against the distribution of randomized values
niveis<-row.names(randomized.patef)
for(k in niveis)
	{
		if(any(is.na(randomized.patef[k,]) == TRUE))
			{
			print(c(k, "metrica tem NA"))
			} else {
	nome.arq<- paste("Patefield_null_",k,".png", sep="")
	png(filename= nome.arq, res= 300, height= 15, width=21, unit="cm")
	plot(density(randomized.patef[k,]), main="Observed vs. randomized",)
	abline(v=orig[k], col="red", lwd=2, xlab="")
	dev.off()
	print(k)
	nome.arq<- paste("Patefield_Null_mean_sd_",k,".txt", sep="")
	write.table(cbind(mean(randomized.patef[k,]),sd(randomized.patef[k,])), file=paste(nome.arq,sep=""), 
      sep=" ",row.names=TRUE,col.names=FALSE)
		}
	}

#Estimate the P-value
significance.patef=matrix(nrow=nrow(randomized.patef),ncol=3)
row.names(significance.patef)=row.names(randomized.patef)
colnames(significance.patef)=c("p (rand <= orig)", "p (rand >= orig)", "p (rand=orig)")

signif.sup=function(x) sum(x>=x[1])/length(x)
signif.inf=function(x) sum(x<=x[1])/length(x)
signif.two=function(x) ifelse(min(x)*2>1,1,min(x)*2)

significance.patef[,1]=apply(randomized.patef,1,signif.inf)
significance.patef[,2]=apply(randomized.patef,1,signif.sup)
significance.patef[,3]=apply(significance.patef[,-3],1,signif.two)

significance.patef



##############################################################


####4. P-value for QuanBiMo modularity####

#Dormann & Strauss (2014)

#Delete all previous objects
rm(list= ls())

#Create the object to be analyzed
data<-read.table("rede1.txt", head=TRUE)
data

#Run the QuanBiMo algorithm
mod <- computeModules(web=data, steps=1E6)
mod

#Show the M-value (modularity)
mod@likelihood

#Export the matrix with modules
file_name= paste("Modulos_Quanbimo.png", sep="")
png(filename= file_name, res= 300, height= 25, width=30, unit="cm")
plotModuleWeb(mod)
dev.off()


#Calculate hierarchical modules
modn <- computeModules(data, steps=1E6, deep=T) 

#Estimate the P-value using the Patefield (1981) algorithm
nulls <- nullmodel(data, N=999, method="r2d")
modules.nulls <- sapply(nulls, computeModules)
like.nulls <- sapply(modules.nulls, function(x) x@likelihood)
(z <- (mod@likelihood - mean(like.nulls))/sd(like.nulls))
czvalues(mod) # for all nodes
czvalues(mod, level="lower") # for rows
czvalues(mod, level="lower", weighted=TRUE) #based on edge weights

#Plot the observed value against the distribution of randomized values
plot(density(like.nulls), xlim=c(min((mod@likelihood), min(like.nulls)), max((mod@likelihood), max(like.nulls))), 
     main="Observed vs. randomized")
abline(v=(mod@likelihood), col="red", lwd=2)    

#Estimate the P-value
mean(like.nulls)
sd(like.nulls)
mod@likelihood
praw <- sum(like.nulls>(mod@likelihood)) / length(like.nulls)
ifelse(praw > 0.5, 1-praw, praw)


##########################################################


####5. Comparing one metric between two networks####

#Patefield (1981) null model
#(It may take long, depending on your computer's processing power)


#Delete all previous objects
rm(list= ls())

#Load the two networks as objects
data<-read.table("rede1.txt", head=TRUE)
data
data2<-read.table("rede2.txt", head=TRUE)
data2
data_neo=as.matrix(data)
data_neo
data2_neo=as.matrix(data2)
data2_neo

#Choose the metric
metrics=c("H2")
metrics

#Calculate the difference in the chosen metrics between pairs of randomized networks
orig=abs(networklevel(data_neo,index=metrics)-networklevel(data2_neo,index=metrics))
orig

#Calculate the difference in the chosen metrics between pairs of randomized networks
Nperm = 999
i=1
randomized.patef=matrix(nrow=length(orig),ncol=Nperm+1)
row.names(randomized.patef)=names(orig)
randomized.patef[,1]=orig

while(i <=Nperm){ 
  
  data_aleat=permatfull(data_neo,fixedmar="both",mtype="count",times=1)
  data_aleat=data_aleat$perm[[1]]
  data2_aleat=permatfull(data2_neo,fixedmar="both",mtype="count",times=1)
  data2_aleat=data2_aleat$perm[[1]]
  linha<-abs(networklevel(data_aleat, index=metrics)-networklevel(data2_aleat, index=metrics))
  randomized.patef[,i+1]=linha
  print(i)
  i=i+1
  
} 

randomized.patef

#Plot the observed difference against the distribution of randomized pairwise differences
niveis<-row.names(randomized.patef)
for(k in niveis)
	{
		if(any(is.na(randomized.patef[k,]) == TRUE))
			{
			print("k tem NA")
			} else {
	nome.arq<- paste("Hist_DIFERENCES_patef_null_",k,".png", sep="")
	png(filename= nome.arq, res= 300, height= 15, width=21, unit="cm")
	plot(density(randomized.patef[k,]), main="Observed vs. randomized",)
	abline(v=orig[k], col="red", lwd=2, xlab="")
	dev.off()
	print(k)
	nome.arq<- paste("DIFERENCES_patef_Null_mean_sd_",k,".txt", sep="")
	write.table(cbind(mean(randomized.patef[k,]),sd(randomized.patef[k,])), file=paste(nome.arq,sep=""), 
      sep=" ",row.names=TRUE,col.names=FALSE)
		}
	}


#Estimate the P-value
significance.patef=matrix(nrow=nrow(randomized.patef),ncol=3)
row.names(significance.patef)=row.names(randomized.patef)
colnames(significance.patef)=c("p (rand <= orig)", "p (rand >= orig)", "p (rand=orig)")

signif.sup=function(x) sum(x>=x[1])/length(x)
signif.inf=function(x) sum(x<=x[1])/length(x)
signif.two=function(x) ifelse(min(x)*2>1,1,min(x)*2)

significance.patef[,1]=apply(randomized.patef,1,signif.inf)
significance.patef[,2]=apply(randomized.patef,1,signif.sup)
significance.patef[,3]=apply(significance.patef[,-3],1,signif.two)

significance.patef



##############################################################


####6. Comparing QuanBiMo between two networks####

#Patefield (1981) null model
#(It may take long, depending on your computer's processing power)

#Delete all previous objects
rm(list= ls())

#Create the objects
data<-read.table("net1.txt", head=TRUE)
data2<-read.table("net2.txt", head=TRUE)
data_neo=as.matrix(data)
data_neo
data2_neo=as.matrix(data2)
data2_neo

#Calculate QuanBiMo with standarized Q-values (z-values)
mod1 <- computeModules(web=data, steps=1E6)
mod2 <- computeModules(web=data2, steps=1E6)

#Standardize the Q-value of net1 and choose the number of permutations
nulls1 <- nullmodel(data, N=999, method="r2d")
modules.nulls1 <- sapply(nulls1, computeModules)
like.nulls1 <- sapply(modules.nulls1, function(x) x@likelihood)
mod_val1<- (mod1@likelihood - mean(like.nulls1))/sd(like.nulls1)
mod_val1

#Standardize the Q-value of net2 and choose the number of permutations
nulls2 <- nullmodel(data2, N=999, method="r2d")
modules.nulls2 <- sapply(nulls2, computeModules)
like.nulls2 <- sapply(modules.nulls2, function(x) x@likelihood)
mod_val2<- (mod2@likelihood - mean(like.nulls2))/sd(like.nulls2)
mod_val2

#As the Q-values are standardized, values higher than 2 are usually
#significant To estimate the P-value, do a simple Monte Carlo count 
#or run a z-test.

 
#Calculate the difference in standardized Q-values between the two networks
orig = abs(mod_val1-mod_val2)
orig

names(orig)<- "Diferenca em Q padronizado"

#Calculate the difference in the chosen metrics between pairs of randomized networks
Nperm = 999 

randomized.patef=matrix(nrow=length(orig),ncol=Nperm+1)
row.names(randomized.patef)=names(orig)
randomized.patef[,1]=orig

randomized.patef<- as.matrix(randomized.patef)
i<-1

while (i<=Nperm){ 
	data_aleat=permatfull(data_neo,fixedmar="both",mtype="count",times=1)
	data_aleat=data_aleat$perm[[1]]
	data2_aleat=permatfull(data2_neo,fixedmar="both",mtype="count",times=1)
	data2_aleat=data2_aleat$perm[[1]]
	#data_aleat<- as.matrix(data_aleat)
	#data2_aleat<- as.matrix(data2_aleat)
	print("passou")
	mod1_aleat <- try(computeModules(web=data_aleat, steps=1E6)) #1E6 no minimo quando for pra valer
	print("mod 1")
	print(i)
	mod2_aleat <- try(computeModules(web=data2_aleat, steps=1E6)) #1E6
	
	print("mod 2 GO GO GO")
	print(i)
	nulls1_aleat <- nullmodel(data_aleat, N=9, method="r2d")
 	modules.nulls1_aleat <- sapply(nulls1_aleat, computeModules)
 	print("passed nulls aleat1")
	like.nulls1_aleat <- sapply(modules.nulls1, function(x) x@likelihood)

	mod_val1_aleat<- (mod1@likelihood - mean(like.nulls1_aleat))/sd(like.nulls1_aleat)
	#Rede 2
 	nulls2 <- nullmodel(data2_aleat, N=9, method="r2d")
 	modules.nulls2_aleat <- sapply(nulls2, computeModules)
	print("passed nulls aleat2") 
	like.nulls2_aleat <- sapply(modules.nulls2_aleat, function(x) x@likelihood)
	mod_val2_aleat<- (mod2@likelihood - mean(like.nulls2_aleat))/sd(like.nulls2_aleat)
	
	print("passed TRY!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
		
	linha<-abs(mod_val1_aleat-mod_val2_aleat)

	randomized.patef[,i+1]=linha
	
	print(i)
	print(Sys.time())
	i<-i+1
	nome<- paste("RANDOMIZED_values_MOD_patef_sink","_",i,".txt")
	sink(nome)
	cat(randomized.patef)
	sink()
	} 

nometxt<- paste("RANDOMIZED_values_MOD_Patefield_TXT","_",i,".txt")
write.table(randomized.patef,nometxt, sep="\t", quote=F)
nometxt


#Plot the observed difference against the distribution of randomized pairwise differences
niveis<-row.names(randomized.patef)
for(k in niveis)
	{
		if(any(is.na(randomized.patef[k,]) == TRUE))
			{
			print("k tem NA")
			} else {
	nome.arq<- paste("Hist_DIFERENCES_Quanbimo_null_",k,".png", sep="")
	png(filename= nome.arq, res= 300, height= 15, width=21, unit="cm")
	plot(density(randomized.patef[k,]), main="Observed vs. randomized",)
	abline(v=orig[k], col="red", lwd=2, xlab="")
	dev.off()
	print(k)
	nome.arq<- paste("DIFERENCES_patef_Null_mean_sd_QuanBiMo",k,".txt", sep="")
	write.table(cbind(mean(randomized.patef[k,]),sd(randomized.patef[k,])), file=paste(nome.arq,sep=""), 
      sep=" ",row.names=TRUE,col.names=FALSE)
		}
	}

#Estimate the P-value
randomized.patef

significance.patef=matrix(nrow=nrow(randomized.patef),ncol=3)
row.names(significance.patef)=row.names(randomized.patef)
colnames(significance.patef)=c("p (rand <= orig)", "p (rand >= orig)", "p (rand=orig)")

signif.sup=function(x) sum(x>=x[1])/length(x) #unicaudal
signif.inf=function(x) sum(x<=x[1])/length(x) #unicaudal
signif.two=function(x) ifelse(min(x)*2>1,1,min(x)*2)#bicaudal

significance.patef[,1]=apply(randomized.patef,1,signif.inf)
significance.patef[,2]=apply(randomized.patef,1,signif.sup)

significance.patef2<-data.frame(significance.patef)
significance.patef2[,3]=apply(significance.patef2[,-3],1,signif.two)
colnames(significance.patef2)=c("p (rand <= orig)", "p (rand >= orig)", "p (rand=orig)")

significance.patef2




#########################################################################
#######################BE HAPPY :)#######################################
	
####7. Suggested readings##################################

#Barabasi, A.-L. (2016) Network Science, 1st ed. Cambridge University
 #Press, Cambridge.

#Bascompte, J. & Jordano, P. (2014) Mutualistic Networks, 1st ed. Princeton
 #University Press, Princeton.

#Bascompte, J., Jordano, P., Melian, C.J. & Olesen, J.M. (2003) The nested assembly
	#of plant-animal mutualistic networks. Proceedings of the National Academy
	#of Sciences of the United States of America, 100, 9383–9387.

#Bezerra, E.L.S., Machado, I.C. & Mello, M.A.R. (2009) Pollination networks of
  #oil-flowers: a tiny world within the smallest of all worlds. Journal of Animal
  #Ecology, 78, 1096–101.
  #Fonte da rede1 do exemplo (Catimbau).

#Blondel, V.D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E. (2008) Fast 
	#unfolding of communities in large networks. Journal of Statistical Mechanics:
	#Theory and Experiment, 2008, P10008.

#Dormann, C.F. & Strauss, R. (2014) A method for detecting modules in quantitative
	#bipartite networks (ed P Peres-Neto). Methods in Ecology and Evolution, 5, 90–98.

#Patefield, W.M. 1981. Algorithm AS159. An efficient method of generating
	#r x c tables with given row and column totals. Applied Statistics 30, 91-97.

#Vazquez, D.P. & Simberloff, D. 2003. Changes in interaction biodiversity
  #induced by an introduced ungulate. Ecology Letters, 6, 1077-1083.
  #Fonte da rede2 do exemplo (Safariland, parcial).

