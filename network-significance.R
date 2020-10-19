################################################################################
##### SCRIPT FOR ESTIMATING P-VALUES FOR NETWORK METRICS
################################################################################


##### Ecological Synthesis Lab (SintECO)
##### https://marcomellolab.wordpress.com
##### Authors: Renata Muylaert, Pavel Dodonov & Marco Mello
##### E-mail: renatamuy@gmail.com
##### See README for further info:
##### https://github.com/marmello77/network-significance/blob/master/README.md


################################################################################
##### Summary
################################################################################


#1. Get ready
#2. P-value of a topology metric
#3. P-value of a modularity metric
#4. Comparing topology between two networks
#5. Comparing modularity between two networks


################################################################################
##### 1. Get ready
################################################################################


#Set the working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Delete all previous objects
rm(list= ls())

#Load the packages
library(bipartite)
library(reshape2)
library(tidyverse)

#Import the data
net1 <- as.matrix(read.delim("net1.txt", row.names=1))
net2 <- as.matrix(read.delim("net2.txt", row.names=1))

#Inspect the networks
class(net1)
net1
dim(net1)
min(net1)
max(net1)

class(net2)
net2
dim(net2)
min(net2)
max(net2)

#Plot the matrices
visweb(net1)
visweb(net2)

#Set the number of permutations to be used in all null model analyses
#Permutations analysis is quite resource-demanding. Therefore, when setting
#this number, take into account your computer's power and memory. 
#When doing this kind of analysis for real in a paper, we use at least 1,000
#permutations. Here we set this number lower just for testing the script.
permutations <- 10

#Set seed to allow reproducible results
set.seed(14)

#Generate randomized versions of each network using a null model
nulls1 <- nullmodel(net1, 
                    N=permutations, #Use the number set before
                    method="vaznull") #You can choose a different null model

nulls2 <- nullmodel(net2, 
                    N=permutations, 
                    method="vaznull")


################################################################################
##### 2. P-value of a topology metric
################################################################################


#Calculate a metric of your choice for the original network
Nest <- networklevel(net1,
                     index="NODF") #Choose any network-level metric

#Check the value
Nest

#Calculate metric for the randomized networks
randomized.Nest <- unlist(sapply(nulls1, networklevel, index="NODF"))
(Nest - mean(randomized.Nest))/sd(randomized.Nest) # Z value
Nest.sig <- sum(randomized.Nest>Nest)/length(randomized.Nest) # p value

#Plot the observed value against the distribution of randomized values
par(mar = c(4,4,5,4))
plot(density(randomized.Nest), main="Observed vs. randomized",
     xlim=c(min((Nest), min(randomized.Nest)), 
            max((Nest), max(randomized.Nest))))
abline(v=Nest, col="red", lwd=2, xlab="")

Nest #observed value
mean(randomized.Nest) #randomized mean
sd(randomized.Nest) #randomized SD
(Nest - mean(randomized.Nest))/sd(randomized.Nest) # Z-value
sum(randomized.Nest>(Nest)) / length(randomized.Nest) #P randomized > observed
sum(randomized.Nest<(Nest)) / length(randomized.Nest) #P randomized < observed


################################################################################
##### 3. P-value of a modularity metric
################################################################################


#Modularity needs a slightly different code because it uses a different function
#from the package bipartite.

#Choose the modularity algorithm
algorithm=c("Beckett") #this is the DIRTLPAwb+ algorithm

#Calculate modularity for the original network
Mod <- computeModules(net1, 
                      method = algorithm) #in this case DIRTLPAwb+

#Check the value
Mod@likelihood

#Extract module membership
Part <- bipartite::module2constraints(Mod)
row.Part <- Part[1:nrow(net1)]
col.Part <- Part[(nrow(net1)+1):(nrow(net1)+ncol(net1))]

#Calculate metric for the randomized networks
nullmod <- sapply(nulls1, computeModules, method = algorithm)
modnull <- sapply(nullmod, function(x) x@likelihood)
(Mod@likelihood - mean(modnull))/sd(modnull) # Z value
Mod.sig <- sum(modnull>(Mod@likelihood)) / length(modnull) # p value

#Plot the observed value against the distribution of randomized values
par(mar = c(4,4,5,4))
plot(density(modnull), main="Observed vs. randomized",
     xlim=c(min((Mod@likelihood), min(modnull)), 
            max((Mod@likelihood), max(modnull))))
abline(v=Mod@likelihood, col="red", lwd=2, xlab="")

Mod@likelihood #observed value
mean(modnull) #randomized mean
sd(modnull) #randomized SD
(Mod@likelihood - mean(modnull))/sd(modnull) # Z-value
sum(modnull>(Mod@likelihood)) / length(modnull) #P randomized > observed
sum(modnull<(Mod@likelihood)) / length(modnull) #P randomized < observed


################################################################################
##### 4. Comparing topology between two networks
################################################################################


#Use the same two networks imported before
net1
net2

#Choose the network-level metric
metrics=c("NODF")

#Calculate the difference in the chosen metric between the observed networks
orig1 = abs(networklevel(net1,index=metrics)
         -networklevel(net2,index=metrics))

#Check the difference
orig1

#Calculate the difference between pairs of randomized networks
randomized1=matrix(nrow=length(orig1),ncol=permutations+1)
row.names(randomized1)=names(orig1)
randomized1[,1]=orig1
randomized1

i=1

while(i <= permutations){ 
  
  net_rand1=permatfull(net1,fixedmar="both",mtype="count",times=1)
  net_rand1=net_rand1$perm[[1]]
  net_rand2=permatfull(net2,fixedmar="both",mtype="count",times=1)
  net_rand2=net_rand2$perm[[1]]
  lines<-abs(networklevel(net_rand1, index=metrics)-networklevel(net_rand2, index=metrics))
  randomized1[,i+1]=lines
  print(i)
  i=i+1
  
} 

randomized1

#Plot the observed difference against the distribution of randomized differences
levels<-row.names(randomized1)
for(k in levels){
		if(any(is.na(randomized1[k,]) == TRUE))
			{
			print("k tem NA")
			} else {
	plot(density(randomized1[k,]), main="Observed vs. randomized",)
	abline(v=orig1[k], col="red", lwd=2, xlab="")
		}
	}

#Estimate the P-value
orig1 #observed difference
mean(randomized1) #mean randomized differences
sd(randomized1) #SD randomized differences
(orig1 - mean(randomized1))/sd(randomized1) # Z-value
sum(randomized1>(orig1)) / length(randomized1) #P randomized > observed
sum(randomized1<(orig1)) / length(randomized1) #P randomized < observed


################################################################################
##### 5. Comparing modularity between two networks
################################################################################


#Use the same two networks imported before
net1
net2

#Choose the modularity algorithm
algorithm=c("Beckett") #this is the DIRTLPAwb+ algorithm

#Calculate modularity for the original networks
mod1 <- computeModules(net1, 
                      method = algorithm) 
mod2 <- computeModules(net2, 
                       method = algorithm)

#Use the same randomized versions created before
nulls1
nulls2

#Standardize the Q-value (modularity)
modules.nulls1 <- sapply(nulls1, computeModules,  method = algorithm)
like.nulls1 <- sapply(modules.nulls1, function(x) x@likelihood)
mod_val1<- (mod1@likelihood - mean(like.nulls1))/sd(like.nulls1)
mod_val1

modules.nulls2 <- sapply(nulls2, computeModules,  method = algorithm)
like.nulls2 <- sapply(modules.nulls2, function(x) x@likelihood)
mod_val2<- (mod2@likelihood - mean(like.nulls2))/sd(like.nulls2)
mod_val2

#As the Q-values are standardized, values higher than 2 are usually
#significant. To estimate the P-value, do a simple Monte Carlo count.

#Calculate the difference in standardized Q-values between the two networks
orig2 = abs(mod_val1-mod_val2)
names(orig2)<- "Difference in standardized Q"
orig2

#Calculate the difference in the chosen metrics between pairs of randomized networks
randomized2=matrix(nrow=length(orig2),ncol=permutations+1)
row.names(randomized2)=names(orig2)
randomized2[,1]=orig2
randomized2<- as.matrix(randomized2)
randomized2

i<-1

while (i<=permutations){ 
	net_rand1=permatfull(net1,fixedmar="both",mtype="count",times=1)
	net_rand1=net_rand1$perm[[1]]
	
	net_rand2=permatfull(net2,fixedmar="both",mtype="count",times=1)
	net_rand2=net_rand2$perm[[1]]

	print("passed")
	mod1_rand <- try(computeModules(net_rand1, method = algorithm))
	print("mod 1")
	print(i)
	
	mod2_rand <- try(computeModules(net_rand2, method = algorithm))
	print("mod 2")
	print(i)
	
	nulls1_rand <- nullmodel(net_rand1, N=permutations, method="vaznull")
 	modules.nulls1_rand <- sapply(nulls1_rand, 
 	                              computeModules, method = algorithm)
 	print("passed nulls rand1")
	like.nulls1_rand <- sapply(modules.nulls1, function(x) x@likelihood)
	mod_val1_rand<- (mod1@likelihood - 
	                     mean(like.nulls1_rand))/sd(like.nulls1_rand)

	nulls2_rand <- nullmodel(net_rand1, N=permutations, method="vaznull")
	modules.nulls2_rand <- sapply(nulls2_rand, 
	                              computeModules, method = algorithm)
	print("passed nulls rand2")
	like.nulls2_rand <- sapply(modules.nulls1, function(x) x@likelihood)
	mod_val1_rand<- (mod1@likelihood - 
	                     mean(like.nulls2_rand))/sd(like.nulls2_rand)
	
	print("passed TRY")
		
	lines<-abs(mod_val1_rand-mod_val2_rand)

	randomized2[,i+1]=lines
	
	print(i)
	print(Sys.time())
	i<-i+1
	name<- paste("RANDOMIZED_values_MOD_patef_sink","_",i,".txt")
	} 

nometxt<- paste("RANDOMIZED_values_MOD_Patefield_TXT","_",i,".txt")
nometxt

#Plot the observed difference against the distribution of randomized pairwise differences
levels<-row.names(randomized2)
for(k in levels)
	{
		if(any(is.na(randomized2[k,]) == TRUE))
			{
			print("k tem NA")
			} else {
	plot(density(randomized2[k,]), main="Observed vs. randomized",)
	abline(v=orig2[k], col="red", lwd=2, xlab="")
		}
	}

#Estimate the P-value
orig2 #observed difference
mean(randomized2) #mean randomized differences
sd(randomized2) #SD randomized differences
(orig2 - mean(randomized2))/sd(randomized2) # Z-value
sum(randomized2>(orig2)) / length(randomized2) #P randomized > observed
sum(randomized2<(orig2)) / length(randomized2) #P randomized < observed


#################################### END #######################################

