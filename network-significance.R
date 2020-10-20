################################################################################
##### SCRIPT FOR ESTIMATING P-VALUES OF NETWORK METRICS
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

#Set the number of permutations to be used in all permutation analyses.
#Permutation is quite resource-demanding. Therefore, when setting
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


#Choose a network-level metric to be passed as an argument to all functions
metric <- c("NODF") #an example using binary NODF

#Calculate the metric for the original network
top <- networklevel(net1,
                     index=metric) 

#Check the value
top

#Calculate metric for the randomized networks
randomized.top <- unlist(sapply(nulls1, networklevel, index = metric))

#Plot the observed value against the distribution of randomized values
par(mar = c(4,4,5,4))
plot(density(randomized.top), main="Observed vs. randomized",
     xlim=c(min((top), min(randomized.top)), 
            max((top), max(randomized.top))))
abline(v=top, col="red", lwd=2, xlab="")

top #observed value
mean(randomized.top) #randomized mean
sd(randomized.top) #randomized SD
(top - mean(randomized.top)) / sd(randomized.top) # Z-value
sum(randomized.top>=(top)) / length(randomized.top) #P randomized > observed
sum(randomized.top<=(top)) / length(randomized.top) #P randomized < observed


################################################################################
##### 3. P-value of a modularity metric
################################################################################


#Modularity needs a slightly different code because it uses a different function
#from the package bipartite.

#Choose a modularity algorithm to be passed as an argument to all functions
algorithm=c("Beckett") #an example using the DIRTLPAwb+ algorithm

#Calculate modularity for the original network
mod <- computeModules(net1, 
                      method = algorithm) 

#Check the value
mod@likelihood

#Extract module membership
part <- bipartite::module2constraints(mod)
row.part <- part[1:nrow(net1)]
col.part <- part[(nrow(net1)+1):(nrow(net1)+ncol(net1))]
length(unique((part))) #number of modules

#Calculate metric for the randomized networks
nullmod <- sapply(nulls1, computeModules, method = algorithm)
modnull <- sapply(nullmod, function(x) x@likelihood)

#Plot the observed value against the distribution of randomized values
par(mar = c(4,4,5,4))
plot(density(modnull), main="Observed vs. randomized",
     xlim=c(min((mod@likelihood), min(modnull)), 
            max((mod@likelihood), max(modnull))))
abline(v=mod@likelihood, col="red", lwd=2, xlab="")

mod@likelihood #observed value
mean(modnull) #randomized mean
sd(modnull) #randomized SD
(mod@likelihood - mean(modnull)) / sd(modnull) # Z-value
sum(modnull>=(mod@likelihood)) / length(modnull) #P randomized > observed
sum(modnull<=(mod@likelihood)) / length(modnull) #P randomized < observed


################################################################################
##### 4. Comparing topology between two networks
################################################################################


#Use the same two networks imported before
net1
net2

#Choose a network-level metric to be passed as an argument to all functions
metric=c("NODF") #an example using binary NODF

#Caculate the same metric for both networks
net1.metric <- networklevel(net1,
                      index = metric) 
net2.metric <- networklevel(net2,
                      index = metric) 

net1.metric
net2.metric

#Calculate the difference in the metric between the two networks
diff <- abs(net1.metric - net2.metric)
diff

#Calculate the same metric for all randomized versions of both networks
randomized.net1.metric <- sapply(nulls1, #the randomized list created before
                                 networklevel, 
                                 index = metric)
randomized.net2.metric <- sapply(nulls2, 
                                 networklevel, 
                                 index = metric)

#Calculate the difference between all pairs of randomized networks
diff.rand <- abs(randomized.net1.metric - randomized.net2.metric)

#Plot the observed value against the distribution of randomized values
par(mar = c(4,4,5,4))
plot(density(diff.rand), main="Observed vs. randomized (pairwise differences)",
     xlim=c(min((diff), min(diff.rand)), 
            max((diff), max(diff.rand))))
abline(v=diff, col="red", lwd=2, xlab="")

#Estimate the P-values
diff #observed difference
mean(diff.rand) #mean randomized differences
sd(diff.rand) #SD randomized differences
(diff - mean(diff.rand)) / sd(diff.rand) # Z-value
sum(diff.rand>=(diff)) / length(diff.rand) #P randomized > observed
sum(diff.rand<=(diff)) / length(diff.rand) #P randomized < observed


################################################################################
##### 5. Comparing modularity between two networks
################################################################################


#Modularity needs a slightly different code because it uses a different function
#from the package bipartite.
                  
#Choose the modularity algorithm to be passed as an argument to all functions
algorithm=c("Beckett") #an example using the DIRTLPAwb+ algorithm

#Calculate modularity for the original networks
mod1 <- computeModules(net1, 
                       method = algorithm)
mod2 <- computeModules(net2, 
                       method = algorithm)

#Check the values
mod1@likelihood #modularity
mod2@likelihood #modularity

#Calculate the difference between both networks
diff <- abs(mod1@likelihood - mod2@likelihood)
diff

#Extract module membership
part1 <- bipartite::module2constraints(mod1)
row.part1 <- part1[1:nrow(net1)]
col.part1 <- part1[(nrow(net1)+1):(nrow(net1)+ncol(net1))]
length(unique((part1))) #number of modules

part2 <- bipartite::module2constraints(mod2)
row.part2 <- part2[1:nrow(net2)]
col.part2 <- part2[(nrow(net2)+1):(nrow(net2)+ncol(net2))]
length(unique((part2))) #number of modules

#Calculate the same metric for all randomized versions of both networks
nullmod1 <- sapply(nulls1, computeModules, method = algorithm)
mod1null <- sapply(nullmod1, function(x) x@likelihood)

nullmod2 <- sapply(nulls2, computeModules, method = algorithm)
mod2null <- sapply(nullmod2, function(x) x@likelihood)

#Calculate the difference between all pairs of randomized networks
diff.rand <- abs(mod1null - mod2null)

#Plot the observed value against the distribution of randomized values
par(mar = c(4,4,5,4))
plot(density(diff.rand), main="Observed vs. randomized",
     xlim=c(min((diff), min(diff.rand)), 
            max((diff), max(diff.rand))))
abline(v=diff, col="red", lwd=2, xlab="")

diff #observed difference
mean(diff.rand) #mean randomized differences
sd(diff.rand) #SD randomized differences
(diff - mean(diff.rand)) / sd(diff.rand) # Z-value
sum(diff.rand>=(diff)) / length(diff.rand) #P randomized > observed
sum(diff.rand<=(diff)) / length(diff.rand) #P randomized < observed


#################################### END #######################################

