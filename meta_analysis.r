## Meta analysis of braak NFT stage in R -  extract z values, SE and apply bacon
## Author: Gemma Shireby
## date: 29.03.2021

#################### 
## load libraries ##
#################### 

library(meta)
library(dplyr)
library(bacon)

###############
## functions ##
###############


extractFromList<-function(mat, row = NULL, col = NULL)
{	if(is.null(row) & !is.null(col)){		return(mat[,col])	}	
if(!is.null(row) & is.null(col)){		
return(mat[row,col])	}	
if(!is.null(row) & !is.null(col)){		
return(mat[row,col])	}	
if(is.null(row) & is.null(col)){		
return(mat)	}}


## lambda genomic inflation factor function
LambdaInf<-function(pvals){ # pvals = vector of p values
  chisq <- qchisq(1-pvals,1)
  lambda = median(chisq)/qchisq(0.5,1)
  print(lambda) 
}

###########
## setwd ##
###########

setwd()

###################################### 
## read in all meta results as list ##
######################################

dataFiles <- lapply(Sys.glob("*.assoc"), read.table, h=T, row.names=1)
Sys.glob("*.assoc") # check order of names for labels below 
names(dataFiles)<-c("AZ1","AZ2","BDR","LBB1","LBB2","MS","ROSMAP")

#######################################################
## minimum cohorts to consider for the meta analysis ##
#######################################################

ncohorts<-length(dataFiles) # number of cohors 
mincohorts<-2 # minimum number of cohorts probe present in to include in analysis (e.g. 2 = at least 2 cohorts have that probe)

#########################################################
## create vector of probes to include in meta-analysis ##
#########################################################

probeCount<-table(unlist(lapply(dataFiles, rownames)))
probes<-names(probeCount[which(probeCount >= mincohorts)])


#######################
## run meta-analysis ##
#######################

require(parallel)

Braak_meta <- function(i){

	probe<-probes[i]
	effects<-unlist(lapply(dataFiles, extractFromList, probe, "BETA"))
	ses<-unlist(lapply(dataFiles, extractFromList, probe, "SE"))
	effects<-effects[!is.na(effects)]
	ses<-ses[!is.na(ses)]
	out<-metagen(effects, ses)
	res<-c(out$TE.fixed,out$pval.fixed, out$statistic.fixed, out$seTE.fixed, out$lower.fixed, out$upper.fixed, out$TE.random, out$pval.random, out$statistic.random, out$seTE.random, out$lower.random, out$upper.random, out$I2, 1-pchisq(out$Q, out$df.Q))
	return(res)

}

###Run EWAS using foreach() and %dopar% to tell R to run the analysis is parallel

Out <- mclapply(1:length(probes), Braak_meta,mc.cores=20) # returns list
res <- do.call(rbind,Out) # make list a matrix

res<-as.data.frame(res)
rownames(res)<-probes
colnames(res)<-c("Effect_fixed", "P_fixed", "Z_fixed", "SE.fixed", "lower.fixed_CI", "upper.fixed_CI", "Effect_random", "P_random", "Z_random","SE.random", "lower.random_CI", "upper.random_CI" ,"I2", "Het P") # must have SE and Z values for bacon

write.csv(res, "")

## Check inflation
LambdaInf(res$P_fixed) #
LambdaInf(res$P_random) # 