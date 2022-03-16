##### script to run EWAS in BDR 
##### AuthorGemma Shireby

args <- commandArgs(trailingOnly = TRUE)
phenotype<-args[1] 

## set working directory 
setwd("")

## load libraries

library(lmerTest)
library(bacon)
library(ggplot2)


## load files 

load("file.rdat") # betas # opheno

## match betas
betas<-betas[,match(pheno$Basename, colnames(betas))]
dim(betas)


## setting up parallel processors

library(doParallel)
cl<-makeCluster(20)
registerDoParallel(cl)
clusterEvalQ(cl, library(lmerTest))

#Create function which performs analysis for a probes

if(class(pheno$phenotype) != "numeric") {
  #Create function which performs analysis for a probes
  
  testCpG<-function(row, pheno){
    
    
    modellmer<-lmer(betas[i,] ~ phenotype + Age + Gender+ NeuN+ Double + Plate + PC1 + (1|Brain_ID), data=pheno)
    Beta<-fixef(modellmer)["phenotype1"]
    SE<-summary(modellmer)$coefficients["phenotype1",2]
    T<-summary(modellmer)$coefficients["phenotype1",4]
    model.null<-lmer(betas[i,] ~ Age + Gender+ NeuN+ Double + Plate + PC1 + (1|Brain_ID), data=pheno)
    P<-anova(modellmer,model.null)["modellmer","Pr(>Chisq)"]
    #P<-summary(modellmer)$coefficients["phenotype1",5]
    return(c(Beta,SE,P,T))
  }
  
  
  ###Run EWAS using foreach() and %dopar% to tell R to run the analysis is parallel
  res<-foreach(i=1:nrow(betas), .combine=rbind) %dopar%{
    testCpG(betas[i,], pheno)	
  }
  
  
} else {
  
  testCpG<-function(row, pheno){
    
    
    modellmer<-lmer(betas[i,] ~ phenotype + Age + Gender + NeuN + Double + Plate + PC1 + (1|Brain_ID), data=pheno)
    Beta<-fixef(modellmer)["phenotype"]
    SE<-summary(modellmer)$coefficients["phenotype",2]
    T<-summary(modellmer)$coefficients["phenotype",4]
    model.null<-lmer(betas[i,] ~ Age + Gender + NeuN + Double + Plate +PC1 + (1|Brain_ID), data=pheno)
    P<-anova(modellmer,model.null)["modellmer","Pr(>Chisq)"]
    #P<-summary(modellmer)$coefficients["phenotype",5]
    return(c(Beta,SE,P,T))
    
  }
  

  ###Run EWAS using foreach() and %dopar% to tell R to run the analysis is parallel
  res<-foreach(i=1:nrow(betas), .combine=rbind) %dopar%{
    testCpG(betas[i,], pheno)	
  }
  
}


colnames(res)<-c("Beta", "SE", "P","T")
rownames(res)<-rownames(betas)
res<-as.data.frame(res)

#Convert 0-1 proportions to percentages as they are easier to understand like this
res[,"Beta"]<-res[,"Beta"]*100
res[,"SE"]<-res[,"SE"]*100


# annotate

epicManifest<-read.csv("MethylationEPIC_v-1-0_B4.csv", skip = 7)
epicManifest<-epicManifest[match(rownames(res), epicManifest$Name),]
res<-cbind(res, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")])


#Sorting results by P.value
res<-res[order(res$P),]
## save

write.csv(res,paste0("",phenotype,".csv", sep=""))


## bacon
res_bacon <- bacon(teststatistics = NULL, effectsizes = res$Beta, standarderrors = res$SE)
res$Pbacon<-pval(res_bacon)
res$SEbacon<-se(res_bacon)
res$Betabacon<-es(res_bacon)

# save

write.csv(res,paste0(",phenotype,".csv", sep=""))


