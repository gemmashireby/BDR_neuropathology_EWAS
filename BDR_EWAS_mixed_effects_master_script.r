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
library(sva)

## load files 

load("file.rdat") # betas # opheno

## match betas
betas<-betas[,match(pheno$Basename, colnames(betas))]
dim(betas)


## run regression 
if(class(pheno[,phenotype]) != "numeric") {
  
  for(i in which(is.na(res$Beta))){
      modellmer<-lmer(betas[i,] ~  pheno[,phenotype] + pheno$Gender +  pheno$NeuN + pheno$Double  + pheno$Plate + pheno$Age +  PC$PC1  + (1|pheno$Brain_ID), REML =F, na.action=na.exclude)
      res[i,1]<-fixef(modellmer)["pheno[, phenotype]1"]
      res[i,2]<-summary(modellmer)$coefficients["pheno[, phenotype]1",2]
      res[i,3]<-summary(modellmer)$coefficients["pheno[, phenotype]1",5]
      res[i,4]<-summary(modellmer)$coefficients["pheno[, phenotype]1",4]

    
  }
  
} else {
  
  for(i in which(is.na(res$Beta))){
      modellmer<-lmer(betas[i,] ~ pheno[,phenotype]  + pheno$Gender +  pheno$NeuN + pheno$Double + pheno$Plate + pheno$Age + PC$PC1 + (1|pheno$Brain_ID), na.action=na.exclude)
      res[i,1]<-fixef(modellmer)["pheno[, phenotype]"]
      res[i,2]<-summary(modellmer)$coefficients["pheno[, phenotype]",2]
      res[i,3]<-summary(modellmer)$coefficients["pheno[, phenotype]",5]
      res[i,4]<-summary(modellmer)$coefficients["pheno[, phenotype]",4]
    
  }
  
}

######################################################################################
######################################################################################

#Convert 0-1 proportions to percentages as they are easier to understand like this
res<-as.data.frame(res)
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


