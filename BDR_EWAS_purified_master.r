##### master script to run EWAS in BDR single cell data
##### Author:Gemma Shireby
##### Date: 13-04-2021

#################################################
###### set up args for command line #############
#################################################

args <- commandArgs(trailingOnly = TRUE)
cell_type<-args[1]
phenotype<-args[2]
filename<-args[3] 

# ensure character structure
cell_type<-as.character(cell_type)
filename<-as.character(filename)
phenotype<-as.character(phenotype)

#################################################
############# set working directory #############
#################################################

setwd("")

##########################################
############# load libraries #############
##########################################

library(qqman)
library(bacon)

##########################################
############# load files #################
##########################################


load()
pheno<-SampleSheet
rm(SampleSheet)

### reformat cell type column (no spaces)
summary(as.factor(pheno$Cell_Type))
pheno$Cell_Type<-gsub(" ","",pheno$Cell_Type)
pheno$Cell_Type<-as.factor(pheno$Cell_Type)

summary(as.factor(pheno$Cell_Type))
#DOUBLENEG     IRF8+     IRF8-     NEUN+    SOX10+     TOTAL
#       21         3         2        27        28        26


## add braak numeric column

pheno$Braak_tangle<-pheno$Phenotype
pheno$Braak_tangle<-as.character(pheno$Braak_tangle)
pheno$Braak_tangle<-gsub("Braak ","",pheno$Braak_tangle)
pheno$Braak_tangle<-as.numeric(as.character(pheno$Braak_tangle))
summary(as.factor(pheno$Braak_tangle))
summary(as.factor(pheno$Phenotype))

## add status columns

pheno$AD<-ifelse(pheno$Braak_tangle < 3, paste(0),paste(1))
pheno$AD<-as.factor(pheno$AD)
summary(as.factor(pheno$AD))
 # 0  1
# 53 54

## ensure correct structures for covariates 

pheno$Plate<-as.factor(pheno$Plate)
pheno$Sentrix_ID <-as.factor(pheno$Sentrix_ID)
pheno$Age <-as.numeric(pheno$Age)
pheno$Sex <-as.factor(pheno$Sex)


## subset cell type and match to betas 
pheno<-pheno[which(pheno$Cell_Type ==cell_type),]
betas<-betas[,match(pheno$Basename, colnames(betas))]

## subset phenotype and combine to pheno - need additional phenotype col as 'pheno[,phenotype]' not working with the clustering function and need it to work with args 
pheno$phenotype<-pheno[,phenotype]

## setting up parallel processors
library(doParallel)
cl<-makeCluster(10)
registerDoParallel(cl)

# Create clustering functions within if statement depending on data structure 

if(class(pheno$phenotype) != "numeric") { 

# Create function which performs analysis for a probe
testCpG<-function(row, pheno){
    model<-lm(betas[i,] ~ phenotype + Age + Sex + Sentrix_ID, data=pheno)
  return(summary(model)$coefficients["phenotype1",c(1,2,4,3)])
}



# Run EWAS using foreach() and %dopar% to tell R to run the analysis is parallel
res<-foreach(i=1:nrow(betas), .combine=rbind) %dopar%{
  testCpG(betas[i,], pheno)
}


} else {
 
# Create function which performs analysis for a probe 
testCpG<-function(row, pheno){
  model<-lm(betas[i,] ~ phenotype + Age + Sex + Sentrix_ID, data=pheno)
  return(summary(model)$coefficients["phenotype",c(1,2,4,3)])
}


# Run EWAS using foreach() and %dopar% to tell R to run the analysis is parallel
res<-foreach(i=1:nrow(betas), .combine=rbind) %dopar%{
  testCpG(betas[i,], pheno)
}


}


### name res
colnames(res)<-c("Beta", "SE", "P","T")
rownames(res)<-rownames(betas)

#Convert 0-1 proportions to percentages as they are easier to understand like this
res[,"Beta"]<-res[,"Beta"]*100
res[,"SE"]<-res[,"SE"]*100


## read annotation file and annotate
epicManifest<-read.csv("MethylationEPIC_v-1-0_B4.csv", skip = 7)
epicManifest<-epicManifest[match(rownames(res), epicManifest$Name),]
res<-cbind(res, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")])

#Sorting results by P.value
res<-res[order(res$P),]

## save

write.csv(res,paste0("purified_cell_normalisedwithincell_BDR_EWAS_Age_Sex_SentrixID_",filename,".csv", sep=""))