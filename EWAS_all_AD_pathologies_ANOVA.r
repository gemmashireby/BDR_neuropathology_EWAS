##### Run ANOVA of all AD pathologies and extract stats when controlling for all pathologies 
##### Author:Gemma Shireby

## setwd
setwd("")

## load libraries
library(lmerTest)

## load files 

load("")

## match to betas

betas<-betas[,match(pheno$Basename, colnames(betas))]
dim(betas) 

###### setting up parallel processors

library(doParallel)
cl<-makeCluster(20)
registerDoParallel(cl)
clusterEvalQ(cl, library(lmerTest))

#Create function which performs analysis for a probes

testCpG<-function(row, pheno){
  
  
  modellmer<-lmer(betas[i,] ~  Braak_tangle + cerad_density + Thal_amyloid + Gender + Age + Double + NeuN  + Plate  + (1|Brain_ID) + PC1, data=pheno, REML =F)
  
  model.null<-lmer(betas[i,] ~ Gender + Age + Double + NeuN  + Plate  + (1|Brain_ID) + PC1, data=pheno, REML =F)
  
  P_all<-anova(modellmer,model.null)["modellmer","Pr(>Chisq)"]
  
  ### Braak
  BetaBraak<-fixef(modellmer)["Braak_tangle"]
  SEBraak<-summary(modellmer)$coefficients["Braak_tangle",2]            
  TBraak<-summary(modellmer)$coefficients["Braak_tangle",4]
  model.nullBraak<-lmer(betas[i,] ~ cerad_density + Thal_amyloid + Gender + Age + Double + NeuN  + Plate  + (1|Brain_ID) + PC1, data=pheno, REML =F)
  PBraak<-anova(modellmer,model.nullBraak)["modellmer","Pr(>Chisq)"]

 ### Cerad

  BetaCerad<-fixef(modellmer)["cerad_density"]
  SECerad<-summary(modellmer)$coefficients["cerad_density",2]
  TCerad<-summary(modellmer)$coefficients["cerad_density",4]
  model.nullCerad<-lmer(betas[i,] ~ Braak_tangle  + Thal_amyloid + Gender + Age + Double + NeuN  + Plate  + (1|Brain_ID) + PC1, data=pheno, REML =F)
  PCerad<-anova(modellmer,model.nullCerad)["modellmer","Pr(>Chisq)"]
  
 ### Thal

  BetaThal<-fixef(modellmer)["Thal_amyloid"]
  SEThal<-summary(modellmer)$coefficients["Thal_amyloid",2]
  TThal<-summary(modellmer)$coefficients["Thal_amyloid",4]
  model.nullThal<-lmer(betas[i,] ~ Braak_tangle + cerad_density + Gender + Age + Double + NeuN  + Plate  + (1|Brain_ID)+ PC1, data=pheno, REML =F)
  PThal<-anova(modellmer,model.nullThal)["modellmer","Pr(>Chisq)"]
 
  return(c(P_all,BetaBraak,SEBraak,PBraak,TBraak,BetaCerad,SECerad,PCerad,TCerad,BetaThal,SEThal,PThal,TThal)) # 
}


###Run EWAS using foreach() and %dopar% to tell R to run the analysis is parallel
res<-foreach(i=1:nrow(betas), .combine=rbind) %dopar%{
  testCpG(betas[i,], pheno)	
}


colnames(res)<-c("P_all","Beta_Braak", "SE_Braak", "P_Braak","T_Braak","Beta_Cerad", "SE_Cerad", "P_Cerad","T_Cerad","Beta_Thal", "SE_Thal","P_Thal","T_Thal")
rownames(res)<-rownames(betas)
res<-as.data.frame(res)

#Convert 0-1 proportions to percentages as they are easier to understand like this
res[,c("Beta_Braak","Beta_Cerad","Beta_Thal")]<-res[,c("Beta_Braak","Beta_Cerad","Beta_Thal")]*100
res[,c("SE_Braak","SE_Cerad","SE_Thal")]<-res[,c("SE_Braak","SE_Cerad","SE_Thal")]*100

## annotate
epicManifest<-read.csv("/gpfs/mrc0/projects/Research_Project-MRC190311/References/EPICArray/MethylationEPIC_v-1-0_B4.csv", skip = 7)
epicManifest<-epicManifest[match(rownames(res), epicManifest$Name),]
res<-cbind(res, epicManifest[,c("CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island", "Regulatory_Feature_Name", "Regulatory_Feature_Group", "DNase_Hypersensitivity_NAME", "OpenChromatin_NAME", "TFBS_Evidence_Count", "Methyl450_Loci", "SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")])

## save results
write.csv(res, "")


