#=====================================================
# script03_Phenodat_Single_Env_Analysis_BLUE
# 2022-09-02
#=====================================================

# Aims: For each Single Environment(i.e. Experiment), compute BLUE

# Clean workspace
rm(list=ls())

library(foreach)
library(reshape2)
library(ggplot2)
library(MegaLMM) #Image()
library(lme4)
library(emmeans)
library(admisc) #capture warnings, errors and messages
library(tibble)

Date = "2023-09-15"
DIR_Input="../02_Phenodat_OutlierRemoval/"
DIR_Output="../03_Phenodat_BLUE/"
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)
dir.exists(DIR_Output)

# define global variables
traitnames = rev(c("Silk_DAP_days","Plant_Height_cm","Yield_Mg_ha"))
FamSize_cutoff=50 #tester family size cutoff

# ------------------------------------------------------------------------------
# Data Loading - Training phenotypic data
# ------------------------------------------------------------------------------
list.files(DIR_Input)
ls=load(
  file = sprintf("%s2023-09-10Pheno_Trn_Clean.RData",DIR_Input)
)
ls
# [1] "Pheno_Trn_Clean" 

# I change original Block ID to Replicate::Block in each Env, which is going to 
# be used to compute overlapping genotypes across Replicate::Block below 
Pheno_Trn_Clean$Block=paste(Pheno_Trn_Clean$Replicate,Pheno_Trn_Clean$Block,sep="::")

dim(Pheno_Trn_Clean) #248475     23
Pheno_Trn_Clean[1:3,]
length(unique(Pheno_Trn_Clean$TraitEnvNew)) #811
length(unique(Pheno_Trn_Clean$EnvNew[Pheno_Trn_Clean$TraitName==traitnames[1]]))#302
length(unique(Pheno_Trn_Clean$EnvNew[Pheno_Trn_Clean$TraitName==traitnames[2]]))#278
length(unique(Pheno_Trn_Clean$EnvNew[Pheno_Trn_Clean$TraitName==traitnames[3]]))#231

# ------------------------------------------------------------------------------
# Take a look at EnvNew  Categories
# ------------------------------------------------------------------------------
ctgy_BlkRandom=NULL
dat_modelfit=Pheno_Trn_Clean
ctgy_BlkRandom=foreach(i = 1:length(unique(dat_modelfit$TraitEnvNew)),.combine = "rbind") %do% {
  # extract data for a specific TraitEnvNew
  TraitEnvNew = unique(dat_modelfit$TraitEnvNew)[i]
  dat = dat_modelfit[dat_modelfit$TraitEnvNew==TraitEnvNew,]
  # categorize EnvNew
  if(length(unique(dat$Replicate))==1){
    if(length(unique(dat$Block))==1){
      ctgy="singleRepSingleBlock"
    }else{
      ctgy="singleRepMultiBlocks"
    }
  }else if(length(unique(dat$Replicate))==2){
    numintersect=length(intersect(dat$Hybrid[dat$Replicate==1],
                                  dat$Hybrid[dat$Replicate==2]))
    if(numintersect<4){
      ctgy="singleRepMultiBlocks"
    }else if(length(unique(dat$Block[dat$Replicate==1]))==1 & 
             length(unique(dat$Block[dat$Replicate==2]))==1){
      ctgy="TwoRepsOneBlockEach"
    }else {
      ctgy="TwoRepsMultiBlocks"
    }
  }else if (length(unique(dat$Replicate))==3) {
    if(length(unique(dat$Block[dat$Replicate==1]))==1 & 
       length(unique(dat$Block[dat$Replicate==2]))==1 &
       length(unique(dat$Block[dat$Replicate==3]))==1){
      ctgy="ThreeRepsOneBlockEach"
    }
  }
  print(paste("#################",i,TraitEnvNew,ctgy,"#######################"))
  data.frame(i=i,TraitEnvNew=TraitEnvNew,Category=ctgy)
}

head(ctgy_BlkRandom)
table(ctgy_BlkRandom$Category)
# singleRepMultiBlocks  singleRepSingleBlock ThreeRepsOneBlockEach    TwoRepsMultiBlocks 
# 80                    10                     2                   518 
# TwoRepsOneBlockEach 
# 201
write.csv(ctgy_BlkRandom,sprintf("%s%sExpDsgnCategory_EnvNew_BlkRandom.csv",
                                 DIR_Output,Date),row.names = F)

# ------------------------------------------------------------------------------
# Model fitting - compute trait BLUEs for each genotype
# ------------------------------------------------------------------------------

# empty variable going to be used below
Pheno_Trn_BLUE=NULL;df_BLUE_notes=NULL;
Pheno_Trn_ModFit_Input=NULL # to output Pheno data for model fitting
dat_modelfit=Pheno_Trn_Clean

# model fitting with for loops
start.time <- Sys.time()
for(i in 1:length(unique(dat_modelfit$TraitEnvNew))) {
  # extract data for a specific TraitEnvNew
  TraitEnvNew = unique(dat_modelfit$TraitEnvNew)[i]
  dat = dat_modelfit[dat_modelfit$TraitEnvNew==TraitEnvNew,]
  # categorize EnvNew
  if(length(unique(dat$Replicate))==1){
    if(length(unique(dat$Block))==1){
      ctgy="singleRepSingleBlock"
      }else{
      ctgy="singleRepMultiBlocks"
      }
    }else if(length(unique(dat$Replicate))==2){
    numintersect=length(intersect(dat$Hybrid[dat$Replicate==1],
                                  dat$Hybrid[dat$Replicate==2]))
    if(numintersect<4){
      ctgy="singleRepMultiBlocks"
      }else if(length(unique(dat$Block[dat$Replicate==1]))==1 & 
             length(unique(dat$Block[dat$Replicate==2]))==1){
        ctgy="TwoRepsOneBlockEach"
        }else {
          ctgy="TwoRepsMultiBlocks"
          }
    }else if (length(unique(dat$Replicate))==3) {
      if(length(unique(dat$Block[dat$Replicate==1]))==1 & 
       length(unique(dat$Block[dat$Replicate==2]))==1 &
       length(unique(dat$Block[dat$Replicate==3]))==1){
      ctgy="ThreeRepsOneBlockEach"
    }
  }
  print(paste("#################",i,TraitEnvNew,ctgy,"#######################"))
  
  # convert to factor
  dat$Block=factor(dat$Block)
  dat$Replicate=factor(dat$Replicate)
  dat$Hybrid=factor(dat$Hybrid)
  levels(dat$Block)

  # set up a df to save hybrid means
  hybrid_means=NULL
  hybrid_means = data.frame(dat[match(unique(dat$Hybrid),dat$Hybrid),],
                            Category=ctgy,
                            Prediction = NA)
  
  if(ctgy=="singleRepSingleBlock"){
    # TraitEnvNew="Yield_Mg_ha::MOH1_1_2018::PHT69"
    mod_formula=TraitVal~Hybrid 
    fm = NULL; #empty before using
    message=NULL; message=tryCatchWEM(fm <- lm(mod_formula,data = dat), capture = FALSE)
    if(is.null(message)){
      Message="NULL"
    }else if(names(message)=='message'){
      Message=gsub("\n","",message$message)
    }else if(names(message)=='warning'){
      Message=gsub("\n","",message$warning)
      if(length(message$warning)>1){
        Message=paste(message$warning,collapse = ";")
      }
    }
    # using predict()
    hybrid_means$Prediction = predict(fm,newdata = data.frame(Hybrid=hybrid_means$Hybrid))
    df_notes=NULL; df_notes=data.frame(i=i,
                                       TraitEnvNew=TraitEnvNew,
                                       Category=ctgy,
                                       MessageName=ifelse(is.null(names(message)),"NULL",
                                                          ifelse(length(message)==1,names(message),paste(names(message),collapse = ";"))),
                                       MessageContent=Message)
    
  }else if(ctgy=="singleRepMultiBlocks"){
    if(all(table(dat$Hybrid)==1)){
      fm=lmer(TraitVal ~ (1|Block),data = dat)
      hybrid_means$Prediction=resid(fm)+fixef(fm)
      df_notes=NULL; df_notes=data.frame(i=i,
                                         TraitEnvNew=TraitEnvNew,
                                         Category=ctgy,
                                         MessageName="NULL",
                                         MessageContent="singleRepMultiBlocks,every hybrid has only one observation"
      )
    }else{
      mod_formula=TraitVal ~ Hybrid+ (1|Block)
      fm = NULL;message=NULL
      message=tryCatchWEM(fm <- lmer(mod_formula,data = dat), capture = FALSE)
      if(is.null(message)){
        Message="NULL"
      }else if(names(message)=='message'){
        Message=gsub("\n","",message$message)
      }else if(names(message)=='warning'){
        Message=gsub("\n","",message$warning)
        if(length(message$warning)>1){
          Message=paste(message$warning,collapse = ";")
        }
      }
      # using predict() function if blocks are random:
      hybrid_means$Prediction = predict(fm,newdata = data.frame(Hybrid = hybrid_means$Hybrid),re.form=~0)
      df_notes=NULL; df_notes=data.frame(i=i,
                                         TraitEnvNew=TraitEnvNew,
                                         Category=ctgy,
                                         MessageName=ifelse(is.null(names(message)),"NULL",
                                                            ifelse(length(message)==1,names(message),paste(names(message),collapse = ";"))),
                                         MessageContent=Message
      )
    }
  }else if(ctgy=="TwoRepsOneBlockEach"|ctgy=="ThreeRepsOneBlockEach"){
    mod_formula=TraitVal ~ Replicate+Hybrid
    fm = NULL; #empty before using
    message=NULL; message=tryCatchWEM(fm <- lm(mod_formula,data = dat), capture = FALSE)
    if(is.null(message)){
      Message="NULL"
    }else if(names(message)=='message'){
      Message=gsub("\n","",message$message)
    }else if(names(message)=='warning'){
      Message=gsub("\n","",message$warning)
      if(length(message$warning)>1){
        Message=paste(message$warning,collapse = ";")
      }
    }
    # using predict() function if Replicates are fixed
    predictions = sapply(unique(dat$Replicate),function(rep)
      predict(fm,newdata = data.frame(Hybrid = hybrid_means$Hybrid,Replicate = rep)))
    hybrid_means$Prediction = rowMeans(predictions)
    
    df_notes=NULL; df_notes=data.frame(i=i,
                                       TraitEnvNew=TraitEnvNew,
                                       Category=ctgy,
                                       MessageName=ifelse(is.null(names(message)),"NULL",
                                                          ifelse(length(message)==1,names(message),paste(names(message),collapse = ";"))),
                                       MessageContent=Message
    )
    
  }else if(ctgy=="TwoRepsMultiBlocks"){
    mod_formula=TraitVal ~ Hybrid+Replicate + (1|Replicate:Block)
    fm = NULL;message=NULL
    message=tryCatchWEM(fm <- lmer(mod_formula,data = dat), capture = FALSE)
    if(is.null(message)){
      Message="NULL"
    }else if(length(names(message))==1){
      if(names(message)=='message'){
        Message=gsub("\n","",message$message)
      }else if(names(message)=='warning'){
        Message=gsub("\n","",message$warning)
        if(length(message$warning)>1){
          Message=paste(message$warning,collapse = ";")
        }
      }else if(names(message)=='error'){
        Message=gsub("\n","",message$error)
      }
    }else if(!is.null(message) & length(names(message)) >1){
      Message=paste(names(unlist(message)),unlist(message),collapse = ";")
    }
    
    if(!is.null(fm)){
      # using prediction()
      predictions = sapply(unique(dat$Replicate),function(rep)
        predict(fm,newdata = data.frame(Hybrid = hybrid_means$Hybrid,Replicate = rep),re.form=~0))
      hybrid_means$Prediction = rowMeans(predictions)
      df_notes=NULL; df_notes=data.frame(i=i,
                                         TraitEnvNew=TraitEnvNew,
                                         Category=ctgy,
                                         MessageName=ifelse(is.null(names(message)),"NULL",
                                                            ifelse(length(message)==1,names(message),
                                                                   paste(names(message),collapse = ";"))),
                                         MessageContent=Message)
      }
    }#end of else if
  
  print(head(hybrid_means))
  print(df_notes)
  # ----------------------------------------------------
  # accumulate results across trials
  # de-factor dat before stacking input phenotypic data; otherwise get a warning
  dat$Replicate = as.integer(as.character(dat$Replicate))
  dat$Block = as.character(dat$Block)
  dat$Hybrid = as.character(dat$Hybrid)
  
  if(TraitEnvNew == unique(dat_modelfit$TraitEnvNew)[1]){
    Pheno_Trn_ModFit_Input=dat
    Pheno_Trn_BLUE=hybrid_means
    df_BLUE_notes=df_notes
  }else{
    if(identical(colnames(Pheno_Trn_ModFit_Input),colnames(dat))){
      Pheno_Trn_ModFit_Input=rbind(Pheno_Trn_ModFit_Input,dat)
    }
    Pheno_Trn_BLUE=rbind(Pheno_Trn_BLUE,hybrid_means)
    df_BLUE_notes=rbind(df_BLUE_notes,df_notes)
  }
  
  # empty local variables within each loop
  dat = NULL
  hybrid_means = NULL
  df_notes = NULL
  
}
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken #Time difference of 3.211563 mins

# --------------------------------------------------------------------
# validating
sapply(list(Pheno_Trn_ModFit_Input,Pheno_Trn_BLUE,df_BLUE_notes), dim)
#        [,1]   [,2] [,3]
# [1,] 248475 182503  811
# [2,]     23     25    5
sapply(list(Pheno_Trn_ModFit_Input,Pheno_Trn_BLUE,df_BLUE_notes), function(x)
  length(unique(x$TraitEnvNew)))
# [1] 811 811 811
#is there any NA generated in Pheno_Trn_BLUE?
sum(is.na(Pheno_Trn_BLUE$Prediction)) #0
# what types of message/warnings we got during model fitting/prediction?
table(df_BLUE_notes$MessageName)
# message    NULL warning 
#   169     614      28 
df_BLUE_notes[df_BLUE_notes$MessageName=="message",]
df_BLUE_notes[df_BLUE_notes$MessageName=="warning",]
head(Pheno_Trn_BLUE)
sum(is.na(Pheno_Trn_BLUE$Prediction))

# is Pheno_Trn_ModFit_Input identical to Pheno_Trn_Clean? And why?
identical(Pheno_Trn_ModFit_Input,Pheno_Trn_Clean) #FALSE
sapply(list(Pheno_Trn_ModFit_Input,Pheno_Trn_Clean), dim)
#        [,1]   [,2]
# [1,] 248475 248475
# [2,]     23     23
lapply(list(Pheno_Trn_ModFit_Input,Pheno_Trn_Clean), head)

identical(setdiff(Pheno_Trn_ModFit_Input$TraitEnvNew,Pheno_Trn_Clean$TraitEnvNew),
          character(0)) #TRUE
all.equal(Pheno_Trn_ModFit_Input$TraitEnvNew,Pheno_Trn_Clean$TraitEnvNew) #"59628 string mismatches"
all.equal(unique(Pheno_Trn_ModFit_Input$TraitEnvNew),unique(Pheno_Trn_Clean$TraitEnvNew)) #TRUE
all.equal(table(Pheno_Trn_ModFit_Input$TraitEnvNew),table(Pheno_Trn_Clean$TraitEnvNew)) #TRUE
for(i in 1: length(unique(Pheno_Trn_Clean$TraitEnvNew))){
  print(paste(i,
              all.equal(
                Pheno_Trn_ModFit_Input$TraitEnvNew[Pheno_Trn_ModFit_Input$TraitEnvNew==unique(Pheno_Trn_ModFit_Input$TraitEnvNew)[i]],
                Pheno_Trn_Clean$TraitEnvNew[Pheno_Trn_Clean$TraitEnvNew==unique(Pheno_Trn_Clean$TraitEnvNew)[i]])
              ,sep = "::"
  ))
}

Pheno_Trn_Clean[1:200,c("TraitEnvNew","Replicate","Block","Plot","Range","Pass")]
Pheno_Trn_ModFit_Input[1:200,c("TraitEnvNew","Replicate","Block","Plot","Range","Pass")]
unique(Pheno_Trn_Clean$TraitEnvNew[1:200]) #"Yield_Mg_ha::DEH1_2014::LH185" "Yield_Mg_ha::DEH1_2014::LH195"
unique(Pheno_Trn_ModFit_Input$TraitEnvNew[1:200]) #"Yield_Mg_ha::DEH1_2014::LH185"
# conclusion:
#  - Pheno_Trn_ModFit_Input and Pheno_Trn_Clean are not identical to each other
#  - This is because within each Env (site:year), TraitEnvNew is in different order
#  - if order Pheno_Trn_Clean by TraitEnvNew, they would be identical to each other
for(i in 1: length(unique(Pheno_Trn_Clean$TraitEnvNew))){
  print(paste(i,
              all.equal(
                Pheno_Trn_ModFit_Input[Pheno_Trn_ModFit_Input$TraitEnvNew==unique(Pheno_Trn_ModFit_Input$TraitEnvNew)[i],],
                Pheno_Trn_Clean[Pheno_Trn_Clean$TraitEnvNew==unique(Pheno_Trn_Clean$TraitEnvNew)[i],])
              ,sep = "::"
  ))
}

# --------------------------------------------------------------------
save(Pheno_Trn_ModFit_Input,Pheno_Trn_BLUE,df_BLUE_notes,
     file=sprintf("%s%sPhenodat_Single_Env_Analysis_BLUE_RepFixed_BlkRandom.RData",DIR_Output,Date))
