# ===================================================
# script09_summarize_prediction_results
# ===================================================

# Clean workspace
rm(list=ls())

library(lineup)
library(MegaLMM)
library(foreach)
library(tibble)
library(ggplot2)
library(ggbeeswarm)
source("./custom_functions.R")
options(max.print = 9999)

Date="2023-10-23"
DIR_Inputs = c("../07_G2F_MegaLMM/","../07b_G2F_MegaLMM/","../07c_G2F_MegaLMM/")
DIR_Output="../09_Prediction_Results_Summary/"
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)
CVScns = c("CV1","CV2","CV3NewTrial","CV3NewState","CV3NewTester","CV4")
minTrnEnvs=1

# ------------------------------------------------------------------------------
# res: Prepare all results except for CV3/CV4 - Xenv=="Zero" 
# ------------------------------------------------------------------------------
# for testing
DIR_Input = DIR_Inputs[1]
i=301
res=NA
res=
foreach(DIR_Input = DIR_Inputs,.combine = "rbind") %do% {
  (files=list.files(DIR_Input,pattern = ".RData"))
  stopifnot(length(files)==1495)
  # (files=files[grepl("Fold00_CV1|Fold00_CV2",files)])
  foreach(i = 1:length(files),.combine = "rbind") %do% {
  # foreach(i = 1:3,.combine = "rbind") %do% {
    # file="G2F_PlantHeightcm_Xenv_Zero_CV2_Fold00_CV1_Fold05.RData"
    file = files[i]
    print(file)
    
    # --------------------------------------------------------------------------
    # extract meta info
    (traitname = sapply(strsplit(file,"_"), "[",2))
    (Xenv = strsplit( sapply(strsplit(file,"_Xenv_"), "[",2), "_") [[1]][1])
    (CVScn = sapply(strsplit(file,"_"), "[",7))
    (ID_Fold = as.numeric(gsub("\\D", "", sapply(strsplit(file,"_"), "[",6))))
    (ID_Fold_Level2 = as.numeric(gsub("\\D", "", sapply(strsplit(file,"_"), "[",8))))
    
    # --------------------------------------------------------------------------
    # load prediction results
    ls = load(file = paste0(DIR_Input,file))
    sapply(results, dim)
    #      Y_Test Y_Train U_Train U_Test  Eta_mean Eta_test rrBLUP
    # [1,]   1700    1700    1700   1700   1700     1700     1700   1700
    # [2,]    302     302     302    302    302      302      302    302
  
    # --------------------------------------------------------------------------
    # extract results matrix
    Y_Test = results$Y_Test
    Y_Train = results$Y_Train
    U_Train = results$U_Train
    U_Test = results$U_Test
    U_Pred = results$U_Pred
    Eta_mean = results$Eta_mean
    Eta_test = results$Eta_test
    rrBLUP = results$rrBLUP
    if(!is.null(Eta_mean)) rownames(Eta_mean)=rownames(Y_Train)
    if(!is.null(Eta_test)) rownames(Eta_test)=rownames(Y_Test)
    
    # get MegaLMM_Pred
    if(CVScn=="CV1"){
      MegaLMM_Pred = U_Pred
    }else if(CVScn=="CV2"){
      MegaLMM_Pred = Eta_mean
    }else if(substr(CVScn,1,3)=="CV3") {
      if(Xenv=="Zero"){
        MegaLMM_Pred = Eta_mean
      }else{
        MegaLMM_Pred = Eta_test
      }
    }else if(CVScn=="CV4"){
      if(Xenv=="Zero"){
        MegaLMM_Pred = U_Pred
      }else{
        MegaLMM_Pred = U_Test
      }
    }
    
    # --------------------------------------------------------------------------
    # Prepare General and Specific prediction matrix
    cor_MegaLMM_Gen=NULL
    cor_MegaLMM_Spec=NULL
    cor_GBLUP_Gen=NULL
    df_GBLUP_Spec_State=NULL
    df_GBLUP_Spec_Tester=NULL
    df_GBLUP_Spec_StatePlusTester=NULL
    
    MegaLMM_Gen=NULL
    MegaLMM_Spec=NULL
    MegaLMM_Spec_State=NULL
    MegaLMM_Spec_Tester=NULL
    MegaLMM_Spec_StatePlusTester=NULL
    
    # prep MegaLMM_Gen
    MegaLMM_Gen = average_by_all(Y_Test,MegaLMM_Pred) #S::O/T::O/ST::O/O::O
    # prep MegaLMM_Spec
    if(Xenv=="Zero"){
      MegaLMM_Spec_State=average_by_state(Y_Test,MegaLMM_Pred,minTrnEnvs) #O::S
      MegaLMM_Spec_Tester=average_by_tester(Y_Test,MegaLMM_Pred,minTrnEnvs) #O::T
      MegaLMM_Spec_StatePlusTester=average_by_State_Plus_Tester(Y_Test,MegaLMM_Pred,minTrnEnvs) #O::ST
      
      GBLUP_Gen = average_by_all(Y_Test, rrBLUP) #G::O
      GBLUP_Spec_State=average_by_state(Y_Test,rrBLUP,minTrnEnvs) #G::S
      GBLUP_Spec_Tester=average_by_tester(Y_Test,rrBLUP,minTrnEnvs) #G::T
      GBLUP_Spec_StatePlusTester=average_by_State_Plus_Tester(Y_Test,rrBLUP,minTrnEnvs) #G::ST
      
      # calculate prediction accuracy
      if(CVScn %in% c("CV1","CV2")){
        cor_GBLUP_Spec_Env = corbetw2mat(Y_Test, rrBLUP)
        cor_MegaLMM_Spec_Env = corbetw2mat(Y_Test, MegaLMM_Pred)
      }
      
      cor_GBLUP_Gen  =corbetw2mat(Y_Test,GBLUP_Gen)
      cor_GBLUP_Spec_State=corbetw2mat(Y_Test,GBLUP_Spec_State)
      cor_GBLUP_Spec_Tester=corbetw2mat(Y_Test,GBLUP_Spec_Tester)
      cor_GBLUP_Spec_StatePlusTester=corbetw2mat(Y_Test,GBLUP_Spec_StatePlusTester)
      
      cor_MegaLMM_Gen=corbetw2mat(Y_Test,MegaLMM_Gen)
      cor_MegaLMM_Spec_State=corbetw2mat(Y_Test,MegaLMM_Spec_State)
      cor_MegaLMM_Spec_Tester=corbetw2mat(Y_Test,MegaLMM_Spec_Tester)
      cor_MegaLMM_Spec_StatePlusTester=corbetw2mat(Y_Test,MegaLMM_Spec_StatePlusTester)
    }else{
      if(CVScn %in% c("CV1","CV2")){
        if(Xenv=="State"){ # S::S
          MegaLMM_Spec=average_by_state(Y_Test,MegaLMM_Pred,minTrnEnvs)
        }else if(Xenv=="Tester"){ # T::T
          MegaLMM_Spec=average_by_tester(Y_Test,MegaLMM_Pred,minTrnEnvs)
        }else if(Xenv=="State+Tester"){ #ST::ST
          MegaLMM_Spec=average_by_State_Plus_Tester(Y_Test,MegaLMM_Pred,minTrnEnvs)
        }
      }else{ #CVScn != c("CV1","CV2")
        MegaLMM_Spec = MegaLMM_Pred
      }
    }
    
    # --------------------------------------------------------------------------    
    # accumulate results
    if(Xenv=="Zero"){
      df_GBLUP_Gen=data.frame(file=gsub(".RData","",file),
                              TraitName=traitname,
                              CVScn=CVScn,
                              ID_Fold=ID_Fold,
                              ID_Fold_Level2=ID_Fold_Level2,
                              EnvName=colnames(results$Y_Test),
                              NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                              Method="GBLUP",
                              PriorType="GBLUP",
                              #GBLUP does not use prior but should put Xenv here to trace its source
                              PriorLevel=Xenv,
                              PredictionType="General",
                              PredictionModel="General",
                              accu_r=cor_GBLUP_Gen,
                              accu_z=R_to_Z(cor_GBLUP_Gen))
      df_GBLUP_Spec_State=data.frame(file=gsub(".RData","",file),
                                     TraitName=traitname,
                                     CVScn=CVScn,
                                     ID_Fold=ID_Fold,
                                     ID_Fold_Level2=ID_Fold_Level2,
                                     EnvName=colnames(results$Y_Test),
                                     NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                                     Method="GBLUP",
                                     PriorType="GBLUP",
                                     #GBLUP does not use prior but should put Xenv here to trace its source
                                     PriorLevel=Xenv, 
                                     PredictionType="Specific",
                                     PredictionModel="State",
                                     accu_r=cor_GBLUP_Spec_State,
                                     accu_z=R_to_Z(cor_GBLUP_Spec_State))
      df_GBLUP_Spec_Tester=data.frame(file=gsub(".RData","",file),
                                      TraitName=traitname,
                                      CVScn=CVScn,
                                      ID_Fold=ID_Fold,
                                      ID_Fold_Level2=ID_Fold_Level2,
                                      EnvName=colnames(results$Y_Test),
                                      NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                                      Method="GBLUP",
                                      PriorType="GBLUP",
                                      #GBLUP does not use prior but should put Xenv here to trace its source
                                      PriorLevel=Xenv, 
                                      PredictionType="Specific",
                                      PredictionModel="Tester",
                                      accu_r=cor_GBLUP_Spec_Tester,
                                      accu_z=R_to_Z(cor_GBLUP_Spec_Tester))
      df_GBLUP_Spec_StatePlusTester=data.frame(file=gsub(".RData","",file),
                                               TraitName=traitname,
                                               CVScn=CVScn,
                                               ID_Fold=ID_Fold,
                                               ID_Fold_Level2=ID_Fold_Level2,
                                               EnvName=colnames(results$Y_Test),
                                               NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                                               Method="GBLUP",
                                               PriorType="GBLUP",
                                               #GBLUP does not use prior but should put Xenv here to trace its source
                                               PriorLevel=Xenv, 
                                               PredictionType="Specific",
                                               PredictionModel="State+Tester",
                                               accu_r=cor_GBLUP_Spec_StatePlusTester,
                                               accu_z=R_to_Z(cor_GBLUP_Spec_StatePlusTester))
      df_MegaLMM_Gen=data.frame(file=gsub(".RData","",file),
                              TraitName=traitname,
                              CVScn=CVScn,
                              ID_Fold=ID_Fold,
                              ID_Fold_Level2=ID_Fold_Level2,
                              EnvName=colnames(results$Y_Test),
                              NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                              Method="MegaLMM",
                              PriorType=ifelse(Xenv=="Zero","No","Yes"),
                              #GBLUP does not use prior but should put Xenv here to trace its source
                              PriorLevel=Xenv,
                              PredictionType="General",
                              PredictionModel="General",
                              accu_r=cor_MegaLMM_Gen,
                              accu_z=R_to_Z(cor_MegaLMM_Gen))
      
      df_MegaLMM_Spec_State=data.frame(file=gsub(".RData","",file),
                                     TraitName=traitname,
                                     CVScn=CVScn,
                                     ID_Fold=ID_Fold,
                                     ID_Fold_Level2=ID_Fold_Level2,
                                     EnvName=colnames(results$Y_Test),
                                     NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                                     Method="MegaLMM",
                                     PriorType=ifelse(Xenv=="Zero","No","Yes"),                                     #GBLUP does not use prior but should put Xenv here to trace its source
                                     PriorLevel=Xenv, 
                                     PredictionType="Specific",
                                     PredictionModel="State",
                                     accu_r=cor_MegaLMM_Spec_State,
                                     accu_z=R_to_Z(cor_MegaLMM_Spec_State))
      df_MegaLMM_Spec_Tester=data.frame(file=gsub(".RData","",file),
                                      TraitName=traitname,
                                      CVScn=CVScn,
                                      ID_Fold=ID_Fold,
                                      ID_Fold_Level2=ID_Fold_Level2,
                                      EnvName=colnames(results$Y_Test),
                                      NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                                      Method="MegaLMM",
                                      PriorType=ifelse(Xenv=="Zero","No","Yes"),                                      #GBLUP does not use prior but should put Xenv here to trace its source
                                      PriorLevel=Xenv, 
                                      PredictionType="Specific",
                                      PredictionModel="Tester",
                                      accu_r=cor_MegaLMM_Spec_Tester,
                                      accu_z=R_to_Z(cor_MegaLMM_Spec_Tester))
      df_MegaLMM_Spec_StatePlusTester=data.frame(file=gsub(".RData","",file),
                                               TraitName=traitname,
                                               CVScn=CVScn,
                                               ID_Fold=ID_Fold,
                                               ID_Fold_Level2=ID_Fold_Level2,
                                               EnvName=colnames(results$Y_Test),
                                               NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                                               Method="MegaLMM",
                                               PriorType=ifelse(Xenv=="Zero","No","Yes"),                                               #GBLUP does not use prior but should put Xenv here to trace its source
                                               PriorLevel=Xenv, 
                                               PredictionType="Specific",
                                               PredictionModel="State+Tester",
                                               accu_r=cor_MegaLMM_Spec_StatePlusTester,
                                               accu_z=R_to_Z(cor_MegaLMM_Spec_StatePlusTester))
      
      if(CVScn %in% c("CV1","CV2")){
        df_GBLUP_Spec_Env=data.frame(file=gsub(".RData","",file),
                                     TraitName=traitname,
                                     CVScn=CVScn,
                                     ID_Fold=ID_Fold,
                                     ID_Fold_Level2=ID_Fold_Level2,
                                     EnvName=colnames(results$Y_Test),
                                     NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                                     Method="GBLUP",
                                     PriorType="GBLUP",
                                     #GBLUP does not use prior but should put Xenv here to trace its source
                                     PriorLevel=Xenv, 
                                     PredictionType="Specific",
                                     PredictionModel="Env",
                                     accu_r= cor_GBLUP_Spec_Env,
                                     accu_z=R_to_Z(cor_GBLUP_Spec_Env))
        
        df_MegaLMM_Spec_Env=data.frame(file=gsub(".RData","",file),
                                       TraitName=traitname,
                                       CVScn=CVScn,
                                       ID_Fold=ID_Fold,
                                       ID_Fold_Level2=ID_Fold_Level2,
                                       EnvName=colnames(results$Y_Test),
                                       NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                                       Method="MegaLMM",
                                       PriorType=ifelse(Xenv=="Zero","No","Yes"),                                   #GBLUP does not use prior but should put Xenv here to trace its source
                                       PriorLevel=Xenv, 
                                       PredictionType="Specific",
                                       PredictionModel="Env",
                                       accu_r= cor_MegaLMM_Spec_Env,
                                       accu_z=R_to_Z(cor_MegaLMM_Spec_Env))
        df  = rbind(df_GBLUP_Gen,df_GBLUP_Spec_State,df_GBLUP_Spec_Tester,df_GBLUP_Spec_StatePlusTester,df_GBLUP_Spec_Env,
                    df_MegaLMM_Gen,df_MegaLMM_Spec_State,df_MegaLMM_Spec_Tester,df_MegaLMM_Spec_StatePlusTester,df_MegaLMM_Spec_Env)
      }else{
        df  = rbind(df_GBLUP_Gen,df_GBLUP_Spec_State,df_GBLUP_Spec_Tester,df_GBLUP_Spec_StatePlusTester,
                    df_MegaLMM_Gen,df_MegaLMM_Spec_State,df_MegaLMM_Spec_Tester,df_MegaLMM_Spec_StatePlusTester)
      }
    }else { ##########---------- Xenv!="Zero" ----------##########
      cor_MegaLMM_Gen  =corbetw2mat(Y_Test,MegaLMM_Gen)
      cor_MegaLMM_Spec  =corbetw2mat(Y_Test,MegaLMM_Spec)
      df_MegaLMM_Gen=data.frame(file=gsub(".RData","",file),
                                TraitName=traitname,
                                CVScn=CVScn,
                                ID_Fold=ID_Fold,
                                ID_Fold_Level2=ID_Fold_Level2,
                                EnvName=colnames(results$Y_Test),
                                NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                                Method="MegaLMM",
                                PriorType=ifelse(Xenv=="Zero","No","Yes"),
                                PriorLevel=Xenv,
                                PredictionType="General",
                                PredictionModel="General",
                                accu_r=cor_MegaLMM_Gen,
                                accu_z=R_to_Z(cor_MegaLMM_Gen)
      )
      df_MegaLMM_Spec=data.frame(file=gsub(".RData","",file),
                                 TraitName=traitname,
                                 CVScn=CVScn,
                                 ID_Fold=ID_Fold,
                                 ID_Fold_Level2=ID_Fold_Level2,
                                 EnvName=colnames(results$Y_Test),
                                 NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                                 Method="MegaLMM",
                                 PriorType=ifelse(Xenv=="Zero","No","Yes"),
                                 PriorLevel=Xenv,
                                 PredictionType="Specific",
                                 PredictionModel=Xenv,
                                 accu_r=cor_MegaLMM_Spec,
                                 accu_z=R_to_Z(cor_MegaLMM_Spec))
      if(CVScn %in% c("CV1","CV2")){
        cor_MegaLMM_Spec_Env  =corbetw2mat(Y_Test,MegaLMM_Pred)
        df_MegaLMM_Spec_Env=data.frame(file=gsub(".RData","",file),
                                       TraitName=traitname,
                                       CVScn=CVScn,
                                       ID_Fold=ID_Fold,
                                       ID_Fold_Level2=ID_Fold_Level2,
                                       EnvName=colnames(results$Y_Test),
                                       NumObs_TSTSet=colSums(!is.na(results$Y_Test)),
                                       Method="MegaLMM",
                                       PriorType=ifelse(Xenv=="Zero","No","Yes"),
                                       PriorLevel=Xenv,
                                       PredictionType="Specific",
                                       PredictionModel="Env",
                                       accu_r=cor_MegaLMM_Spec_Env,
                                       accu_z=R_to_Z(cor_MegaLMM_Spec_Env))
        df=rbind(df_MegaLMM_Gen,df_MegaLMM_Spec,df_MegaLMM_Spec_Env)
      }else{
        df=rbind(df_MegaLMM_Gen,df_MegaLMM_Spec)
      }
    }
    # print out to show computation progress
    print(paste("##########",traitname,"::",CVScn,":: fileID=", i, "##########"))
    # str(df)
    # return df
    return(df)  
    }
}
dim(res) #460660     14
head(res)

# ------------------------------------------------------------------------------
# exclude accu_r is NA
res_NAs = res[is.na(res$accu_r),]
dim(res_NAs) #21398    14
head(res_NAs)
table(res_NAs$CVScn)
# CV1          CV2  CV3NewState CV3NewTester  CV3NewTrial          CV4 
# 19           19         8110         8230           40         4980 
unique(res_NAs$EnvName[res_NAs$NumObs_TSTSet==1]) #"NEH4_2016::PHZ51"

# excluding NAs
res = res[!is.na(res$accu_r),]
dim(res) # 439262     14

# ------------------------------------------------------------------------------
# shorten PriorLevel/PredictionModel names
res0=res

# shorten PriorLevel names
res$PriorLevel[res$PriorLevel=="State"]="S"
res$PriorLevel[res$PriorLevel=="State+Tester"]="S+T"
res$PriorLevel[res$PriorLevel=="Tester"]="T"
res$PriorLevel[res$PriorLevel=="TesterGRM"]="K"
res$PriorLevel[res$PriorLevel=="Wea"]="W"
res$PriorLevel[res$PriorLevel=="Zero"]="O"
res$PriorLevel[res$PriorLevel=="Wea+Tester"]="W+T"
res$PriorLevel[res$PriorLevel=="State+TesterGRM"]="S+K"
res$PriorLevel[res$PriorLevel=="Wea+TesterGRM"]="W+K"
unique(res$PriorLevel)
# "S"   "S+T" "S+K" "T"   "K"   "W"   "W+T" "W+K" "O" 

# shorten PredictionModel names
res$PredictionModel[res$PredictionModel=="State"]="S"
res$PredictionModel[res$PredictionModel=="State+Tester"]="S+T"
res$PredictionModel[res$PredictionModel=="Tester"]="T"
res$PredictionModel[res$PredictionModel=="TesterGRM"]="K"
res$PredictionModel[res$PredictionModel=="Wea"]="W"
res$PredictionModel[res$PredictionModel=="Zero"]="O"
res$PredictionModel[res$PredictionModel=="Wea+Tester"]="W+T"
res$PredictionModel[res$PredictionModel=="State+TesterGRM"]="S+K"
res$PredictionModel[res$PredictionModel=="Wea+TesterGRM"]="W+K"
res$PredictionModel[res$PredictionModel=="General"]="O"
unique(res$PredictionModel)
# "O"   "S"   "Env" "S+T" "S+K" "T"   "K"   "W"   "W+T" "W+K"

# modify TraitName
unique(res$TraitName)
res$TraitName[res$TraitName=="SilkDAPdays"]="Silk Days"
res$TraitName[res$TraitName=="PlantHeightcm"]="Plant Height"
res$TraitName[res$TraitName=="YieldMgha"]="Grain Yield"

# ------------------------------------------------------------------------------
# remove duplicated PriorLevel::PredictionModel combinations for some CVScns
for(CVScn in unique(res$CVScn)){
  print(unique(paste(CVScn,res$PriorLevel[res$CVScn==CVScn],
                     res$PredictionModel[res$CVScn==CVScn],sep = "::")))
}
unique(res$CVScn)
# "CV1" "CV2" "CV3NewTester" "CV3NewTrial"  "CV4" "CV3NewState" 
sapply(list(CV1=res[res$CVScn=="CV1",],CV2=res[res$CVScn=="CV2",],
            CV3NewTrial=res[res$CVScn=="CV3NewTrial",], 
            CV3NewState=res[res$CVScn=="CV3NewState",],
            CV3NewTester=res[res$CVScn=="CV3NewTester",],
            CV4=res[res$CVScn=="CV4",]), 
       function(x) sort(unique(paste(x$Method,x$PriorLevel,x$PredictionModel,sep = "::"))) )


# remove duplicated/unwanted PriorLevel::PredictionModel combinations
res$accu_r[(res$CVScn=="CV3NewState" & paste(res$PriorLevel,res$PredictionModel,sep = "::")=="O::S+T")|
             (res$CVScn=="CV3NewTester" & paste(res$PriorLevel,res$PredictionModel,sep = "::")=="O::S+T")|
             (res$CVScn=="CV4" & paste(res$PriorLevel,res$PredictionModel,sep = "::") %in% c("O::S+T","O::T"))
           ]=NA
res=res[!is.na(res$accu_r),]
sum(is.na(res0$accu_r))
sapply(list(res0,res), dim)
#        [,1]   [,2]
# [1,] 439262 411762
# [2,]     14     14
# checking again 
sapply(list(CV1=res[res$CVScn=="CV1",],CV2=res[res$CVScn=="CV2",],
            CV3NewTrial=res[res$CVScn=="CV3NewTrial",], 
            CV3NewState=res[res$CVScn=="CV3NewState",],
            CV3NewTester=res[res$CVScn=="CV3NewTester",],
            CV4=res[res$CVScn=="CV4",]), 
       function(x) sort(unique(paste(x$Method,x$PriorLevel,x$PredictionModel,sep = "::"))) )

# ------------------------------------------------------------------------------
# Average Prediction Results for each Env/Trial across Five CV Folds
# ------------------------------------------------------------------------------
head(res)
factlevels(res)

table(res$CVScn)
#   CV1          CV2  CV3NewState CV3NewTester  CV3NewTrial          CV4 
# 77026        77026        72300        64360        72720        48330 
# ------------------------------------------------------------------------------
# examining whether each Env has been predicted with 5 CV2 reps for CV3/CV4
# holding out ID_Fold (CV2 Fold), I expect 
#  - each Env for CV3/4 has 5 accu values
#  - each Env for CV1/2 has 1 accu value

tmp1=aggregate(accu_r ~ TraitName + CVScn + ID_Fold_Level2 + 
                 EnvName + Method + PriorLevel + PredictionModel,
               data = res,
               FUN = length)
sapply(list(CV1=tmp1$accu_r[tmp1$CVScn=="CV1"],
            CV2=tmp1$accu_r[tmp1$CVScn=="CV2"],
            CV3NewTrial=tmp1$accu_r[tmp1$CVScn=="CV3NewTrial"],
            CV3NewState=tmp1$accu_r[tmp1$CVScn=="CV3NewState"],
            CV3NewTester=tmp1$accu_r[tmp1$CVScn=="CV3NewTester"],
            CV4=tmp1$accu_r[tmp1$CVScn=="CV4"]), 
       function(x) data.frame(table(x)))
#      CV1   CV2   CV3NewTrial CV3NewState CV3NewTester CV4 
# x    1     1     5           5           5            5   
# Freq 77026 77026 14544       14460       12872        9666

# ------------------------------------------------------------------------------
# holding out the ID_Fold_Level2, I expect:
#  - each Env for CV1/2 has 5 accu values
#  - each Env for CV3/CV4 has 1 accu value

tmp2=aggregate(accu_r ~ TraitName + CVScn + ID_Fold + 
                 EnvName + Method + PriorLevel + PredictionModel,
               data = res,
               FUN = length)
data.frame(lapply(list(CV1=tmp2$accu_r[tmp2$CVScn=="CV1"],
                       CV2=tmp2$accu_r[tmp2$CVScn=="CV2"],
                       CV3NewTrial=tmp2$accu_r[tmp2$CVScn=="CV3NewTrial"],
                       CV3NewState=tmp2$accu_r[tmp2$CVScn=="CV3NewState"],
                       CV3NewTester=tmp2$accu_r[tmp2$CVScn=="CV3NewTester"],
                       CV4=tmp2$accu_r[tmp2$CVScn=="CV4"]), 
                  function(x) t(data.frame(table(x))) ))
#      CV1.1 CV1.2 CV2.1 CV2.2 CV3NewTrial CV3NewState CV3NewTester   CV4
# x        4     5     4     5           1           1            1     1
# Freq    19 15390    19 15390       72300        64360       48370 48330

# ------------------------------------------------------------------------------
# For each Env, average prediction accuracy across Folds
#  - CV1/CV2: For each Env, average accu_r/accu_z across five Folds of ID_Fold_Level2
#  - CV3/CV4: For each Env in each TestSet, accu_r/accu_z across five Folds of ID_Fold 

head(res)
res_CV1CV2 = res[res$CVScn %in% c("CV1","CV2"),]
res_CV3CV4 = res[!(res$CVScn %in% c("CV1","CV2")),]
sapply(list(res_CV1CV2,res_CV3CV4), dim)
#        [,1]   [,2]
# [1,] 154052 257710
# [2,]     14     14

# average CV1/CV2
# note: for CV1/CV2, NAs values generated by masking are different among the five folds
#       so NumObs_TSTSet is also different between different testset fold for the same Env
res_CV1CV2_EnvMean =  
  aggregate(cbind(NumObs_TSTSet,accu_r,accu_z) ~ TraitName + CVScn + ID_Fold + 
              EnvName + Method + PriorType + PriorLevel + PredictionType + PredictionModel,
            data = res_CV1CV2,FUN = mean, na.rm=T)
dim(res_CV1CV2_EnvMean) #30818    12
811*19*2==dim(res_CV1CV2_EnvMean)[1]
head(res_CV1CV2_EnvMean)

# average CV3/CV4
# note: for a specific ID_Fold_Level2, it was based on five-folds of CV2 TRN
#       for CV2, NAs values generated by masking are different among the five folds
#       so NumObs_TSTSet is also different for the same Env across five folds of CV2 TRN
res_CV3CV4_EnvMean=
  aggregate(cbind(NumObs_TSTSet,accu_r,accu_z) ~ TraitName + CVScn + 
              ID_Fold_Level2 + EnvName + Method + PriorType +PriorLevel + 
              PredictionType + PredictionModel, 
            data = res_CV3CV4,FUN = mean, na.rm=T)
dim(res_CV3CV4_EnvMean) #51542    12
head(res_CV3CV4_EnvMean)

# merge res_CV1CV2_EnvMean and res_CV3CV4_EnvMean
colnames(res_CV1CV2_EnvMean)[colnames(res_CV1CV2_EnvMean)=="ID_Fold"]="TestSetID"
colnames(res_CV3CV4_EnvMean)[colnames(res_CV3CV4_EnvMean)=="ID_Fold_Level2"]="TestSetID"
if(identical(colnames(res_CV1CV2_EnvMean),colnames(res_CV3CV4_EnvMean))){
  res_EnvMean = rbind(res_CV1CV2_EnvMean, res_CV3CV4_EnvMean)
}

# verify res_EnvMean
dim(res_EnvMean) #82360    12
head(res_EnvMean)
table(res_EnvMean$CVScn)
#   CV1          CV2  CV3NewState CV3NewTester  CV3NewTrial          CV4 
# 15409        15409        14460        12872        14544         9666 

# number of PriorLevel::PredictionModel
sapply(list(CV1=res_EnvMean[res_EnvMean$CVScn=="CV1",],
            CV2=res_EnvMean[res_EnvMean$CVScn=="CV2",],
            CV3NewState=res_EnvMean[res_EnvMean$CVScn=="CV3NewState",],
            CV3NewTester=res_EnvMean[res_EnvMean$CVScn=="CV3NewTester",],
            CV3NewTrial=res_EnvMean[res_EnvMean$CVScn=="CV3NewTrial",],
            CV4=res_EnvMean[res_EnvMean$CVScn=="CV4",]), 
       function(x) length(unique(paste(x$Method, x$PriorLevel, x$PredictionModel,sep = "::") ))
)
# CV1          CV2  CV3NewState CV3NewTester  CV3NewTrial          CV4 
#  19           19           18           16           18           12 
# number of Envs in each CV Scenario
traitnames=unique(res$TraitName) #"Grain Yield"  "Silk Days"    "Plant Height"
numEnvs=foreach(i = 1:3,.combine = "rbind") %do% {
  sapply(list(CV1=res$EnvName[res$CVScn=="CV1"&res$TraitName==traitnames[i]],
              CV2=res$EnvName[res$CVScn=="CV2"&res$TraitName==traitnames[i]],
              CV3NewState=res$EnvName[res$CVScn=="CV3NewState"&res$TraitName==traitnames[i]],
              CV3NewTester=res$EnvName[res$CVScn=="CV3NewTester"&res$TraitName==traitnames[i]],
              CV3NewTrial=res$EnvName[res$CVScn=="CV3NewTrial"&res$TraitName==traitnames[i]],
              CV4=res$EnvName[res$CVScn=="CV4"&res$TraitName==traitnames[i]]), 
         function(x) length(unique(x)))
}
rownames(numEnvs)=traitnames
print(numEnvs)
#              CV1 CV2 CV3NewState CV3NewTester CV3NewTrial CV4
# Grain Yield  302 302         302          302         302 302
# Silk Days    231 231         231          231         231 231
# Plant Height 278 278         278          278         278 278

# output results
res_EnvMean = res_EnvMean
res_Original= res
sapply(list(res_Original,res_EnvMean), dim)
#        [,1]  [,2]
# [1,] 411762 82360
# [2,]     14    12
save(res_Original,res_EnvMean,
     file = sprintf("%s%sG2F_Prediction_Results_Formatted.RData",DIR_Output,Date) )

# ------------------------------------------------------------------------------
# ggplot
# ------------------------------------------------------------------------------
# DIR_Output="../09_Prediction_Results_Summary/"
# ls=load(file=sprintf("%s%sG2F_Prediction_Results_Formatted.RData",DIR_Output,"2023-10-23"))
# ls #"res_Original" "res_EnvMean" 

resPlot=NULL
resPlot = res_EnvMean
head(resPlot)

resPlot = add_column(.data = resPlot,
                     .after = "PredictionModel",
                     Meth_PriorLvl_PredMod=paste(resPlot$Method,
                                                 resPlot$PriorLevel,
                                                 resPlot$PredictionModel,
                                                 sep = "::"))
# (Meth_PriorLvl_PredMod_Order=rev(sort(unique(resPlot$Meth_PriorLvl_PredMod))))
# (Meth_PriorLvl_PredMod_Order=sort(unique(resPlot$Meth_PriorLvl_PredMod)))
sort(unique(resPlot$Meth_PriorLvl_PredMod))
Meth_PriorLvl_PredMod_Order=
  c("GBLUP::O::O","GBLUP::O::S","GBLUP::O::S+T","GBLUP::O::T","GBLUP::O::Env",
    "MegaLMM::O::O", "MegaLMM::O::S", "MegaLMM::O::S+T", "MegaLMM::O::T", "MegaLMM::O::Env",
    "MegaLMM::S::O","MegaLMM::S::S","MegaLMM::S::Env",
    "MegaLMM::W::O","MegaLMM::W::W",
    "MegaLMM::T::O","MegaLMM::T::T","MegaLMM::T::Env",
    "MegaLMM::K::K", "MegaLMM::K::O",
    "MegaLMM::S+T::O","MegaLMM::S+T::S+T","MegaLMM::S+T::Env",
    "MegaLMM::S+K::O","MegaLMM::S+K::S+K",
    "MegaLMM::W+T::O",   "MegaLMM::W+T::W+T",
    "MegaLMM::W+K::O",   "MegaLMM::W+K::W+K")
stopifnot(identical(sort(Meth_PriorLvl_PredMod_Order),sort(unique(resPlot$Meth_PriorLvl_PredMod))))
# geom_text
accu_basicstat=NULL
accu_basicstat = aggregate(accu_r ~ CVScn + TraitName + Method +
                             PriorLevel + PredictionModel,
                           data = resPlot, FUN = summary)
accu_basicstat = data.frame(accu_basicstat[,1:5],data.frame(accu_basicstat$accu_r))
colnames(accu_basicstat)

df_geom_text_accu_stat = 
  add_column(.data = accu_basicstat,
             CV_Meth_PriorLvl_PredMod=paste(accu_basicstat$CVScn,
                                            accu_basicstat$Method, 
                                            accu_basicstat$PriorLevel,
                                            accu_basicstat$PredictionModel,
                                            sep = "::"),
             Meth_PriorLvl_PredMod=paste(accu_basicstat$Method, 
                                         accu_basicstat$PriorLevel,
                                         accu_basicstat$PredictionModel,
                                         sep = "::"))

head(df_geom_text_accu_stat)

# set up outliers and plotting boundaries
max(resPlot$accu_r) # 0.9203672
accu_uplimit = 1; accu_lowlimit=-0.1
resPlot$accu_r[resPlot$accu_r<accu_lowlimit]=accu_lowlimit
resPlot_outliers = resPlot[resPlot$accu_r==accu_lowlimit,]
resPlot$accu_r[resPlot$accu_r==accu_lowlimit]=NA
head(resPlot_outliers)
range(resPlot$accu_r,na.rm = T)

# factorize
resPlot$Meth_PriorLvl_PredMod=factor(resPlot$Meth_PriorLvl_PredMod,
                                     levels = Meth_PriorLvl_PredMod_Order)
str(resPlot_outliers)

# output plot in pdf
head(resPlot)
i=3
for(i in 1:length(unique(resPlot$CVScn))){
  CVScn = unique(resPlot$CVScn)[i]
  (figure_caption=sprintf('%s%s_G2F_MegaLMM_General_vs_Specific_vs_EnvSpec_horizontal_%s',
                          DIR_Output,Date,CVScn))
  pdf(paste0(figure_caption,".pdf"),width = 12,height = 5) #10.98,5.74
  head(resPlot)
  
  # extrat plot data for CV_i
  resPlot_CV_i = resPlot[resPlot$CVScn==CVScn,]
  resPlot_outliers_CV_i=resPlot_outliers[resPlot_outliers$CVScn==CVScn,]
  
  # add empty Meth_PriorLvl_PredMod for each CVScn to make every CVScn has the 
  # same set of Meth_PriorLvl_PredMod 
  if(!(all(Meth_PriorLvl_PredMod_Order %in% resPlot_CV_i$Meth_PriorLvl_PredMod))){
    # resPlot_CV_i
    head(resPlot_CV_i)
    Meth_PriorLvl_PredMod=
    Meth_PriorLvl_PredMod_Order[!(Meth_PriorLvl_PredMod_Order %in% resPlot_CV_i$Meth_PriorLvl_PredMod)]
    df_GrainYield=data.frame(TraitName="Grain Yield",
               CVScn=CVScn,
               TestSetID=NA,
               EnvName=NA,
               Method=NA,
               PriorType=NA,
               PriorLevel=NA,
               PredictionType=NA,
               PredictionModel=NA, 
               Meth_PriorLvl_PredMod=Meth_PriorLvl_PredMod,
               NumObs_TSTSet =NA,
               accu_r=NA,
               accu_z=NA
               )
    df_SilkDays=df_GrainYield; df_SilkDays$TraitName=="Silk Days"
    df_PlantHeight=df_GrainYield; df_PlantHeight$TraitName=="Plant Height"
    resPlot_CV_i=rbind(resPlot_CV_i,
                       df_GrainYield,
                       df_SilkDays,
                       df_PlantHeight
                       )
  }
  unique(unique(resPlot_CV_i$Meth_PriorLvl_PredMod))
  
  # ggplot
  plot=ggplot(data = resPlot_CV_i, 
    # data = resPlot[resPlot$CVScn==unique(resPlot$CVScn)[i],], 
              aes(x = Meth_PriorLvl_PredMod,
                  y = accu_r,
                  fill=Meth_PriorLvl_PredMod)
  ) +
    labs(y = expression(paste("cor(y, ", hat(g), ")")))+
    facet_wrap(~TraitName)+
    geom_boxplot(outlier.size=0.5,
                 lwd = 0.4,#to make the lines thinner
                 fatten = 1 #make median line thinner
                 )+
    geom_text(                            ###plot means on top
      data = df_geom_text_accu_stat[df_geom_text_accu_stat$CVScn==unique(resPlot$CVScn)[i],],
      aes(y = max(resPlot$accu_r[resPlot$CVScn==unique(resPlot$CVScn)[i]], na.rm = T)-max(Mean)+Mean,
          # y=0.4+Mean,
          label = sprintf("%.3f", Mean)
          # color=factor(TraitName,levels = TraitSel )
      ),
      angle = 0,
      hjust =0.45
      # inherit.aes = FALSE,
      # position = position_dodge(width = 0.9),  # Apply dodge here
      # # vjust = -0.5  # Adjust vjust as needed
    ) +
    # guides(fill = guide_legend(nrow = 1))+
    theme_bw()+
    theme(#axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size=12),
          axis.text.x = element_text(angle=0),#vjust = 0.5,hjust=1),
          # axis.text.y = element_text(angle=90),
          strip.text.x = element_text(size = 12),
          # panel.grid = element_blank(),
          # panel.grid.major = element_blank(),
          legend.position="none",
          # legend.title = element_blank(),
          # legend.direction ="horizontal",
          # # legend.position = c(legend.position.x,legend.position.y),
          # legend.position="bottom",
          # # legend.justification = c("right","top")
          # legend.justification = "right",
          # legend.margin = margin(t = 0, unit = "cm")  # Adjust the top margin here
    )+
    geom_beeswarm(data = resPlot_outliers_CV_i, 
                  aes(x=Meth_PriorLvl_PredMod,y = accu_r),
                  method = "swarm",#must be one of: swarm, compactswarm, hex, square, center, or centre
                  corral="random",#"gutter","none","wrap","random", omit"
                  # # "random" places runaway points randomly in the region
                  # # "gutter" collects runaway points along the boundary between groups
                  corral.width=0.6,
                  color="red",show.legend = FALSE,
                  alpha = 0.5, size = 0.5)+
    coord_flip()+
    # geom_vline(xintercept = 20.5,lty=2,lwd=0.4)+
    geom_vline(xintercept = 5.5,lty=2,lwd=0.36)+
    ggtitle(unique(resPlot$CVScn)[i])
  print(plot)
  dev.off()
}



