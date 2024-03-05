# =================================================================
# script07_G2F_MegaLMM
# =================================================================

# Clean workspace
rm(list=ls())

library(rrBLUP)
library(MegaLMM); K=NULL
library(data.table)
library(lineup)
library(cvTools)
library(foreach)
library(doParallel)
library(rrBLUP)
library(dplyr)
# library(data.tree) #list str

#-------------------------------------------------------------------------------
# SECTION 01: set up DIR_Output
#-------------------------------------------------------------------------------
Date="2023-10-22"
DIR_Pheno_TRN_TST="../06_MegaLMM_PhenoDat_TRN_TST/"
DIR_Geno="../../P5G2F/02_GRMs/DER/"
DIR_EnvDat="../../P5G2F/04_Env_Dat/2022-12-26/"
DIR_Pheno_FMT2="../05_MegaLMM_PhenoDat_Format2/"
DIR_Output="../07c_G2F_MegaLMM/"
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)

#-------------------------------------------------------------------------------
# SECTION 02: set up global parameters
#-------------------------------------------------------------------------------
traitname = commandArgs(t=T)[1]
X_Env_Model = commandArgs(t=T)[2]
CVScn = commandArgs(t=T)[3]
ID_Fold = as.numeric(commandArgs(t=T)[4])
ID_Fold_Level2= as.numeric(commandArgs(t=T)[5])
if(is.na(traitname)) traitname="Plant_Height_cm"#
if(is.na(X_Env_Model)) X_Env_Model="Wea+TesterGRM" #"Zero","State","Tester" "State:Tester" "State+Tester"
if(is.na(CVScn)) CVScn="CV3NewState" #CVScn="CV1"#CVScn="CV3NewState"; CVScn="CV3NewTrial"
if(is.na(ID_Fold)) ID_Fold=1  #ID_fold for Y_Train_CV2, is 1-5 for CV3/4; 0 for CV1/2
if(is.na(ID_Fold_Level2)) ID_Fold_Level2=20 #actual Fold ID for CVScn under running

# set up MegaLMM running parameters 
centering_trn_set=T
scale_Y = F # Step1_intcpt=F; XF_intcpt=T

# get run_ID
run_ID = sprintf('G2F_%s_Xenv_%s_CV2_Fold%02d_%s_Fold%02d',
                 gsub("_","",traitname),X_Env_Model,ID_Fold,CVScn,ID_Fold_Level2)
if(dir.exists(run_ID)) unlink(run_ID,recursive = T)
print(paste0("####################",run_ID,"####################"))

#-------------------------------------------------------------------------------
# SECTION 03: Load data and set up MegaLMM Input
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# load Pheno_TRN_TST
list.files(DIR_Pheno_TRN_TST)
ls=load(sprintf("%s%sMegaLMM_TRN_TST.RData",DIR_Pheno_TRN_TST,"2023-10-20"))
ls
# [1] "Y_Train_CV1"          "Y_Test_CV1"           "Y_Train_CV2"          "Y_Test_CV2"          
# [5] "Y_Train_CV3NewTrial"   "Y_Test_CV3NewTrial"    "Y_Train_CV3NewState"  "Y_Test_CV3NewState"  
# [9] "Y_Train_CV3NewTester" "Y_Test_CV3NewTester"  "Y_Train_CV4"          "Y_Test_CV4"   
# sapply(Y_Train_CV2,function(x) sapply(x,dim ))
# sapply(Y_Train_CV3NewState,function(x) sapply(x,dim ))

# load Pheno_FMT2
list.files(DIR_Pheno_FMT2)
ls2=load(sprintf("%s%sMegaLMM_PhenoDat_Format2.RData",DIR_Pheno_FMT2,"2023-09-15"))
ls2 #"Phe_Trn_Fmt1Wide" "Phe_Trn_Fmt2Wide" "Phe_Trn_Fmt2Long" "K_info" 
sapply(Phe_Trn_Fmt2Wide, dim)

# load K matrix for P1 (created by DER)
K0 = readRDS(paste0(DIR_Geno,'A_P1.rds'))

# load K matrix for P2 (created by DER)
K_Tester0 = readRDS(paste0(DIR_Geno,'A_P2_new.rds'))

# load env data
list.files(DIR_EnvDat)
ls2=load(paste0(DIR_EnvDat,"G2F_EnvDat_QC.RData"))
ls2 #"Wea_Dat"     "Soil_Dat"    "ECs_Dat"     "metadat_Wea"

#-------------------------------------------------------------------------------
# get Pheno matrix Y
if(CVScn=="CV1"){
  Y_Train = Y_Train_CV1[[traitname]][[paste0("Fold_",ID_Fold_Level2)]]
  Y_Test = Y_Test_CV1[[traitname]][[paste0("Fold_",ID_Fold_Level2)]]
}else if(CVScn=="CV2"){
  Y_Train = Y_Train_CV2[[traitname]][[paste0("Fold_",ID_Fold_Level2)]]
  Y_Test = Y_Test_CV2[[traitname]][[paste0("Fold_",ID_Fold_Level2)]]
}else if(CVScn=="CV3NewTrial"){
  Y_Train = Y_Train_CV3NewTrial[[traitname]][[paste0("CV2_FoldID_",ID_Fold)]][[ID_Fold_Level2]]
  Y_Test = Y_Test_CV3NewTrial[[traitname]][[paste0("CV2_FoldID_",ID_Fold)]][[ID_Fold_Level2]]
}else if(CVScn=="CV3NewState"){
  Y_Train = Y_Train_CV3NewState[[traitname]][[paste0("CV2_FoldID_",ID_Fold)]][[ID_Fold_Level2]]
  Y_Test = Y_Test_CV3NewState[[traitname]][[paste0("CV2_FoldID_",ID_Fold)]][[ID_Fold_Level2]]
  stopifnot(length(Y_Train_CV3NewState[[traitname]][[paste0("CV2_FoldID_",ID_Fold)]])==20)
  stopifnot(length(Y_Test_CV3NewState[[traitname]][[paste0("CV2_FoldID_",ID_Fold)]])==20)
}else if(CVScn=="CV3NewTester"){
  Y_Train = Y_Train_CV3NewTester[[traitname]][[paste0("CV2_FoldID_",ID_Fold)]][[ID_Fold_Level2]]
  Y_Test = Y_Test_CV3NewTester[[traitname]][[paste0("CV2_FoldID_",ID_Fold)]][[ID_Fold_Level2]]
}else if(CVScn=="CV4"){
  Y_Train = Y_Train_CV4[[traitname]][[paste0("CV2_FoldID_",ID_Fold)]][[ID_Fold_Level2]]
  Y_Test = Y_Test_CV4[[traitname]][[paste0("CV2_FoldID_",ID_Fold)]][[ID_Fold_Level2]]
}
Y_Train0=Y_Train #keep an original version of Y_Train
sapply(list(Y_Train,Y_Test), dim)
# lapply(list(Y_Train,Y_Test), function(x) x[1:3,1:4])
# colnames(Y_Test)

# check whether the Y_Train and Y_Test are as expected
if(CVScn=="CV1"|CVScn=="CV2"){
  stopifnot(identical(rownames(Y_Train),rownames(Y_Test)))
  stopifnot(identical(colnames(Y_Train),colnames(Y_Test)))
}else if (CVScn=="CV4"){
  stopifnot(identical(intersect(rownames(Y_Train),rownames(Y_Test)),character(0)))
  stopifnot(identical(intersect(colnames(Y_Train),colnames(Y_Test)),character(0)))
}else{ #For CV3s
  stopifnot(identical(rownames(Y_Train),rownames(Y_Test)))
  stopifnot(identical(intersect(colnames(Y_Train),colnames(Y_Test)),character(0)))
}
sapply(list(Y_Train,Y_Test), dim)
# lapply(list(Y_Train,Y_Test), function(x) x[1:3,1:4])

# centering training set by global mean of training set
if(centering_trn_set){
  globalmean = mean(as.matrix(Y_Train),na.rm = T)
  Y_Train = Y_Train - globalmean
}
# check to ensure that the subtraction is only performed once
stopifnot(all.equal( (Y_Train+globalmean)[!is.na(Y_Train+globalmean)],
                     Y_Train0[!is.na(Y_Train0)] ) )

#-------------------------------------------------------------------------------
# set up Y which contains all lines in row and all EnvNew in column
# used for CV3/CV4 specifically 
if(CVScn=="CV1" | CVScn=="CV2"){
  Y=Y_Train
}else {
  Y = Y_Train_CV2[[traitname]][[1]]
  # for CV3/CV4, we only use rownames/colnames of Y, not use its data for computation
}
# get lines_train and lines_test (used for CV3/CV4)
lines_train = rownames(Y) %in% rownames(Y_Train)
lines_test  = rownames(Y) %in% rownames(Y_Test)
envs_train  = colnames(Y) %in% colnames(Y_Train)
envs_test   = colnames(Y) %in% colnames(Y_Test)

# get line_data and X_F_train that are going to be used for MegaLMM modeling
line_data = data.frame(Line = rownames(Y))
line_data_train = line_data[line_data$Line %in% rownames(Y_Train),,drop=F]
line_data_test  = line_data[line_data$Line %in% rownames(Y_Test),,drop=F]
sapply(list(line_data,line_data_train,line_data_test), dim)

# X_F_train will be used for extra_regressions in setup_model_MegaLMM
X_F_train = model.matrix(~1,line_data_train)
X_F_test = model.matrix(~1,line_data_test)
# head(line_data)
# head(X_F_train)
# head(X_F_test)

#-------------------------------------------------------------------------------
# Prepare K for Parent1 (note: K should contain both TRN and TST)
# subset K0 for those Hybrids existing in Pheno data
if(all(K_info$Hybrid %in% rownames(K0))){
  K=K0[K_info$Hybrid,K_info$Hybrid]
}
dim(K0) #4918 4918
dim(K)  #4057 4057
stopifnot(identical(rownames(K),K_info$Hybrid))
K[1:6,1:6]

# extract K for P1
if(identical(rownames(K),K_info$Hybrid) &
   all(rownames(Y) %in% K_info$Hybrid_Parent1)){
  i = match(rownames(Y),K_info$Hybrid_Parent1)
  # subset K_P1
  K_P1 = K[i,i] 
  # update row/colnames for K_P1
  rownames(K_P1) = colnames(K_P1) = K_info$Hybrid_Parent1[i] 
}
dim(K_P1) #1700 1700
K_P1[1:6,1:6]
# checking whether Parent1 in the same order between Pheno and Geno data
stopifnot(identical(rownames(Y),rownames(K_P1)))
K = K_P1 # to be convenient, I assign K_P1 to K, and use K for K_P1 hereafter

#-------------------------------------------------------------------------------
# Prepare K for Parent2
testers_test_lst=unique(sapply(strsplit(colnames(Y),"::"), "[",2))
K_Tester = K_Tester0[testers_test_lst,testers_test_lst]
dim(K_Tester0); dim(K_Tester)
# K_Tester0[1:4,1:4]
K_Tester

# ------------------------------------------------------------------------------
# Prepare Xenv model for predicting lambda
# take a look at the Env data
ls2 #"Wea_Dat"     "Soil_Dat"    "ECs_Dat"     "metadat_Wea"
sapply(list(Wea=Wea_Dat,Soil=Soil_Dat,ECs=ECs_Dat,metaWea=metadat_Wea), dim)
#      Wea Soil ECs metaWea
# [1,] 235  162 189     238
# [2,] 281   26 764      18
# lapply(list(Wea=Wea_Dat,Soil=Soil_Dat,ECs=ECs_Dat,metaWea=metadat_Wea), function(x)
#   x[1:4,1:6])
# Wea_Dat[1:4,1:6]
# unique(Wea_Dat$Year)

# prepare Z_Tester for using Tester GRM as predictor
Envs = unique(sapply(strsplit(colnames(Y),"::"), "[",1))
M_Env= matrix(1,nrow=length(Envs),ncol = length(Envs))
rownames(M_Env)=colnames(M_Env)=Envs
Testers = sapply(strsplit(colnames(Y),"::"), "[",2)
M_Tester = model.matrix(~0+Testers)
colnames(M_Tester) = sub('Testers','',colnames(M_Tester))
K_Tester = K_Tester[colnames(M_Tester),colnames(M_Tester)]
L_Tester = t(chol(K_Tester))
Z_Tester = M_Tester %*% L_Tester
rownames(Z_Tester)=colnames(Y)

# setup a df for categorical predictors
Y_info = data.frame(EnvNew=colnames(Y),
                    State=substr(colnames(Y), start = 1, stop = 2),
                    Tester=sapply(strsplit(colnames(Y),"::"), "[",2)
                    )
head(Y_info)
df=NULL; df=data.frame(EnvNew=Y_info$EnvNew,
              State=factor(Y_info$State),
              Tester=factor(Y_info$Tester),
              State_Tester = as.factor(paste(Y_info$State,Y_info$Tester,sep = ":")))

# set up Xenv design matrix and X_Env_Model_groups
if(X_Env_Model=="State"){
  X_Env = model.matrix(~State,
                       data = df,
                       contrasts.arg = list(State = contrasts(df$State, contrasts = FALSE))
  )
  rownames(X_Env) = Y_info$EnvNew
  X_Env_groups = c(1, rep(2,length(unique(Y_info$State)) ) )
}else if(X_Env_Model=="State:Tester"){
  X_Env = model.matrix(~State_Tester,
                       data = df,
                       contrasts.arg = list(State_Tester = contrasts(df$State_Tester, contrasts = FALSE))
  ) 
  rownames(X_Env) = Y_info$EnvNew
  X_Env_groups = c(1, rep(2,length(unique(df$State_Tester)) ) )
}else if(X_Env_Model=="Tester"){
  X_Env = model.matrix(~Tester,
                       data = df,
                       contrasts.arg = list(Tester = contrasts(df$Tester, contrasts = FALSE)))
  rownames(X_Env) = Y_info$EnvNew
  X_Env_groups = c(1,rep(2,length(unique(Y_info$Tester)) ) )
}else if(X_Env_Model=="State+Tester"){
  X_Env = model.matrix(~State+Tester,
                       data = df,
                       contrasts.arg = lapply(df[, sapply(df, is.factor), drop = FALSE],
                                              contrasts, contrasts = FALSE))
  rownames(X_Env) = Y_info$EnvNew
  X_Env_groups = c(1,
                   rep(2,length(unique(Y_info$State)) ),
                   rep(3,length(unique(Y_info$Tester)) ) )
}else if(X_Env_Model=="Wea"){ #Intcp+Wea
  # step 1: define env data
  env_data = Wea_Dat
  # step 2: harmonize env_data with Y
  Env_Y = data.frame(Env=sapply(strsplit(colnames(Y),"::"), "[",1), EnvNew=colnames(Y))
  env_data_Y = inner_join(Env_Y,env_data,by="Env")
  dim(Env_Y); dim(env_data_Y)
  Env_Y = Env_Y[Env_Y$EnvNew %in% env_data_Y$EnvNew,] #update Env_Y
  Y = Y[,colnames(Y) %in% env_data_Y$EnvNew]          #updated Y
  stopifnot(identical(env_data_Y$EnvNew,colnames(Y)))
  dim(Y)
  # step 3: Select and scale env covariates
  rownames(env_data_Y)=env_data_Y$EnvNew
  env_covariates = env_data_Y[,-c(1:4)]
  env_data = env_data_Y[,1:4]
  env_covariates = env_covariates[,apply(env_covariates,2,var)>1e-10]
  env_covariates = scale(env_covariates)
  # step 4a: prepare X_Env intercept
  X_Env_1 = model.matrix(~1,env_data)
  # step 4b: prepare X_Env_wea
  s_env_covariates = svd(env_covariates)
  X_Env_wea = t(sqrt(s_env_covariates$d) * t(s_env_covariates$u))
  rownames(X_Env_wea)=env_data$EnvNew
  # X_Env_wea[1:4,1:4]
  # step 4c: prepare X_Env and X_Env_groups
  X_Env = cbind(X_Env_1, X_Env_wea)
  X_Env_groups = c(rep(1,ncol(X_Env_1)),rep(2,ncol(X_Env_wea)) )
  # length(X_Env_groups) == ncol(X_Env) #186
  # dim(X_Env)# 185 186
  # X_Env[1:4,1:6]
}else if(X_Env_Model=="Soil"){ #Intcp+Soil
  # step 1: define env data
  env_data = Soil_Dat
  # step 2: harmonize env_data with Y
  Env_Y = data.frame(Env=sapply(strsplit(colnames(Y),"::"), "[",1), EnvNew=colnames(Y))
  env_data_Y = inner_join(Env_Y,env_data,by="Env")
  dim(Env_Y); dim(env_data_Y)
  Env_Y = Env_Y[Env_Y$EnvNew %in% env_data_Y$EnvNew,] #update Env_Y
  Y = Y[,colnames(Y) %in% env_data_Y$EnvNew]          #updated Y
  stopifnot(identical(env_data_Y$EnvNew,colnames(Y)))
  # step 3: Select and scale env covariates
  rownames(env_data_Y)=env_data_Y$EnvNew
  env_covariates = env_data_Y[,-c(1:4)]
  env_data = env_data_Y[,1:4]
  env_covariates = env_covariates[,apply(env_covariates,2,var)>1e-10]
  env_covariates = scale(env_covariates)
  # step 4a: prepare X_Env intercept
  X_Env_1 = model.matrix(~1,env_data)
  # step 4b: prepare X_Env_wea
  s_env_covariates = svd(env_covariates)
  X_Env_wea = t(sqrt(s_env_covariates$d) * t(s_env_covariates$u))
  rownames(X_Env_wea)=env_data$EnvNew
  X_Env_wea[1:4,1:4]
  # step 4c: prepare X_Env and X_Env_groups
  X_Env = cbind(X_Env_1, X_Env_wea)
  X_Env_groups = c(rep(1,ncol(X_Env_1)),rep(2,ncol(X_Env_wea)) )
}else if(X_Env_Model=="ECs"){ #Intcp+ECs
  # step 1: define env data
  env_data = ECs_Dat
  # step 2: harmonize env_data with Y
  Env_Y = data.frame(Env=sapply(strsplit(colnames(Y),"::"), "[",1), EnvNew=colnames(Y))
  env_data_Y = inner_join(Env_Y,env_data,by="Env")
  dim(Env_Y); dim(env_data_Y)
  Env_Y = Env_Y[Env_Y$EnvNew %in% env_data_Y$EnvNew,] #update Env_Y
  Y = Y[,colnames(Y) %in% env_data_Y$EnvNew]          #updated Y
  stopifnot(identical(env_data_Y$EnvNew,colnames(Y)))
  # step 3: Select and scale env covariates
  rownames(env_data_Y)=env_data_Y$EnvNew
  env_covariates = env_data_Y[,-c(1:4)]
  env_data = env_data_Y[,1:4]
  env_covariates = env_covariates[,apply(env_covariates,2,var)>1e-10]
  env_covariates = scale(env_covariates)
  # step 4a: prepare X_Env intercept
  X_Env_1 = model.matrix(~1,env_data)
  # step 4b: prepare X_Env_wea
  s_env_covariates = svd(env_covariates)
  X_Env_wea = t(sqrt(s_env_covariates$d) * t(s_env_covariates$u))
  rownames(X_Env_wea)=env_data$EnvNew
  X_Env_wea[1:4,1:4]
  # step 4c: prepare X_Env and X_Env_groups
  X_Env = cbind(X_Env_1, X_Env_wea)
  X_Env_groups = c(rep(1,ncol(X_Env_1)),rep(2,ncol(X_Env_wea)) )
  # dim(X_Env) #144 145
}else if(X_Env_Model=="TesterGRM"){ #Intcpt + TesterGRM
  # Prepare X_Env intercept
  X_Env_1 = model.matrix(~1,Y_info)
  rownames(X_Env_1)=Y_info$EnvNew
  # prepare X_Env_Tester
  X_Env_Tester = Z_Tester
  # prepare X_Env
  X_Env = cbind(X_Env_1, X_Env_Tester)
  X_Env_groups = c(rep(1,ncol(X_Env_1)),rep(2,ncol(X_Env_Tester)))
}else if(X_Env_Model=="Wea+Tester"){
  # step 1: define env data
  env_data = Wea_Dat
  env_data[1:4,1:6]
  sum(is.na(env_data))
  
  # step 2: harmonize env_data with Y
  Env_Y = data.frame(Env=sapply(strsplit(colnames(Y),"::"), "[",1), EnvNew=colnames(Y))
  env_data_Y = inner_join(Env_Y,env_data,by="Env")
  Env_Y = Env_Y[Env_Y$EnvNew %in% env_data_Y$EnvNew,] #update Env_Y
  Y = Y[,colnames(Y) %in% env_data_Y$EnvNew]          #updated Y
  stopifnot(identical(env_data_Y$EnvNew,colnames(Y)))
  
  # step 3: Select and scale env covariates
  rownames(env_data_Y)=env_data_Y$EnvNew
  env_covariates = env_data_Y[,-c(1:4)]
  env_data = env_data_Y[,1:4]
  env_covariates = env_covariates[,apply(env_covariates,2,var)>1e-10]
  env_covariates = scale(env_covariates)
  
  # step 4a: prepare X_Env intercept
  X_Env_1 = model.matrix(~1,env_data)
  
  # step 4b: prepare X_Env_wea
  s_env_covariates = svd(env_covariates)
  X_Env_wea = t(sqrt(s_env_covariates$d) * t(s_env_covariates$u))
  rownames(X_Env_wea)=env_data$EnvNew
  X_Env_wea[1:4,1:4]
  
  # step 4c: prepare X_Env_tester based on updated Y
  Y_info = data.frame(EnvNew=colnames(Y),
                      State=substr(colnames(Y), start = 1, stop = 2),
                      Tester=sapply(strsplit(colnames(Y),"::"), "[",2)
  )
  df=data.frame(State=factor(Y_info$State),Tester=factor(Y_info$Tester))
  X_Env_tester = model.matrix(~0+Tester,
                              data = df,
                              contrasts.arg = list(Tester = contrasts(df$Tester, contrasts = FALSE)))
  sapply(list(X_Env_1, X_Env_wea, X_Env_tester), dim)

  # step 4d: prepare X_Env_groups
  X_Env = cbind(X_Env_1, X_Env_wea, X_Env_tester)
  X_Env_groups = c(rep(1,ncol(X_Env_1)),rep(2,ncol(X_Env_wea)),rep(3,ncol(X_Env_tester)))
  
}else if(X_Env_Model=="State+TesterGRM"){
  # Prepare X_Env intercept
  X_Env_1 = model.matrix(~1,Y_info)
  rownames(X_Env_1)=Y_info$EnvNew
  # prepare X_Env_state
  X_Env_state = model.matrix(~0+State,
                             data = df)
  rownames(X_Env_state)=df$EnvNew
  # Prepare X_Env_Tester
  X_Env_Tester = Z_Tester
  
  # prepare X_Env
  if(all(rownames(X_Env_state)==rownames(X_Env_1) & rownames(X_Env_Tester)==rownames(X_Env_state))){
    X_Env = cbind(X_Env_1, X_Env_state, X_Env_Tester)
    X_Env_groups = c(rep(1,ncol(X_Env_1)),rep(2,ncol(X_Env_state)),rep(3,ncol(X_Env_Tester)))
  }
}else if(X_Env_Model=="Wea+TesterGRM"){ #W+K
  # step 1: define env data
  env_data = Wea_Dat
  env_data[1:4,1:6]
  dim(env_data)
  sum(is.na(env_data))
  colnames(Y_Test)
  unique(substr(colnames(Y_Test),1,2)) %in% unique(substr(env_data$Env,1,2))
  
  # step 2: harmonize env_data with Y
  Env_Y = data.frame(Env=sapply(strsplit(colnames(Y),"::"), "[",1), EnvNew=colnames(Y))
  env_data_Y = inner_join(Env_Y,env_data,by="Env")
  Env_Y = Env_Y[Env_Y$EnvNew %in% env_data_Y$EnvNew,] #update Env_Y
  Y = Y[,colnames(Y) %in% env_data_Y$EnvNew]          #updated Y
  stopifnot(identical(env_data_Y$EnvNew,colnames(Y)))
  dim(Y)
  
  # step 3: Select and scale env covariates
  rownames(env_data_Y)=env_data_Y$EnvNew
  env_covariates = env_data_Y[,-c(1:4)]
  env_data = env_data_Y[,1:4]
  env_covariates = env_covariates[,apply(env_covariates,2,var)>1e-10]
  env_covariates = scale(env_covariates)
  
  # step 4a: prepare X_Env intercept
  X_Env_1 = model.matrix(~1,env_data)
  
  # step 4b: prepare X_Env_wea
  s_env_covariates = svd(env_covariates)
  X_Env_wea = t(sqrt(s_env_covariates$d) * t(s_env_covariates$u))
  rownames(X_Env_wea)=env_data$EnvNew
  X_Env_wea[1:4,1:4]
  
  # step 4c: Prepare X_Env_Tester based on updated Y
  Envs = unique(sapply(strsplit(colnames(Y),"::"), "[",1))
  M_Env= matrix(1,nrow=length(Envs),ncol = length(Envs))
  rownames(M_Env)=colnames(M_Env)=Envs
  Testers = sapply(strsplit(colnames(Y),"::"), "[",2)
  M_Tester = model.matrix(~0+Testers)
  colnames(M_Tester) = sub('Testers','',colnames(M_Tester))
  K_Tester = K_Tester[colnames(M_Tester),colnames(M_Tester)]
  L_Tester = t(chol(K_Tester))
  Z_Tester = M_Tester %*% L_Tester
  X_Env_Tester = Z_Tester
  # sapply(list(X_Env_1, X_Env_wea, X_Env_Tester,Y), dim)
  # X_Env_Tester[1:3,1:4]
  
  # step 4d: prepare X_Env_groups
  if(nrow(X_Env_1)==nrow(X_Env_wea) & nrow(X_Env_1)==nrow(X_Env_Tester)){
    X_Env = cbind(X_Env_1, X_Env_wea, X_Env_Tester)
    X_Env_groups = c(rep(1,ncol(X_Env_1)),rep(2,ncol(X_Env_wea)),rep(3,ncol(X_Env_Tester)))
  }
}

# checking X_Env and update Y_Train and Y_Test
if(X_Env_Model!="Zero"){
  stopifnot(length(X_Env_groups) == ncol(X_Env))
  # get X_Env_Train and X_Env_Test
  X_Env_Train = X_Env[rownames(X_Env) %in% colnames(Y_Train),]
  X_Env_Test  = X_Env[rownames(X_Env) %in% colnames(Y_Test),,drop=F]
  
  # Updating Y_Train and Y_Test for EnvNew since not all EnvNeW have Env variables
  Y_Train=Y_Train[,colnames(Y_Train) %in% colnames(Y)]
  Y_Test = Y_Test[, colnames(Y_Test) %in% colnames(Y),drop=F]
  
  # sort X_Env_Train and X_Env_Test 
  X_Env_Train=X_Env_Train[colnames(Y_Train),]
  X_Env_Test =X_Env_Test[colnames(Y_Test),,drop=F]
  
  # checking to make sure Y_Train/X_Env_Train and Y_Test/X_Env_Test identical for EnvNew
  stopifnot(identical(colnames(Y_Train),rownames(X_Env_Train)))
  stopifnot(identical(colnames(Y_Test),rownames(X_Env_Test)))
}

#-------------------------------------------------------------------------------
# checking before running MegaLMM
if(CVScn=="CV1"|CVScn=="CV2"){
  # check rows
  stopifnot(all(sapply(list(Y_Train=Y_Train, Y_Test=Y_Test, K=K_P1),function(x) 
    identical(rownames(x),line_data$Line)) ))
  # checking columns
  stopifnot(identical(colnames(Y_Train),colnames(Y_Test)))
  if(X_Env_Model!="Zero"){
    stopifnot(all(sapply(list(Y_Train, Y_Test), function(x) 
      identical(colnames(x),rownames(X_Env)))))
  }
}else if(CVScn=="CV4"){
  # checking rows
  stopifnot(identical(intersect(rownames(Y_Train), rownames(Y_Test)),character(0)))
  # checking columns
  stopifnot(identical(intersect(colnames(Y_Train), colnames(Y_Test)),character(0)))
  if(X_Env_Model!="Zero"){
    stopifnot(identical(colnames(Y_Train),rownames(X_Env_Train)))
    stopifnot(identical(colnames(Y_Test),rownames(X_Env_Test)))
  }
}else{ #CV3
  # checking rows
  stopifnot( all(sapply(list(Y_Train=Y_Train,Y_Test=Y_Test,K=K_P1),function(x) 
    identical(rownames(x),line_data$Line)) ) )
  # checking columns
  stopifnot(identical(intersect(colnames(Y_Train), colnames(Y_Test)),character(0)))
  if(X_Env_Model!="Zero"){
    stopifnot(identical(colnames(Y_Train),rownames(X_Env_Train)))
    stopifnot(identical(colnames(Y_Test),rownames(X_Env_Test)))
  }
}
dim(Y_Train) #1700  294
dim(Y_Test)  #1700  294
dim(K) # 1700 1700

stopifnot(!all(is.na(Y_Train))) #if all NAs in Y_Train, stop running
stopifnot(!all(is.na(Y_Test)))  #if all NAs in Y_Test,  stop running

#-------------------------------------------------------------------------------
# SECTION 04: run rrBLUP
#-------------------------------------------------------------------------------
# rrBLUP/GBLUP can only predict within MET, cannot predict TPE

# cl <- makeCluster(RcppParallel::defaultNumThreads()-1)
cl <- makeCluster(5)
registerDoParallel(cl)
i=1

rrBLUP = foreach(i = 1:ncol(Y_Train),.combine = cbind) %dopar% {
  print(i)
  library(rrBLUP)
  y_train = Y_Train[,i]
  j = !is.na(y_train)
  Kj = K[lines_train,lines_train][j,j]
  while(T) {
    offset = 0
    # r = try(mixed.solve(y_train[j],X = X_F_train[j,],K = Kj+diag(offset,nrow(Kj))))
    r = try(mixed.solve(y_train[j],K = Kj+diag(offset,nrow(Kj))))
    if(is.character(r)){
      offset = offset + 0.01
    } else if(offset > .1){
      break
    } else{
      break
    }
  }
  a2=X_F_test %*% r$beta + as.matrix(K[lines_test,rownames(Y_Train)[j]] %*% MASS::ginv(Kj)) %*% r$u
  # a2= as.matrix(K[lines_test,rownames(Y_Train)[j]] %*% MASS::ginv(Kj)) %*% r$u
  a2
}
rownames(rrBLUP)=rownames(Y_Test)
colnames(rrBLUP)=colnames(Y_Train)
stopCluster(cl)
# rrBLUP[1:3,1:4]

#-------------------------------------------------------------------------------
# SECTION 05: run MegaLMM
#-------------------------------------------------------------------------------

# add a small number in the diag of K to avoid SVD decompose problem
# offset=runif(n=nrow(K),min=0,max=0.00001)
# offset=runif(n=nrow(K),min=0,max=0.0001)
# K = K + diag(offset,nrow(K))

#-------------------------------------------------------------------------------
### Set control parameters
run_parameters = MegaLMM_control(
  h2_divisions = 20, 
  # Each variance component is allowed to explain between 0% and 100% of the
  # total variation. How many segments should the range [0,100) be divided 
  # into for each random effect?
  K_eigen_tol = 1e-1, #new argument added by DER on 01/23/23
  scale_Y = scale_Y,
  drop0_tol = 1e-10,
  burn = 0,  
  # number of burn in samples before saving posterior samples. I set this to 
  # zero and instead run the chain in small chunks, doing the burning manually, 
  # as described below.
  thin = 2,
  # during sampling, we'll save every 2nd sample to the posterior database.
  K = 50 # number of factors
)

#-------------------------------------------------------------------------------
### Create the model object
# The function `setup_model_MegaLMM` parses the model formulas, links the GRM to
# the random effects, and creates an object to store all components of the model.
# Y_Train[1:3,1:4]

MegaLMM_state = setup_model_MegaLMM(
  Y = Y_Train,  
  # The n x p trait matrix
  formula = ~ 0 + (1|Line),
  # formula = ~ 1 + (1|Line),  
  # This is syntax like lme4 for mixed effect models. 
  # We specify a fixed effect of population and a random effect for genotype (Line)
  data = line_data_train,         
  # the data.frame with information for constructing the model matrices
  extra_regressions = list(X=X_F_train,factors=T), #factors=T is required otherwise get error
  relmat = list(Line = K), 
  # A list of covariance matrices to link to the random effects in formula.
  # each grouping variable in formula can be linked to a covariance matrix.
  # If so, every level of the grouping variable must be in the rownames of K.
  # additional rows of K not present in data will still be predicted 
  # (and therefore will use memory and computational time!)
  run_parameters=run_parameters,
  # This list of control parameters created above
  run_ID = run_ID
  # A run identifier. The function will create a folder with this name 
  # and store lots of useful data inside it
)

#-------------------------------------------------------------------------------
# get K_train_inv from ZL for U_cond (newly added by DER on 01/23/23)
L_train = as.matrix(MegaLMM_state$data_matrices$ZL)
sK_train_sqrt = svd(L_train)
d = sK_train_sqrt$d>1e-10
K_train_inv = sK_train_sqrt$u[,d] %*% (1/sK_train_sqrt$d[d]^2 * t(sK_train_sqrt$u[,d]))
rownames(K_train_inv) = line_data_train$Line

#-------------------------------------------------------------------------------
### Set priors
# We need to set priors for the variance components ($\sigma^2$ and $h^2$ for `Y_R` and `F`,
# and for the parameters of the factor loadings `Lambda`.
if(X_Env_Model=="Zero"){
  Lambda_prior = list(
    sampler = sample_Lambda_prec_ARD,
    Lambda_df = 3,
    delta_1   = list(shape = 2,  rate = 1),
    delta_2   = list(shape = 3, rate = 1),  # parameters of the gamma distribution giving the expected change in proportion of non-zero loadings in each consecutive factor
    # X = X_Env_Train,
    # X_group = X_Env_groups,
    fit_X = F,  # we start by letting Lambda converge without X, but then turn it on during burnins.
    Lambda_beta_var_shape = 3,
    Lambda_beta_var_rate = 1,
    delta_iterations_factor = 100   # parameter that affects mixing of the MCMC sampler. This value is generally fine.
  )
}else{
  Lambda_prior = list(
    sampler = sample_Lambda_prec_ARD,
    Lambda_df = 3,
    delta_1   = list(shape = 2,  rate = 1),
    delta_2   = list(shape = 3, rate = 1),  # parameters of the gamma distribution giving the expected change in proportion of non-zero loadings in each consecutive factor
    X = X_Env_Train,
    X_group = X_Env_groups,
    fit_X = F,  # we start by letting Lambda converge without X, but then turn it on during burnins.
    Lambda_beta_var_shape = 3,
    Lambda_beta_var_rate = 1,
    delta_iterations_factor = 100   # parameter that affects mixing of the MCMC sampler. This value is generally fine.
  )
  
}

# For the remaining priors we use `MegaLMM_priors`
priors = MegaLMM_priors(
  tot_Y_var = list(V = 0.5,   nu = 5),      
  # Prior variance of trait residuals after accounting for fixed effects and factors
  # See MCMCglmm for meaning of V and nu
  tot_F_var = list(V = 18/20, nu = 20),     
  # Prior variance of factor traits. This is included to improve MCMC mixing, 
  # but can be turned off by setting nu very large
  h2_priors_resids_fun = function(h2s,n)  1,  
  # Function that returns the prior density for any value of the h2s vector 
  # (ie the vector of random effect proportional variances across all random effects. 
  # 1 means constant prior. 
  # n is the number of h2 divisions above (here=20)
  # 1-n*sum(h2s)/n linearly interpolates between 1 and 0, 
  # giving more weight to lower values
  h2_priors_factors_fun = function(h2s,n) 1, 
  # See above. 
  # sum(h2s) linearly interpolates between 0 and 1,
  # giving more weight to higher values
  # Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
  Lambda_prior = Lambda_prior
  # from above
)

# We then assign them to the `MegaLMM_state` object:
MegaLMM_state = set_priors_MegaLMM(MegaLMM_state,priors)

# Deal with the missing data
maps = make_Missing_data_map(MegaLMM_state,max_NA_groups = ncol(Y_Train)+1,verbose=F)

# maps$map_results
all_nas = lapply(maps$Missing_data_map_list,function(map) {
  nas=sapply(1:ncol(Y_Train),function(col) {
    map_id = which(sapply(map,function(x) col %in% x$Y_cols))
    c(sum(!is.na(Y_Train[,col])),sum(is.na(Y_Train[map[[map_id]]$Y_obs,col])))
  })
})
max_per_NAs = sapply(all_nas,function(nas) max(nas[2,]/colSums(nas)))
maps$map_results[which(max_per_NAs<0.1)[1],]
MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map_list[[which(max_per_NAs<0.1)[1]]])

#-------------------------------------------------------------------------------
### Initialize the model object
# Next, we create random starting values for all parameters:
MegaLMM_state$current_state = NULL
MegaLMM_state$current_state$delta = NULL
MegaLMM_state$priors$Lambda_prior$fit_X = F
MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
MegaLMM_state$run_parameters$burn = 0

# Estimate Memory usage
# estimate_memory_initialization_MegaLMM(MegaLMM_state)

# run the initial calculations
MegaLMM_state = initialize_MegaLMM(MegaLMM_state,verbose = T)

#-------------------------------------------------------------------------------
### Prepare the Posterior datasets
# specify which posterior parameters to save
MegaLMM_state$Posterior$posteriorSample_params = c('Lambda','F_h2','resid_h2',
                                                   'tot_Eta_prec','B1','U_F','F',
                                                   'B2_F','Lambda_beta','Lambda_beta_var')
MegaLMM_state$Posterior$posteriorMean_params = 'Eta_mean'

#-------------------------------------------------------------------------------
### prepare prediction matrices
U_cond = K[lines_test,lines_train] %*% K_train_inv %*% MegaLMM_state$data_matrices$Z
lines_pred = match(rownames(K)[lines_test],line_data_train$Line)
# it's ok to leave U_cond and lines_pred here, although not all CVs use them
if(X_Env_Model=="Zero"){
  MegaLMM_state$Posterior$posteriorFunctions = list(
    U_Train = '(U_F) %*% Lambda + U_R', # The predictions! 12/13/22 +X1%*%B1
    U_Pred  = '(U_cond %*% U_F) %*% Lambda + U_cond %*% U_R', # The predictions! 12/13/22 +X_F_test %*%B1
    # U_Test  = '(U_cond %*% U_F) %*% Lambda_beta %*% t(X_Env_Test)', # The predictions! updated 01/02/23
    # Eta_test = 'F[lines_pred,] %*% Lambda_beta %*% t(X_Env_Test)',
    h2 = '(colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,])/(colSums(Lambda^2)+1/tot_Eta_prec[1,])'
  )
  
}else{
  MegaLMM_state$Posterior$posteriorFunctions = list(
    U_Train = '(U_F) %*% Lambda + U_R', # The predictions! 12/13/22 +X1%*%B1
    U_Pred  = '(U_cond %*% U_F) %*% Lambda + U_cond %*% U_R', # The predictions! 12/13/22 +X_F_test %*%B1
    U_Test  = '(U_cond %*% U_F) %*% Lambda_beta %*% t(X_Env_Test)', # The predictions! updated 01/02/23
    Eta_test = 'F[lines_pred,] %*% Lambda_beta %*% t(X_Env_Test)',
    h2 = '(colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,])/(colSums(Lambda^2)+1/tot_Eta_prec[1,])'
  )
}

# initialize the posterior database
MegaLMM_state = clear_Posterior(MegaLMM_state)

#-------------------------------------------------------------------------------
# Run the MCMC and diagnose convergence
#-------------------------------------------------------------------------------
n_iter = 100
MegaLMM_state$run_parameters$which_sampler$Y = 2
for(i in 1:10) {
  if(i > 3) MegaLMM_state$priors$Lambda_prior$fit_X = T
  print(sprintf('Burnin run %d',i))
  # Factor order doesn't "mix" well in the MCMC.
  # We can help it by manually re-ordering from biggest to smallest
  MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6)
  # clear any previous collected samples because we've re-started the chain 
  MegaLMM_state = clear_Posterior(MegaLMM_state)
  # Draw n_iter new samples, storing the chain
  MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)
  # # make diagnostic plots
  # traceplot_array(MegaLMM_state$Posterior$Lambda,name = file.path(run_ID,'Lambda.pdf'))
  # # traceplot_array(MegaLMM_state$Posterior$U_Test,name = 'U_Test.pdf',
  # #                 facet_dim = 3)
  # traceplot_array(MegaLMM_state$Posterior$U_Test,name = file.path(run_ID,'U_Test.pdf'))
  # traceplot_array(MegaLMM_state$Posterior$B2_F,facet_dim=3,name = file.path(run_ID,'B2_F.pdf'))
  # traceplot_array(MegaLMM_state$Posterior$Lambda_beta,facet_dim=2,name = file.path(run_ID,'Lambda_beta.pdf'))
  # traceplot_array(MegaLMM_state$Posterior$Lambda_beta_var,facet_dim=2,name = file.path(run_ID,'Lambda_beta_var.pdf'))
  print(sprintf('Completed %d burnin samples', MegaLMM_state$current_state$nrun))
  }
MegaLMM_state = clear_Posterior(MegaLMM_state)

#-------------------------------------------------------------------------------
#### Collect posterior samples
# If we think the model is reasonable converged to stationary, we can now collect
# posterior samples.
n_iter = 250
for(i in 1:4) {
  print(sprintf('Sampling run %d',i))
  MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter) 
  MegaLMM_state = save_posterior_chunk(MegaLMM_state)
  print(MegaLMM_state)
}

#-------------------------------------------------------------------------------
### MegaLMM Predictions
MegaLMM_state$Posterior = reload_Posterior(MegaLMM_state)
if(X_Env_Model=="Zero"){
  U_Train = get_posterior_mean(MegaLMM_state$Posterior$U_Train)
  U_Pred =  get_posterior_mean(MegaLMM_state$Posterior$U_Pred)
  Eta_mean = load_posterior_param(MegaLMM_state,'Eta_mean')
}else{
  U_Train = get_posterior_mean(MegaLMM_state$Posterior$U_Train)
  U_Pred =  get_posterior_mean(MegaLMM_state$Posterior$U_Pred)
  U_Test =  get_posterior_mean(MegaLMM_state$Posterior$U_Test)
  Eta_mean = load_posterior_param(MegaLMM_state,'Eta_mean')
  Eta_test = get_posterior_mean(MegaLMM_state$Posterior$Eta_test)
}

#-------------------------------------------------------------------------------
# SECTION 06: save results
#-------------------------------------------------------------------------------

# save results
if(X_Env_Model=="Zero"){
  results=list(Y_Test=Y_Test,
               Y_Train=Y_Train,
               U_Train=U_Train,
               U_Pred=U_Pred,
               Eta_mean=Eta_mean,
               rrBLUP=rrBLUP)
}else{
  results=list(Y_Test=Y_Test,
               Y_Train=Y_Train,
               U_Train=U_Train,
               U_Test=U_Test,
               U_Pred=U_Pred,
               Eta_mean=Eta_mean,
               Eta_test=Eta_test,
               rrBLUP=rrBLUP)
}
save(results,file=paste0(DIR_Output,run_ID,".RData"))
