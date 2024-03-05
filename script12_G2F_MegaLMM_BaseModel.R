# =================================================================
# script12_G2F_MegaLMM_BaseModel
# =================================================================

# Clean workspace
rm(list=ls())

# devtools::install_github('deruncie/MegaLMM@restructure')
library(rrBLUP)
library(MegaLMM); K=NULL
library(data.table)
library(lineup)
library(cvTools)
library(foreach)
library(doParallel)
library(rrBLUP)
library(dplyr)
library(plyr)
# library(data.tree) #list str

#-------------------------------------------------------------------------------
# SECTION 01: set up DIR_Output
#-------------------------------------------------------------------------------
Date="2023-01-05"
DIR_Geno="../../P5G2F/02_GRMs/DER/"
DIR_Pheno_FMT2="../05_MegaLMM_PhenoDat_Format2/"
DIR_Output="../12_G2F_MegaLMM_BaseModel/"
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)

#-------------------------------------------------------------------------------
# SECTION 02: set up global parameters
#-------------------------------------------------------------------------------
traitname = commandArgs(t=T)[1]
X_Env_Model = commandArgs(t=T)[2]
CVScn = commandArgs(t=T)[3]
ID_Rep = as.numeric(commandArgs(t=T)[4])
# ID_Fold_Level2= as.numeric(commandArgs(t=T)[5])
if(is.na(traitname)) traitname="Yield_Mg_ha"#
if(is.na(X_Env_Model)) X_Env_Model="Zero" #"Zero","State","Tester" "State:Tester" "State+Tester"
if(is.na(CVScn)) CVScn="CV2" #CVScn="CV1"#CVScn="CV3NewState"; CVScn="CV3NewTrial"
if(is.na(ID_Rep)) ID_Rep=1
# if(is.na(ID_Fold)) ID_Fold=1  #ID_fold for Y_Train_CV2, is 1-5 for CV3/4; 0 for CV1/2
# if(is.na(ID_Fold_Level2)) ID_Fold_Level2=20 #actual Fold ID for CVScn under running

# set up MegaLMM running parameters 
centering_trn_set=F # 12/15/23 DER said: if scale Y, it does not matter whether substract global mean
scale_Y = T # Step1_intcpt=F; XF_intcpt=T
# get run_ID
run_ID = sprintf('G2F_%s_Xenv_%s_%s_Run%02d',
                 gsub("_","",traitname),X_Env_Model,CVScn,ID_Rep)
# if(dir.exists(run_ID)) unlink(run_ID,recursive = T)
print(paste0("####################",run_ID,"####################"))

#-------------------------------------------------------------------------------
# SECTION 03: Load data and set up MegaLMM Input
#-------------------------------------------------------------------------------

# load Pheno_FMT2
list.files(DIR_Pheno_FMT2)
ls2=load(sprintf("%s%sMegaLMM_PhenoDat_Format2.RData",DIR_Pheno_FMT2,"2023-09-15"))
ls2 #"Phe_Trn_Fmt1Wide" "Phe_Trn_Fmt2Wide" "Phe_Trn_Fmt2Long" "K_info"
sapply(Phe_Trn_Fmt2Wide, dim)
Y = Phe_Trn_Fmt2Wide[[traitname]]
names(Phe_Trn_Fmt2Wide)

# centering training set by global mean of training set
if(centering_trn_set){
  globalmean = mean(as.matrix(Y),na.rm = T)
  Y = Y - globalmean
}

# load K matrix for P1 (created by DER)
K0 = readRDS(paste0(DIR_Geno,'A_P1.rds'))

# load K matrix for P2 (created by DER)
K_Tester = readRDS(paste0(DIR_Geno,'A_P2_new.rds'))

dim(Y) #1702  302
dim(K)

#-------------------------------------------------------------------------------
# Prepare K for Parent1

# subset K0 for those Hybrids existing in Pheno data
if(all(K_info$Hybrid %in% rownames(K0))){
  K=K0[K_info$Hybrid,K_info$Hybrid]
}

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
# K_P1[1:6,1:6]
# checking whether Parent1 in the same order between Pheno and Geno data
stopifnot(identical(rownames(Y),rownames(K_P1)))
K = K_P1 # to be convenient, I assign K_P1 to K, and use K for K_P1 hereafter

#-------------------------------------------------------------------------------
# checking before running MegaLMM
line_data = data.frame(Line = rownames(Y))
stopifnot(all(sapply(list(Y=Y, K=K_P1),function(x)
  identical(rownames(x),line_data$Line)) ))
dim(Y) #1702  302
dim(K) # 1700 1700

#-------------------------------------------------------------------------------
# SECTION 04: run rrBLUP
#-------------------------------------------------------------------------------
# rrBLUP/GBLUP can only predict within MET, cannot predict TPE

# cl <- makeCluster(RcppParallel::defaultNumThreads()-1)
cl <- makeCluster(5)
registerDoParallel(cl)
i=1
X_F = model.matrix(~1,line_data)
rrBLUP = foreach(i = 1:ncol(Y),.combine = cbind) %dopar% {
  print(i)
  library(rrBLUP)
  y = Y[,i]
  j = !is.na(y)
  Kj = K[j,j]
  while(T) {
    offset = 0
    # r = try(mixed.solve(y_train[j],X = X_F_train[j,],K = Kj+diag(offset,nrow(Kj))))
    r = try(mixed.solve(y[j],K = Kj+diag(offset,nrow(Kj))))
    if(is.character(r)){
      offset = offset + 0.01
    } else if(offset > .1){
      break
    } else{
      break
    }
  }
  a2=X_F %*% r$beta + as.matrix(K[,rownames(Y)[j]] %*% MASS::ginv(Kj)) %*% r$u
  # a2= as.matrix(K[lines_test,rownames(Y_Train)[j]] %*% MASS::ginv(Kj)) %*% r$u
  a2
}
rownames(rrBLUP)=rownames(Y)
colnames(rrBLUP)=colnames(Y)
stopCluster(cl)
# dim(rrBLUP)
# rrBLUP[1:3,1:4]

#-------------------------------------------------------------------------------
# SECTION 05: run MegaLMM
#-------------------------------------------------------------------------------

library(MegaLMM)
# default_values <- formals(MegaLMM_control)
# print(default_values)
Y = scale(Y)

#-------------------------------------------------------------------------------
### Set control parameters
run_parameters = MegaLMM_control(
  h2_divisions = 20,
  # Each variance component is allowed to explain between 0% and 100% of the
  # total variation. How many segments should the range [0,100) be divided
  # into for each random effect?
  K_eigen_tol = 1e-1, #new argument added by DER on 01/23/23
  # scale_Y = c(T,T),
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
MegaLMM_state = setup_model_MegaLMM(
  Y = Y,
  # The n x p trait matrix
  # formula = ~ 0 + (1|Line),
  formula = ~ 1 + (1|Line),
  # This is syntax like lme4 for mixed effect models.
  # We specify a fixed effect of population and a random effect for genotype (Line)
  data = line_data,
  # the data.frame with information for constructing the model matrices
  #extra_regressions = list(X=X_F,factors=T), #factors=T is required otherwise get error
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
### MegaLMM Predictions
MegaLMM_state$Posterior = reload_Posterior(MegaLMM_state)
U_Train = get_posterior_mean(MegaLMM_state$Posterior$U_Train)
Eta_mean = load_posterior_param(MegaLMM_state,'Eta_mean')
Lambda_samples = load_posterior_param(MegaLMM_state,'Lambda')
F_samples = load_posterior_param(MegaLMM_state,'F')
# F_posterior_mean = get_posterior_mean(MegaLMM_state$Posterior$F)
G = MegaLMM_state$Posterior$G
R = MegaLMM_state$Posterior$R
Pcov=MegaLMM_state$Posterior$G+MegaLMM_state$Posterior$R
# aaply: For each slice of an array, apply function, keeping results as an array.
Pcor=aaply(Pcov, 1, cov2cor)

#-------------------------------------------------------------------------------
# SECTION 06: save results
#-------------------------------------------------------------------------------

# save results
results=list(Y=Y,
             U_Train=U_Train,
             G=G,
             R=R,
             Pcov=Pcov,
             Pcor=Pcor,
             F_samples=F_samples,
             Lambda_samples=Lambda_samples,
             Eta_mean=Eta_mean,
             rrBLUP=rrBLUP
             )
sapply(results, dim)
# sapply(results, function(x) x[1:3,1:4])
save(results,file=paste0(DIR_Output,run_ID,".RData"))
