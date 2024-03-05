#====================================================================
# script01_Obtain_Tester_FamSize
#=====================================================================

# Aims: Obtain_FamSize_and_Evaluate_influence_of_Parent2_name_swap


# Clean workspace
rm(list=ls())

library(tibble)

Date = "2023-09-07"
DIR_TRN="../../P5G2F/Maize_GxE_Competition_Data/Training_Data/"
DIR_Output="../01_Tester_FameSize_Obtained_From_DiffParent2/"
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)
dir.exists(DIR_Output)

# define global variables
traitnames = rev(c("Silk_DAP_days","Plant_Height_cm","Yield_Mg_ha"))
FamSize_cutoff=50 #tester family size cutoff

# ------------------------------------------------------------------------------
# Data Loading - Training phenotypic data
# ------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Training data - pheno
phe_trn=read.csv(paste0(DIR_TRN,"1_Training_Trait_Data_2014_2021.csv"),
                 stringsAsFactors = F)
phe_trn=add_column(.data = phe_trn,
                   .after = "Env",
                   State=substr(phe_trn$Env,1,2))
dim(phe_trn) #136012     26
head(phe_trn)
sapply(list(phe_trn$Env, phe_trn$Year, phe_trn$Field_Location, phe_trn$Experiment,
            phe_trn$Hybrid, phe_trn$Hybrid_orig_name), 
       function(x) length(unique(x)))
# [1] 217    8   45   11 4683 4959

# summary:
# => 217  Environments
# => 8    years
# => 45   Field Locations
# => 4683 Unique curated Hybrid names
# => 4959 Unique original Hybrid names
unique(phe_trn$Env)

# From How many Stats
length(unique(substr(phe_trn$Env,1,2)))#21
library(usdata)
df_State=data.frame(StateNameAbbr=unique(substr(phe_trn$Env,1,2)),
           StateNameFull=abbr2state(unique(substr(phe_trn$Env,1,2)))
           )
# ON (Ontario)
# GE=Gottingen


# DesCat$Env[DesCat$DesCat=="SingleRepTrial"]
# "MOH1_1_2018"     "MOH1_2_2018"     "TXH1-Early_2018" "MOH1_1_2020"     "MOH1_2_2020"
Envs_SingleRepTrial=c("MOH1_1_2018","MOH1_2_2018","TXH1-Early_2018","MOH1_1_2020","MOH1_2_2020")
for(traitname in traitnames){
  for(Env in Envs_SingleRepTrial){
    tmp=phe_trn[phe_trn$Env==Env,traitname]
    print(paste(traitname,Env,sum(!is.na(phe_trn[phe_trn$Env==Env,traitname]))))
  }
}

# ------------------------------------------------------------------------------
# filter out local check/commercial hybrids of which trait values are NAs
phe_trn=phe_trn[!is.na(phe_trn$Yield_Mg_ha),]
dim(phe_trn) #128167     26
sapply(list(phe_trn$Env, phe_trn$Year, phe_trn$Field_Location, phe_trn$Experiment,
            phe_trn$Hybrid, phe_trn$Hybrid_orig_name), 
       function(x) length(unique(x)))
# [1]  217    8   45   11 4513 4753

for(traitname in traitnames){
  print(paste(traitname,sum(is.na(phe_trn[,traitname])),sep="::"))
}
# [1] "Yield_Mg_ha::0"
# [1] "Plant_Height_cm::9716"
# [1] "Silk_DAP_days::26804"

# ------------------------------------------------------------------------------
# number of Envs with missing Range and Pass
length(unique(phe_trn$Env[is.na(phe_trn$Range)|is.na(phe_trn$Pass)])) #43
length(unique(phe_trn$Env[!(is.na(phe_trn$Range)) | (!is.na(phe_trn$Pass))])) #189

# keep raw phe_trn
phe_trn_raw = phe_trn
# I need to assign phe_trn to phe_trn_raw because I used phe_trn as local variable
# in the below command lines
dim(phe_trn_raw) #128167     27

# ------------------------------------------------------------------------------
# get fam size
# ------------------------------------------------------------------------------

dat = phe_trn_raw

traitname = traitnames[1]
head(dat)

i=1
FamSize_dat=NULL;FamSize_geno=NULL
for(traitname in traitnames){
  # extract clean phenotypic data for trait i
  phe_trn = dat[!is.na(dat[,traitname]),]
  head(phe_trn)
  # dim(phe_trn)
  # dim(dat_modelfit)
  # length(unique(phe_trn$Env)) #179
  
  for(i in 1: length(unique(phe_trn$Env)) ) {
    print(paste("############",traitname,"::",i,"####################"))
    # for(i in 1: 50) {
    Env=unique(phe_trn$Env)[i]
    # Env="DEH1_2014" #LMM sigular fit
    
    # extract phe_trn_envi
    phe_trn_envi=phe_trn[phe_trn$Env==Env,]
    head(phe_trn_envi)
    sum(is.na(phe_trn_envi[,traitname]))
    dim(phe_trn_envi)
    
    # add a column of Parent2
    phe_trn_envi = add_column(.data = phe_trn_envi,
                              .after = "Hybrid",
                              Parent2=sapply(strsplit(phe_trn_envi$Hybrid,"/"), "[",2)
                              )
    
    # compute FamSize
    FamSize_dat_envi=aggregate(Hybrid~Hybrid_Parent2,data = phe_trn_envi,FUN = function(x) length(unique(x)))
    FamSize_geno_envi=aggregate(Hybrid~Parent2,data = phe_trn_envi,FUN = function(x) length(unique(x)))
    sum(FamSize_dat_envi$Hybrid)
    sum(FamSize_geno_envi$Hybrid)
    colnames(FamSize_dat_envi)=c("Hybrid_Parent2","FamSizedat")
    colnames(FamSize_geno_envi)=c("Parent2","FamSizegeno")
    FamSize_dat_envi = data.frame(traitname=traitname,Env=Env,FamSize_dat_envi)
    FamSize_geno_envi = data.frame(traitname=traitname,Env=Env,FamSize_geno_envi)
    
    # select Parent2 according to family size
    FamSize_dat_envi = FamSize_dat_envi[FamSize_dat_envi$FamSizedat>=FamSize_cutoff,]
    FamSize_geno_envi = FamSize_geno_envi[FamSize_geno_envi$FamSizegeno>=FamSize_cutoff,]
    
    if(traitname==traitnames[1] & i==1){
      FamSize_dat = FamSize_dat_envi
      FamSize_geno = FamSize_geno_envi
    }else{
      FamSize_dat = rbind(FamSize_dat, FamSize_dat_envi)
      FamSize_geno = rbind(FamSize_geno, FamSize_geno_envi)
    }
    FamSize_dat_envi=NULL
    FamSize_geno_envi=NULL
  }
}
sapply(list(FamSize_dat=FamSize_dat,FamSize_geno=FamSize_geno), dim)
#      FamSize_dat FamSize_geno
# [1,]         820          816
# [2,]           4            4
lapply(list(FamSize_dat=FamSize_dat,FamSize_geno=FamSize_geno), head)

# --------------------------------------------------------------------
# add a column of EnvNew
FamSize_dat = add_column(.data = FamSize_dat,
                         .after = "Hybrid_Parent2",
                         EnvNew = paste(FamSize_dat$Env,
                                        FamSize_dat$Hybrid_Parent2,sep = "::")
)
FamSize_geno = add_column(.data = FamSize_geno,
                          .after = "Parent2",
                          EnvNew = paste(FamSize_geno$Env,
                                         FamSize_geno$Parent2,sep = "::")
)
lapply(list(FamSize_dat,FamSize_geno), head)
sapply(list(FamSize_dat,FamSize_geno), dim)
# [,1] [,2]
# [1,]  820  816
# [2,]    5    5

# --------------------------------------------------------------------
# output FamSize
if(all(paste(FamSize_geno$traitname,FamSize_geno$EnvNew,sep = "::") %in% 
       paste(FamSize_dat$traitname,FamSize_dat$EnvNew,sep = "::"))){
  save(FamSize_dat,FamSize_geno,file = sprintf("%s%sFamSize_based_on_different_Parent2.RData",DIR_Output,Date))
}


table(FamSize_dat$traitname)
# Plant_Height_cm   Silk_DAP_days     Yield_Mg_ha 
#              281             235             304 
lapply(list(FamSize_dat,FamSize_geno), head)
sum(duplicated(paste(FamSize_dat$traitname,FamSize_dat$Env,FamSize_dat$Hybrid_Parent2)))#0
sum(duplicated(paste(FamSize_geno$traitname,FamSize_geno$Env,FamSize_geno$Parent2)))#0

all(paste(FamSize_geno$traitname,FamSize_geno$Env,FamSize_geno$Parent2,sep = "::") %in%
      paste(FamSize_dat$traitname,FamSize_dat$Env,FamSize_dat$Hybrid_Parent2,sep = "::")
    )#TRUE


# ------------------------------------------------------------------------------
# Conclusions
# ------------------------------------------------------------------------------

# 1) with the FamSize_cutoff=50,
#   - the number of EnvNew from the column of "Hybrid_Parent2" = 820
#   - the number of EnvNew from Parent2 from splitting the column of "Hybrid"=816
#   - all the 816 EnvNew are included in the 820 EnvNew
# 2) Therefore, we probably can use the column of "Hybrid_Parent2" as Tester
