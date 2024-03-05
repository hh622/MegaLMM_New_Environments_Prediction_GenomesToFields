#=====================================================
# script02_Phenodat_Remove_Outliers_in_each_TesterFam
#=====================================================

#Aims: Phenodat_Remove_Outliers_in_each_TesterFam

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

Date = "2023-09-10"
DIR_TRN="../../P5G2F/Maize_GxE_Competition_Data/Training_Data/"
DIR_Output="../02_Phenodat_OutlierRemoval/"
DIR_FamSize="../01_Tester_FameSize_Obtained_From_DiffParent2/"
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)
dir.exists(DIR_Output)

# define global variables
traitnames = rev(c("Silk_DAP_days","Plant_Height_cm","Yield_Mg_ha"))
FamSize_cutoff=50 #tester family size cutoff
Pheno_Trn_Clean=NULL

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
# keep raw phe_trn
phe_trn_raw = phe_trn
# I need to assign phe_trn to phe_trn_raw because I used phe_trn as local variable
# in the below command lines
dim(phe_trn_raw) #128167     27
head(phe_trn_raw)

# ------------------------------------------------------------------------------
# loading FamSize
list.files(DIR_FamSize)
ls = load(file=paste0(DIR_FamSize,"2023-09-07FamSize_based_on_different_Parent2.RData"))
ls #"FamSize_dat"  "FamSize_geno"
head(FamSize_dat)
dim(FamSize_dat) #820   5
table(FamSize_dat$traitname)
# Plant_Height_cm   Silk_DAP_days     Yield_Mg_ha 
#             281             235             304 
range(FamSize_dat$FamSizedat)
# ------------------------------------------------------------------------------
# Step 1: remove outlier data points based on joint distribution across Envs
# ------------------------------------------------------------------------------
traitname = traitnames[1]
Pheno_Trn_Clean=foreach(traitname = traitnames,.combine = "rbind") %do% {
  dat=phe_trn_raw[,c(colnames(phe_trn_raw)[1:14],traitname)]
  dim(dat) #128167     15
  head(dat)
  # extract for selected Tester families only
  dat = dat[paste(dat$Env,dat$Hybrid_Parent2,sep = "::") %in% 
              FamSize_dat$EnvNew[FamSize_dat$traitname==traitname],]
  dim(dat) #92113    15
  
  # NumEnvs=length(unique(phe_trn$Env[!is.na(phe_trn[,traitname])]))
  
  # Env=unique(dat$Env)
  # calculate mean and sd
  Mean = mean(dat[,traitname],na.rm=T)
  SD=sd(dat[,traitname],na.rm=T)
  # calculate pr for each data point
  Pr=sapply(dat[,traitname], function(x){
    if(is.na(x)){
      return(NA)
    }else if(!is.na(x) & x<Mean) {
      pnorm(x,mean = Mean,sd=SD)
    }else if(!is.na(x) & x>Mean) {
      pnorm(x,mean = Mean,sd=SD,lower.tail=F)
    }
  }
  )
  # given total number of data points for each trait across trails and Pr. of observing
  # a specific data point, we compute how many times that each data point is expected to
  # be observed
  NumDatPoints=sum( !is.na(dat[,traitname]) )
  tmp=data.frame(dat[1:14],TraitName=traitname,TraitVal=dat[,traitname],Mean,SD,
                 NumNADatPoints=NumDatPoints,Pr,NumExpected=NumDatPoints*Pr)
  head(tmp)
  dim(tmp)
  
  # if time of observing a data <1, that data point is considered as outlier,
  # and set as NA
  tmp$TraitVal[tmp$NumExpected<1 & !is.na(tmp$NumExpected)]=NA
  print(tmp[tmp$NumExpected<1 & !is.na(tmp$NumExpected),])
  return(tmp)
}
head(Pheno_Trn_Clean)
stopifnot(sum(!is.na(phe_trn_raw[,traitnames])) > sum(!is.na(Pheno_Trn_Clean$TraitVal)))

# summary: how many outlier data points for each trait going to be removed
for(traitname in traitnames){
  print(paste(traitname,sum(is.na(Pheno_Trn_Clean$TraitVal[Pheno_Trn_Clean$TraitName==traitname]))))
}
# [1] "Yield_Mg_ha 3"
# [1] "Plant_Height_cm 692"
# [1] "Silk_DAP_days 956"

# update Pheno_Trn_Clean
Pheno_Trn_Clean = Pheno_Trn_Clean[!is.na(Pheno_Trn_Clean$TraitVal),]

dim(phe_trn_raw)#128167     27
sum(!is.na(phe_trn_raw[,traitnames[1]]))+sum(!is.na(phe_trn_raw[,traitnames[2]]))+
  sum(!is.na(phe_trn_raw[,traitnames[3]]))#347981
dim(Pheno_Trn_Clean)#250818     21

# ------------------------------------------------------------------------------
# Step 2: remove outlier trials based on distribution of population means
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# calculate population mean and its Pr in each EnvNew
traitname = traitnames[1]
res_EnvNewMeans=NULL
res_EnvNewMeans=foreach (traitname = traitnames, .combine = "rbind") %do% {
  mydat = Pheno_Trn_Clean[Pheno_Trn_Clean$TraitName==traitname,]
  head(mydat)
  # add a column of EnvNew
  mydat = add_column(.data = mydat,.after = "Hybrid_Parent2",
                   EnvNew=paste(mydat$Env,mydat$Hybrid_Parent2,sep = "::"))
  # EnvNew_mean 
  EnvNew_mean = aggregate(data=mydat,TraitVal~EnvNew,FUN = mean,na.rm=T)
  dim(EnvNew_mean)
  head(EnvNew_mean)
  
  # colnames(Env_mean)=c("Env","EnvNewMean")
  data.frame(TraitName=traitname,EnvNew=EnvNew_mean[,1],EnvNewMean=EnvNew_mean[,2])
}
head(res_EnvNewMeans)
# calculate mean and sd of population means
df_mean_sd = foreach (traitname = traitnames, .combine = "rbind") %do% {
  data.frame(TraitName=traitname,Mean=mean(res_EnvNewMeans$EnvNewMean[res_EnvNewMeans$TraitName==traitname]),
             Sd=sd(res_EnvNewMeans$EnvNewMean[res_EnvNewMeans$Trait==traitname]))
}

# calculate Pr/pnorm for Env mean given Norm Distribution
# traitname = traitnames[1]
res_EnvNewMeans_Pr=NULL
res_EnvNewMeans_Pr=foreach(traitname = traitnames,.combine = "rbind") %do% {
  mydat = res_EnvNewMeans[res_EnvNewMeans$Trait==traitname,]
  NumEnvNew_Obs=length(unique(mydat$EnvNew))
  
  head(mydat)
  # compute mean and sd for 
  Mean = mean(mydat$EnvNewMean)
  SD=sd(mydat$EnvNewMean)
  
  Pr=sapply(mydat$EnvNewMean, function(x){
    if(x<Mean) {
      pnorm(x,mean = Mean,sd=SD)
    }else{
      pnorm(x,mean = Mean,sd=SD,lower.tail=F)
    }
  }
  )
  data.frame(mydat,Mean,SD,Pr,
             NumEnvNew_Obs,NumExpected_NumEnvsObs=round(NumEnvNew_Obs*Pr,3))
}
head(res_EnvNewMeans_Pr)
res_EnvNewMeans_Pr=res_EnvNewMeans_Pr[order(res_EnvNewMeans_Pr$TraitName,
                                      res_EnvNewMeans_Pr$NumExpected_NumEnvsObs),]

# ------------------------------------------------------------------------------
# output df_OutlierEnvs
df_OutlierEnvs=res_EnvNewMeans_Pr[res_EnvNewMeans_Pr$NumExpected_NumEnvsObs<1,]
dim(df_OutlierEnvs)
head(df_OutlierEnvs)
# output df_OutlierEnvs
if(!file.exists(paste0(DIR_Output,"DF_OutlierEnvNews_based_on_PopulationMeanDistrib.csv"))){
  write.csv(df_OutlierEnvs,
            sprintf("%s%sDF_OutlierEnvNews_based_on_PopulationMeanDistrib.csv",DIR_Output,Date),
            row.names = F)
}

# ------------------------------------------------------------------------------
# exclude outlier trials to get clean data
head(df_OutlierEnvs)
head(Pheno_Trn_Clean)

# get sel_idx
sel_idx=!(paste(Pheno_Trn_Clean$Env,Pheno_Trn_Clean$Hybrid_Parent2,
                Pheno_Trn_Clean$TraitName,sep="::") %in%
  paste(df_OutlierEnvs$EnvNew,df_OutlierEnvs$TraitName,sep = "::"))
table(sel_idx)

# update Pheno_Trn_Clean
Pheno_Trn_Clean=Pheno_Trn_Clean[sel_idx,]
dim(Pheno_Trn_Clean) # 248475     21
head(Pheno_Trn_Clean)

# add a column of EnvNew to clean data
Pheno_Trn_Clean = add_column(.data = Pheno_Trn_Clean,
                             .after = "Hybrid_Parent2",
                             EnvNew = paste(Pheno_Trn_Clean$Env,Pheno_Trn_Clean$Hybrid_Parent2,sep = "::")
                             )

# validating
length(unique(Pheno_Trn_Clean$EnvNew[Pheno_Trn_Clean$TraitName==traitnames[1]]))#302
length(unique(Pheno_Trn_Clean$EnvNew[Pheno_Trn_Clean$TraitName==traitnames[2]]))#278
length(unique(Pheno_Trn_Clean$EnvNew[Pheno_Trn_Clean$TraitName==traitnames[3]]))#231
length(unique(paste(Pheno_Trn_Clean$TraitName,Pheno_Trn_Clean$EnvNew,sep = "::")))+
  dim(df_OutlierEnvs)[1]==dim(FamSize_dat)[1] #TRUE
df_OutlierEnvs[,1:4]
sum(is.na(Pheno_Trn_Clean$TraitVal))==0  #TRUE
dim(Pheno_Trn_Clean) # 248475     22

# add a column of TraitEnvNew
Pheno_Trn_Clean = add_column(.data = Pheno_Trn_Clean,
                             .after = "TraitName",
                             TraitEnvNew=paste(Pheno_Trn_Clean$TraitName,
                                               Pheno_Trn_Clean$EnvNew,
                                               sep = "::")
)
dim(Pheno_Trn_Clean) #248475     23

# ------------------------------------------------------------------------------
# output Clean Pheno data and EnvNew Category
# ------------------------------------------------------------------------------

dim(Pheno_Trn_Clean)#248475     23
length(unique(paste(Pheno_Trn_Clean$TraitName,Pheno_Trn_Clean$EnvNew)))#811
length(unique(Pheno_Trn_Clean$EnvNew[Pheno_Trn_Clean$TraitName==traitnames[1]]))#302
length(unique(Pheno_Trn_Clean$EnvNew[Pheno_Trn_Clean$TraitName==traitnames[2]]))#278
length(unique(Pheno_Trn_Clean$EnvNew[Pheno_Trn_Clean$TraitName==traitnames[3]]))#231
head(Pheno_Trn_Clean)
head(phe_trn_raw)
sum(!is.na(Pheno_Trn_Clean$TraitVal))/(sum(!is.na(phe_trn_raw[,traitnames[1]]))+
                                         sum(!is.na(phe_trn_raw[,traitnames[2]]))+
                                         sum(!is.na(phe_trn_raw[,traitnames[3]])))
# 0.7140476
length(unique(Pheno_Trn_Clean$Hybrid_Parent2))#12
length(unique(phe_trn_raw$Hybrid_Parent2))#77

save(Pheno_Trn_Clean,file=sprintf("%s%sPheno_Trn_Clean.RData",DIR_Output,Date))

