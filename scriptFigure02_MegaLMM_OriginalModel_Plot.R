# ===================================================
# scriptFigure02_MegaLMM_Original_Model_Plot
# ===================================================

# Clean workspace
rm(list=ls())

library(reshape2)
library(ggplot2)
library(foreach)
library(tibble)
library(abind)
library(broom)
library(dplyr)
library(ggpubr)
library(scales)
library(ggbeeswarm)

# library(ggbeeswarm) #plot outliers in a line
source("./custom_functions.R")
source("./meta-analyzing_correlations.R")

DIR_Input_MegaLMM_Results = "../12_G2F_MegaLMM_BaseModel/"
DIR_Input="../09_Prediction_Results_Summary/"
DIR_Output="../Figure02_MegaLMM_OriginalModel/"
DIR_Input_Phedat_Format2="../05_MegaLMM_PhenoDat_Format2/"
DIR_EnvDat="../../P5G2F/04_Env_Dat/2022-12-26/"
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)

Date="2024-02-15"
TraitSel=rev(c("Silk Days","Plant Height","Grain Yield"))
CVScns = c("CV1","CV2","CV3NewTrial","CV3NewState","CV3NewTester","CV4")

# ------------------------------------------------------------------------------
# loading data
# ------------------------------------------------------------------------------

# Loading Formatted Prediction Results 
list.files(DIR_Input)
ls=load(file=sprintf("%s%sG2F_Prediction_Results_Formatted.RData",DIR_Input,"2023-10-16"))
ls #"res_Original" "res_EnvMean" 
head(res_EnvMean)
unique(paste(res_EnvMean$PriorLevel,res_EnvMean$PredictionModel,sep = "::"))
unique(paste(res_EnvMean$PriorLevel[res_EnvMean$CVScn=="CV1"],
             res_EnvMean$PredictionModel[res_EnvMean$CVScn=="CV1"],sep = "::"))

# load Pheno dat Format2
list.files(DIR_Input_Phedat_Format2)
ls=load(sprintf("%s%sMegaLMM_PhenoDat_Format2.RData",
                DIR_Input_Phedat_Format2,"2023-09-15"))
ls #"Phe_Trn_Fmt1Wide" "Phe_Trn_Fmt2Wide" "Phe_Trn_Fmt2Long" "K_info"
lapply(Phe_Trn_Fmt2Wide, function(x) x[1:4,1:6])

sapply(Phe_Trn_Fmt2Wide, function(x) rowSums(!is.na(x)))
range(rowSums(!is.na(Phe_Trn_Fmt2Wide$Yield_Mg_ha))) #4 114
range(rowSums(!is.na(Phe_Trn_Fmt2Wide$Plant_Height_cm)))#4 100
range(rowSums(!is.na(Phe_Trn_Fmt2Wide$Silk_DAP_days)))   #4 81
sapply(Phe_Trn_Fmt2Wide, dim)
#     Yield_Mg_ha Plant_Height_cm Silk_DAP_days
# [1,]        1702            1699          1663
# [2,]         302             278           231
length(unique(sapply(strsplit(colnames(Phe_Trn_Fmt2Wide$Yield_Mg_ha),"::"), "[", 1)))
#195
length(unique(substr(colnames(Phe_Trn_Fmt2Wide$Yield_Mg_ha),1,4)))
# 37
length(unique(sapply(strsplit(colnames(Phe_Trn_Fmt2Wide$Plant_Height_cm),"::"), "[", 1)))
#183
length(unique(sapply(strsplit(colnames(Phe_Trn_Fmt2Wide$Silk_DAP_days),"::"), "[", 1)))
#160
sum(is.na(Phe_Trn_Fmt2Wide[[1]]))/(nrow(Phe_Trn_Fmt2Wide[[1]])*
                                            ncol(Phe_Trn_Fmt2Wide[[1]]))
#0.8710302
sum(is.na(Phe_Trn_Fmt2Wide[[2]]))/(nrow(Phe_Trn_Fmt2Wide[[2]])*
                                            ncol(Phe_Trn_Fmt2Wide[[2]]))
# 0.8703469
sum(is.na(Phe_Trn_Fmt2Wide[[3]]))/(nrow(Phe_Trn_Fmt2Wide[[3]])*
                                            ncol(Phe_Trn_Fmt2Wide[[3]]))
#0.8637626
sort(table(substr(colnames(Phe_Trn_Fmt2Wide$`Grain Yield`),1,2)))
sort(table(substr(colnames(Phe_Trn_Fmt2Wide$`Grain Yield`),1,2)))
tmp_trials = unique(sapply(strsplit(colnames(Phe_Trn_Fmt2Wide[[1]]),"::"), "[", 1))
# how many trials per State
sort(table(substr(tmp_trials,1,2)))
# SD GE KS SC AR MI MN OH DE IL NC ON IN GA MO WI NE NY TX IA 
# 1  3  3  4  5  6  6  6  7  7  7  7  8 11 11 17 18 19 19 30 
# how many experiments per State
sort(table(substr(colnames(Phe_Trn_Fmt2Wide[[1]]),1,2)))
# SD GE SC AR KS MI IL OH IN MN ON GA NC DE MO TX NE NY WI IA 
# 1  3  4  5  6  6  9  9 10 11 11 15 15 17 24 24 26 26 28 52 
mean(table(substr(colnames(Phe_Trn_Fmt2Wide[[1]]),1,2))) #15
# how many experiments per tester
sort(table(sapply(strsplit(colnames(Phe_Trn_Fmt2Wide[[1]]),"::"), "[", 2)))
# CG102    PB80   LH185   LH198   PHK76   PHP02   PHB47   PHT69    LH82 DK3IIH6   PHZ51   LH195 
#   4       6      12      12      18      22      24      25      26      31      56      66 
mean(table(sapply(strsplit(colnames(Phe_Trn_Fmt2Wide[[1]]),"::"), "[", 2)))

# get number of non-NAs in each Env
names(Phe_Trn_Fmt2Wide)=TraitSel
TraitName=TraitSel[1]
df_NumObs_Env=foreach(TraitName=TraitSel, .combine = "rbind") %do% {
  tempdat = Phe_Trn_Fmt2Wide[[TraitName]]
  dim(tempdat)
  data.frame(TraitName=TraitName,
             EnvName=colnames(tempdat),
             NumObs_Env=colSums(!is.na(tempdat)))
}
dim(df_NumObs_Env)
head(df_NumObs_Env)

# read in and prepare lambda
list.files(DIR_Input_MegaLMM_Results)
file=list.files(DIR_Input_MegaLMM_Results,pattern = ".RData")[1]
TraitNames=c("YieldMgha","PlantHeightcm","SilkDAPdays")

# load env data
list.files(DIR_EnvDat)
ls2=load(paste0(DIR_EnvDat,"G2F_EnvDat_QC.RData"))
ls2 #"Wea_Dat"     "Soil_Dat"    "ECs_Dat"     "metadat_Wea"
ECs_Dat[1:3,1:10]
ECs_Dat[1:3,grepl("^TT_",colnames(ECs_Dat))]
sapply(list(Wea=Wea_Dat,Soil=Soil_Dat,ECs=ECs_Dat,metaWea=metadat_Wea), dim)
#      Wea Soil ECs metaWea
# [1,] 235  162 189     238
# [2,] 281   26 764      18
lapply(list(Wea=Wea_Dat,Soil=Soil_Dat,ECs=ECs_Dat,metaWea=metadat_Wea), function(x)
  x[1:4,1:6])
head(metadat_Wea)

# ------------------------------------------------------------------------------
# preparing data for plotting
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# compute Square Root of the Mean Squares (RMS)
# RMS = sqrt(mean(x^2))
list.files(DIR_Input_MegaLMM_Results,pattern = "RData")
res_Lambda_Squared=NULL
i=11
res_Lambda_Squared=foreach(i = 1:length(list.files(DIR_Input_MegaLMM_Results,pattern = "RData")),.combine = "rbind") %do% {
  print(paste(i, "of", length(list.files(DIR_Input_MegaLMM_Results,pattern = "RData"))))
  (file=list.files(DIR_Input_MegaLMM_Results,pattern = "RData")[i])
  # load data
  ls = load(file = paste0(DIR_Input_MegaLMM_Results,file))
  sapply(results, dim)
  
  # get original F and Lambda
  F_samples = results$F_samples
  Lambda_samples = results$Lambda_samples
  sapply(list(F_samples,Lambda_samples), dim)
  
  # get sd of each factor
  F_samples_sd <- apply(F_samples, c(1, 3), sd)
  
  # rescale Lambda: multiply Lambda_samples by sd of F
  Lambda_samples_rescaled = array(0, dim = dim(Lambda_samples))
  for(i in 1:500){
    sd_i = apply(F_samples[i,,],2,sd)
    Lambda_samples_rescaled[i,,] = sweep(Lambda_samples[i,,],1,sd_i,'*')
    # for(j in 1:50){
    #   # Lambda_samples_rescaled[i, j, ]=F_samples_sd[i,j] * Lambda_samples[i,j,]
    #   Lambda_samples_rescaled[i, j, ]= Lambda_samples[i,j,]*F_samples_sd[i,j]
    # }
  }
  dim(Lambda_samples_rescaled)
  # 500  50 278
  # compute Lambda_SS
  Lambda_SS=apply(Lambda_samples_rescaled, c(1, 3), function(x) sum(x^2) )
  dim(Lambda_SS)
  heatmap(Lambda_SS)
  # Let's square Lambda_samples_rescaled
  Lambda_samples_squared <- Lambda_samples_rescaled^2
  # Let's average across 500 samples
  Lambda_samples_squared_mean = apply(Lambda_samples_squared,c(2, 3), mean)
  dim(Lambda_samples_squared_mean)
  Lambda_samples_squared_mean[1:10,1:10]

  # plot(lambda_RMS)
  # accumulate results
  (TraitName = sapply(strsplit(file,"_"), "[",2))
  # modify TraitName
  (TraitName=ifelse(TraitName=="YieldMgha","Grain Yield",
          ifelse(TraitName=="PlantHeightcm","Plant Height",
                 ifelse(TraitName=="SilkDAPdays","Silk Days","Wrong Trait Name!"))))
  (RunID = as.numeric(gsub("\\D", "", sapply(strsplit(file,"_"), "[",6))))
  data.frame(TraitName=TraitName,
             RunID=RunID,
             Lambda_ID=1:nrow(Lambda_samples_squared_mean),
             Lambda_samples_squared_mean_max = apply(Lambda_samples_squared_mean,1,max)
             )
}
dim(res_Lambda_Squared)#1500    4
head(res_Lambda_Squared)

# ------------------------------------------------------------------------------
# compute Lambda sampe mean
# notes: 
#  - cor among 500 samples within each run are all positively correlated
#  - however, sample means between runs are not always positively correlated
#    Therefore, cannot simply average across runs

res_Lambda_SampleMean=list()
RunID=1
TraitName = TraitNames[1]
for(TraitName in TraitNames) {
  # Lambda_SampleMean_Lst=list()
  res_Lambda_SampleMean[[TraitName]]=list()
  for(RunID in 1:10) {
    print(paste(TraitName,RunID,sep = "::"))
    list.files(DIR_Input_MegaLMM_Results)
    (file=list.files(DIR_Input_MegaLMM_Results,
                     pattern = sprintf("%s_Xenv_Zero_CV2_Run%02d",TraitName,RunID)))
    
    # load data
    ls = load(file = paste0(DIR_Input_MegaLMM_Results,file))
    sapply(results, dim)
    # extract raw lambda - 3d
    lambda_raw = results$Lambda_samples
    dim(lambda_raw) #500  50 278
    # average lambda across 500 samplings
    lambda_samplemean <- apply(lambda_raw, c(2, 3), mean)
    dim(lambda_samplemean) #50 278
    # modify TraitName
    # (TraitName=ifelse(TraitName=="YieldMgha","Grain Yield",
    #         ifelse(TraitName=="PlantHeightcm","Plant Height",
    #                ifelse(TraitName=="SilkDAPdays","Silk Days","Wrong Trait Name!"))))
    # Lambda_SampleMean_Lst[[RunID]]=lambda_samplemean
    res_Lambda_SampleMean[[TraitName]][[RunID]]=lambda_samplemean
    # dim(Lambda_SampleMean_Lst[[RunID]])
    # print(head(apply(lambda_samplemean, 1, function(x) sqrt(sum(x^2)) )))
  }
}
library(data.tree)
FromListSimple(res_Lambda_SampleMean)
sapply(res_Lambda_SampleMean[["YieldMgha"]], dim)
dim(res_Lambda_SampleMean[["YieldMgha"]][[1]])

# ------------------------------------------------------------------------------
# regressions of factor loadings (x-axis) on latitude, longitude, state, tester
TraitName = names(res_Lambda_SampleMean)[1]
RunID = 1

# add environmental variables to factor loading matrix
res_lm=NULL
res_lm=foreach(TraitName = names(res_Lambda_SampleMean), .combine = "rbind") %do% {
  foreach(RunID = 1:10, .combine = "rbind") %do% {
    tmp = res_Lambda_SampleMean[[TraitName]][[RunID]]
    
    # add latitude, longitude, state, tester
    tmp = t(tmp)
    tmp[1:3,1:4]
    tmp=rownames_to_column(.data = data.frame(tmp),var = "Loc_Year_Tester")
    tmp=add_column(.data = tmp,
                   .after = "Loc_Year_Tester",
                   Env=sapply(strsplit(tmp$Loc_Year_Tester,"::"),"[",1),
                   State=substr(tmp$Loc_Year_Tester,1,2),
                   Tester=sapply(strsplit(tmp$Loc_Year_Tester,"::"),"[",2)
                   )
    tmp=cbind(left_join(tmp[,1:4],metadat_Wea,by="Env"),
               tmp[,-1*1:4])
    head(tmp)
    ID_lambda = 1
    # model fitting
    foreach(ID_lambda = 1:50, .combine = "rbind") %do% {
      print(paste(TraitName,RunID,ID_lambda,sep = "::"))
      # extract response variable
      y = tmp[,grepl("X",colnames(tmp))][,ID_lambda]
      # fit State
      fit=NULL;fit=lm(y ~ State,data = tmp)
      pvalue = as.numeric(glance(fit)$p.value)
      (df_State=data.frame(TraitName=TraitName,RunID=RunID,ID_lambda=ID_lambda,
                           Predictor="State",Pvalue=pvalue))
      # fit Tester
      fit=NULL;fit=lm(y ~ Tester,data = tmp)
      pvalue = as.numeric(glance(fit)$p.value)
      (df_Tester=data.frame(TraitName=TraitName,RunID=RunID,ID_lambda=ID_lambda,
                           Predictor="Tester",Pvalue=pvalue))
      
      # fit latitude
      fit=NULL;fit=lm(y ~ Weather_Station_Latitude,data = tmp)
      summary(fit)
      pvalue = as.numeric(glance(fit)$p.value)
      df_Latitude=data.frame(TraitName=TraitName,RunID=RunID,ID_lambda=ID_lambda,
                          Predictor="Latitude",Pvalue=pvalue)
      
      # fit State
      fit=NULL;fit=lm(y ~ Weather_Station_Longitude,data = tmp)
      pvalue = as.numeric(glance(fit)$p.value)
      df_Longitude=data.frame(TraitName=TraitName,RunID=RunID,ID_lambda=ID_lambda,
                          Predictor="Longitude",Pvalue=pvalue)
      df=rbind(df_Latitude,df_Longitude,df_State,df_Tester)
      return(df)
    }
  }
}
dim(res_lm)
head(res_lm)

# ------------------------------------------------------------------------------
# MegaLMMInputFormat2 -Matrix plot - axis = Env::Tester (Site_Year_Tester)
# ------------------------------------------------------------------------------

# ----------------------------------------------------------------------
# get dat_Plot
FamSize_cutoff=50 #tester family size cutoff
names(Phe_Trn_Fmt2Long) #"Yield_Mg_ha"     "Plant_Height_cm" "Silk_DAP_days"
TraitName="Yield_Mg_ha"
dat_Plot = Phe_Trn_Fmt2Long[[TraitName]]
dat_Plot=dat_Plot[order(dat_Plot$Year),]
dat_Plot$Year=factor(dat_Plot$Year)
dat_Plot$EnvNew=factor(dat_Plot$EnvNew,levels = unique(dat_Plot$EnvNew))
dat_Plot$Hybrid_Parent1=factor(dat_Plot$Hybrid_Parent1,
                               levels = rev(unique(dat_Plot$Hybrid_Parent1)))

# ----------------------------------------------------------------------
# Prepare geom_text
# Convert the "Env" variable to numeric
head(dat_Plot)
dat_Plot$EnvNew <- as.numeric(dat_Plot$EnvNew)

# Compute boundaries for each fill color
color_boundaries <- dat_Plot %>%
  group_by(Year) %>%
  summarize(x_min = min(EnvNew), x_max = max(EnvNew)) %>%
  ungroup()
color_boundaries$x_bound_right = color_boundaries$x_max+0.5
# Prepare y position for geom_text
dat_Plot$Year = as.integer(as.character(dat_Plot$Year))
Hybrid_Parent1_CumSums=foreach(Year = unique(dat_Plot$Year),.combine = "c") %do% {
  length(unique(dat_Plot$Hybrid_Parent1[dat_Plot$Year %in% c(min(dat_Plot$Year):Year)]))
}
color_boundaries$yposi = length(unique(dat_Plot$Hybrid_Parent1))-Hybrid_Parent1_CumSums
color_boundaries$yposi[color_boundaries$Year==2020]=color_boundaries$yposi[color_boundaries$Year==2018]
color_boundaries$yposi[color_boundaries$Year==2021]=color_boundaries$yposi[color_boundaries$Year==2018]

# ----------------------------------------------------------------------
# ggplot
dat_Plot$Year=as.factor(dat_Plot$Year)
dat_Plot$EnvNew <- as.factor(dat_Plot$EnvNew)
head(dat_Plot)
head(color_boundaries)
length(unique(dat_Plot$EnvNew))#302
length(unique(dat_Plot$Env))#195
length(unique(dat_Plot$Hybrid))#4002
length(unique(dat_Plot$Hybrid_Parent1))#1702
DatStr_Plot =
  ggplot(dat_Plot, aes(
    y=Hybrid_Parent1,
    x=EnvNew,
    fill=Year)) + 
  geom_raster() + 
  # ylim(0,length(levels(dat_Plot$Hybrid_Parent1))+10)+
  # scale_fill_gradient(low="white", high="blue",na.value = 'white') +
  # labs(y="Hybrid Parent1", x="Location::Year::Tester") +
    labs(y="P1", x="Experiment") +
  theme_classic() + theme(axis.text=element_blank(),
                          axis.ticks = element_blank(),
                          axis.title = element_text(size=12), #, face="bold"
                          plot.title = element_text(size=12),
                          legend.position = 'none' #c(0.1,0.05),
                          # legend.position = 'bottom',
                          # legend.direction = "horizontal",
                          # legend.justification='right'
  )+
  geom_text(data = color_boundaries[color_boundaries$Year%in%c(2014,2015,2018:2019),],
            aes(x = (x_min+x_max)/2,
                y=yposi,
                label=Year),inherit.aes = F,vjust = 1
            
  )+
  geom_text(data = color_boundaries[color_boundaries$Year==2016,],
            aes(x = (x_min+x_max)/2,
                y=yposi,
                label=Year),
            inherit.aes = F,vjust = 0.5,hjust=1,angle=90
            
  )+
  geom_text(data = color_boundaries[color_boundaries$Year==2017,],
            aes(x = (x_min+x_max)/2,
                y=yposi,
                label=Year),
            inherit.aes = F,vjust = 0.5,hjust=1,angle=90
            
  )+
  geom_text(data = color_boundaries[color_boundaries$Year==2020,],
            aes(x = x_min*0.98,
                y=yposi,
                label=Year),inherit.aes = F,hjust = 0
            
  )+
    geom_text(data = color_boundaries[color_boundaries$Year==2021,],
              aes(x=x_max*1.01,
                  y=yposi,
                  label=Year),inherit.aes = F,hjust = 1
              
    )
  
# ------------------------------------------------------------------------------
# ggplot
# ------------------------------------------------------------------------------
#  a)	Boxplots: CV2,CV1 using old MegaLMM
#  b)	Scatterplot: CV2 vs GBLUP using old MegaLMM for 1 trait

# exclude sample size too small for plot a) and b)
res = res_EnvMean 
sum(res$NumObs_TSTSet<11)/length(res$NumObs_TSTSet) #0.02296325
res = res_EnvMean[res$NumObs_TSTSet>=11, ]
resPlot = res[res$CVScn %in% c("CV1","CV2") &
                        res$PriorLevel=="O" &
                        res$PredictionModel=="Env",]
dim(resPlot)
head(resPlot)
unique(resPlot$Method)
table(paste(resPlot$CVScn,resPlot$TraitName,resPlot$Method,sep = "::"))

# -------------------------------------------------------------
# Boxplots - suppl: CV2,CV1 using old MegaLMM
# prepare data
resPlot_bxp_CV1=
dcast(resPlot[resPlot$CVScn=="CV1",], TraitName+CVScn+TestSetID+EnvName~Method, 
      value.var ="accu_r")
resPlot_bxp_CV2=
  dcast(resPlot[resPlot$CVScn=="CV2",], TraitName+CVScn+TestSetID+EnvName~Method, 
        value.var ="accu_r")
resPlot_bxp_CV1$accu_diff=resPlot_bxp_CV1$MegaLMM-resPlot_bxp_CV1$GBLUP
resPlot_bxp_CV2$accu_diff=resPlot_bxp_CV2$MegaLMM-resPlot_bxp_CV2$GBLUP
resPlot_bxp=rbind(resPlot_bxp_CV1,resPlot_bxp_CV2)
head(resPlot_bxp)
df_mean = aggregate(accu_diff~TraitName+CVScn,data = resPlot_bxp,FUN = mean)

axis.title.size=12
plot_bxp_suppl=
  ggplot(resPlot_bxp, 
       aes(x = TraitName, y = accu_diff,fill=CVScn)) +
  # labs(y = expression(paste("cor(y, ", hat(g), ")")))+
  labs(y=expression(r[MegaLMM] - r[GBLUP]))+
  geom_boxplot(outlier.size = 0.7)+
  theme_bw()+
  theme(axis.text.y = element_text(angle =90,hjust=0.5),
        axis.title.y= element_text(size=axis.title.size),
        legend.title = element_blank(),
        legend.background=element_blank(),
        legend.direction = "horizontal",
        legend.position = c(1,1),
        legend.justification = c("right","top")
  )+
  labs(x ="Agronomic traits")+
  geom_text(data=df_mean,
            aes(x=TraitName,y=0.45,
                label=sprintf("%.2f",accu_diff)),
            vjust=1,
            position = position_dodge(width = 0.8))
  
dim(resPlot_bxp)
length(unique(resPlot_bxp$EnvName[resPlot_bxp$TraitName==unique(resPlot_bxp$TraitName)[1]]))
length(unique(resPlot_bxp$EnvName[resPlot_bxp$TraitName==unique(resPlot_bxp$TraitName)[2]]))
length(unique(resPlot_bxp$EnvName[resPlot_bxp$TraitName==unique(resPlot_bxp$TraitName)[3]]))

(figure_caption=sprintf('%s%s_Figure02Suppl_MegaLMM_OriginalModel_boxplot_accudiff',DIR_Output,Date))
pdf(paste0(figure_caption,".pdf"),width = 9,height = 6)
print(plot_bxp_suppl)
dev.off()

# -------------------------------------------------------------
# b)	Boxplots - suppl: CV2,CV1 using old MegaLMM
plot_bxp=
  ggplot(resPlot[resPlot$Method=="MegaLMM",], 
         aes(x = TraitName, y = accu_r,fill=CVScn)) +
  labs(y = expression(paste("cor(y, ", hat(g), ")")))+
  # labs(y=expression(r[MegaLMM] - r[GBLUP]))+
  geom_boxplot(outlier.size = 0.7)+
  theme_bw()
df_CV1CV2_accu_mean = aggregate(accu_r~TraitName+CVScn,
                                data = resPlot[resPlot$Method=="MegaLMM",],
                                FUN = mean)


# -------------------------------------------------------------
#  c)	Scatterplot: CV2 vs GBLUP using old MegaLMM for 1 trait
res_scatterplot = reshape2::dcast(resPlot[resPlot$CVScn=="CV2"&
                                                    resPlot$TraitName=="Grain Yield",], 
                                  TraitName + CVScn + EnvName  ~ Method, value.var ="accu_r")
res_scatterplot=left_join(res_scatterplot,df_NumObs_Env[df_NumObs_Env$TraitName=="Grain Yield",],
          by=c("TraitName","EnvName"))
colnames(res_scatterplot)[colnames(res_scatterplot)=="NumObs_Env"]="SampleSize"

plot_sct=
ggplot(res_scatterplot, aes(x = GBLUP, y = MegaLMM, color = SampleSize)) +
  geom_point() +
  geom_abline(slope=1, intercept=0) +
  scale_color_gradient(low = "blue", high = "red") +  # Adjust colors as needed
  # labs(x = "GBLUP",y = "MegaLMM",color = "Sample Size") +
  theme_bw()
  # guides(fill=guide_legend(title="New Legend Title"))
  # theme(#axis.title.x = element_blank(),
  #       # axis.title.y = element_text(size=axis.title.size),
  #       # axis.title= element_text(size=axis.title.size),
  #       legend.title = element_text(angle = 90,size=axis.title.size-2),#element_blank(),
  #       legend.background=element_blank(),
  #       # legend.position = c(legend_posi_x,legend_posi_y),
  #       legend.justification = c("right","bottom"))+
  # guides(fill = guide_colorbar(title.position = "left",name=""))
plot_sct

# -------------------------------------------------------------
#  d)	Line plot: sd(factor_loadings)~factor index for each model 
#     - showing ~# important factor per model
#     - one line = one run of one trait
head(res_Lambda_Squared)
tail(res_Lambda_Squared)
dim(res_Lambda_Squared)
unique(res_Lambda_Squared$TraitName)
res_Lambda_Squared[res_Lambda_Squared$TraitName=="Grain Yield",]
# Define colors for each group
group_colors <- hue_pal()(3)
group_colors <- c("blue", "red", "green")
group_colors <- c("blue", alpha("red", 0.5),"green")
legend_line_width=1
# legend_posi_x=0.98; legend_posi_y=0.98; 
plot_line=
  ggplot(res_Lambda_Squared,
         aes(x=Lambda_ID,
             y=Lambda_samples_squared_mean_max,
             color=TraitName,
             group=paste(TraitName,RunID,sep = "::")))+
  # facet_wrap(~TraitName)+
  # labs(x = "Lambda ID",
  #      y = "RMS")+
  geom_line(lwd=0.3)+
  theme_bw()
  # theme(#axis.title= element_text(size=axis.title.size),
  #       legend.title = element_blank(),
  #       legend.background=element_blank(),
  #       # legend.position = c(legend_posi_x,legend_posi_y),
  #       legend.justification = c("right","top"),
  #       panel.spacing = unit(0, "points"), strip.text = element_blank())+
  #   scale_color_manual(values = group_colors)+
  # guides(color = guide_legend(override.aes = list(lwd = legend_line_width)))
plot_line

# How many factors with lambda^2>0.01 
res_Lambda_Squared_filtered = res_Lambda_Squared[res_Lambda_Squared$Lambda_samples_squared_mean_max>0.01,]
length(unique(res_Lambda_Squared_filtered$Lambda_ID[res_Lambda_Squared_filtered$TraitName=="Grain Yield"]))
# 29
length(unique(res_Lambda_Squared_filtered$Lambda_ID[res_Lambda_Squared_filtered$TraitName=="Plant Height"]))
# 19
length(unique(res_Lambda_Squared_filtered$Lambda_ID[res_Lambda_Squared_filtered$TraitName=="Silk Days"]))
# 20

length(res_Lambda_Squared_filtered$Lambda_ID[res_Lambda_Squared_filtered$TraitName=="Grain Yield"])/10
# 20.3
length(res_Lambda_Squared_filtered$Lambda_ID[res_Lambda_Squared_filtered$TraitName=="Plant Height"])/10
# 13.2
length(res_Lambda_Squared_filtered$Lambda_ID[res_Lambda_Squared_filtered$TraitName=="Silk Days"])/10
# 13

table(res_Lambda_Squared_filtered$Lambda_ID[res_Lambda_Squared_filtered$TraitName=="Grain Yield"])
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 29 31 39 
# 10 10 10 10 10 10 10 10 10 10 10 10 10  7  9  8  9  8  8  2  5  5  4  1  2  1  1  2  1 
table(res_Lambda_Squared_filtered$Lambda_ID[res_Lambda_Squared_filtered$TraitName=="Plant Height"])
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 22 24 
# 10 10 10 10 10 10 10 10 10 10  5  6  4  6  4  3  2  1  1 
table(res_Lambda_Squared_filtered$Lambda_ID[res_Lambda_Squared_filtered$TraitName=="Silk Days"])
# 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 19 23 28 
# 10 10 10 10 10 10 10  9  8  7  8  5  7  4  4  1  3  3  1  1 

# ------------------------------------------------------------------------------
#  e)	-log10(p-values) (y-axis) of regressions of factor loadings (x-axis) on 
#      latitude, longitude, state, tester
# DER: 
#   - don’t average them across runs
#   - show one run only
#   - show one trait only

# convert original pvalue to -log10(p-values) score
res_lm$Score=-log10(res_lm$Pvalue)
dim(res_lm)
head(res_lm)
str(res_lm)
tmp=res_lm[res_lm$TraitName=="YieldMgha"&res_lm$RunID==1,]
dim(tmp)
head(tmp)
table(tmp$Predictor[tmp$Score> (-log10(0.05/50))])
table(tmp$Predictor[tmp$Score> (-log10(0.05))])
# table(tmp$Predictor[tmp$Score> (-log10(0.1/50))])

plot_regpval_GrainYield=
ggplot(data=res_lm[res_lm$TraitName=="YieldMgha"&res_lm$RunID==1,],
       aes(x=as.factor(ID_lambda),y=Score))+
  geom_point(aes(color=Predictor),cex=0.8)+
  # labs(x = expression(paste(Lambda," Index")),
  #      y = expression("-log"[10]*"(p-values)")) +
  scale_x_discrete(breaks = c("1", "10", "20", "30", "40", "50")) +
  theme_bw()
  # theme(#axis.title= element_text(size=axis.title.size),
  #       legend.title =element_blank(),
  #       legend.key=element_blank(),
  #       legend.background=element_blank(),
  #       # legend.position = c(legend_posi_x,legend_posi_y),
  #       legend.justification = c("right","top")) +
  # guides(color = guide_legend(override.aes = list(cex=1)))

plot_regpval_PlantHeight=
ggplot(data=res_lm[res_lm$TraitName=="PlantHeightcm"&res_lm$RunID==1,],
       aes(x=as.factor(ID_lambda),y=Score))+
  geom_point(aes(color=Predictor),cex=0.8)+
  # labs(x = expression(paste(Lambda," Index")),
  #      y = expression("-log"[10]*"(p-values)")) + 
  scale_x_discrete(breaks = c("1", "10", "20", "30", "40", "50")) +
  theme_bw()+
  theme(#axis.title= element_text(size=axis.title.size),
        legend.title =element_blank(),
        legend.key=element_blank(),
        legend.background=element_blank(),
        # legend.position = c(legend_posi_x,legend_posi_y),
        legend.justification = c("right","top")) 

plot_regpval_SilkDays=
ggplot(data=res_lm[res_lm$TraitName=="SilkDAPdays"&res_lm$RunID==1,],
       aes(x=as.factor(ID_lambda),y=Score))+
  geom_point(aes(color=Predictor),cex=0.8)+
  labs(x = "Factor Loading Index",
       y = expression("-log"[10]*"(p-values)")) + 
  scale_x_discrete(breaks = c("1", "10", "20", "30", "40", "50")) +
  theme_bw()+
  theme(#axis.title= element_text(size=axis.title.size),
        legend.title =element_blank(),
        legend.key=element_blank(),
        legend.background=element_blank(),
        # legend.position = c(legend_posi_x,legend_posi_y),
        legend.justification = c("right","top")) 

plot_regpval=plot_regpval_GrainYield

# ------------------------------------------------------------------------------
#  f)	Boxplots of factor loadings for 1-4 example factors grouped by state or tester
# DER: 1) show one run only; 2) show one trait only

length(res_Lambda_SampleMean[["YieldMgha"]])
sapply(res_Lambda_SampleMean[["YieldMgha"]], dim)
res_Lambda_SampleMean[["YieldMgha"]]

# 
tmp = t(res_Lambda_SampleMean[["YieldMgha"]][[1]])
tmp[1:3,1:6]
tmp=rownames_to_column(.data = data.frame(tmp),var = "Loc_Year_Tester")
tmp=add_column(.data = tmp,
               .after = "Loc_Year_Tester",
               Env=sapply(strsplit(tmp$Loc_Year_Tester,"::"),"[",1),
               State=substr(tmp$Loc_Year_Tester,1,2),
               Tester=sapply(strsplit(tmp$Loc_Year_Tester,"::"),"[",2)
)
head(tmp)
tmp_long = reshape2::melt(tmp,id.vars = c("Loc_Year_Tester","Env","State","Tester"),
                          variable.name = "Lambda_ID",value.name = "Lambda")
tmp_long$Lambda_ID=gsub("X","",tmp_long$Lambda_ID)
head(tmp_long)

# set plot pars
legend_posi_x=0.98; legend_posi_y=0.98

# grouped by state
# plot_group_by_state=
  ggplot(data=tmp_long[tmp_long$Lambda_ID%in%as.character(1),],
       aes(x=State,y=Lambda,fill=Lambda_ID))+
  geom_boxplot(outlier.size = 0.6)+
  labs(
       x = "")+
  theme_bw()
  # theme(#axis.title= element_text(size=axis.title.size),
  #       legend.background=element_blank(),
  #       legend.position = c(legend_posi_x,legend_posi_y),
  #       legend.direction = "horizontal",
  #       legend.justification = c("right","top")
  #       ) 

# grouped by stester
plot_group_by_tester=
  ggplot(data=tmp_long[tmp_long$Lambda_ID%in%as.character(1),],
       aes(x=Tester,y=Lambda))+ #,fill=Lambda_ID
  geom_boxplot(outlier.shape = NA)+
    # geom_dotplot(aes(fill=Tester,color=Tester),binaxis='y', stackdir='center', 
    #              dotsize=0.4,show.legend = F)+
    # geom_jitter(shape=16, position=position_jitter(0.2))
    geom_quasirandom(aes(color=Tester),size = 0.7,show.legend = F)+
    # geom_beeswarm(aes(color=Tester),size = 0.7,priority = "descending",show.legend = F)+
  labs(x = "",y=expression(Lambda[1]))+
  theme_bw()
plot_bxplambda=plot_group_by_tester

# ------------------------------------------------------------------------------
# merge
library(ggpubr)
library(gridExtra)
legend_reduction_factor <- 0.5
axis.title.size=12
(figure_caption=sprintf('%s%s_Figure02_MegaLMM_OriginalModel',DIR_Output,Date))
pdf(paste0(figure_caption,".pdf"),width = 10,height = 6)
ggarrange(
  ggarrange(DatStr_Plot,
            plot_bxp+theme(axis.text.y = element_text(angle =90,hjust=0.5),
                           axis.title.y= element_text(size=axis.title.size),
                           legend.title = element_blank(),
                           legend.background=element_blank(),
                           legend.direction = "horizontal",
                           legend.position = c(1,1),
                           legend.justification = c("right","top")
            )+
              labs(x ="Agronomic traits")+
              geom_text(data=df_CV1CV2_accu_mean,
                        aes(x=TraitName,
                            y=c(rep(-0.10,3),rep(0.03,3)),
                            label=sprintf("%.3f",accu_r)),
                        vjust=0,
                        position = position_dodge(width = 0.7)),
            plot_sct+labs(x = "GBLUP",y = "MegaLMM",color = "Sample Size") +
              theme(axis.text.y = element_text(angle =90,hjust=0.5),
                    legend.title = element_text(size=10,angle = 90),
                    legend.position = c(1,0.02),
                    legend.background=element_blank(),
                    legend.justification = c("right","bottom"),
                    legend.key.width = unit(legend_reduction_factor, "cm"),
                    legend.key.height = unit(0.4, "cm"))+
              guides(color = guide_colorbar(title.position = "left",ticks.colour = NA)),
            labels = toupper(letters[1:3]),nrow=1,
            align = "hv"
  ),
  ggarrange(plot_line+theme(axis.text.y = element_text(angle =90,hjust=0.5),
                         legend.title = element_blank(),
                         legend.background=element_blank(),
                         legend.position = c(1,1),
                         legend.justification = c("right","top"),
                         panel.spacing = unit(0, "points"), strip.text = element_blank())+
              labs(x = expression(paste(Lambda," Index")),y=expression(Lambda^2))+
              scale_color_manual(values = group_colors)+
              guides(color = guide_legend(override.aes = list(lwd = legend_line_width)))
            ,
            plot_regpval + labs(x = expression(paste(Lambda," Index")),
                         y = expression("-log"[10]*"(p-values)")) +
              theme(#axis.title= element_text(size=axis.title.size),
                legend.title =element_blank(),
                legend.key=element_blank(),
                legend.background=element_blank(),
                legend.position = c(1,1),
                legend.justification = c("right","top")) +
              geom_hline(yintercept = -log10(0.05/50),lty=2,lwd=0.4)+
              guides(color = guide_legend(override.aes = list(cex=1)))
            ,
            plot_bxplambda+theme(axis.text.y = element_text(angle =90,hjust=0.5),
                         axis.text.x = element_text(angle =45,hjust = 1),
                         axis.title.x = element_blank()),
            labels = toupper(letters[4:6]),nrow=1
  ),
  nrow = 2,
  align = "v"
)+
theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm"))# margin(t, r, b, l)
dev.off()

