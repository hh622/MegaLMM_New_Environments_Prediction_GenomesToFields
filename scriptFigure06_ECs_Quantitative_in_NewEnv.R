# =============================================
# scriptFigure06_ECs_Quantitative_in_NewEnv
# =============================================

# Clean workspace
rm(list=ls())

library(dplyr)
library(tibble)
library(reshape2)
library(ggplot2)
library(foreach)
library(ggbeeswarm)
source("./custom_functions.R")
source("./meta-analyzing_correlations.R")

DIR_Input="../09_Prediction_Results_Summary/"
DIR_Input_Phedat_Format2="../05_MegaLMM_PhenoDat_Format2/"
DIR_Output="../Figure06_QuantECs_Substitute_QualECs/"
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)

Date="2024-01-30"
TraitNames=rev(c("Silk Days","Plant Height","Grain Yield"))
# CVScns = c("CV3NewTrial","CV3NewState","CV3NewTester","CV4")
CVScns = c("NewState","NewTester")

# ------------------------------------------------------------------------------
# loading and prepare data
# ------------------------------------------------------------------------------

# loading
list.files(DIR_Input)
ls=load(file=sprintf("%s%sG2F_Prediction_Results_Formatted_for_MetaAnalysis_ECs_Improves_GP.RData",
                     DIR_Input,"2024-01-08"))
ls # "res_MegaLMM_cor_within_Scn" "res_MegaLMM_cor_btw_Scn" 
res_within_Scn=res_MegaLMM_cor_within_Scn
head(res_within_Scn)
dim(res_within_Scn) #40205     8

# load Pheno dat Format2
list.files(DIR_Input_Phedat_Format2)
ls=load(sprintf("%s%sMegaLMM_PhenoDat_Format2.RData",
                DIR_Input_Phedat_Format2,"2023-09-15"))
ls #"Phe_Trn_Fmt1Wide" "Phe_Trn_Fmt2Wide" "Phe_Trn_Fmt2Long" "K_info"
lapply(Phe_Trn_Fmt2Wide, function(x) x[1:4,1:6])
sapply(Phe_Trn_Fmt2Wide, dim)
#     Yield_Mg_ha Plant_Height_cm Silk_DAP_days
# [1,]        1702            1699          1663
# [2,]         302             278           231

# get number of non-NAs in each Env
names(Phe_Trn_Fmt2Wide)=TraitNames
TraitName=TraitNames[1]
df_NumObs_Env=foreach(TraitName=TraitNames, .combine = "rbind") %do% {
  tempdat = Phe_Trn_Fmt2Wide[[TraitName]]
  dim(tempdat)
  data.frame(TraitName=TraitName,
             EnvName=colnames(tempdat),
             NumObs_Env=colSums(!is.na(tempdat)))
}
dim(df_NumObs_Env)
head(df_NumObs_Env)

# ------------------------------------------------------------------------------
# Meta analysis within CVScn
# ------------------------------------------------------------------------------
# data preparation and selection
res = res_MegaLMM_cor_within_Scn
head(res)
dim(res) #40205     8

TraitName = unique(res$TraitName)[1]
CVScn = unique(res$CVScn)[1]
PriorPredModel = unique(res$PriorPredModel[res$CVScn==CVScn])[2]
head(res[res$TraitName==TraitName&res$CVScn==CVScn&res$PriorPredModel==PriorPredModel,])

z=1.96
res_meta_within_Scn=NA
res_meta_within_Scn=foreach(TraitName = unique(res$TraitName),.combine = "rbind") %do% {
  foreach(CVScn = unique(res$CVScn), .combine = "rbind") %do%{
    foreach(PriorPredModel = unique(res$PriorPredModel[res$CVScn==CVScn]), .combine = "rbind") %do%{
      res_tmp = res[res$TraitName==TraitName & res$CVScn==CVScn & 
                      res$PriorPredModel==PriorPredModel,]
      #compute average across CV2 Folds
      res_tmp = aggregate(cor_MegaLMM~TraitName+CVScn+PriorPredModel+PriorLvl+PredictionModel+EnvName,
                          data = res_tmp,FUN = mean, na.rm=T)
      dim(res_tmp)
      head(res_tmp)
      tail(df_NumObs_Env)
      # TraitName translation
      res_tmp$TraitName[res_tmp$TraitName=="YieldMgha"]="Grain Yield"
      res_tmp$TraitName[res_tmp$TraitName=="PlantHeightcm"]="Plant Height"
      res_tmp$TraitName[res_tmp$TraitName=="SilkDAPdays"]="Silk Days"
      # add NumObs_Env
      res_tmp = left_join(res_tmp,df_NumObs_Env,by=c("TraitName","EnvName"))
      
      # meta_correlations (note: n should be # Observed values in each Env)
      cor_meta=meta_correlations(r=res_tmp$cor_MegaLMM,n=res_tmp$NumObs_Env)
      cor_mean_meta=as.vector(cor_meta$mean_cor)
      cor_sd_meta=cor_meta$sd
      NumEnv=length(res_tmp$NumObs_Env)
      data.frame(TraitName=unique(res_tmp$TraitName),
                 CVScn=unique(res_tmp$CVScn),
                 PriorLevel=unique(res_tmp$PriorLvl),
                 PredictionModel=unique(res_tmp$PredictionModel),
                 PriorPredModel=unique(res_tmp$PriorPredModel),
                 NumEnv=NumEnv,
                 cor_mean_raw=mean(res_tmp$cor_MegaLMM,na.rm=T),
                 cor_mean_meta=cor_mean_meta,
                 cor_sd_raw=sd(res_tmp$cor_MegaLMM,na.rm=T),
                 cor_sd_meta=cor_sd_meta,
                 # CI = xba +/- z*sd/sqrt(n) # 95% CI for sample mean
                 lower.CL=cor_mean_meta - z*(cor_sd_meta/sqrt(NumEnv)),
                 upper.CL=cor_mean_meta + z*(cor_sd_meta/sqrt(NumEnv))
      )
      
    }# end of foreach(PriorPredModel
  }# end of foreach(CVScn)
}# end of foreach(TraitName)
dim(res_meta_within_Scn) #30 12
head(res_meta_within_Scn)

# ------------------------------------------------------------------------------
# Meta analysis between CVScn
# ------------------------------------------------------------------------------

# --------------------------------------------------------
# create a custom function to summarize meta comp results
# meta_correlations 
# - Let 1 be the observed values, 2 be MegaLMM's estimate and 3 be GBLUP's estimate. 
# - r12 is the correlation between MegaLMM and observed, 
# - r13 is the correlation between GBLUP and observed, and 
# - r23 is the correlation between MegaLMM and GBLUP. 
# - n is the number of observations(note: n should be # Observed values in each Env)
# - output: - the average (mu) 
#           - standard deviation (sd) of the difference r12-r13, 
#           - p-value of the comparison p. 
meta_comp_res = function(r12,r13,r23,n,PriorPredMode1,PriorPredMode2){
  # perform meta analysis
  cor_meta=meta_compare_correlations(r12,r13,r23,n) 
  cor_mean_meta=as.vector(cor_meta$mu)
  cor_sd_meta=cor_meta$sd
  cor_pval_meta=cor_meta$p
  data.frame(TraitName=TraitName,
             CVScn=CVScn,
             NumEnv=length(n),
             PriorPredMode1 = PriorPredMode1,
             PriorPredMode2 = PriorPredMode2,
             cor_mean_raw=mean(r12-r13,na.rm=T),
             cor_mean_meta=cor_mean_meta,
             cor_sd_raw=sd(r12-r13,na.rm=T),
             cor_sd_meta=cor_sd_meta,
             cor_pval_meta=cor_pval_meta,
             # CI = xba +/- z*sd/sqrt(n) # 95% CI for sample mean
             lower.CL=cor_mean_meta - z*(cor_sd_meta/sqrt(NumEnv)),
             upper.CL=cor_mean_meta + z*(cor_sd_meta/sqrt(NumEnv)))
  
}

# perform meta analysis
res_meta_btw_Scn=NA
res=res_MegaLMM_cor_btw_Scn
sapply(res, dim)
(TraitName = unique(res[[1]]$TraitName)[1])
(CVScn = names(res)[2])
res_meta_btw_Scn=foreach(TraitName = unique(res[[1]]$TraitName),.combine = "rbind") %do% {
  foreach(CVScn = names(res), .combine = "rbind") %do%{
    res_tmp = res[[CVScn]]
    res_tmp = res_tmp[res_tmp$TraitName==TraitName & res_tmp$CVScn==CVScn,]
    dim(res_tmp)
    head(res_tmp)
    colnames(res_tmp)
    # TraitName translation
    res_tmp$TraitName[res_tmp$TraitName=="YieldMgha"]="Grain Yield"
    res_tmp$TraitName[res_tmp$TraitName=="PlantHeightcm"]="Plant Height"
    res_tmp$TraitName[res_tmp$TraitName=="SilkDAPdays"]="Silk Days"
    
    if(CVScn=="CV3NewState"){
      # average cor across five CV2 Folds
      res_tmp = aggregate(cbind(CV3NewState_O_O_cor,CV3NewState_ST_O_cor,CV3NewState_O_O_vs_ST_O_cor,
                                CV3NewState_ST_T_cor,CV3NewState_ST_O_vs_ST_T_cor,
                                CV3NewState_ST_T_vs_WT_WT_cor,CV3NewState_WT_WT_cor,
                                CV3NewState_WK_WK_cor,CV3NewState_WT_WT_vs_WK_WK_cor)~TraitName+
                            CVScn+EnvName,
                          data = res_tmp,FUN = mean, na.rm=T)
      # add NumObs_Env
      sapply(list(res_tmp,df_NumObs_Env), dim)
      res_tmp = left_join(res_tmp,df_NumObs_Env,by=c("TraitName","EnvName"))
      
      # - CV3NewState meta analysis between scenarios: 
      #   - compare O::O vs S+T::O
      #   - compare S+T::O vs S+T::T
      #   - compare W+T::W+T vs W+K::W+K
      df_CV3NewState_ST_O_vs_O_O=
      meta_comp_res(r12=res_tmp$CV3NewState_ST_O_cor,
                    r13=res_tmp$CV3NewState_O_O_cor,
                    r23=res_tmp$CV3NewState_O_O_vs_ST_O_cor,
                    n=res_tmp$NumObs_Env,
                    PriorPredMode1="O::O",
                    PriorPredMode2="S+T::O"
                    )
      df_CV3NewState_ST_T_vs_ST_O=
        meta_comp_res(r12=res_tmp$CV3NewState_ST_T_cor,
                      r13=res_tmp$CV3NewState_ST_O_cor,
                      r23=res_tmp$CV3NewState_ST_O_vs_ST_T_cor,
                      n=res_tmp$NumObs_Env,
                      PriorPredMode1="S+T::O",
                      PriorPredMode2="S+T::T")
      df_CV3NewState_WT_WT_vs_ST_T=
        meta_comp_res(r12=res_tmp$CV3NewState_WT_WT_cor,
                      r13=res_tmp$CV3NewState_ST_T_cor,
                      r23=res_tmp$CV3NewState_ST_T_vs_WT_WT_cor,
                      n=res_tmp$NumObs_Env,
                      PriorPredMode1="S+T::T",
                      PriorPredMode2="W+T::W+T")
      df_CV3NewState_WT_WT_vs_WK_WK=
        meta_comp_res(r12=res_tmp$CV3NewState_WT_WT_cor,
                      r13=res_tmp$CV3NewState_WK_WK_cor,
                      r23=res_tmp$CV3NewState_WT_WT_vs_WK_WK_cor,
                      n=res_tmp$NumObs_Env,
                      PriorPredMode1="W+T::W+T",
                      PriorPredMode2="W+K::W+K")
      df = rbind(df_CV3NewState_ST_O_vs_O_O,df_CV3NewState_ST_T_vs_ST_O,
                 df_CV3NewState_WT_WT_vs_ST_T,df_CV3NewState_WT_WT_vs_WK_WK)
    }else if(CVScn=="CV3NewTester"){
      # average cor across CV2 Folds
      res_tmp = aggregate(cbind(CV3NewTester_O_O_cor,CV3NewTester_ST_O_cor,CV3NewTester_O_O_vs_ST_O_cor,
                                CV3NewTester_ST_S_cor,CV3NewTester_ST_O_vs_ST_S_cor,
                                CV3NewTester_SK_SK_cor,CV3NewTester_ST_S_vs_SK_SK_cor,
                                CV3NewTester_WK_WK_cor,CV3NewTester_SK_SK_vs_WK_WK_cor)~TraitName+
                            CVScn+EnvName,
                          data = res_tmp,FUN = mean, na.rm=T)
      # add NumObs_Env
      sapply(list(res_tmp,df_NumObs_Env), dim)
      res_tmp = left_join(res_tmp,df_NumObs_Env,by=c("TraitName","EnvName"))
      
      # - CV3NewTester meta analysis between scenarios: 
      #   - compare O::O vs S+T::O
      #   - compare S+T::O vs S+T::S
      #   - compare S+K::S+K vs W+K::W+K
      df_CV3NewTester_ST_O_vs_O_O=
        meta_comp_res(r12=res_tmp$CV3NewTester_ST_O_cor,
                      r13=res_tmp$CV3NewTester_O_O_cor,
                      r23=res_tmp$CV3NewTester_O_O_vs_ST_O_cor,
                      n=res_tmp$NumObs_Env,
                      PriorPredMode1="O::O",
                      PriorPredMode2="S+T::O")
      df_CV3NewTester_ST_S_vs_ST_O=
        meta_comp_res(r12=res_tmp$CV3NewTester_ST_S_cor,
                      r13=res_tmp$CV3NewTester_ST_O_cor,
                      r23=res_tmp$CV3NewTester_ST_O_vs_ST_S_cor,
                      n=res_tmp$NumObs_Env,
                      PriorPredMode1="S+T::O",
                      PriorPredMode2="S+T::S")
      
      df_CV3NewTester_SK_SK_vs_ST_S=
        meta_comp_res(r12=res_tmp$CV3NewTester_SK_SK_cor,
                      r13=res_tmp$CV3NewTester_ST_S_cor,
                      r23=res_tmp$CV3NewTester_ST_S_vs_SK_SK_cor,
                      n=res_tmp$NumObs_Env,
                      PriorPredMode1="S+T::S",
                      PriorPredMode2="S+K::S+K")
      
      df_CV3NewTester_SK_SK_vs_WK_WK=
        meta_comp_res(r12=res_tmp$CV3NewTester_SK_SK_cor,
                      r13=res_tmp$CV3NewTester_WK_WK_cor,
                      r23=res_tmp$CV3NewTester_SK_SK_vs_WK_WK_cor,
                      n=res_tmp$NumObs_Env,
                      PriorPredMode1="S+K::S+K",
                      PriorPredMode2="W+K::W+K")
      df = rbind(df_CV3NewTester_ST_O_vs_O_O,df_CV3NewTester_ST_S_vs_ST_O,
                 df_CV3NewTester_SK_SK_vs_ST_S,df_CV3NewTester_SK_SK_vs_WK_WK)
    }
    return(df)
  }
}
dim(res_meta_btw_Scn)
head(res_meta_btw_Scn)

# ------------------------------------------------------------------------------
# ggplot
# ------------------------------------------------------------------------------
resPlot = res_meta_within_Scn
head(resPlot)

# update CVScn
resPlot$CVScn = gsub("CV3","",resPlot$CVScn)
# update PriorPredModel
unique(resPlot$PriorPredModel)
resPlot$PriorPredModel[resPlot$PriorPredModel=="Zero::Zero"]="O::O"
resPlot$PriorPredModel[resPlot$PriorPredModel=="State+Tester::Zero"]="S+T::O"
resPlot$PriorPredModel[resPlot$PriorPredModel=="State+Tester::Tester"]="S+T::T"
resPlot$PriorPredModel[resPlot$PriorPredModel=="Wea+Tester::Wea+Tester"]="W+T::W+T"
resPlot$PriorPredModel[resPlot$PriorPredModel=="Wea+TesterGRM::Wea+TesterGRM"]="W+K::W+K"
resPlot$PriorPredModel[resPlot$PriorPredModel=="State+Tester::State"]="S+T::S"
resPlot$PriorPredModel[resPlot$PriorPredModel=="State+TesterGRM::State+TesterGRM"]="S+K::S+K"

# differen colors for different models
library(scales)
show_col(hue_pal()(4))
PriorPredModel_All=c("O::O","S+T::T","W+T::W+T","W+K::W+K","S+T::S","S+K::S+K")
# ------------------------------------------------------------------------------
# plot resPlot_NewState
PriorPredModel_Sel=c("O::O","S+T::T","W+T::W+T","W+K::W+K")
resPlot_NewState = resPlot[resPlot$CVScn=="NewState" &
                             resPlot$PriorPredModel %in% PriorPredModel_Sel,]
# factorization
resPlot_NewState$PriorPredModel=factor(resPlot_NewState$PriorPredModel,
                                       levels = PriorPredModel_All)

# get plot data
plot=ggplot(resPlot_NewState, 
            aes(x = TraitName, y = cor_mean_meta,fill=PriorPredModel)) +
  geom_bar(position="dodge", stat="identity")
(plot_data <- data.frame(ggplot_build(plot)$data))
resPlot_NewState = cbind(resPlot_NewState,plot_data)
head(resPlot_NewState)
plot_data

# Create a df for Adding segments between the first and second boxplots within each group
segments_vertical_len=0.01
segments_padding=c(0.02,0.035, # adjust y posi of individual segments
                   0.02,0.02,
                   0.022,0.02)-0.005

pval_padding=0.01
idx_xstart=c(2:3,2:3+4,2:3+8)
idx_xend  =c(2:3,2:3+4,2:3+8)+1

# segments_vertical_len=0.01
# segments_padding=c(0.02,0.02,0.02, # adjust y posi of individual segments
#                    0.02,0.035,0.02,
#                    0.02,0.022,0.02)-0.005
# 
# pval_padding=0.01
# idx_xstart=c(1:3,1:3+4,1:3+8)
# idx_xend  =c(1:3,1:3+4,1:3+8)+1


df_segments_horiz = data.frame(idx_xstart=idx_xstart,
                               idx_xend=idx_xend,
                               xstart=resPlot_NewState$x[idx_xstart],
                               xend=resPlot_NewState$x[idx_xend],
                               ystart=apply(cbind(resPlot_NewState$upper.CL[idx_xstart],
                                                  resPlot_NewState$upper.CL[idx_xend]),1,max)+
                                 segments_padding, 
                               yend=apply(cbind(resPlot_NewState$upper.CL[idx_xstart],
                                                resPlot_NewState$upper.CL[idx_xend]),1,max)+
                                 segments_padding)
df_segments_left_tip = data.frame(idx_xstart=idx_xstart,
                                  idx_xend=idx_xend,
                                  xstart=resPlot_NewState$x[idx_xstart],
                                  xend=resPlot_NewState$x[idx_xstart],
                                  ystart=apply(cbind(resPlot_NewState$upper.CL[idx_xstart],
                                                     resPlot_NewState$upper.CL[idx_xend]),1,max)+
                                    segments_padding, 
                                  yend=apply(cbind(resPlot_NewState$upper.CL[idx_xstart],
                                                   resPlot_NewState$upper.CL[idx_xend]),1,max)+
                                    segments_padding-segments_vertical_len)
df_segments_right_tip = data.frame(idx_xstart=idx_xstart,
                                   idx_xend=idx_xend,
                                   xstart=resPlot_NewState$x[idx_xend],
                                   xend=resPlot_NewState$x[idx_xend],
                                   ystart=apply(cbind(resPlot_NewState$upper.CL[idx_xstart],
                                                      resPlot_NewState$upper.CL[idx_xend]),1,max)+
                                     segments_padding, 
                                   yend=apply(cbind(resPlot_NewState$upper.CL[idx_xstart],
                                                    resPlot_NewState$upper.CL[idx_xend]),1,max)+
                                     segments_padding-segments_vertical_len)

# Create a df for Adding p value symbols
head(resPlot_NewState)
head(res_meta_btw_Scn)
head(df_segments_horiz)
df_pval = res_meta_btw_Scn[res_meta_btw_Scn$CVScn=="CV3NewState"&
                             res_meta_btw_Scn$PriorPredMode1 %in%PriorPredModel_Sel[-1],]
df_pval$label=ifelse(df_pval$cor_pval_meta  < 0.01,"**",
                     ifelse(df_pval$cor_pval_meta< 0.05, "*","ns"))
df_pval$x=(df_segments_horiz$xstart+df_segments_horiz$xend)/2
df_pval$y= df_segments_horiz$yend + pval_padding
head(df_pval)
# df_pval$color = ifelse(is.na(df_pval$pval_z),"white",ifelse(df_pval$pval_z < 0.05,"blue","black"))

# plot pars
errorbarwidth=as.numeric(resPlot_NewState$xmax-resPlot_NewState$xmin)[1]/5
padding=0.02
head(resPlot_NewState)
library(scales)
plot_NewState=
ggplot(resPlot_NewState, 
       aes(x = TraitName, y = cor_mean_meta,
           fill=factor(PriorPredModel,levels = PriorPredModel_All))) +
  labs(y = "cor")+
  geom_bar(position="dodge2", stat="identity")+
  # scale_fill_manual(values = hue_pal()(6)[1:3])+
  scale_fill_manual(#name = 'PriorPredModel', 
    values = c("O::O"=hue_pal()(6)[1],
               "S+T::T"=hue_pal()(6)[2],
               "W+T::W+T"=hue_pal()(6)[3],
               "W+K::W+K"=hue_pal()(6)[4]
               # "S+T::S"=hue_pal()(6)[5],
               # "S+K::S+K"=hue_pal()(6)[6]
               ),
    drop = FALSE)+
  labs(y = expression(paste("cor(y, ", hat(g), ") ")))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(angle=0),
        # panel.grid = element_blank(),
        legend.title = element_blank(),
        # legend.direction ="horizontal",
        legend.direction = "vertical",
        # legend.box = "horizontal",
        # legend.position = c(legend.position.x,legend.position.y),
        legend.position = "right",
        legend.justification = "left"
        # legend.justification = c("right","bottom")
  )+
  geom_segment(#data = plot_data,
    aes(x=x,
        xend=x,
        y=lower.CL,
        yend=upper.CL
        # y=cor_mean_meta-cor_sd_meta,
        # yend=cor_mean_meta+cor_sd_meta
    )
  )+
  geom_segment(
    aes(x=x-errorbarwidth/2,
        xend=x+errorbarwidth/2,
        y=lower.CL,
        yend=lower.CL
        # y=cor_mean_meta-cor_sd_meta,
        # yend=cor_mean_meta-cor_sd_meta
    )
  )+
  geom_segment(
    aes(x=x-errorbarwidth/2,
        xend=x+errorbarwidth/2,
        y=upper.CL,
        yend=upper.CL
        # y=cor_mean_meta+cor_sd_meta,
        # yend=cor_mean_meta+cor_sd_meta
    )
  )+
  annotate(geom="text", x=resPlot_NewState$x,
           y=0-padding,
           label=sprintf("%.3f",resPlot_NewState$cor_mean_meta))+
  geom_segment(data = df_segments_horiz,
               aes(x=xstart,
                   xend=xend,
                   y=ystart,
                   yend=yend),
               inherit.aes = F
  )+
  geom_segment(data = df_segments_left_tip,
               aes(x=xstart,
                   xend=xend,
                   y=ystart,
                   yend=yend),
               inherit.aes = F
  )+
  geom_segment(data = df_segments_right_tip,
               aes(x=xstart,
                   xend=xend,
                   y=ystart,
                   yend=yend),
               inherit.aes = F
  )+
  geom_text(data = df_pval,
            aes(x=x,
                y=y,
                label = label
            ),inherit.aes = F,vjust=0.2
            )

# ------------------------------------------------------------------------------
# plot resPlot_NewTester
PriorPredModel_Sel=c("O::O","S+T::S","S+K::S+K","W+K::W+K")
resPlot_NewTester = resPlot[resPlot$CVScn=="NewTester" &
                              resPlot$PriorPredModel %in% PriorPredModel_Sel,]

# factorization
resPlot_NewTester$PriorPredModel=factor(resPlot_NewTester$PriorPredModel,
                                        levels = c("O::O","S+T::T","W+T::W+T",
                                                   "S+T::S","S+K::S+K","W+K::W+K"))
# get plot data
plot=ggplot(resPlot_NewTester, 
            aes(x = TraitName, y = cor_mean_meta,fill=PriorPredModel)) +
  geom_bar(position="dodge", stat="identity")
(plot_data <- data.frame(ggplot_build(plot)$data))
resPlot_NewTester = cbind(resPlot_NewTester,plot_data)
head(resPlot_NewTester)
plot_data

# Create a df for Adding segments between the first and second boxplots within each group
# segments_padding=c(0.02,0.03,0.035,0.02,
#                    0.02,0.04,0.035,0.02,
#                    0.035,0.02,0.04,0.02)-0.005
segments_vertical_len=0.01
pval_padding=0.01
# idx_xstart=c(1:4,1:4+5,1:4+10)
# idx_xend  =c(1:4,1:4+5,1:4+10)+1
# idx_ystart=c(1,2,4,c(1,2,4)+5,c(1,2,4)+10)
segments_padding=c(0.02,0.035, # adjust y posi of individual segments
                   0.02,0.035,
                   0.02,0.037)-0.005
pval_padding=0.01
# idx_xstart=c(1:2,1:2+3,1:2+6)
# idx_xend  =c(1:2,1:2+3,1:2+6)+1
idx_xstart=c(2:3,2:3+4,2:3+8)
idx_xend  =c(2:3,2:3+4,2:3+8)+1


head(resPlot_NewTester)
df_segments_horiz = data.frame(xstart=resPlot_NewTester$x[idx_xstart],
                         xend=resPlot_NewTester$x[idx_xend],
                         ystart=apply(cbind(resPlot_NewTester$upper.CL[idx_xstart],
                                            resPlot_NewTester$upper.CL[idx_xend]),1,max)+
                           segments_padding, 
                         yend=apply(cbind(resPlot_NewTester$upper.CL[idx_xstart],
                                          resPlot_NewTester$upper.CL[idx_xend]),1,max)+
                           segments_padding
)


df_segments_left_tip = data.frame(xstart=resPlot_NewTester$x[idx_xstart],
                                  xend=resPlot_NewTester$x[idx_xstart],
                                  ystart=apply(cbind(resPlot_NewTester$upper.CL[idx_xstart],
                                                     resPlot_NewTester$upper.CL[idx_xend]),1,max)+
                                    segments_padding, 
                                  yend=apply(cbind(resPlot_NewTester$upper.CL[idx_xstart],
                                                   resPlot_NewTester$upper.CL[idx_xend]),1,max)+
                                    segments_padding-segments_vertical_len
                                  )
df_segments_right_tip = data.frame(xstart=resPlot_NewTester$x[idx_xend],
                                   xend=resPlot_NewTester$x[idx_xend],
                                   ystart=apply(cbind(resPlot_NewTester$upper.CL[idx_xstart],
                                                      resPlot_NewTester$upper.CL[idx_xend]),1,max)+
                                     segments_padding, 
                                   yend=apply(cbind(resPlot_NewTester$upper.CL[idx_xstart],
                                                    resPlot_NewTester$upper.CL[idx_xend]),1,max)+
                                     segments_padding-segments_vertical_len
                                   )
# Create a df for Adding p value symbols
head(resPlot_NewTester)
df_pval = res_meta_btw_Scn[res_meta_btw_Scn$CVScn=="CV3NewTester" &
                             res_meta_btw_Scn$PriorPredMode1 %in% PriorPredModel_Sel[-1],]
df_pval$label=ifelse(df_pval$cor_pval_meta  < 0.01,"**",
                     ifelse(df_pval$cor_pval_meta< 0.05, "*","ns"))
df_pval$x=(df_segments_horiz$xstart+df_segments_horiz$xend)/2
df_pval$y= df_segments_horiz$yend + pval_padding
head(df_pval)

# plot pars
errorbarwidth=as.numeric(resPlot_NewTester$xmax-resPlot_NewTester$xmin)[1]/5
padding=0.02
head(resPlot_NewTester)
plot_NewTester=
ggplot(resPlot_NewTester, 
       aes(x = TraitName, y = cor_mean_meta,fill=PriorPredModel)) +
  labs(y = "cor")+
  geom_bar(position="dodge2", stat="identity")+
  # scale_fill_manual(values = hue_pal()(6)[c(4,5,3)])+
  scale_fill_manual(#name = 'PriorPredModel', 
    values = c("O::O"=hue_pal()(6)[1],
               # "S+T::T"=hue_pal()(6)[2],
               # "W+T::W+T"=hue_pal()(6)[3],
               "W+K::W+K"=hue_pal()(6)[4],
               "S+T::S"=hue_pal()(6)[5],
               "S+K::S+K"=hue_pal()(6)[6]
    ),
    drop = FALSE)+
  labs(y = expression(paste("cor(y, ", hat(g), ") ")))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.x = element_text(angle=0),
        # panel.grid = element_blank(),
        legend.title = element_blank(),
        # legend.direction ="horizontal",
        legend.direction = "vertical",
        # legend.box = "horizontal",
        # legend.position = c(legend.position.x,legend.position.y),
        legend.position = "right",
        legend.justification = "left"
        # legend.justification = c("right","bottom")
  )+
  geom_segment(#data = plot_data,
    aes(x=x,
        xend=x,
        y=lower.CL,
        yend=upper.CL
        # y=cor_mean_meta-cor_sd_meta,
        # yend=cor_mean_meta+cor_sd_meta
    )
  )+
  geom_segment(
    aes(x=x-errorbarwidth/2,
        xend=x+errorbarwidth/2,
        y=lower.CL,
        yend=lower.CL
        # y=cor_mean_meta-cor_sd_meta,
        # yend=cor_mean_meta-cor_sd_meta
    )
  )+
  geom_segment(
    aes(x=x-errorbarwidth/2,
        xend=x+errorbarwidth/2,
        y=upper.CL,
        yend=upper.CL
        # y=cor_mean_meta+cor_sd_meta,
        # yend=cor_mean_meta+cor_sd_meta
    )
  )+
  annotate(geom="text", x=resPlot_NewTester$x,
           y=0-padding,
           label=sprintf("%.3f",resPlot_NewTester$cor_mean_meta))+
  geom_segment(data = df_segments_horiz,
               aes(x=xstart,
                   xend=xend,
                   y=ystart,
                   yend=yend),
               inherit.aes = F
  )+
  geom_segment(data = df_segments_left_tip,
               aes(x=xstart,
                   xend=xend,
                   y=ystart,
                   yend=yend),
               inherit.aes = F
  )+
  geom_segment(data = df_segments_right_tip,
               aes(x=xstart,
                   xend=xend,
                   y=ystart,
                   yend=yend),
               inherit.aes = F
  )+
  geom_text(data = df_pval,
            aes(x=x,
                y=y,
                label = label
            ),inherit.aes = F,vjust=0.2
  )

library(ggpubr)

# output in pdf
figure=ggpubr::ggarrange(plot_NewState+rremove("y.title")+
                           annotate(geom="text", x=as.numeric(plot_data$xmin[1]), 
                                    y=0.60, label="NewState",color="black",cex=4.5,
                                    hjust=0,vjust=0)+
                           theme(legend.justification = "bottom",
                                 # legend.direction="horizontal"
                                 ),
                           # guides(fill = guide_legend(nrow = 1)),
                         plot_NewTester+rremove("y.title")+
                           annotate(geom="text", x=as.numeric(plot_data$xmin[1]), y=0.48, 
                                    label="NewTester",color="black",cex=4.5,
                                    hjust=0,vjust=0)+
                           theme(legend.justification = "bottom",
                                 # legend.direction="horizontal"
                           ),
                         # guides(fill = guide_legend(nrow = 1)),
                         # common.legend = TRUE,legend = c("bottom"),
                         label.x = -0.00, #the number represents percentage of total x length(+move right; -, left)
                         hjust=0.7,
                         labels = c("A", "B"),
                         widths = c(1, 1),
                         nrow = 2,ncol = 1)+
  theme(plot.margin = margin(0.1,0.1,0.1,0.2, "cm"))# margin(t, r, b, l)

(figure_caption=sprintf("%s%s_Figure06_ECs_Quantitative_in_NewEnv.pdf",DIR_Output,Date))
pdf(figure_caption,width = 8, height = 6)
print(
  annotate_figure(figure,
                  # left = expression(paste("cor(y, ", hat(g), ")"))
                  left = text_grob(expression(paste("cor(y, ", hat(g), ")")),
                                   size=12,face="bold",rot = 90)
                  # top = trait
  ))

# ggpubr ggarrange legend on right side legend width difference results in plot width difference

dev.off()

