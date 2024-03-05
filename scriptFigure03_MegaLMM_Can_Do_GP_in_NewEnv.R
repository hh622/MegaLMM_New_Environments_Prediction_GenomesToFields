# =====================================================
# scriptFigure03_MegaLMM_Can_Do_GP_in_NewEnv
# =====================================================

# Clean workspace
rm(list=ls())

library(reshape2)
library(ggplot2)
library(foreach)
library(tibble)
library(ggbeeswarm) #plot outliers in a line
source("./custom_functions.R")
source("./meta-analyzing_correlations.R")

DIR_Input="../09_Prediction_Results_Summary/"
DIR_Output="../Figure03_MegaLMM_Works_NewEnv/"
DIR_Input_Phedat_Format2="../05_MegaLMM_PhenoDat_Format2/"
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)

Date="2023-12-31"
TraitSel=rev(c("Silk Days","Plant Height","Grain Yield"))
# CVScns = c("CV1","CV2","CV3NewTrial","CV3NewState","CV3NewTester","CV4")
CVScns = c("CV1","CV2","NewTrial","NewState","NewTester","NewGenoNewYear")
# Ref="CV2Original"
Ref="CV2S+T"
# ------------------------------------------------------------------------------
# loading and prepare data
# ------------------------------------------------------------------------------

# Loading Formatted Prediction Results 
list.files(DIR_Input,".RData")
ls=load(file=sprintf("%s%sG2F_Prediction_Results_Formatted.RData",DIR_Input,"2023-10-23"))
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
sapply(Phe_Trn_Fmt2Wide, dim)
#     Yield_Mg_ha Plant_Height_cm Silk_DAP_days
# [1,]        1702            1699          1663
# [2,]         302             278           231

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

# ------------------------------------------------------------------------------
# ggplot - accu_r based on trail/env values - Prior + specific prediction
# ------------------------------------------------------------------------------

# data preparation and selection
res = res_EnvMean
head(res)
dim(res) #72812    12

# add NumObs_Env
library(dplyr)
res=left_join(res,df_NumObs_Env,by=c("TraitName","EnvName"))
res=res[,colnames(res)[c(1:(ncol(res)-3),ncol(res),ncol(res)-2, ncol(res)-1)]]

# add other columns
res=add_column(.data = res,
                  .after = "PriorLevel",
                  CV_Prior_Prediction=paste(res$CVScn,res$PriorLevel,res$PredictionModel,sep = "::"))
res=add_column(.data = res,
                  .after = "CV_Prior_Prediction",
                  Trait_Method_CV_Prior_Prediction=paste(
                    res$TraitName,res$Method,
                    res$CVScn,res$PriorLevel,res$PredictionModel,sep = "::"))

label = unique(res$Trait_Method_CV_Prior_Prediction)[1]
z=1.96
res_meta=NA
res_meta=foreach(label = unique(res$Trait_Method_CV_Prior_Prediction),.combine = "rbind") %do% {
  res_tmp = res[res$Trait_Method_CV_Prior_Prediction==label,]
  # meta_correlations (note: n should be # Observed values in each Env)
  r=res_tmp$accu_r; NumObs_Env=res_tmp$NumObs_Env
  cor_meta=meta_correlations(r=r,n=NumObs_Env)
  cor_mean_meta=as.vector(cor_meta$mean_cor)
  cor_sd_meta=cor_meta$sd
  NumEnv=length(NumObs_Env)
  data.frame(TraitName=unique(res_tmp$TraitName),
             CVScn=unique(res_tmp$CVScn),
             Method=unique(res_tmp$Method),
             PriorType=unique(res_tmp$PriorType),
             PriorLevel=unique(res_tmp$PriorLevel),
             PredictionType=unique(res_tmp$PredictionType),
             PredictionModel=unique(res_tmp$PredictionModel),
             CV_Prior_Prediction=unique(res_tmp$CV_Prior_Prediction),
             NumEnv=NumEnv,
             cor_mean_raw=mean(res_tmp$accu_r,na.rm=T),
             cor_mean_meta=cor_mean_meta,
             cor_sd_raw=sd(res_tmp$accu_r,na.rm=T),
             cor_sd_meta=cor_sd_meta,
             # CI = xba +/- z*sd/sqrt(n) # 95% CI for sample mean
             lower.CL=cor_mean_meta - z*(cor_sd_meta/sqrt(NumEnv)),
             upper.CL=cor_mean_meta + z*(cor_sd_meta/sqrt(NumEnv))
             )
  
  
}
dim(res_meta) #270  12
head(res_meta)
unique(res_meta$CV_Prior_Prediction[res_meta$CVScn=="CV2"])
plot(res_meta$cor_mean_meta,res_meta$cor_mean_raw);abline(0,1)
plot(res_meta$cor_sd_meta,res_meta$cor_sd_raw);abline(0,1)

# -------------------------------------------------------------------
# Figure plot with ggplot

Figure="Main"
# Figure="Suppl"
resPlot = res_meta
# subset scenarios for res_meta
if(Figure=="Main"){
  resPlot= resPlot[resPlot$CVScn!="CV1" &
                     resPlot$CVScn!="CV2" &
                   resPlot$Method=="MegaLMM" &
                   resPlot$PriorType=="Yes" &
                   resPlot$PredictionType=="Specific" &
                   resPlot$PriorLevel=="S+T"&
                   resPlot$PredictionModel=="S+T",]
}else if(Figure=="Suppl"){
  resPlot= resPlot[resPlot$CVScn!="CV1" &
                     resPlot$CVScn!="CV2" &
                   resPlot$Method=="MegaLMM" &
                   ( 
                     # (resPlot$CVScn =="CV2" & resPlot$PriorLevel=="O" &
                     #    resPlot$PredictionModel=="Env")|
                       (resPlot$CVScn =="CV3NewTrial" & resPlot$PriorLevel=="S+T" &
                          resPlot$PredictionModel=="S+T") |
                       (resPlot$CVScn=="CV3NewState" & resPlot$PriorLevel=="S+T" &
                          resPlot$PredictionModel=="S+T") |
                       (resPlot$CVScn=="CV3NewTester" & resPlot$PriorLevel=="S+T"&
                          resPlot$PredictionModel=="S+T") |
                       (resPlot$CVScn=="CV4" & resPlot$PriorLevel=="S+T" &
                          resPlot$PredictionModel=="S+T")
                   ),]
}

stopifnot(ncol(resPlot)==15)
dim(resPlot)
head(resPlot)
unique(resPlot$Method)
unique(resPlot$CV_Prior_Prediction)

# modify CVScn names
resPlot$CVScn = gsub("CV3","",resPlot$CVScn)
resPlot$CVScn[resPlot$CVScn=="CV4"] = "NewGenoNewYear"
CVScns = c("NewTrial","NewState","NewTester","NewGenoNewYear")
unique(resPlot$CVScn)

# ggplot
# define orders
resPlot$CVScn = factor(resPlot$CVScn,levels=CVScns)

# get plot data
plot=ggplot(resPlot, 
            aes(x = TraitName, y = cor_mean_meta,fill=CVScn)) +
  geom_bar(position="dodge", stat="identity")
(plot_data <- data.frame(ggplot_build(plot)$data))
resPlot = cbind(resPlot,plot_data)
head(resPlot)
unique(resPlot$CV_Prior_Prediction)

# plot pars
errorbarwidth=as.numeric(plot_data$xmax-plot_data$xmin)[1]/4
padding=0.02
Ref="NULL"

(figure_caption=sprintf('%s%s_G2F_Figure03_MegaLMM_Works_NewEnv_%s_Ref%s',
                        DIR_Output,Date,Figure,Ref))
pdf(paste0(figure_caption,".pdf"),width = 10,height = 5)
ggplot(resPlot, 
       aes(x = TraitName, y = cor_mean_meta,fill=CVScn)) +
       # aes(x = TraitName, y = cor_mean_meta,color=CVScn)) +
  labs(y = expression(paste("cor(y, ", hat(g), ")")))+
  geom_bar(position="dodge2", stat="identity")+
  # geom_point(position=position_dodge(width = 0.9))+
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
  annotate(geom="text", x=resPlot$x,
           # y=max(resPlot$upper.CL)+padding,
           y=0-padding,
           label=sprintf("%.3f",resPlot$cor_mean_meta))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        legend.justification = "right") +
  guides(fill = guide_legend(nrow = 1))
  # guides(color = guide_legend(nrow = 1))
dev.off()
