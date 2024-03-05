# ==============================================================================
# scriptFigure04_GBLUP_Specific_vs_General
# ==============================================================================

# Clean workspace
rm(list=ls())

library(reshape2)
library(ggplot2)
library(tibble)
library(foreach)
library(ggbeeswarm)
library(dplyr)

DIR_Input="../09_Prediction_Results_Summary/"
DIR_Output="../Figure04_MegaLMM_Versus_GBLUP/"
DIR_Input_Phedat_Format2="../05_MegaLMM_PhenoDat_Format2/"
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)

Date="2024-01-08"
source("./custom_functions.R")
source("./meta-analyzing_correlations.R")
TraitNames = rev(c("Silk Days","Plant Height","Grain Yield"))
CVScns = c("CV1","CV2","CV3NewTrial","CV3NewState","CV3NewTester","CV4")
# MegaLMM_Models_Sel="S+T"
# GBLUP_Models_Sel="S+T"
MegaLMM_Models_Sel="Best"
GBLUP_Models_Sel="Best"

# ------------------------------------------------------------------------------
# loading and prepare data
# ------------------------------------------------------------------------------

# Loading Formatted Prediction Results 
list.files(DIR_Input,pattern = "RData")
ls=load(file=sprintf("%s%sG2F_Prediction_Results_Formatted_MegaLMM_%s_vs_GBLUP_%s.RData",
                     DIR_Input,"2024-01-08",MegaLMM_Models_Sel,GBLUP_Models_Sel))
ls #"res_Original" "res_EnvMean" 
res=res_EnvMean
head(res)
dim(res)
(302+278+231)

# checking with original results
tmp = res[res$CVScn=="CV4",]
tmp = melt(data = tmp,id.vars = c("TraitName", "CVScn", "TestSetID","EnvName"),
           measure.vars = c("cor_MegaLMM","cor_GBLUP","cor_MegaLMMGBLUP"),
           variable.name = "Category",
           value.name = "accu")
head(tmp)
ggplot(data=tmp,
       aes(x=TraitName,y=accu,fill=Category))+
  geom_boxplot()+
  stat_summary(
    fun = mean, 
    geom = "point", 
    shape = 4,
    size = 3,
    position = position_dodge(width = 0.75),
    color = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "text", 
    aes(label = round(..y.., digits = 3)),  # Use round(..y.., digits = 2) to round mean values
    vjust = -0.5,  # Adjust vertical position of text
    position = position_dodge(width = 0.75),
    color = "black"
  )
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
# perform meta analysis
head(res)
dim(res) #4850   10

TraitName = unique(res$TraitName)[1]; CVScn = unique(res$CVScn)[1]
TraitName="Grain Yield";CVScn="CV4"
z=1.96
CVScn = "CV3NewTester"
res_meta=NA
res_meta=foreach(TraitName = unique(res$TraitName),.combine = "rbind") %do% {
  foreach(CVScn = unique(res$CVScn), .combine = "rbind") %do%{
    res_tmp = res[res$TraitName==TraitName & res$CVScn==CVScn,]
    dim(res_tmp)
    head(res_tmp)
    head(df_NumObs_Env)
    unique(res_tmp$TraitName)
    unique(res_tmp$CVScn)
    unique(res_tmp$TestSetID)
    
    # checking data size: for each scenario should with less than 302 EnvName
    stopifnot(nrow(res_tmp)<=302)
    
    res_tmp = left_join(res_tmp,df_NumObs_Env,by=c("TraitName","EnvName"))
    
    # meta_correlations 
    # - Let 1 be the observed values, 2 be MegaLMM's estimate and 3 be GBLUP's estimate. 
    # - r12 is the correlation between MegaLMM and observed, 
    # - r13 is the correlation between GBLUP and observed, and 
    # - r23 is the correlation between MegaLMM and GBLUP. 
    # - n is the number of observations(note: n should be # Observed values in each Env)
    # - output: - the average (mu) 
    #           - standard deviation (sd) of the difference r12-r13, 
    #           - p-value of the comparison p. 
    r12=res_tmp$cor_MegaLMM
    r13=res_tmp$cor_GBLUP
    r23=res_tmp$cor_MegaLMMGBLUP
    NumObs_Env=res_tmp$NumObs_Env
    sapply(list(r12,r13,r23), function(x) sum(is.na(x)))
    sapply(list(r12,r13,r23),length)
    
    # excluding outliers
    t = r.test(n=NumObs_Env,r12=r12,r13=r13,r23=r23)$t
    if(sum(is.nan(t))>0 & length(t)==nrow(res_tmp)){
      res_tmp=res_tmp[!is.nan(t),]
      r12=res_tmp$cor_MegaLMM
      r13=res_tmp$cor_GBLUP
      r23=res_tmp$cor_MegaLMMGBLUP
      NumObs_Env=res_tmp$NumObs_Env
    }
    
    # perform meta analysis
    cor_meta=meta_compare_correlations(r12,r13,r23,n=NumObs_Env) 
    cor_mean_meta=as.vector(cor_meta$mu)
    cor_sd_meta=cor_meta$sd
    cor_pval_meta=cor_meta$p
    NumEnv=length(NumObs_Env)
    head(res_tmp)
    data.frame(TraitName=TraitName,
               CVScn=CVScn,
               Model_MegaLMM=unique(res_tmp$Model_MegaLMM),
               Model_GBLUP=unique(res_tmp$Model_GBLUP),
               NumEnv=NumEnv,
               cor_mean_raw=mean(res_tmp$cor_MegaLMM-res_tmp$cor_GBLUP,na.rm=T),
               cor_mean_meta=cor_mean_meta,
               cor_sd_raw=sd(res_tmp$cor_MegaLMM-res_tmp$cor_GBLUP,na.rm=T),
               cor_sd_meta=cor_sd_meta,
               cor_pval_meta=cor_pval_meta,
               # CI = xba +/- z*sd/sqrt(n) # 95% CI for sample mean
               lower.CL=cor_mean_meta - z*(cor_sd_meta/sqrt(NumEnv)),
               upper.CL=cor_mean_meta + z*(cor_sd_meta/sqrt(NumEnv))
    )
  }
}
dim(res_meta)
head(res_meta)
head(res)
res_CV4_GY=res[res$CVScn=="CV4"&res$TraitName=="Grain Yield",]
dim(res_CV4_GY)
head(res_CV4_GY)
res_CV4_GY = reshape2::melt(data = res_CV4_GY, id.vars = "EnvName",
                            measure.vars = c("cor_MegaLMM","cor_GBLUP","cor_MegaLMMGBLUP"))
ggplot(res_CV4_GY[res_CV4_GY$variable %in% c("cor_MegaLMM","cor_GBLUP"),], 
       aes(x=variable, y = value)) +
  geom_boxplot(aes(fill = variable)) + #,alpha = 0.5
  geom_point(position=position_dodge(width=0.75),cex=0.5)+
  ylab("Prediction Accuracy")+
  geom_line(aes(group = EnvName),
            alpha = 0.4, colour = "grey") 

ggplot(res_CV4_GY, 
       aes(x=variable, y = value)) +
  geom_boxplot(aes(fill = variable))
head(res_CV4_GY)

# ------------------------------------------------------------------------------
unique(res_meta$CVScn)
# "CV1"          "CV2"          "CV3NewTester" "CV4"          "CV3NewTrial"  "CV3NewState" 

resPlot = res_meta

resPlot$CVScn=gsub("CV3","",resPlot$CVScn)
resPlot$CVScn[resPlot$CVScn=="CV4"]="NewGenoNewYear"
unique(resPlot$CVScn)
CVScns_Sel=c("NewTrial","NewState","NewTester","NewGenoNewYear")
resPlot = resPlot[resPlot$CVScn %in% CVScns_Sel,]


# define orders
resPlot$CVScn = factor(resPlot$CVScn,levels=CVScns_Sel)

# get plot data
plot=ggplot(resPlot, 
            aes(x = TraitName, y = cor_mean_meta,fill=CVScn)) +
  geom_bar(position="dodge", stat="identity")
(plot_data <- data.frame(ggplot_build(plot)$data))
resPlot = cbind(resPlot,plot_data)

(figure_caption=sprintf('%s%s_G2F_MegaLMM_%s_vs_GBLUP_%s_barplot.pdf',
                        DIR_Output,Date,MegaLMM_Models_Sel,GBLUP_Models_Sel))
pdf(paste0(figure_caption,".pdf"),width = 10,height = 5)
padding=0.02
errorbarwidth=as.numeric(plot_data$xmax-plot_data$xmin)[1]/4
ggplot(resPlot, 
       aes(x = TraitName, y = cor_mean_meta,fill=CVScn)) +
  # aes(x = TraitName, y = cor_mean_meta,color=CVScn)) +
  labs(y = expression(r[MegaLMM] - r[GBLUP]))+
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
  geom_text(aes(y=max(upper.CL)*1.01,
                label=ifelse(cor_pval_meta>0.05,"ns",
                             ifelse(cor_pval_meta<0.01,"**","*"))
                # position="dodge2",
                ),
            cex=5,
            # vjust=1,
            position = position_dodge(width = 0.9)
            # position="dodge2"
            )+
  annotate(geom="text", x=resPlot$x,
           y=0-0.002, #-padding,
           vjust=1,
           label=sprintf("%.3f",resPlot$cor_mean_meta))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        legend.title = element_blank(),
        # legend.direction = "horizontal",
        legend.position = "right",
        legend.justification = "bottom") 
  # guides(fill = guide_legend(nrow = 1))
dev.off()

# ---------------------------------------------
# based on raw data
head(resPlot)
(figure_caption=sprintf('%s%s_G2F_MegaLMM_%s_vs_GBLUPBest_barplot_rawdata.pdf',
                        DIR_Output,Date,MegaLMM_Models_Sel))
pdf(paste0(figure_caption,".pdf"),width = 10,height = 5)
ggplot(resPlot, 
       aes(x = TraitName, y = cor_mean_raw,fill=CVScn)) +
  labs(y = expression(r[MegaLMM] - r[GBLUP]))+
  geom_bar(position="dodge2", stat="identity")+
  annotate(geom="text", x=resPlot$x,
           # y=max(resPlot$upper.CL)+padding,
           y=0,
           angle=90,
           label=sprintf("%.3f",resPlot$cor_mean_raw))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=14),
        legend.title = element_blank(),
        # legend.direction = "horizontal",
        legend.position = "right",
        legend.justification = "bottom"
        ) 
dev.off()

