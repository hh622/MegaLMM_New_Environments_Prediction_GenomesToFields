# ============================================
# script05_Prepare_MegaLMM_Y_Matrix
# ============================================

# Prepare_MegaLMM_Y_Matrix
# P1 in rows
# Year-Location-Tester Combination in Column

# Clean workspace
rm(list=ls())

library(reshape2)
library(tibble)
library(cvTools)
library(ggplot2)
library(ggpubr)
library(dplyr)

Date="2023-09-15"
DIR_Input=paste0("../03_Phenodat_BLUE/")
DIR_Output = paste0("../05_MegaLMM_PhenoDat_Format2/")
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)

# define global variables
traitnames = rev(c("Silk_DAP_days","Plant_Height_cm","Yield_Mg_ha"))
FamSize_cutoff=50 #tester family size cutoff
numParent1PerRow_cutoff=4

# ------------------------------------------------------------------------------
# loading data
# ------------------------------------------------------------------------------
list.files(DIR_Input)

# ------------------------------------------------------------------------------
# loading phenotypic data
ls=load(file = sprintf("%s%sPhenodat_Single_Env_Analysis_BLUE_RepFixed_BlkRandom.RData",
                       DIR_Input,"2023-09-15"))
ls #"Pheno_Trn_ModFit_Input" "Pheno_Trn_BLUE" "df_BLUE_notes"
Pheno_Trn_BLUE0=Pheno_Trn_BLUE
dim(Pheno_Trn_BLUE0) #182503     25
head(Pheno_Trn_BLUE0)
sum(is.na(Pheno_Trn_BLUE0$Prediction))/dim(Pheno_Trn_BLUE0)[1] #0

length(unique(Pheno_Trn_BLUE0$EnvNew[Pheno_Trn_BLUE0$TraitName==traitnames[1]]))#302
length(unique(Pheno_Trn_BLUE0$EnvNew[Pheno_Trn_BLUE0$TraitName==traitnames[2]]))#278
length(unique(Pheno_Trn_BLUE0$EnvNew[Pheno_Trn_BLUE0$TraitName==traitnames[3]]))#231

# ------------------------------------------------------------------------------
# load K matrix for P1 (created by DER on 02/03/23)
list.files('../../P5G2F/02_GRMs/DER/') 
#"A_P1.rds"     "A_P2_new.rds" "A_P2.rds"     "README.txt"  
K0 = readRDS('../../P5G2F/02_GRMs/DER/A_P1.rds')
K_info = data.frame(Hybrid = rownames(K0))
K_info$Parent1 = sapply(rownames(K0),function(x) strsplit(x,'/')[[1]][1])
K_info$Parent2 = sapply(rownames(K0),function(x) strsplit(x,'/')[[1]][2])
head(K_info)

dim(K0) #4918 4918
identical(rownames(K0),colnames(K0))
K0[1:6,1:6]

# ------------------------------------------------------------------------------
# load K matrix for P2 (created by DER)
K_Tester0 = readRDS('../../P5G2F/02_GRMs/DER/A_P2_new.rds')
dim(K_Tester0) #57 57
K_Tester0[1:4,1:4]

# ------------------------------------------------------------------------------
# filtering for Hybrids with P1/P2 genotype data available
Pheno_Trn_BLUE = Pheno_Trn_BLUE0[Pheno_Trn_BLUE0$Hybrid %in% rownames(K0),]
dim(Pheno_Trn_BLUE) #180225     25
head(Pheno_Trn_BLUE)
if(all(Pheno_Trn_BLUE$Hybrid %in% K_info$Hybrid)){
  print(
    all(unique(Pheno_Trn_BLUE$Hybrid_Parent1) %in% K_info$Parent1|
          unique(Pheno_Trn_BLUE$Hybrid_Parent1) %in% K_info$Parent2)
  ) #TRUE
  print(sum(!(unique(Pheno_Trn_BLUE$Hybrid_Parent1) %in% K_info$Parent1)))#2
}
length(unique(Pheno_Trn_BLUE$Hybrid_Parent1)) #1749
length(unique(K_info$Parent1)) #2189
dim(K_info) #4918    3
# ------------------------------------------------------------------------------
# harmonize inconsistency of Parent1 and Parent2 between Pheno and K_info
K_info = K_info[K_info$Hybrid %in% Pheno_Trn_BLUE$Hybrid,]
Pheno_Trn_BLUE_uniqHybrid = Pheno_Trn_BLUE[!duplicated(Pheno_Trn_BLUE$Hybrid),]

dim(K_info) #4057    3
dim(Pheno_Trn_BLUE_uniqHybrid) #4057   25

if(all(Pheno_Trn_BLUE_uniqHybrid$Hybrid %in% K_info$Hybrid)){
  K_info$Hybrid_Parent1 = Pheno_Trn_BLUE_uniqHybrid$Hybrid_Parent1[
    match(K_info$Hybrid,Pheno_Trn_BLUE_uniqHybrid$Hybrid)]
  K_info$Hybrid_Parent2 = Pheno_Trn_BLUE_uniqHybrid$Hybrid_Parent2[
    match(K_info$Hybrid,Pheno_Trn_BLUE_uniqHybrid$Hybrid)]
}
# checking
head(K_info)
all(unique(Pheno_Trn_BLUE$Hybrid_Parent1) %in% K_info$Hybrid_Parent1)
all(unique(Pheno_Trn_BLUE$Hybrid_Parent2) %in% K_info$Hybrid_Parent2)
sum(K_info$Parent1!=K_info$Hybrid_Parent1)

sum(!(K_info$Parent1%in%K_info$Hybrid_Parent1))
sum(!(K_info$Hybrid_Parent1 %in% K_info$Parent1))
length(unique(K_info$Hybrid_Parent1))#1749
length(unique(K_info$Parent1)) #1746
length(unique(Pheno_Trn_BLUE$Hybrid_Parent1)) #1749

stopifnot(all(Pheno_Trn_BLUE$Hybrid_Parent1 %in% K_info$Hybrid_Parent1))
stopifnot(all(Pheno_Trn_BLUE$Hybrid_Parent2 %in% rownames(K_Tester0)))

# ------------------------------------------------------------------------------
# reshape Pheno_Trn_BLUE to MegaLMM Input Format 1 amd Format 2
# ------------------------------------------------------------------------------

# create empty list()
Phe_Trn_Fmt1Wide=list()
Phe_Trn_Fmt2Wide=list()
Phe_Trn_Fmt2Long=list()
traitname = traitnames[1]
for(traitname in traitnames){
  # ----------------------------------------------------------------------------
  # Step 1: select phenotypic data
  phe_trn_sel=Pheno_Trn_BLUE[Pheno_Trn_BLUE$TraitName==traitname,]
  dim(phe_trn_sel) #67184    25; 66376    25
  phe_trn_sel[1:4,]
  
  # -------------------------------------------------------
  # convert to EnvNew format 1
   phe_trn_fmt1=
     reshape2::dcast(data = phe_trn_sel,
                     Hybrid~EnvNew,
                     value.var = "Prediction"
     )
   phe_trn_fmt1=column_to_rownames(.data = phe_trn_fmt1, var = "Hybrid")
   phe_trn_fmt1[1:10,1:4]
   dim(phe_trn_fmt1)#4125  302; 4057  302
   sum(duplicated(rownames(phe_trn_fmt1)))
   # filtering columns in which all values are NAs
   sort(colSums(!is.na(phe_trn_fmt1)))[1:20]
   phe_trn_fmt1=phe_trn_fmt1[,colSums(!is.na(phe_trn_fmt1))!=0]
   dim(phe_trn_fmt1) #4125  302; 4057  302
   
   # -------------------------------------------------------
   # convert to EnvNew format 2
   phe_trn_fmt2=
     reshape2::dcast(data = phe_trn_sel,
                     Hybrid_Parent1~EnvNew,
                     value.var = "Prediction"
     )
   stopifnot(all(phe_trn_fmt2$Hybrid_Parent1 %in% K_info$Hybrid_Parent1))
   dim(phe_trn_fmt2)#1787  303; 1749  303
   sum(duplicated(phe_trn_fmt2$Hybrid_Parent1)) #0
   
   # put Hybrid_Parent1 as rownames
   phe_trn_fmt2=column_to_rownames(.data = phe_trn_fmt2, var = "Hybrid_Parent1")
   dim(phe_trn_fmt2)#1787  302; 1749  302
   phe_trn_fmt2[1:10,1:4]
   
   # filtering out columns in which all values are NAs
   sort(colSums(!is.na(phe_trn_fmt2)))[1:20]
   phe_trn_fmt2=phe_trn_fmt2[,colSums(!is.na(phe_trn_fmt2))!=0]
   dim(phe_trn_fmt2) #1787  302; 1749  302
   
   # filtering out row in which non NAs fewer than 5 columns
   sort(rowSums(!is.na(phe_trn_fmt2)))[1:100]
   phe_trn_fmt2=phe_trn_fmt2[rowSums(!is.na(phe_trn_fmt2))>=numParent1PerRow_cutoff,]
   dim(phe_trn_fmt2) #1706  302;1702  302
   
   # checking
   sort(colSums(!is.na(phe_trn_fmt2)))[1:20]
   sort(rowSums(!is.na(phe_trn_fmt2)))[1:20]
   
   # collect results
   print(traitname)
   phe_trn_fmt1[1:3,1:4]
   phe_trn_fmt2[1:3,1:4]
   
   Phe_Trn_Fmt1Wide[[traitname]]=phe_trn_fmt1
   Phe_Trn_Fmt2Wide[[traitname]]=phe_trn_fmt2
   phe_trn_sel2=phe_trn_sel[phe_trn_sel$Hybrid_Parent1 %in% rownames(phe_trn_fmt2),]
   Phe_Trn_Fmt2Long[[traitname]]=phe_trn_sel2
}
sapply(list(Phe_Trn_Fmt1Wide,Phe_Trn_Fmt2Wide), names)
Phe_Trn_Fmt2Long[[traitname]][1:3,1:4]

# output MegaLMM_Input_Format2
sapply(Phe_Trn_Fmt2Wide, dim)
#      Yield_Mg_ha Plant_Height_cm Silk_DAP_days
# [1,]        1702            1699          1663
# [2,]         302             278           231
if(all(rownames(Phe_Trn_Fmt2Wide[[traitnames[1]]]) %in% K_info$Hybrid_Parent1) &
   all(rownames(Phe_Trn_Fmt2Wide[[traitnames[2]]]) %in% K_info$Hybrid_Parent1) &
   all(rownames(Phe_Trn_Fmt2Wide[[traitnames[3]]]) %in% K_info$Hybrid_Parent1) &
   all(Phe_Trn_Fmt2Long[[traitnames[1]]]$Hybrid %in% K_info$Hybrid) &
   all(Phe_Trn_Fmt2Long[[traitnames[2]]]$Hybrid %in% K_info$Hybrid) &
   all(Phe_Trn_Fmt2Long[[traitnames[3]]]$Hybrid %in% K_info$Hybrid)
){
  save(Phe_Trn_Fmt1Wide,
       Phe_Trn_Fmt2Wide,
       Phe_Trn_Fmt2Long,
       K_info,
       file = sprintf("%s%sMegaLMM_PhenoDat_Format2.RData",DIR_Output,Date)
  )
}

# ------------------------------------------------------------------------------
# MegaLMMInputFormat2 -Matrix plot - axis = Env::Tester (Site_Year_Tester)
# ------------------------------------------------------------------------------

(figure_caption=sprintf('%s%s_Phenotye_Distri_%s_xaxis%s',
                        DIR_Output,Date,"PheTrnBLUE_MegaLMMInputFormat2_Matrixplot","EnvNew"))
traitname = traitnames[1]
pdf(paste0(figure_caption,".pdf"),width = 12,height = 6)
for(traitname in traitnames){
  phe_trn_sel = Phe_Trn_Fmt2Long[[traitname]]
  phe_trn_sel=phe_trn_sel[order(phe_trn_sel$Year),]
  phe_trn_sel$Year=factor(phe_trn_sel$Year)
  phe_trn_sel$EnvNew=factor(phe_trn_sel$EnvNew,levels = unique(phe_trn_sel$EnvNew))
  phe_trn_sel$Hybrid_Parent1=factor(phe_trn_sel$Hybrid_Parent1,
                                    levels = rev(unique(phe_trn_sel$Hybrid_Parent1)))
  dat_Plot=phe_trn_sel
  dat_Plot[1:3,]
  str(dat_Plot)
  gg =
    ggplot(dat_Plot, aes(
      y=Hybrid_Parent1,
      x=EnvNew,
      fill=Year)) + 
    geom_raster() + 
    # ylim(0,length(levels(dat_Plot$Hybrid_Parent1))+10)+
    # scale_fill_gradient(low="white", high="blue",na.value = 'white') +
    labs(y="Hybrid_Parent1", x="Site::Year::Hybrid_Parent2", title="") +
    theme_classic() + theme(axis.text=element_blank(),
                            axis.ticks = element_blank(),
                            axis.title =element_text(size=12, face="bold"),
                            plot.title=element_text(size=14),
                            legend.position = 'none' #c(0.1,0.05),
                            # legend.position = 'bottom',
                            # legend.direction = "horizontal",
                            # legend.justification='right'
    )+# guides(fill = guide_legend(nrow = 1))
    ggtitle( 
      sprintf("%s (#Site::Year::Hybrid_Parent2 = %i), FamSizeCutoff = %i",
              traitname,
              length(unique(dat_Plot$EnvNew[!is.na(dat_Plot$Prediction)&
                                              dat_Plot$Prediction!=0])),
              FamSize_cutoff) )
  
  
  
  # Convert the "Env" variable to numeric
  head(dat_Plot)
  dat_Plot$EnvNew <- as.numeric(dat_Plot$EnvNew)
  
  # Compute boundaries for each fill color
  color_boundaries <- dat_Plot %>%
    group_by(Year) %>%
    summarize(x_min = min(EnvNew), x_max = max(EnvNew)) %>%
    ungroup()
  color_boundaries$x_bound_right = color_boundaries$x_max+0.5
  length(unique(dat_Plot$EnvNew))
  head(color_boundaries)
  
  # add a single boundary
  plot=
    gg + 
    # geom_vline(data = color_boundaries[-nrow(color_boundaries),], #most right section does not need a bound
    #            aes(xintercept = x_bound_right), color = "black",
    #            linetype = "dashed",lwd=0.3)+ #the line width will default to 0.5
    geom_text(data = color_boundaries,
              aes(x = (x_min+x_max)/2,
                  y = length(levels(dat_Plot$Hybrid_Parent1)),
                  label=Year),inherit.aes = F,vjust = 1
              
    )
  print(plot)
  
  
}
dev.off()

# ------------------------------------------------------------------------------
# MegaLMMInputFormat2 -Boxplot - axis = Env::Tester (Site_Year_Tester)
# ------------------------------------------------------------------------------

(figure_caption=sprintf('%s%s_Phenotye_Distri_%s_xaxis%s',
                        DIR_Output,Date,"PheTrnBLUE_MegaLMMInputFormat2_Boxplot","EnvNew"))
traitname = traitnames[1]
pdf(paste0(figure_caption,".pdf"),width = 12,height = 6)
for(traitname in traitnames){
  #must reassign for each trait since Env is converted to numeric below
  # dat_Plot=Pheno_Trn_BLUE[Pheno_Trn_BLUE$Traitname==traitname &
  #                         !is.na(Pheno_Trn_BLUE$Prediction),]
  dat_Plot = Phe_Trn_Fmt2Long[[traitname]]
  dat_Plot = dat_Plot[order(dat_Plot$Year),]
  dat_Plot$Year=factor(dat_Plot$Year)
  dat_Plot$EnvNew = factor(dat_Plot$EnvNew, levels = unique(dat_Plot$EnvNew))
  dat_Plot=dat_Plot[!is.na(dat_Plot$Prediction),]
  head(dat_Plot)
  sum(is.na(dat_Plot$Prediction))
  gg <-
    ggplot(data = dat_Plot, 
           aes(x = EnvNew, 
               y = Prediction, 
               fill = Year)) +
    geom_boxplot(outlier.size = 0.4, outlier.colour = "darkgray", lwd = 0.3) +
    xlab("Site::Year::Tester") +
    ylab(traitname) +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 90, size = 4),
      axis.ticks.x = element_line(linewidth = 0.3),
      legend.position="none"
      # legend.position = c(0.99, 0.99),
      # legend.justification = c("right", "top"),
      # legend.direction = "horizontal"
    ) +
    ggtitle( 
      sprintf("Model Fit Output - BLUE (#EnvNew(Site::Year::Tester) = %i), FamSizeCutoff = %i", 
              length(unique(dat_Plot$EnvNew[!is.na(dat_Plot$Prediction)&
                                              dat_Plot$Prediction!=0])),50 ) )
    
  # Convert the "Env" variable to numeric
  head(dat_Plot)
  dat_Plot$EnvNew <- as.numeric(dat_Plot$EnvNew)
  
  # Compute boundaries for each fill color
  color_boundaries <- dat_Plot %>%
    group_by(Year) %>%
    summarize(x_min = min(EnvNew), x_max = max(EnvNew)) %>%
    ungroup()
  color_boundaries$x_bound_right = color_boundaries$x_max+0.5
  length(unique(dat_Plot$EnvNew))
  
  # add a single boundary
  plot=gg + 
    geom_vline(data = color_boundaries[-nrow(color_boundaries),], #most right section does not need a bound
               aes(xintercept = x_bound_right), color = "black",
               linetype = "dashed",lwd=0.3)+ #the line width will default to 0.5
    geom_text(data = color_boundaries,
              aes(x = (x_min+x_max)/2,
                  y=max(dat_Plot$Prediction,na.rm = T),
                  label=Year),inherit.aes = F
    )
  print(plot)
  
}
dev.off()

  