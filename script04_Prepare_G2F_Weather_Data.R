#========================================================
# script04_Prepare_G2F_Weather_Data
#========================================================

# 2022-12-26
# Prepare_G2F_Weather_Data

# Clean workspace
rm(list=ls())

library(tibble)
library(weathermetrics)
library(dplyr)
library(foreach)
library(ComplexHeatmap) #Heatmap

Date="2022-12-26"
DIR_TRN="../Maize_GxE_Competition_Data/Training_Data/"
DIR_TST="../Maize_GxE_Competition_Data/Testing_Data/"
DIR_Output=paste0("../04_Env_Dat/",Date,"/")
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)

#------------------------------------
# Training data
#------------------------------------

#------------------------------------
# Training data - meta
meta_trn=read.csv(paste0(DIR_TRN,"2_Training_Meta_Data_2014_2021.csv"),
                  check.names = F,
                  stringsAsFactors = F)
meta_trn$Experiment_Code=gsub(" ","",meta_trn$Experiment_Code)
colnames(meta_trn)[colnames(meta_trn)=="Experiment_Code"]="Field_Location"
dim(meta_trn) #217  36
meta_trn[1:3,1:10]
meta_trn[,1:4]

sapply(list(meta_trn$Year,meta_trn$Env,meta_trn$Field_Location), 
       function(x) length(unique(x)))
# 8 217  45

# summary (TRN::meta):
# => 217  Environments
# => 8    years
# => 45   Field Locations

#------------------------------------
# Training data - pheno
phe_trn=read.csv(paste0(DIR_TRN,"1_Training_Trait_Data_2014_2021.csv"),
                 stringsAsFactors = F)
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

#------------------------------------
# Training data - Weather
wea_trn=read.csv(paste0(DIR_TRN,"4_Training_Weather_Data_2014_2021.csv"),
                 stringsAsFactors = F)
wea_trn=add_column(.data = wea_trn,.before = "Env",
                   Year=substr(wea_trn$Date,1,4) )
wea_trn=add_column(.data = wea_trn,.after = "Env",
                   Field_Location=sapply(strsplit(wea_trn$Env,"_"),"[",1) )
wea_trn$Year=as.integer(wea_trn$Year)

dim(wea_trn)#77431    20
wea_trn[1:10,]
sapply(list(wea_trn$Year,wea_trn$Env), function(x) length(unique(x)))
# 8 212
sort(unique(wea_trn$Year))
sum(is.na(wea_trn))

# summary(TRN::wea)
# =>   8  years (2014:no soil data)
# => 212  Environments
# => each Env has one observation

#------------------------------------
# Testing data
#------------------------------------

list.files(DIR_TST)
# [1] "1_Submission_Template_2022.csv"  "2_Testing_Meta_Data_2022.csv"   
# [3] "3_Testing_Soil_Data_2022.csv"    "4_Testing_Weather_Data_2022.csv"
# [5] "6_Testing_EC_Data_2022.csv" 

#------------------------------------------------
# Testing data - meta
meta_tst=read.csv(paste0(DIR_TST,"2_Testing_Meta_Data_2022.csv"),
                  check.names = F, stringsAsFactors = F)
colnames(meta_tst)[colnames(meta_tst)=="Experiment_Code"]="Field_Location"
meta_tst[1:3,1:4]
sapply(list(meta_tst$Year,meta_tst$Env,meta_tst$Field_Location), 
       function(x) length(unique(x)))
# 1 26 26

# summary (TRN::meta):
# => 26  Environments
# => 1    years
# => 26   Field Locations

#------------------------------------
# Training data - Weather
wea_tst=read.csv(paste0(DIR_TST,"4_Testing_Weather_Data_2022.csv"),
                 stringsAsFactors = F)
wea_tst=add_column(.data = wea_tst,.before = "Env",
                   Year=substr(wea_tst$Date,1,4) )
wea_tst=add_column(.data = wea_tst,.after = "Env",
                   Field_Location=sapply(strsplit(wea_tst$Env,"_"),"[",1) )
wea_tst$Year = as.integer(wea_tst$Year)
dim(wea_tst)#8164   20
wea_tst[1:3,]
sapply(list(wea_tst$Year,wea_tst$Env), function(x) length(unique(x)))
# 1 26
# summary(TST::wea)
# =>   1  years
# =>  26  Environments
# => each Env has one observation

# ------------------------------------------------------------------------------
# Merge wea_trn and wea_tst 
# ------------------------------------------------------------------------------

# ----------------------------------------------------
# is colnames(wea_trn) identical to colnames(wea_tst)?
sapply(list(wea_trn,wea_tst), dim)
#   [,1] [,2]
# [1,] 77431 8164
# [2,]    20   20

identical(colnames(wea_trn),colnames(wea_tst))#FALSE
identical(setdiff(colnames(wea_trn),colnames(wea_tst)),character(0))#TRUE
# colnames(wea_trn) are identical to colnames(wea_tst)
# but in different order

# ---------------------------------------------------------------
# what is the proportion of missing values in wea_trn and wea_tst
apply(wea_trn, 2, function(x) sum(is.na(x)))/nrow(wea_trn)
# Year                Env     Field_Location               Date               QV2M 
# 0                  0                  0                  0                  0 
# T2MDEW                 PS               RH2M               WS2M            GWETTOP 
# 0                  0                  0                  0                  0 
# ALLSKY_SFC_SW_DWN ALLSKY_SFC_PAR_TOT            T2M_MAX            T2M_MIN             T2MWET 
# 0                  0                  0                  0                  0 
# GWETROOT                T2M           GWETPROF  ALLSKY_SFC_SW_DNI        PRECTOTCORR 
# 0                  0                  0                  0                  0 

apply(wea_tst, 2, function(x) sum(is.na(x)))/nrow(wea_tst)
#       Year                Env     Field_Location               Date               QV2M 
# 0.00000000         0.00000000         0.00000000         0.00000000         0.01592357 
#         PS               WS2M  ALLSKY_SFC_SW_DWN               RH2M            GWETTOP 
# 0.01592357         0.01592357         0.03466438         0.01592357         0.80891720 
#     T2MWET            T2M_MAX             T2MDEW           GWETROOT           GWETPROF 
# 0.01592357         0.01592357         0.01592357         0.80891720         0.80891720 
#        T2M            T2M_MIN ALLSKY_SFC_PAR_TOT        PRECTOTCORR  ALLSKY_SFC_SW_DNI 
# 0.01592357         0.01592357         0.80891720         0.01592357         0.80891720 

colnames(wea_tst)[apply(wea_tst, 2, function(x) sum(is.na(x)))/nrow(wea_tst)>0.8]
# "GWETTOP" "GWETROOT" "GWETPROF" "ALLSKY_SFC_PAR_TOT" "ALLSKY_SFC_SW_DNI"

# summary:
# - wea_trn has no missing values
# - wea_tst has >80% missing values for five weather variables:
#   "GWETTOP" "GWETROOT" "GWETPROF" "ALLSKY_SFC_PAR_TOT" "ALLSKY_SFC_SW_DNI"

# ---------------------------------------------------------------
# merge wea_trn and wea_tst
colnames(wea_trn)
colnames(wea_tst)

wea_tst=wea_tst[,colnames(wea_trn)] #sort wea_tst
sapply(list(wea_trn,wea_tst), function(x) length(unique(x$Env)))
# 212  26

if(identical(colnames(wea_trn),colnames(wea_tst))){
  Weadat = rbind(wea_trn,wea_tst)
}
dim(Weadat)#85595    20
Weadat[1:4,]
sapply(list(Weadat$Year, Weadat$Env, Weadat$Field_Location), function(x) length(unique(x)))
# 9 238  45
# note: TRN has 217 Envs, TST has 26 Envs; total = 243

# --------------------------------------------------------------------
# Convert from Celsius to Fahrenheit.
# celsius.to.fahrenheit(c(10,NA,30)) #e.g.
colnames(Weadat)[colnames(Weadat)=="T2M_MAX"]="T2M_MAX_C"
colnames(Weadat)[colnames(Weadat)=="T2M_MIN"]="T2M_MIN_C"
Weadat=add_column(.data = Weadat,
                  .after = "T2M_MIN_C",
                  T2M_MAX_F=celsius.to.fahrenheit(Weadat$T2M_MAX_C),
                  T2M_MIN_F=celsius.to.fahrenheit(Weadat$T2M_MIN_C))

# add a column of EnvDate
Weadat=add_column(.data = Weadat,
                  .after = "Date",
                  EnvDate=paste(Weadat$Env,Weadat$Date,sep = "::")
)
dim(Weadat) #85595    23

# exclude rows where all climate data are NAs
colnames(Weadat)
idx=apply(Weadat[,!(colnames(Weadat)%in%c("Year","Env","Field_Location","Date","EnvDate"))],
          1, function(x) all(is.na(x)))
sum(idx) #130
Weadat=Weadat[!idx,]
dim(Weadat) #85465    23
length(unique(Weadat$Env))

# ------------------------------------------------------------------------------
# Impute missing Date_Planted and Date_Harvested
# ------------------------------------------------------------------------------
# Methods:
# - step 1: impute Date_Planted using mean value of all years at a specific site
# - step 2: impute Date_Harvested = Date_Planted + average growing days

Date_Plt_Harv_TRN=phe_trn[,c("Env","Year","Field_Location","Date_Planted","Date_Harvested")][!duplicated(phe_trn$Env),]
dim(Date_Plt_Harv_TRN)#217   5
Date_Plt_Harv_TST=meta_tst[,c("Env","Year","Field_Location","Date_Planted")]
Date_Plt_Harv_TST$Date_Harvested=NA
dim(Date_Plt_Harv_TST)#26  5
Date_Plt_Harv = rbind(Date_Plt_Harv_TRN,Date_Plt_Harv_TST)
dim(Date_Plt_Harv)#243   5
head(Date_Plt_Harv)
# to be simple, format Date to YYYY-MM-DD with Excel, then read in again
if(!file.exists(paste0(DIR_Output,"Date_Plt_Harv_raw.csv"))){
  write.csv(Date_Plt_Harv,paste0(DIR_Output,"Date_Plt_Harv_raw.csv"),row.names = F)
}
# read in formmated date with YYYY-MM-DD format
Date_Plt_Harv=read.csv(paste0(DIR_Output,"Date_Plt_Harv_formatted.csv"),stringsAsFactors = F)
head(Date_Plt_Harv)
unique(Date_Plt_Harv$Date_Planted)
unique(Date_Plt_Harv$Date_Harvested)

# convert calendar date to Julian date
require(lubridate)
Date_Plt_Harv$Date_Planted_Juli=yday(Date_Plt_Harv$Date_Planted)
Date_Plt_Harv$Date_Harvested_Juli=yday(Date_Plt_Harv$Date_Harvested)

# calculate Dates_Growth=Date_Harvested-Date_Planted
Date_Plt_Harv$Dates_Growth=Date_Plt_Harv$Date_Harvested_Juli-Date_Plt_Harv$Date_Planted_Juli
# notes: here Dates_Growth is based on raw data, so some are missings

# impute missing planted and harvested dates
Date_Plt_Harv$Date_Planted_Imputed=Date_Plt_Harv$Date_Planted_Juli
Date_Plt_Harv$Date_Harvested_Imputed=Date_Plt_Harv$Date_Harvested_Juli
setdiff(Date_Plt_Harv$Env[is.na(Date_Plt_Harv$Date_Planted)],
        Date_Plt_Harv$Env[is.na(Date_Plt_Harv$Date_Harvested)])#character(0)
length(setdiff(Date_Plt_Harv$Env[is.na(Date_Plt_Harv$Date_Harvested)],
               Date_Plt_Harv$Env[is.na(Date_Plt_Harv$Date_Planted)]))#18
# => all Envs with missing Date_Planted also have missing Date_Harvested
# => but not vice versa

# impute Date_Planted first
# fldloc = Date_Plt_Harv$Field_Location[is.na(Date_Plt_Harv$Date_Planted)][1]
for(fldloc in Date_Plt_Harv$Field_Location[is.na(Date_Plt_Harv$Date_Planted)]){
  # extract Date_Plt_Harv for a specific Field_Location
  tmp=Date_Plt_Harv[Date_Plt_Harv$Field_Location==fldloc,]
  # impute with mean of the same Field_Location across years
  tmp$Date_Planted_Imputed[is.na(tmp$Date_Planted)]=round(mean(tmp$Date_Planted_Juli,na.rm = T),0)
  # update including Date_Planted_Imputed
  Date_Plt_Harv[Date_Plt_Harv$Field_Location==fldloc,]=tmp
}

# Then, impute Date_Harvested
fldloc = Date_Plt_Harv$Field_Location[is.na(Date_Plt_Harv$Date_Harvested)][1]
fldloc = "NEH3"
for(fldloc in Date_Plt_Harv$Field_Location[is.na(Date_Plt_Harv$Date_Harvested)]){
  # extract Date_Plt_Harv for a specific Field_Location
  tmp=Date_Plt_Harv[Date_Plt_Harv$Field_Location==fldloc,]
  
  # calclualte averege datas of growth
  growthdates_ave=round(mean(tmp$Dates_Growth,na.rm = T),0)
  
  # Date_Harvested_Imputed = Date_Planted_Imputed + growthdates_ave
  # impute with mean of the same Field_Location across years
  tmp$Date_Harvested_Imputed[is.na(tmp$Date_Harvested)]=
    tmp$Date_Planted_Imputed[is.na(tmp$Date_Harvested)]+
    growthdates_ave
  
  # update including Date_Harvested_Imputed
  Date_Plt_Harv[Date_Plt_Harv$Field_Location==fldloc,]=tmp
}

sum(is.na(Date_Plt_Harv$Date_Planted_Imputed))==0 #TRUE
sum(is.na(Date_Plt_Harv$Date_Harvested_Imputed))==0 #TRUE

# ------------------------------------------------------------------------------
# add Date_Planted and Date_Harvested to Weadat
lapply(list(Weadat,Date_Plt_Harv), head)
intersect(colnames(Weadat),colnames(Date_Plt_Harv))
Weadat_tmp=left_join(Weadat,Date_Plt_Harv,by=c("Year","Env","Field_Location"))
head(Weadat_tmp)
Weadat = add_column(.data = Weadat, .after = "Date",
                    # Date_Planted_Raw=Weadat_tmp$Date_Planted,
                    Date_Planted_Imputed=Weadat_tmp$Date_Planted_Imputed,
                    # Date_Harvested_Raw=Weadat_tmp$Date_Harvested,
                    Date_Harvested_Imputed=Weadat_tmp$Date_Harvested_Imputed)
# add Dates_Growth based on Date_Planted_Imputed and Date_Harvested_Imputed,
# so should be no NAs
Weadat = add_column(.data = Weadat, 
                    .after = "Date_Harvested_Imputed",
                    Dates_Growth=Weadat$Date_Harvested_Imputed-Weadat$Date_Planted_Imputed)
head(Weadat)
sum(is.na(Weadat$Dates_Growth))==0 #TRUE

# --------------------------------------------------------------------
# Estimate Growth Stages and Average original weather data
# --------------------------------------------------------------------
# reference: https://ndawn.ndsu.nodak.edu/help-corn-growing-degree-days.html
# Methods:
# Step 1: calculate Daily Corn GDD (growing degree days)
# Corn Growing Degree Day Calculation
# The Daily Average Temp (°F) = (Daily Max Temp °F + Daily Min Temp °F) / 2
# Daily Corn GDD (°F) = Daily Average Temperature °F - 50 °F
# If the daily Max and/or Min Temp < 50 °F (10 °C), it's set equal to 50 °F (10 °C)
# If the daily Max Temperature > 86 °F (30 °C), it's set equal to 86 °F (30 °C)

# Step 2: calculate AGDD and estimate growth stages
# Corn generally requires about 82 to 85 GDD from to complete a leaf collar emergence, 
# up to growth stage V10, later vegetative stages only require about 50 GDD for collar emergence. 
# It usually takes 115 to 120 GDD for corn to emerge after planting. 
# Therefore, if a field has accumulated 380 GDD from date of planting, 
# subtracting 115 GDD from 380 GDD equals 265 GDD, so the growth stage estimate for the field would be early V3.
# Ref1: https://www.dekalbasgrowdeltapine.com/en-us/agronomy/corn-growth-stages-and-gdu-requirements.html#:~:text=Corn%20generally%20requires%20about%2082,corn%20to%20emerge%20after%20planting
# Ref2: https://www.agry.purdue.edu/ext/corn/news/timeless/VStagePrediction.html
# Ref3: https://mygeohub.org/resources/879/download/Corn-growth-stage-day-and-GDU-calendar10.pdf

# Step 3: average weather data according to estimated growth stages

# ------------------------------------------------------------------------------
# Step 1: calculate Daily Corn GDD
EnvDate = unique(Weadat$EnvDate)[1]
options(max.print = 99990)
Weadat_GDD=NULL
Weadat_GDD=foreach(EnvDate = unique(Weadat$EnvDate),.combine = "rbind") %do% {
  tmp=Weadat[Weadat$EnvDate==EnvDate,]
  DailyAverageTemperatureF=(tmp$T2M_MAX_F+tmp$T2M_MIN_F)/2
  if(DailyAverageTemperatureF<50) DailyAverageTemperatureF=50
  if(DailyAverageTemperatureF>86) DailyAverageTemperatureF=86
  GDD=DailyAverageTemperatureF-50
  tmp=add_column(.data = tmp,
                 .after = "T2M_MIN_F",
                 GDD=GDD,
                 GrowthStage=NA)
  return(tmp)
}
dim(Weadat)
dim(Weadat_GDD)
Weadat_GDD[1:3,]
length(unique(Weadat_GDD$Env))#238

# convert Weadat_GDD Calendar Data to Julian Date
Weadat_GDD=add_column(.data = Weadat_GDD,
                      .after = "Date",
                      MM=substr(Weadat_GDD$Date,5,6),
                      DD=substr(Weadat_GDD$Date,7,8)
)
Weadat_GDD=add_column(.data = Weadat_GDD,
                      .after = "Date",
                      Date_Cale=paste(Weadat_GDD$Year,Weadat_GDD$MM,Weadat_GDD$DD,sep = "-")
)
Weadat_GDD=add_column(.data = Weadat_GDD,
                      .after = "Date_Cale",
                      Date_Juli=yday(Weadat_GDD$Date_Cale)
)
Weadat_GDD=Weadat_GDD[,!(colnames(Weadat_GDD) %in% c("MM","DD"))]
dim(Weadat_GDD)
# ------------------------------------------------------------------------------
# calculate AGDD and APRE(accumulated precipitation)
Env = "GAH2_2016"
Weadat_AGDD=foreach(Env = unique(Weadat_GDD$Env),.combine = "rbind") %do% {
  # extract weather data for a specific Env
  tmp=Weadat_GDD[Weadat_GDD$Env==Env,]
  dim(tmp) #366  29
  tmp[1:4,]

  # extract weather data for dates between Planted - Harvested
  tmp=tmp[tmp$Date_Juli>=tmp$Date_Planted_Imputed &
            tmp$Date_Juli<=tmp$Date_Harvested_Imputed,]
  dim(tmp)#146  29
  
  # calculate AGDD
  tmp=add_column(.data = tmp,
             .after = "GDD",
             AGDD=cumsum(tmp$GDD),
             APRE=cumsum(tmp$PRECTOTCORR))
  
  # tmp[,c("Env","Date_Cale","Date_Planted_Imputed","Date_Harvested_Imputed",
  #        "Dates_Growth","GDD","AGDD")]
  # dim(tmp)
  tmp
}
dim(Weadat_AGDD)
head(Weadat_AGDD)

# G2F_correlation_between_weather_variables
wea_var_names=colnames(wea_trn)[-1*1:4]
length(wea_var_names)
sum(colnames(Weadat_AGDD) %in% wea_var_names)
head(Weadat_AGDD[,colnames(Weadat_AGDD) %in% c(wea_var_names,"T2M_MAX_C","T2M_MIN_C","GDD","AGDD")])
Heatmap(cor(Weadat_AGDD[,colnames(Weadat_AGDD) %in% c(wea_var_names,"T2M_MAX_C","T2M_MIN_C","GDD","AGDD","APRE")],
            use="c"))
# note: add GDD and AGDD do not make much sense

# extract growth period AGDD for each Env
Weadat_AGDD_eachEnv=aggregate(AGDD~Env,data=Weadat_AGDD,FUN = max)

# distribution of AGDD
hist(Weadat_AGDD_eachEnv$AGDD,50,
     main=paste("Histogram of Weadat_AGDD (N=",length(unique(Weadat_AGDD_eachEnv$Env)),")"))
abline(v=2850,col="red",lwd=2,lty=2)
table(Weadat_AGDD_eachEnv$AGDD>2850)
# FALSE  TRUE 
# 90   148 

# add US state to Weadat_AGDD
head(Weadat_AGDD_eachEnv)
Weadat_AGDD_eachEnv=add_column(.data = Weadat_AGDD_eachEnv,
                               .after = "Env",
                               Country="US",
                               State_abb=substr(Weadat_AGDD_eachEnv$Env,1,2))
Weadat_AGDD_eachEnv=add_column(.data = Weadat_AGDD_eachEnv,
                               .after = "State_abb",
                               State_name=as.character(sapply(Weadat_AGDD_eachEnv$State_abb, function(x) state.name[grep(x, state.abb)])) )
# why cannot assign US states to some Envs (N=12)
sum(Weadat_AGDD_eachEnv$State_name=="character(0)")#12
Weadat_AGDD_nonUSEnv=Weadat_AGDD_eachEnv[Weadat_AGDD_eachEnv$State_name=="character(0)",]
meta_trn[meta_trn$Env %in% Weadat_AGDD_nonUSEnv$Env,c("Env",
                                     "Field_Location",
                                     "City",
                                     "Weather_Station_Latitude (in decimal numbers NOT DMS)",
                                     "Weather_Station_Longitude (in decimal numbers NOT DMS)")]
# Env Field_Location      City
# 36  GEH1_2019           GEH1 Gottingen
# 37  GEH1_2020           GEH1 Gottingen
# 38  GEH1_2021           GEH1 Gottingen
# 166 ONH1_2014           ONH1  Waterloo
# 167 ONH1_2015           ONH1  Waterloo
# 168 ONH1_2016           ONH1  Waterloo
# 169 ONH1_2017           ONH1  Waterloo
# 170 ONH2_2014           ONH2 Ridgetown
# 171 ONH2_2015           ONH2 Ridgetown
# 172 ONH2_2016           ONH2 Ridgetown
# 173 ONH2_2017           ONH2 Ridgetown
# 174 ONH2_2019           ONH2 Ridgetown
# ONH1  Waterloo (City in Ontario, Canada)
# GEH1 Gottingen (City in Germany)
# ONH2 Ridgetown (Town in Ontario, Canada)
# => summary: 3 Envs from Germany; 9 Envs from Canada

# add latitude and logitue to Weadat_AGDD
Weadat_AGDD_eachEnv[order(Weadat_AGDD_eachEnv$AGDD),]
colnames(meta_trn)
meta_dat=rbind(
  meta_trn[,c("Env",
              "Field_Location",
              "City",
              "Weather_Station_Latitude (in decimal numbers NOT DMS)",
              "Weather_Station_Longitude (in decimal numbers NOT DMS)",
              "Latitude_of_Field_Corner_#1 (lower left)",                          
              "Longitude_of_Field_Corner_#1 (lower left)",                        
              "Latitude_of_Field_Corner_#2 (lower right)",                        
              "Longitude_of_Field_Corner_#2 (lower right)",                        
              "Latitude_of_Field_Corner_#3 (upper right)",                         
              "Longitude_of_Field_Corner_#3 (upper right)",                        
              "Latitude_of_Field_Corner_#4 (upper left)",                          
              "Longitude_of_Field_Corner_#4 (upper left)")],
  meta_tst[,c("Env",
              "Field_Location",
              "City",
              "Weather_Station_Latitude (in decimal numbers NOT DMS)",
              "Weather_Station_Longitude (in decimal numbers NOT DMS)",
              "Latitude_of_Field_Corner_#1 (lower left)",                          
              "Longitude_of_Field_Corner_#1 (lower left)",                        
              "Latitude_of_Field_Corner_#2 (lower right)",                        
              "Longitude_of_Field_Corner_#2 (lower right)",                        
              "Latitude_of_Field_Corner_#3 (upper right)",                         
              "Longitude_of_Field_Corner_#3 (upper right)",                        
              "Latitude_of_Field_Corner_#4 (upper left)",                          
              "Longitude_of_Field_Corner_#4 (upper left)")])
colnames(meta_dat)=c("Env","Field_Location",
                     "City","Weather_Station_Latitude",
                     "Weather_Station_Longitude",
                     "Latitude_of_Field_Corner_lowerleft",                          
                     "Longitude_of_Field_Corner_lowerleft",                        
                     "Latitude_of_Field_Corner_lowerright",                        
                     "Longitude_of_Field_Corner_lowerright",                        
                     "Latitude_of_Field_Corner_upperright",                         
                     "Longitude_of_Field_Corner_upperright",                        
                     "Latitude_of_Field_Corner_upperleft",                          
                     "Longitude_of_Field_Corner_upperleft")
head(meta_dat)
apply(meta_dat, 2, function(x) sum(is.na(x)))

# add Latitude Lat/Long to Weadat_AGDD_eachEnv
Weadat_AGDD_eachEnv=left_join(Weadat_AGDD_eachEnv,
                              meta_dat,
                              by="Env")
# update Country/State for Canadian and German Envs
Weadat_AGDD_eachEnv$Country[Weadat_AGDD_eachEnv$Field_Location=="GEH1"]="Germany"
Weadat_AGDD_eachEnv$State_abb[Weadat_AGDD_eachEnv$Field_Location=="GEH1"]="GE"
Weadat_AGDD_eachEnv$State_name[Weadat_AGDD_eachEnv$Field_Location=="GEH1"]="Gottingen"

Weadat_AGDD_eachEnv$Country[Weadat_AGDD_eachEnv$Field_Location=="ONH1"]="Canada"
Weadat_AGDD_eachEnv$State_abb[Weadat_AGDD_eachEnv$Field_Location=="ONH1"]="ON"
Weadat_AGDD_eachEnv$State_name[Weadat_AGDD_eachEnv$Field_Location=="ONH1"]="Ontario"

Weadat_AGDD_eachEnv$Country[Weadat_AGDD_eachEnv$Field_Location=="ONH2"]="Canada"
Weadat_AGDD_eachEnv$State_abb[Weadat_AGDD_eachEnv$Field_Location=="ONH2"]="ON"
Weadat_AGDD_eachEnv$State_name[Weadat_AGDD_eachEnv$Field_Location=="ONH2"]="Ontario"

head(Weadat_AGDD_eachEnv)
colnames(Weadat_AGDD_eachEnv)
Weadat_AGDD_eachEnv[Weadat_AGDD_eachEnv$Env %in% Weadat_AGDD_nonUSEnv$Env,
                    c("Env","Field_Location","City","AGDD","Weather_Station_Latitude")]
#           Env Field_Location      City     AGDD Weather_Station_Latitude
# 34  GEH1_2019           GEH1 Gottingen 1885.345                 51.50218
# 35  GEH1_2020           GEH1 Gottingen 1787.475                 51.49566
# 36  GEH1_2021           GEH1 Gottingen 1608.525                 51.47064
# 180 ONH1_2014           ONH1  Waterloo 1895.865                 43.49703
# 181 ONH1_2015           ONH1  Waterloo 2273.415                 43.49845
# 182 ONH1_2016           ONH1  Waterloo 2649.840                 43.49779
# 183 ONH1_2017           ONH1  Waterloo 2170.345                 43.49718
# 184 ONH2_2014           ONH2 Ridgetown 2305.275                 42.45420
# 185 ONH2_2015           ONH2 Ridgetown 2608.125                 42.45346
# 186 ONH2_2016           ONH2 Ridgetown 3212.220                 42.45282
# 187 ONH2_2017           ONH2 Ridgetown 2700.615                 43.07471
# 188 ONH2_2019           ONH2 Ridgetown 2338.550                 42.45209

# impute missing Weather_Station_Latitude
# Method: 
# - same Field_Location but from different city quite apart from each other, e.g.TXH2
# - so better infer Latitude based on city (step 1)
# - then infer those all field locations with NAs based on Field_Location mean (step 2)

# step 1: infer Latitude based on city
for(i in 1:length(unique(Weadat_AGDD_eachEnv$City))){
  city=unique(Weadat_AGDD_eachEnv$City)[i]
  print(paste("###############",i,"::",city,"#####################"))
  tmp = Weadat_AGDD_eachEnv[Weadat_AGDD_eachEnv$City==city,]
  print(tmp[,c("Env","Field_Location","City","AGDD","Weather_Station_Latitude",
         "Weather_Station_Longitude")])
  if(all(is.na(tmp$Weather_Station_Latitude)) |
     all(!is.na(tmp$Weather_Station_Latitude))) {
    next
  } else{
    tmp$Weather_Station_Latitude=mean(tmp$Weather_Station_Latitude,na.rm = T)
    tmp$Weather_Station_Longitude=mean(tmp$Weather_Station_Latitude,na.rm = T)
    # due to five digits in latitude, those numbers are slightly different btw years
    # so unique() does not work, use mean() is just ok!
    Weadat_AGDD_eachEnv[Weadat_AGDD_eachEnv$City==city,] = tmp
  }
}
# step 2: infer those all field locations with NAs based on Field_Location mean
sum(is.na(Weadat_AGDD_eachEnv$Weather_Station_Latitude))#3
Weadat_AGDD_eachEnv[is.na(Weadat_AGDD_eachEnv$Weather_Station_Latitude),
                    c("Env","State_name","Field_Location","City","AGDD","Weather_Station_Latitude")]
# Env       State_name Field_Location    City     AGDD Weather_Station_Latitude
# NEH2_2018   Nebraska           NEH2   Wahoo 3296.385                       NA
# TXH2_2018      Texas           TXH2 Lubbock 4357.225                       NA
# TXH4_2019      Texas           TXH4 Lubbock 4034.520                       NA
fldloc="NEH2"
for(fldloc in unique(Weadat_AGDD_eachEnv$Field_Location[is.na(Weadat_AGDD_eachEnv$Weather_Station_Latitude)]) ){
  print(paste("###############",fldloc,"#####################"))
  tmp = Weadat_AGDD_eachEnv[Weadat_AGDD_eachEnv$Field_Location==fldloc,]
  print(tmp[,c("Env","Field_Location","City","AGDD","Weather_Station_Latitude",
               "Weather_Station_Longitude")])
  if(all(is.na(tmp$Weather_Station_Latitude)) |
     all(!is.na(tmp$Weather_Station_Latitude))) {
    next
  } else{
    tmp$Weather_Station_Latitude=mean(tmp$Weather_Station_Latitude,na.rm = T)
    tmp$Weather_Station_Longitude=mean(tmp$Weather_Station_Latitude,na.rm = T)
    Weadat_AGDD_eachEnv[Weadat_AGDD_eachEnv$Field_Location==fldloc,] = tmp
  }
 
}

# step 3: for the Env=TXH4_2019
sum(is.na(Weadat_AGDD_eachEnv$Weather_Station_Latitude))#1
Weadat_AGDD_eachEnv[is.na(Weadat_AGDD_eachEnv$Weather_Station_Latitude),
                    c("Env","State_name","Field_Location","City","AGDD","Weather_Station_Latitude")]
# Env State_name Field_Location    City    AGDD Weather_Station_Latitude
# TXH4_2019      Texas           TXH4 Lubbock 4034.52                      NaN
Weadat_AGDD_eachEnv$Weather_Station_Latitude[Weadat_AGDD_eachEnv$Field_Location=='TXH4']=
  mean(Weadat_AGDD_eachEnv$Weather_Station_Latitude[Weadat_AGDD_eachEnv$City=='Lubbock'],na.rm = T)
Weadat_AGDD_eachEnv$Weather_Station_Longitude[Weadat_AGDD_eachEnv$Field_Location=='TXH4']=
  mean(Weadat_AGDD_eachEnv$Weather_Station_Longitude[Weadat_AGDD_eachEnv$City=='Lubbock'],na.rm = T)
sum(is.na(Weadat_AGDD_eachEnv$Weather_Station_Latitude)) #0

# G2F_correlation_between_Latitude_and_AGDD
plot(Weadat_AGDD_eachEnv$Weather_Station_Latitude,
     Weadat_AGDD_eachEnv$AGDD,
     xlab = "Weather_Station_Latitude",ylab = "AGDD (plant-to-harvest)",
     main=paste("N=",length(unique(Weadat_AGDD_eachEnv$Env))))
abline(h=2865,lty=2,lwd=2,col=2)
abline(v=37,lty=2,lwd=2,col=2)
cor(Weadat_AGDD_eachEnv$Weather_Station_Latitude,
    Weadat_AGDD_eachEnv$AGDD,use = "c") # -0.8418439

if(!file.exists(paste0(DIR_Output,"G2F_Weadat_AGDD_LatLong.csv"))){
  write.csv(Weadat_AGDD_eachEnv,paste0(DIR_Output,"G2F_Weadat_AGDD_LatLong.csv"),row.names = F)
}
head(Weadat_AGDD_eachEnv)

# It seems that there are two groups of Envs
table(Weadat_AGDD_eachEnv$AGDD>2850)
# FALSE  TRUE 
# 90   148 
Weadat_AGDD_eachEnv_grp1=
  Weadat_AGDD_eachEnv[Weadat_AGDD_eachEnv$Weather_Station_Latitude<37,
                      c("Env","State_name","Field_Location","City","AGDD","Weather_Station_Latitude")]
Weadat_AGDD_eachEnv_grp1[order(Weadat_AGDD_eachEnv_grp1$AGDD,decreasing = T),]
unique(Weadat_AGDD_eachEnv_grp1$State_name)
# "Arkansas"       "Georgia"        "North Carolina" "South Carolina" "Texas"

Weadat_AGDD_eachEnv_grp2=
  Weadat_AGDD_eachEnv[Weadat_AGDD_eachEnv$Weather_Station_Latitude>37 &
                        Weadat_AGDD_eachEnv$Weather_Station_Latitude<45,
                      c("Env","State_name","Field_Location","City","AGDD","Weather_Station_Latitude")]
Weadat_AGDD_eachEnv_grp2[order(Weadat_AGDD_eachEnv_grp2$AGDD,decreasing = T),]
unique(Weadat_AGDD_eachEnv_grp2$State_name)
# [1] "Colorado"     "Delaware"     "Iowa"         "Illinois"     "Indiana"      "Kansas"      
# [7] "Michigan"     "Minnesota"    "Missouri"     "Nebraska"     "New York"     "Ohio"        
# [13] "Ontario"      "South Dakota" "Wisconsin" 

# ---------------------------------------------
# G2F_correlation_between_Latitude_and_APRE
head(Weadat_AGDD)
Weadat_APRE_eachEnv=aggregate(APRE~Env,data=Weadat_AGDD,FUN = max)
Weadat_APRE_eachEnv=left_join(Weadat_APRE_eachEnv,
                              Weadat_AGDD_eachEnv,
                              by="Env")
plot(Weadat_APRE_eachEnv$Weather_Station_Latitude,
     Weadat_APRE_eachEnv$APRE,
     xlab = "Weather_Station_Latitude",ylab = "APRE (plant-to-harvest)",
     main=paste("N=",length(unique(Weadat_AGDD_eachEnv$Env))))
plot(Weadat_APRE_eachEnv$APRE,Weadat_APRE_eachEnv$AGDD)

Weadat_APRE_eachEnv[order( Weadat_APRE_eachEnv$APRE),
                    c("Env","Field_Location","City","APRE","AGDD","Weather_Station_Latitude",
                      "Weather_Station_Longitude")]
# => accumulated precipitation has no correlation with lat, long and AGDD

if(!file.exists(paste0(DIR_Output,"G2F_Weadat_APRE_LatLong.csv"))){
  write.csv(Weadat_APRE_eachEnv,paste0(DIR_Output,"G2F_Weadat_APRE_LatLong.csv"),row.names = F)
}
head(Weadat_APRE_eachEnv)

# ------------------------------------------------------------------------------
# Step 2: calculate AGDD and estimate growth stages
head(Weadat)
Env = unique(Weadat_GDD$Env)[1]
Env = "GAH2_2016"

# GDD parameters for each corn developmental stage 
# Ref2: https://www.agry.purdue.edu/ext/corn/news/timeless/VStagePrediction.html
# Ref3: https://mygeohub.org/resources/879/download/Corn-growth-stage-day-and-GDU-calendar10.pdf

# - Corn emergence typically requires from 100 to 120 GDDs from planting
GDDs_Per_Stage_Emerge=120
# - From VE to V10, leaf collar emergence occurs at about one leaf every 82 GDDs
GDDs_Per_Stage_VE_to_V10= 82
# - From V10 to the final leaf, leaf collar emergence occurs approximately one leaf every 50 GDDs.
GDDs_Per_Stage_V10_to_V18=50
# - tassel and reproductive stages
GDDs_Per_Stage_VT=100 #tassel
GDDs_Per_Stage_R1=225
GDDs_Per_Stage_R2=200
GDDs_Per_Stage_R3=300
GDDs_Per_Stage_R4=200
GDDs_Per_Stage_R5=250
GDDs_Per_Stage_R6=250
# AGDD till end of V10
(AGDD_V10=GDDs_Per_Stage_Emerge+GDDs_Per_Stage_VE_to_V10*10)# 940
# AGDD till end of V18
(AGDD_V18=GDDs_Per_Stage_Emerge+GDDs_Per_Stage_VE_to_V10*10 + 
  GDDs_Per_Stage_V10_to_V18*8) #1340
(AGDD_VT=AGDD_V18 + 100) #1440
(AGDD_R1=AGDD_VT + 225)  #1665
(AGDD_R2=AGDD_R1 + 200)  #1865
(AGDD_R3=AGDD_R2 + 300)  #2165
(AGDD_R4=AGDD_R3 + 200)  #2365
(AGDD_R5=AGDD_R4 + 250)  #2615
(AGDD_R6=AGDD_R5 + 250)  #2865

colnames(Weadat_AGDD_eachEnv)
Weadat_AGDD_eachEnv[order(Weadat_AGDD_eachEnv$AGDD),1:8][1:20,]

Env="GEH1_2021"
Weadat_GrowthStages=foreach(Env = unique(Weadat_GDD$Env),.combine = "rbind") %do% {
  # extract weather data for a specific Env
  tmp=Weadat_GDD[Weadat_GDD$Env==Env,]
  dim(tmp) #366  29
  tmp[1:4,]
  tail(tmp)

  # extract weather data for dates between Planted - Harvested
  tmp=tmp[tmp$Date_Juli>=tmp$Date_Planted_Imputed &
            tmp$Date_Juli<=tmp$Date_Harvested_Imputed,]
  dim(tmp)#146  29
  
  # calculate AGDD
  tmp=add_column(.data = tmp,
                 .after = "GDD",
                 AGDD=cumsum(tmp$GDD))
  
  # Estimating Corn Growth Stage with Growing Degree Days
  # dim(tmp)
  tmp[1:10,c("Env","Date_Cale","Date_Planted_Imputed","Date_Harvested_Imputed",
         "Dates_Growth","GDD","AGDD")]
  
  # VE
  tmp$GrowthStage[tmp$AGDD>0 & tmp$AGDD<=GDDs_Per_Stage_Emerge]="VE"
  
  # Determine V1-V10
  tmp$GrowthStage[tmp$AGDD>GDDs_Per_Stage_Emerge & tmp$AGDD<=AGDD_V10]=
    paste0("V",as.integer((tmp$AGDD[tmp$AGDD>GDDs_Per_Stage_Emerge & tmp$AGDD<=AGDD_V10]-GDDs_Per_Stage_Emerge)/GDDs_Per_Stage_VE_to_V10)+1
    )
  
  # Determine V11-V18
  tmp$GrowthStage[tmp$AGDD>AGDD_V10 & tmp$AGDD<=AGDD_V18]=
    paste0("V",as.integer((tmp$AGDD[tmp$AGDD>AGDD_V10 & tmp$AGDD<=AGDD_V18]-AGDD_V10)/GDDs_Per_Stage_V10_to_V18)+11
    )
  
  # VT and Reproductive phase
  tmp$GrowthStage[tmp$AGDD>AGDD_V18 & tmp$AGDD<=AGDD_VT]="VT"
  tmp$GrowthStage[tmp$AGDD>AGDD_VT & tmp$AGDD<=AGDD_R1]="R1"
  tmp$GrowthStage[tmp$AGDD>AGDD_R1 & tmp$AGDD<=AGDD_R2]="R2"
  tmp$GrowthStage[tmp$AGDD>AGDD_R2 & tmp$AGDD<=AGDD_R3]="R3"
  tmp$GrowthStage[tmp$AGDD>AGDD_R3 & tmp$AGDD<=AGDD_R4]="R4"
  tmp$GrowthStage[tmp$AGDD>AGDD_R4 & tmp$AGDD<=AGDD_R5]="R5"
  tmp$GrowthStage[tmp$AGDD>AGDD_R5 & tmp$AGDD<=AGDD_R6]="R6"
  
  # after R6
  if(sum(is.na(tmp$GrowthStage))>0){
    tmp$GrowthStage[is.na(tmp$GrowthStage)]='R6ToHarvest'
  }
  
  # colnames(tmp)
  tmp[,c("Env","Date_Cale","Date_Juli","Date_Planted_Imputed","Date_Harvested_Imputed",
         "Dates_Growth","GDD","AGDD","GrowthStage")]
  
  # remove GrowthStage=NA
  tmp=tmp[!is.na(tmp$GrowthStage),]
  
  # return
  return(tmp)
}

dim(Weadat_GrowthStages)#34953    31
length(unique(Weadat_GrowthStages$Env)) #238
sum(is.na(Weadat_GrowthStages$GrowthStage))
unique(Weadat_GrowthStages$GrowthStage)
head(Weadat_GrowthStages)
tail(Weadat_GrowthStages)
length(unique(Weadat_GDD$Env))#238
sort(unique(Weadat_GrowthStages$Dates_Growth))

# add another GrowthStages to avoid NAs for some Envs
Weadat_GrowthStages=add_column(.data = Weadat_GrowthStages,
                               .after = "GrowthStage",
                               GrowthStage2=Weadat_GrowthStages$GrowthStage)
Weadat_GrowthStages$GrowthStage2[Weadat_GrowthStages$GrowthStage %in% c("R3","R4","R5","R6")]="R3ToR6"
unique(Weadat_GrowthStages$GrowthStage2)

if(!file.exists(paste0(DIR_Output,"Weadat_TRN_TST_with_estimated_GrowthStages.csv"))){
  write.csv(Weadat_GrowthStages,paste0(DIR_Output,"Weadat_TRN_TST_with_estimated_GrowthStages.csv"),
            row.names = F)
}

# how many growth stages for each Env
df_growthstage=NULL; df_envi=NULL
Env = unique(Weadat_GrowthStages$Env)[1]
for(Env in unique(Weadat_GrowthStages$Env)){
  tmp=Weadat_GrowthStages[Weadat_GrowthStages$Env==Env,]
  dim(tmp)
  df_envi=data.frame(GrowthStage=unique(tmp$GrowthStage),Envi=unique(tmp$GrowthStage))
  colnames(df_envi)=c("GrowthStage",Env)
  if(Env == unique(Weadat_GrowthStages$Env)[1]){
    df_growthstage=df_envi
  }else{
    df_growthstage=left_join(df_growthstage,df_envi,by="GrowthStage")
  }
}
df_growthstage[,1:10]
# notes: not every Env has all growth stages
#        it means the model does not fit all Envs well
#        but we could not do a better job

# -------------------------------------------------------------------------------
# Step 3: average weather data according to estimated growth stages (GrowthStage2)
# Convert Weadat to long data
library(reshape2)
colnames(Weadat_GrowthStages)
Weadat_long=NA
Weadat_long = melt(data = Weadat_GrowthStages,
                   id.vars = c("Year","Env","Field_Location","GrowthStage2"),
                   measure.vars = c("QV2M","T2MDEW","PS","RH2M","WS2M","GWETTOP",
                                    "ALLSKY_SFC_SW_DWN","ALLSKY_SFC_PAR_TOT",
                                    "T2M_MAX_C","T2M_MIN_C","GDD","T2MWET","GWETROOT",
                                    "T2M","GWETPROF","ALLSKY_SFC_SW_DNI","PRECTOTCORR") ,
                   variable.name = "Variable",
                   value.name = "Value")
# dropped "T2M_MAX_F","T2M_MIN_F" due to same as "T2M_MAX_C","T2M_MIN_C"
head(Weadat_long)
unique(Weadat_long$GrowthStage)
length(unique(Weadat_long$Env))

# within each Env/Trial, for each climate par, compute mean across daily values for each growth stage
Weadat_long=add_column(.data = Weadat_long,
                       .after = "Variable", 
                       Variable_GrowthStage=paste(Weadat_long$Variable,Weadat_long$GrowthStage2,sep = "::"))
Weadat_StageMean_long=NA
Weadat_StageMean_long=aggregate(Value ~ Year+Env+Field_Location+GrowthStage2+Variable+Variable_GrowthStage,
                                data = Weadat_long, mean)
head(Weadat_StageMean_long)
dim(Weadat_long)          #648992      7
dim(Weadat_StageMean_long)#102976      7(12/25); 92669     7(updated 12/26)

unique(Weadat_StageMean_long$GrowthStage2)
length(unique(Weadat_StageMean_long$Env)) #238

# change long data to wide data (i.e. each Env in one row, all featrues in columns)
Weadat_StageMean_wide = dcast(Weadat_StageMean_long,
                              Year+Env+Field_Location~Variable_GrowthStage,
                              value.var = "Value")

head(Weadat_StageMean_wide)
Weadat_StageMean_wide[1:3,1:10]
dim(Weadat_StageMean_wide) #238 462; 238 411
head(Weadat_StageMean_wide)
if(!file.exists(paste0(DIR_Output,"G2F_Weadat_Stage2Mean_Wide_AllEnvs_AllWeaVars.csv"))){
  write.csv(Weadat_StageMean_wide,paste0(DIR_Output,"G2F_Weadat_Stage2Mean_Wide_AllEnvs_AllWeaVars.csv"),
            row.names = F)
}

# boxplot of Weadat_StageMean_wide
boxplot(Weadat_StageMean_wide[,-1*1:3],las=2,outcex=0.2,outpch=20)
# => it shows diff weather variables have different scales
varnames=as.character(unique(Weadat_StageMean_long$Variable))
boxplot(Weadat_StageMean_wide[,-1*1:3][,grepl(varnames[1],colnames(Weadat_StageMean_wide[,-1*1:3]))],
        las=2,outcex=0.2,outpch=20)

# ------------------------------------------------------------------------------
# Quality filtering
# list.files(DIR_Output)
# Weadat_StageMean_wide=read.csv(paste0(DIR_Output,"G2F_Weadat_Stage2Mean_Wide_AllEnvs_AllWeaVars.csv"),
#                                header = T,check.names = F,stringsAsFactors = F)
dim(Weadat_StageMean_wide) #238 411
Weadat_StageMean_wide[1:3,1:10]
length(unique(sapply(strsplit(colnames(Weadat_StageMean_wide)[-1*1:3],"::"), "[",2)))#24
table(sapply(strsplit(colnames(Weadat_StageMean_wide)[-1*1:3],"::"), "[",2))
unique(sapply(strsplit(colnames(Weadat_StageMean_wide)[-1*1:3],"::"), "[",1))

# Step 1: remove Envs from Germany (since no Envs from Germany in G2F competition test set)
German_Envs=Weadat_StageMean_wide$Env[grepl("GE",Weadat_StageMean_wide$Env)]#"GEH1_2019" "GEH1_2020" "GEH1_2021"
Weadat_StageMean_wide_QC=Weadat_StageMean_wide[!(Weadat_StageMean_wide$Env %in% German_Envs),]
dim(Weadat_StageMean_wide_QC) #235 411


# Step 2: remove the five weather variables that wea_tst has >80% missing data points
colnames(wea_tst)[apply(wea_tst, 2, function(x) sum(is.na(x))/nrow(wea_tst))>0.8]
# "GWETTOP" "ALLSKY_SFC_PAR_TOT" "GWETROOT" "GWETPROF" "ALLSKY_SFC_SW_DNI"
Weadat_StageMean_wide_QC=
  Weadat_StageMean_wide_QC[,!grepl("GWETTOP|ALLSKY_SFC_PAR_TOT|GWETROOT|GWETPROF|ALLSKY_SFC_SW_DNI",
                               colnames(Weadat_StageMean_wide))]
dim(Weadat_StageMean_wide_QC)#235 291

# Step 3: remove GrowthStage='R6ToHarvest'
Weadat_StageMean_wide_QC=Weadat_StageMean_wide_QC[,!grepl("R6ToHarvest",colnames(Weadat_StageMean_wide_QC))]
dim(Weadat_StageMean_wide_QC)#235 279
Weadat_StageMean_wide_QC[1:3,1:10]
sort(unique(sapply(strsplit(colnames(Weadat_StageMean_wide_QC)[-1*1:3],"::",2), "[", 2)))
all(apply(Weadat_StageMean_wide_QC, 2, function(x) sum(is.na(x)))==0) #TRUE

# Stpe 4: add AGDD and APRE to Wea_Dat
Wea_Dat = add_column(.data = Weadat_StageMean_wide_QC,
                     .after = "Field_Location",
                     AGDD=Weadat_APRE_eachEnv$AGDD[match(Weadat_StageMean_wide_QC$Env,Weadat_APRE_eachEnv$Env)],
                     APRE=Weadat_APRE_eachEnv$APRE[match(Weadat_StageMean_wide_QC$Env,Weadat_APRE_eachEnv$Env)]
)

# output Weadat_StageMean_wide_QC
dim(Weadat_StageMean_wide_QC) #235 279
if(!file.exists(paste0(DIR_Output,"G2F_Weadat_Stage2Mean_Wide_AllEnvs_QC.csv"))){
  write.csv(Weadat_StageMean_wide_QC,paste0(DIR_Output,"G2F_Weadat_Stage2Mean_Wide_AllEnvs_QC.csv"),
            row.names = F)
}
