# ==============================================================
# script11_Prepare_MegaLMM_BaseModel_Input
# ==============================================================

# 2023-12-08
# Prepare_MegaLMM_BaseModel_Input

# Clean workspace
rm(list=ls())

library(rrBLUP)
library(MegaLMM); K=NULL
library(data.table)
library(lineup)
library(cvTools)
library(foreach)
library(rrBLUP)
library(dplyr)
library(tibble)


# set up DIR_Output
Date="2023-12-08"
# DIR_Input="../../P5G2Fv2/04_TRN_TST_Dat_EnvNewV2/2023-03-27/"
DIR_Input="../05_MegaLMM_PhenoDat_Format2/"
DIR_EnvDat="../../P5G2F/04_Env_Dat/2022-12-26/"
DIR_Output=paste0("../06_MegaLMM_PhenoDat_TRN_TST/")
if(!dir.exists(DIR_Output)) dir.create(DIR_Output)

# define global variables
traitnames = rev(c("Silk_DAP_days","Plant_Height_cm","Yield_Mg_ha"))
numEnvNew_PerState_cutoff=9 #min number of trails must be included in a testing State
YearGroups_CV4 = c("2014,2015","2016,2017","2018,2019","2020,2021")
numFolds=5

#-------------------------------------------------------------------------------
# load Pheno dat
list.files(DIR_Input)
ls=load(sprintf("%s%sMegaLMM_PhenoDat_Format2.RData",DIR_Input,"2023-09-15"))
ls #"Phe_Trn_Fmt1Wide" "Phe_Trn_Fmt2Wide" "Phe_Trn_Fmt2Long" "K_info"
lapply(Phe_Trn_Fmt2Wide, function(x) x[1:4,1:6])
sapply(Phe_Trn_Fmt2Wide, dim)
#     Yield_Mg_ha Plant_Height_cm Silk_DAP_days
# [1,]        1702            1699          1663
# [2,]         302             278           231
# assign phe_trn_wide_fmt2 to Y

#-------------------------------------------------------------------------------
# Create data partition for G2F different CV Scenarios
#-------------------------------------------------------------------------------

# Create empty list to save Traing and Testing Pheno data
Y_Train_CV1=list(); Y_Test_CV1=list()
Y_Train_CV2=list(); Y_Test_CV2=list()
Y_Train_CV3NewTrial=list(); Y_Test_CV3NewTrial=list()
Y_Train_CV3NewState=list(); Y_Test_CV3NewState=list()
Y_Train_CV3NewTester=list(); Y_Test_CV3NewTester=list()
Y_Train_CV4=list(); Y_Test_CV4=list()

# for script testing purpose only
traitname = traitnames[1]
ID_fold = 1

# TRN and TST partition with for loops
for(traitname in traitnames){
  # Create empty list for each trait
  Y_Train_CV1[[traitname]]=list();  Y_Train_CV2[[traitname]]=list()
  Y_Test_CV1[[traitname]]=list();  Y_Test_CV2[[traitname]]=list()
  Y_Train_CV3NewTrial[[traitname]]=list(); Y_Test_CV3NewTrial[[traitname]]=list()
  Y_Train_CV3NewState[[traitname]]=list(); Y_Test_CV3NewState[[traitname]]=list()
  Y_Train_CV3NewTester[[traitname]]=list(); Y_Test_CV3NewTester[[traitname]]=list()
  Y_Train_CV4[[traitname]]=list(); Y_Test_CV4[[traitname]]=list()
  
  # Run five-fold partition for genotypes within each EnvNew
  for(ID_fold in 1:numFolds){
    #-------------------------------------------------------------------------------
    # get Y matrix for a specific trait
    Y = Phe_Trn_Fmt2Wide[[traitname]]
    dim(Y) #1702  302
    Y[1:4,1:6]
    sort(rowSums(!is.na(Y)))[1:10]
    sum(duplicated(rownames(Y)))
    
    #-------------------------------------------------------------------------------
    # step 0: make total number of Parent1 lines (in rows) of Y dividable by 5
    #  - if the total number of lines not dividable by 5 
    #  - remove extra lines with highest missing rate
    #  - this makes shuffle genotype groups for each EnvNew much esier in CV2
    if(nrow(Y) %% 5 != 0){
      num_extra_lines = nrow(Y) %% 5
      
      # calculate missing rate
      missrate=rowMeans(is.na(Y))
      # table(rowSums(is.na(Y)))
      # note: there might be hundreds of lines with the highest missing rate
      #       in this case, will randomly remove extra lines by R
      # select top num_extra_lines lines with highest missing rate
      id_sel = names(missrate) %in% names(sort(missrate,decreasing=T)[1:num_extra_lines])
      rowMeans(is.na(Y[id_sel,]))
      # delete extra lines
      Y = Y[!id_sel,]
    }
    dim(Y)#1700  302
    # checking
    Y[1:4,1:6]
    stopifnot(nrow(Y) %% 5 == 0)
    
    #---------------------------------------------------------------------------
    # Step 1: create CV1 Five-Folds
    # each fold contains the same set of Parent1 genotypes across EnvNew
    set.seed(12386) #set seeds to make results repeatable for later on running
    CV1Folds=cvFolds(n=nrow(Y), K=numFolds, R=1)
    # - n: an integer giving the number of observations to be split into groups
    # - K: an integer giving the number of groups
    # - R: an integer giving the number of replications for repeated KK-fold cross-validation. 
    sapply(CV1Folds, length)
    sapply(CV1Folds, head)
    head(CV1Folds$subsets[,1][CV1Folds$which==1]) #rowID for each CV1 fold
    
    # get CV1 train/test
    Y_train_CV1=Y_test_CV1=Y
    Y_train_CV1[CV1Folds$subsets[,1][CV1Folds$which==ID_fold],]=NA
    Y_test_CV1[CV1Folds$subsets[,1][CV1Folds$which!=ID_fold],]=NA
    stopifnot(all(colSums(!is.na(Y_train_CV1))+colSums(!is.na(Y_test_CV1))==colSums(!is.na(Y))))
    
    # save CV1 train/test
    Y_Train_CV1[[traitname]][[paste0("Fold_",ID_fold)]]=Y_train_CV1
    Y_Test_CV1[[traitname]][[paste0("Fold_",ID_fold)]]=Y_test_CV1
    
    #---------------------------------------------------------------------------
    # Step 2: sort Y's rows according to CV1 Fold ID
    #  - I extract data blocks for each CV1 Fold (i.e. Parent1 groups)
    #  - then stack them one by one according to CV1 Fold ID (i.e. 1 to 5)
    #  - so the line order in CV2b correspond to CV1 fold 111...,222...,333...,444...,555...,
    #  - however, the line order in CV1 is according to CV1Folds$which, which is random
    Y_Sorted_by_CV1Fold=foreach(i = 1:numFolds,.combine = "rbind") %do% {
      print(paste("############",i,"###################"))
      # get rowID for each CV1 fold
      rowID_fold_i=CV1Folds$subsets[,1][CV1Folds$which==i]
      # return Y matrix for fold i
      return(Y[rowID_fold_i,])
    }
    dim(Y_Sorted_by_CV1Fold)==dim(Y)
    Y_Sorted_by_CV1Fold[1:3,1:4]
    stopifnot(identical(match(rownames(Y_Sorted_by_CV1Fold),rownames(Y)),
                        CV1Folds$subsets[,1][order(CV1Folds$which)]))
    # note: the original rowID in CV1Folds$subsets[,1] is not sorted by Fold ID 
    #       but in random order
    
    #---------------------------------------------------------------------------
    # step 3: create ID_fold_matrix for CV2b
    # - because all Parent1s in each EnvNew have been sorted by CV1 Fold ID (increasing=T)
    # - Then, within each EnvNew, I only need to randomize all Parant1s in Fold level (CV2b),
    #   rather than in the indivial level (which is CV2)
    # - For CV, within each EnvNew, each of the numFolds (here 5) fold will have a chance
    #   to be used as testing set, therefore, we going to get a prediction accuracy value
    #   for each fold in each EnvNew
    # - if we put the prediction accuracy into a matrix with Folds in row and EnvNew in Column
    #   for both CV1, CV2b and CV3,
    # - they are going to be comparable, because they all based on same genotype partition
    #   in each EnvNew, i.e. in ech EnvNew, Fold1, Fold2, Fold3, Fold4, and Fold5 contain 
    #   the same set of Parent1s across different cross-validtion scenarios, here CV1,CV2b,
    #   and CV3
    # create an empty CV2b ID_fold_matrix
    ID_fold_matrix = matrix(NA,
                            nrow = nrow(Y),
                            ncol = ncol(Y),
                            dimnames = dimnames(Y_Sorted_by_CV1Fold))
    for(i in 1:ncol(ID_fold_matrix)) {
      if(all(table(CV1Folds$which)==nrow(Y)/numFolds)){ #nrow of Y should be dividable by numFolds
        # keep fold ID consistent with Y_Sorted_by_CV1Fold
        vector = sort(CV1Folds$which) 
        table(sort(CV1Folds$which) )
        #we need to set known seed for each column to make generated results repeatable
        set.seed(i) 
        # randomize sample CV1 FoldID for each EnvNew
        (shuffled_group <- sample(unique(vector)) )
        # get randomized sample CV1 FoldID for each env
        ID_fold_matrix[,i]=rep(shuffled_group, each = nrow(Y)/numFolds)
      }
    }
    # taking a look at the ID_fold_matrix
    ID_fold_matrix[1:3,1:4]
    dim(ID_fold_matrix)
    Image(ID_fold_matrix)
    
    # checking: I expect there is no NA in the ID_fold_matrix
    stopifnot(all(!is.na(ID_fold_matrix)))
    
    # # output the CV2b ID_fold_matrix
    # write.csv(ID_fold_matrix,paste0(DIR_Output,"plots/","CV2b_ID_fold_matrix.csv"))
    
    # testing: to see whether there is all NAs in either TRN or TST for any EnvNew
    for(i in 1:numFolds){
      Y_train = Y_test = Y_Sorted_by_CV1Fold
      Y_train[ID_fold_matrix == i] = NA
      Y_test[ID_fold_matrix != i | is.na(ID_fold_matrix)] = NA
      print(paste("##############i=",i,"##############"))
      print(which(colSums(!is.na(Y_train))==0))
      which(colSums(!is.na(Y_test))==0)
      which(colSums(!is.na(Y))==0)
      print(all(colSums(!is.na(Y_train))+colSums(!is.na(Y_test))==colSums(!is.na(Y))))
      stopifnot(any(colSums(!is.na(Y_train))!=0))
    }
    
    #---------------------------------------------------------------------------
    # step 4: create CV2b TRN amd TST sets and save
    Y_train_CV2 = Y_test_CV2 = Y_Sorted_by_CV1Fold
    Y_train_CV2[ID_fold_matrix == ID_fold] = NA
    Y_test_CV2[ID_fold_matrix != ID_fold | is.na(ID_fold_matrix)] = NA
    stopifnot(all(colSums(!is.na(Y_train_CV2))+colSums(!is.na(Y_test_CV2))==colSums(!is.na(Y))))
    
    # save CV2b TRN amd TST sets
    Y_Train_CV2[[traitname]][[paste0("Fold_",ID_fold)]]=Y_train_CV2
    Y_Test_CV2[[traitname]][[paste0("Fold_",ID_fold)]]=Y_test_CV2
    
    #---------------------------------------------------------------------------
    # Step 5: get EnvNew_info for State, Trial and Hybrid_Parent2
    EnvNew_info = data.frame(EnvNew=colnames(Y_train_CV2))
    EnvNew_info$Year=paste0("20",sapply(strsplit(sapply(strsplit(EnvNew_info$EnvNew,"::"),"[",1),"_20"),
                                        "[",2))
    # 2023-10-05: use 3 chars to avoid issues between "MO" and "OH" when using grepl()
    EnvNew_info$State=substr(EnvNew_info$EnvNew, 1, 2) 
    EnvNew_info$Trial=sapply(strsplit(EnvNew_info$EnvNew,"::"),"[",1)
    EnvNew_info$Hybrid_Parent2=sapply(strsplit(EnvNew_info$EnvNew,"::"),"[",2)
    head(EnvNew_info)
    length(unique(EnvNew_info$Trial)) #195
    length(unique(EnvNew_info$State)) #20
    length(unique(EnvNew_info$Hybrid_Parent2)) #12
    table(EnvNew_info$State)
    # "Yield_Mg_ha" (14)
    # AR DE GA GE IA IL IN KS MI MN MO NC NE NY OH ON SC SD TX WI 
    #  5 17 15  3 52  9 10  6  6 11 24 15 26 26  9 11  4  1 24 28 
    myvect=c(5,17,15,3,52,9,10,6,6,11,24,15,26,26, 9,11,4,1,24,28)
    sum(myvect)#302
    sum(myvect[myvect>=numEnvNew_PerState_cutoff])#277 (277-9=268)
    # "Plant_Height_cm" (13)
    # AR DE GA GE IA IL IN KS MI MN MO NC NE NY OH ON SC SD TX WI 
    #  6 17 15  3 43  9 10  1  5 11 24 15 20 26  7 11  2  1 24 28 
    myvect=c(6,17,15,3,43, 9,10, 1, 5,11,24,15,20,26, 7,11, 2, 1,24,28)
    sum(myvect[myvect>=numEnvNew_PerState_cutoff])#253
    # "Silk_DAP_days" (12) 
    # AR DE GA GE IA IL IN KS MI MN MO NC NE NY OH ON SC SD TX WI 
    #  6 17 15  2 17  9 10  5  6 11 24  4 14 22  7 11  4  1 24 22 
    myvect=c(6,17,15, 2,17, 9,10, 5, 6,11,24, 4,14,22, 7,11, 4, 1,24,22 )
    sum(myvect[myvect>=numEnvNew_PerState_cutoff])
    sum(table(EnvNew_info$State)>=numEnvNew_PerState_cutoff)#196
    length(table(EnvNew_info$State))
    
    table(EnvNew_info$Hybrid_Parent2)
    # "Yield_Mg_ha"
    # CG102 DK3IIH6   LH185   LH195   LH198    LH82    PB80   PHB47   PHK76   PHP02   PHT69   PHZ51 
    #    4      31      12      66      12      26       6      24      18      22      25      56 
    # "Plant_Height_cm"
    # CG102 DK3IIH6   LH185   LH195   LH198    LH82    PB80   PHB47   PHK76   PHP02   PHT69   PHZ51 
    #    4      28      12      64      12      22       6      21      15      20      25      49 
    # "Silk_DAP_days"
    # CG102 DK3IIH6   LH185   LH195   LH198    LH82    PB80   PHB47   PHK76   PHP02   PHT69   PHZ51 
    #    3      26       8      54       8      17       5      15      12      19      25      39 
    #---------------------------------------------------------------------------
    # Step 6: Create TRN and TST for CV3NewTrial
    set.seed(12306) #set seeds to make results repeatable for later on running
    CV3NewTrialFolds=cvFolds(n= length(unique(EnvNew_info$Trial)), K=numFolds, R=1)
    # CV3NewTrialFolds=cvFolds(n= ncol(Y_train_CV2), K=numFolds, R=1)
    sapply(CV3NewTrialFolds, length)
    sapply(CV3NewTrialFolds, head)
    table(CV3NewTrialFolds$which)
    # 1  2  3  4  5 #Trials
    # 39 39 39 39 39
    # 1  2  3  4  5 #EnvNew
    # 61 61 60 60 60
    head(CV3NewTrialFolds$subsets[,1][CV3NewTrialFolds$which==1]) #rowID for each CV1 fold
    
    # FoldID_CV3NewTrial = 1
    for(FoldID_CV3NewTrial in 1:numFolds){
      # empty Train/Test tmp
      Y_Test_CV3NewTrial_tmp = NULL; Y_Train_CV3NewTrial_tmp = NULL
      # get Trial Names for Fold i
      if(CV3NewTrialFolds$n==length(unique(EnvNew_info$Trial))) {
        # TrialNames_tmp = unique(EnvNew_info$Trial)[CV3NewTrialFolds$subsets[,1][CV3NewTrialFolds$which==ID_fold]]
        TrialNames_tmp = unique(EnvNew_info$Trial)[CV3NewTrialFolds$subsets[,1][CV3NewTrialFolds$which==FoldID_CV3NewTrial]]
      }
      # get Train/Test
      if(identical(EnvNew_info$EnvNew,colnames(Y_train_CV2))){
        Y_Test_CV3NewTrial_tmp=Y_train_CV2[,colnames(Y_train_CV2) %in% EnvNew_info$EnvNew[EnvNew_info$Trial %in% TrialNames_tmp]]
        Y_Train_CV3NewTrial_tmp=Y_train_CV2[,!(colnames(Y_train_CV2) %in% EnvNew_info$EnvNew[EnvNew_info$Trial %in% TrialNames_tmp])]
      }
      # save Train/Test
      Y_Test_CV3NewTrial[[traitname]][[paste0("CV2_FoldID_",ID_fold)]][[FoldID_CV3NewTrial]]=
        Y_Test_CV3NewTrial_tmp
      Y_Train_CV3NewTrial[[traitname]][[paste0("CV2_FoldID_",ID_fold)]][[FoldID_CV3NewTrial]]=
        Y_Train_CV3NewTrial_tmp
      # validating Train/Test
      stopifnot(identical(setdiff(c(colnames(Y_Test_CV3NewTrial_tmp),colnames(Y_Train_CV3NewTrial_tmp)),
                                  colnames(Y_train_CV2)),character(0)))
      # all Testing set TrialNames for Fold i should be included in the colnames of Y_Test 
      stopifnot(all(TrialNames_tmp %in% sapply(strsplit(colnames(Y_Test_CV3NewTrial_tmp),"::"), "[",1)))
      # None of Testing set TrialNames for Fold i should be included in the colnames of Y_Test 
      stopifnot(all(!(TrialNames_tmp %in% sapply(strsplit(colnames(Y_Train_CV3NewTrial_tmp),"::"), "[",1))))
    }
    
    #---------------------------------------------------------------------------
    # Step 7a: Create TRN and TST for CV3NewState - PART 1
    Y_Test_CV3NewState_P1 = list()
    Y_Train_CV3NewState_P1 = list()
    (StateSel_P1=names(table(EnvNew_info$State)[table(EnvNew_info$State)>=numEnvNew_PerState_cutoff]))
    # FoldID_CV3NewState = 1
    for(FoldID_CV3NewState in 1:length(StateSel_P1)){
      # empty Train/Test tmp
      Y_Test_CV3NewState_tmp_P1=NULL; Y_Train_CV3NewState_tmp_P1=NULL
      
      # get Train/Test P1
      Y_Test_CV3NewState_tmp_P1=Y_train_CV2[,substr(colnames(Y_train_CV2),1,2)==StateSel_P1[FoldID_CV3NewState],drop=F]
      Y_Train_CV3NewState_tmp_P1=Y_train_CV2[,substr(colnames(Y_train_CV2),1,2)!=StateSel_P1[FoldID_CV3NewState]]
      
      # save Train/Test
      Y_Test_CV3NewState_P1[[StateSel_P1[FoldID_CV3NewState]]] = Y_Test_CV3NewState_tmp_P1
      Y_Train_CV3NewState_P1[[StateSel_P1[FoldID_CV3NewState]]] = Y_Train_CV3NewState_tmp_P1
      
      # Y_Test_CV3NewState[[traitname]][[paste0("CV2_FoldID_",ID_fold)]][[StateSel[FoldID_CV3NewState]]]=
      #   Y_Test_CV3NewState_tmp
      # Y_Train_CV3NewState[[traitname]][[paste0("CV2_FoldID_",ID_fold)]][[StateSel[FoldID_CV3NewState]]]=
      #   Y_Train_CV3NewState_tmp
      # validating Train/Test
      stopifnot(identical(setdiff(c(colnames(Y_Test_CV3NewState_tmp_P1),colnames(Y_Train_CV3NewState_tmp_P1)),
                                  colnames(Y_train_CV2)),character(0)))
      stopifnot(all(substr(colnames(Y_Test_CV3NewState_tmp_P1),1,2) == StateSel_P1[FoldID_CV3NewState]))
      
    }
    
    #---------------------------------------------------------------------------
    # Step 7b: Create TRN and TST for CV3NewState - PART 2
    Y_Test_CV3NewState_P2 = list()
    Y_Train_CV3NewState_P2 = list()
    (StateSel_P2=names(table(EnvNew_info$State)[table(EnvNew_info$State)<numEnvNew_PerState_cutoff]))
    # FoldID_CV3NewState = 1
    for(FoldID_CV3NewState in 1:length(StateSel_P2)){
      # empty Train/Test tmp
      Y_Test_CV3NewState_tmp_P2=NULL; Y_Train_CV3NewState_tmp_P2=NULL
      # get Train/Test P2
      Y_Test_CV3NewState_tmp_P2=Y_train_CV2[,substr(colnames(Y_train_CV2),1,2)==StateSel_P2[FoldID_CV3NewState],drop=F]
      Y_Train_CV3NewState_tmp_P2=Y_train_CV2[,substr(colnames(Y_train_CV2),1,2)!=StateSel_P2[FoldID_CV3NewState]]
      # save Train/Test
      # save Train/Test
      Y_Test_CV3NewState_P2[[StateSel_P2[FoldID_CV3NewState]]] = Y_Test_CV3NewState_tmp_P2
      Y_Train_CV3NewState_P2[[StateSel_P2[FoldID_CV3NewState]]] = Y_Train_CV3NewState_tmp_P2
      
      # Y_Test_CV3NewState[[traitname]][[paste0("CV2_FoldID_",ID_fold)]][[StateSel[FoldID_CV3NewState]]]=
      #   Y_Test_CV3NewState_tmp
      # Y_Train_CV3NewState[[traitname]][[paste0("CV2_FoldID_",ID_fold)]][[StateSel[FoldID_CV3NewState]]]=
      #   Y_Train_CV3NewState_tmp
      # validating Train/Test
      stopifnot(identical(setdiff(c(colnames(Y_Test_CV3NewState_tmp_P2),colnames(Y_Train_CV3NewState_tmp_P2)),
                                  colnames(Y_train_CV2)),character(0)))
      stopifnot(all(substr(colnames(Y_Test_CV3NewState_tmp_P2),1,2) == StateSel_P2[FoldID_CV3NewState]))
      
    }
    
    #---------------------------------------------------------------------------
    # Step 7c: CV3NewState P1 and P2
    sapply(list(Y_Test_CV3NewState_P1,Y_Train_CV3NewState_P1, 
                Y_Test_CV3NewState_P2,Y_Train_CV3NewState_P2), length)
    sapply(list(Y_Test_CV3NewState_P1,Y_Train_CV3NewState_P1, 
                Y_Test_CV3NewState_P2,Y_Train_CV3NewState_P2), names)
    sapply(list(Y_Test_CV3NewState_P1,Y_Train_CV3NewState_P1, 
                Y_Test_CV3NewState_P2,Y_Train_CV3NewState_P2),
           function(x) dim(x[[1]]))
    stopifnot(identical(setdiff(c(names(Y_Test_CV3NewState_P1),names(Y_Test_CV3NewState_P2)),
                                unique(EnvNew_info$State)),character(0)))
    stopifnot(identical(setdiff(c(names(Y_Train_CV3NewState_P1),names(Y_Train_CV3NewState_P2)),
                                unique(EnvNew_info$State)),character(0)))
    # merge P1 and P2
    Y_Test_CV3NewState_ALL=append(Y_Test_CV3NewState_P1,Y_Test_CV3NewState_P2)
    Y_Train_CV3NewState_ALL=append(Y_Train_CV3NewState_P1,Y_Train_CV3NewState_P2)
    sapply(list(Y_Test_CV3NewState_ALL,Y_Train_CV3NewState_ALL), length)
    sapply(list(Y_Test_CV3NewState_ALL,Y_Train_CV3NewState_ALL), names)
    sapply(list(Y_Test_CV3NewState_ALL,Y_Train_CV3NewState_ALL), function(x) dim(x[[1]]))
    sapply(list(Y_Test_CV3NewState_ALL,Y_Train_CV3NewState_ALL), function(x) dim(x[[15]]))
    # assignment
    Y_Test_CV3NewState[[traitname]][[paste0("CV2_FoldID_",ID_fold)]]=Y_Test_CV3NewState_ALL
    Y_Train_CV3NewState[[traitname]][[paste0("CV2_FoldID_",ID_fold)]]=Y_Train_CV3NewState_ALL
 
    #---------------------------------------------------------------------------
    # Step 8: Create TRN and TST for CV3NewTester
    TesterSel=unique(EnvNew_info$Hybrid_Parent2)
    for(FoldID_CV3NewTester in 1:length(TesterSel)){
      # empty Train/Test tmp
      Y_Test_CV3NewTester_tmp=NULL; Y_Train_CV3NewTester_tmp=NULL
      # get Train/Test
      Y_Test_CV3NewTester_tmp=Y_train_CV2[,grepl(TesterSel[FoldID_CV3NewTester],colnames(Y_train_CV2))]
      Y_Train_CV3NewTester_tmp=Y_train_CV2[,!grepl(TesterSel[FoldID_CV3NewTester],colnames(Y_train_CV2))]
      # save Train/Test
      Y_Test_CV3NewTester[[traitname]][[paste0("CV2_FoldID_",ID_fold)]][[TesterSel[FoldID_CV3NewTester]]]=
        Y_Test_CV3NewTester_tmp
      Y_Train_CV3NewTester[[traitname]][[paste0("CV2_FoldID_",ID_fold)]][[TesterSel[FoldID_CV3NewTester]]]=
        Y_Train_CV3NewTester_tmp
      # validating Train/Test
      stopifnot(identical(setdiff(c(colnames(Y_Test_CV3NewTester_tmp),colnames(Y_Train_CV3NewTester_tmp)),
                                  colnames(Y_train_CV2)),character(0)))
    }
    
    #---------------------------------------------------------------------------
    # Step 9: get Line_info
    # the Line_info should be put somewhere after Step 0 where Y has been updated
    Line_info = NULL
    Line_info = data.frame(LineName = rownames(Y))
    Line_info$YR2014=FALSE
    Line_info$YR2014[rowSums(!is.na(Y[grepl("2014",colnames(Y))]))!=0]=TRUE
    Line_info$YR2015=FALSE
    Line_info$YR2015[rowSums(!is.na(Y[grepl("2015",colnames(Y))]))!=0]=TRUE
    Line_info$YR2016=FALSE
    Line_info$YR2016[rowSums(!is.na(Y[grepl("2016",colnames(Y))]))!=0]=TRUE
    Line_info$YR2017=FALSE
    Line_info$YR2017[rowSums(!is.na(Y[grepl("2017",colnames(Y))]))!=0]=TRUE
    Line_info$YR2018=FALSE
    Line_info$YR2018[rowSums(!is.na(Y[grepl("2018",colnames(Y))]))!=0]=TRUE
    Line_info$YR2019=FALSE
    Line_info$YR2019[rowSums(!is.na(Y[grepl("2019",colnames(Y))]))!=0]=TRUE
    Line_info$YR2020=FALSE
    Line_info$YR2020[rowSums(!is.na(Y[grepl("2020",colnames(Y))]))!=0]=TRUE
    Line_info$YR2021=FALSE
    Line_info$YR2021[rowSums(!is.na(Y[grepl("2021",colnames(Y))]))!=0]=TRUE
    #modify Line_info
    Line_info=column_to_rownames(.data = Line_info,var = "LineName")
    colnames(Line_info)=gsub("YR","",colnames(Line_info))
    # sort Line_info
    Line_info = Line_info[rownames(Y_train_CV2),]
    head(Line_info)
    dim(Line_info)
    #---------------------------------------------------------------------------
    # Step 10: Create TRN and TST for CV4
    head(EnvNew_info)
    
    # # get CV3 train/test
    # Y_Train_CV3NewTrial[[traitname]]=list(); Y_Test_CV3NewTrial[[traitname]]=list()
    # Y_Train_CV3NewState[[traitname]]=list(); Y_Test_CV3NewState[[traitname]]=list()
    # Y_Train_CV3NewTester[[traitname]]=list(); Y_Test_CV3NewTester[[traitname]]=list()
    # Y_Train_CV4[[traitname]]=list(); Y_Test_CV4[[traitname]]=list()
    
    # YearGroups_CV4 = c("2014,2015","2016,2017","2018,2019","2020,2021")
    # FoldID_CV4 = 1
    # head(Line_info)
    
    for(FoldID_CV4 in 1:length(YearGroups_CV4)){
      # empty Train/Test tmp
      Y_Test_CV4_tmp=NULL; Y_Train_CV4_tmp=NULL
      # extract years for a test set
      years_test = YearGroups_CV4[FoldID_CV4]
      years_test = unlist(strsplit(years_test,","))
      # subset EnvNew and Lines for Train/Test
      EnvNew_Test  = EnvNew_info$Year %in% years_test
      LineNames_testyears=rownames(Line_info)[rowSums(Line_info[,years_test])>0]
      LineNames_nontestyears=rownames(Line_info)[rowSums(Line_info[,!(colnames(Line_info) %in% years_test)])>0]
      
      # get Train/Test
      if(identical(rownames(Line_info),rownames(Y_train_CV2))){
        Y_Test_CV4_tmp=Y_train_CV2[(rownames(Y_train_CV2) %in% LineNames_testyears) &
                                     (!(rownames(Y_train_CV2) %in% LineNames_nontestyears)),
                                   EnvNew_Test]
        Y_Train_CV4_tmp=Y_train_CV2[(rownames(Y_train_CV2) %in% LineNames_nontestyears),
                                    !EnvNew_Test]
      }
      
      # save Train/Test
      Y_Test_CV4[[traitname]][[paste0("CV2_FoldID_",ID_fold)]][[paste(years_test,collapse = "and")]]=
        Y_Test_CV4_tmp
      Y_Train_CV4[[traitname]][[paste0("CV2_FoldID_",ID_fold)]][[paste(years_test,collapse = "and")]]=
        Y_Train_CV4_tmp
      # validating Train/Test
      stopifnot(identical(setdiff(c(colnames(Y_Test_CV4_tmp),colnames(Y_Train_CV4_tmp)),
                                  colnames(Y_train_CV2)),character(0)))
      stopifnot(identical(setdiff(c(rownames(Y_Test_CV4_tmp),rownames(Y_Train_CV4_tmp)),
                                  rownames(Y_train_CV2)),character(0)))
      stopifnot(identical(
        setdiff(c(rownames(Y_Test_CV4_tmp),rownames(Y_Train_CV4_tmp)),rownames(Y_train_CV2)),
        character(0) ))
      stopifnot(identical(
        setdiff(c(colnames(Y_Test_CV4_tmp),colnames(Y_Train_CV4_tmp)),colnames(Y_train_CV2)),
        character(0) ))
      } #end of for(FoldID_CV4)
    } # end of for(ID_fold)
} # end of for(traitname)

sapply(list(Y_Train_CV1[[traitnames[1]]],
            Y_Test_CV1[[traitnames[1]]],
            Y_Train_CV2[[traitnames[1]]],
            Y_Test_CV2[[traitnames[1]]]
            ), function(x) dim(x[[1]]))

sapply(list(Y_Train_CV1,
            Y_Test_CV1,
            Y_Train_CV2,
            Y_Test_CV2), names )
sapply(list(Y_Train_CV1[[traitnames[1]]],
            Y_Test_CV1[[traitnames[1]]],
            Y_Train_CV2[[traitnames[1]]],
            Y_Test_CV2[[traitnames[1]]]), names )

# save MegaLMM_TRN_TST
save(Y_Train_CV1,Y_Test_CV1,Y_Train_CV2,Y_Test_CV2,
     Y_Train_CV3NewTrial, Y_Test_CV3NewTrial,
     Y_Train_CV3NewState, Y_Test_CV3NewState,
     Y_Train_CV3NewTester, Y_Test_CV3NewTester,
     Y_Train_CV4, Y_Test_CV4,
     file = sprintf("%s%sMegaLMM_TRN_TST.RData",DIR_Output,Date)
     )
