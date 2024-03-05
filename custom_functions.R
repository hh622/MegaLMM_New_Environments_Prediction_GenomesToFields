

# create function for Fisher's Z transformation
R_to_Z = function(r) log((1+r)/(1-r))/2

# get factor levels from a data frame
factlevels=function(x){
  return(sapply(list(TraitName=x$TraitName,
                     CVScn=x$CVScn,
                     ID_Fold=x$ID_Fold,
                     ID_Fold_Level2=x$ID_Fold_Level2,
                     Method=x$Method,
                     PriorType=x$PriorType,
                     PriorLevel=x$PriorLevel,
                     PredictionType=x$PredictionType,
                     PredictionModel=x$PredictionModel),
                    function(x) unique(x))
  )
}

# create function - average by state
average_by_state = function(Y_test,prediction_matrix,minTrnEnvs=1) {
  res=sapply(colnames(Y_test),function(s) {
    (state = substr(s,start = 1, stop = 2))
    # if(length(grep(state,colnames(prediction_matrix)) )>=minTrnEnvs)
    if( sum(substr(colnames(prediction_matrix),start=1,stop=2)==state) >= minTrnEnvs)
      # return(rowMeans(prediction_matrix[,grep(state,colnames(prediction_matrix)),drop=FALSE],na.rm=T))
      return(rowMeans(prediction_matrix[,substr(colnames(prediction_matrix),start=1,stop=2)==state,drop=FALSE],na.rm=T))
    return(rep(NA,nrow(prediction_matrix)))
  })
  return(res)
}


# create function -average by tester
average_by_tester = function(Y_test,prediction_matrix,minTrnEnvs=1){
  res = sapply(colnames(Y_test),function(st) {
    (tester = strsplit(st,'::',fixed=T)[[1]][2])
    idx = grepl(tester,colnames(prediction_matrix))
    if(sum(idx)>=minTrnEnvs)
      return(rowMeans(prediction_matrix[,idx,drop=FALSE],na.rm=T))
    return(rep(NA,nrow(prediction_matrix)))
  })
}

# create function - average by state+tester
average_by_State_Plus_Tester = function(Y_test,prediction_matrix,minTrnEnvs=1){
  res=sapply(colnames(Y_test),function(s) {
    (state = substr(s,start = 1, stop = 2))
    (tester = strsplit(s,'::',fixed=T)[[1]][2])
    state_plus_tester = grepl(state,colnames(prediction_matrix)) | grepl(tester,colnames(prediction_matrix))
    if(sum(state_plus_tester)>=minTrnEnvs)
      return(rowMeans(prediction_matrix[,state_plus_tester,drop=FALSE],na.rm=T))
    return(rep(NA,nrow(prediction_matrix)))
  })
  return(res)
}

# create function -average by state:tester
average_by_State_x_Tester= function(Y_test,prediction_matrix){
  res=sapply(colnames(Y_test),function(s) {
    (state = substr(s,start = 1, stop = 2))
    (tester = strsplit(s,'::',fixed=T)[[1]][2])
    
    if(sum(grepl(state,colnames(prediction_matrix))&grepl(tester,colnames(prediction_matrix)))>0)
      return(rowMeans(prediction_matrix[,grepl(state,colnames(prediction_matrix))&grepl(tester,colnames(prediction_matrix)),
                                        drop=FALSE],na.rm=T))
    return(rep(NA,nrow(prediction_matrix)))
  })
  return(res)
}

# average_by_all
average_by_all = function(Y_test,prediction_matrix){
  res = sapply(colnames(Y_test),function(s) {
    rowMeans(prediction_matrix,na.rm=T)
  })
  return(res)
}
