# =======================================================
# Two functions for meta-analyzing correlations
# DER created on 2023-12-05
# =======================================================

# Here is a script with two functions for meta-analyzing correlations:
# meta_correlations(r,n) 
# - is used to summarize the population of correlations, accounting for measurement error. 
# - It takes a vector of correlation estimates r and a corresponding vector of sample sizes n 
# - returns the estimated mean (mean_cor) and sd (sd) of the population of correlations.

# meta_compare_correlations(r12,r13,r23,n)  
# - is used to compare two different correlations measured on the same set of observations. 
# - Let 1 be the observed values, 2 be MegaLMM's estimate and 3 be GBLUP's estimate. 
# - r12 is the correlation between MegaLMM and observed, 
# - r13 is the correlation between GBLUP and observed, and 
# - r23 is the correlation between MegaLMM and GBLUP. 
# - n is the number of observations. 
# - Each of these is a vector across site:year:testers. 
# - It returns the average (mu) and standard deviation (sd) of the difference r12-r13, 
#   as well as a p-value of the comparison p. 
# - It also returns the meta analysis model which may be useful for somethings model.

meta_correlations = function(r,n){
  # follows instructions here: https://www.metafor-project.org/doku.php/tips:hunter_schmidt_method
  # variances are based on the weighted average correlation: (1-r^2)^2/(n-1)
  # this is recommended, even with variation in true correlation
  # this method works on raw correlations, not z-scores
  require(metafor)

  d = escalc(measure = 'COR',ri = r, ni = n, vtype = 'AV')
  m = metafor::rma(yi = d$yi, vi = d$vi,ni = n,method = 'HS')
  mu=m$beta
  sd = sqrt(m$tau2)

  # result
  list(mean_cor = mu,
       sd = sd
  )
}

meta_compare_correlations = function(r12,r13,r23,n) {
  # uses method here: https://personality-project.org/r/psych/help/r.test.html to compute SE of difference between dependent correlations
  # then does standard meta-analysis of these
  # uses raw r values
  require(metafor)
  require(psych)

  # notes: pooled=T and paired=T do not matter
  # estimates r12-r13
  diff = r12-r13
  t = r.test(n,r12=r12,r13=r13,r23=r23)$t
  se = diff/t

  # fit mean and SD of diff
  m = metafor::rma(diff,vi = se^2,weights=n) # Note: weighted by n
  mu=m$beta
  sd = sqrt(m$tau2)

  # return estimates, p-value, and full model

  list(mu = mu,
       sd = sd,
       p = m$pval,
       model = m)
}


# meta_correlations_old = function(r,n, weight_by_ni = F,quantiles = NULL,n_samp = 1e4){
#   require(metafor)
#   require(psych)
#   # convert r to fz
#   fz = psych::fisherz(r)
#
#   # fit mean and SD of fz values
#   if(weight_by_ni) {
#     m = metafor::rma(yi = fz,vi = 1/(n-3),weights = n) # from wiki page, SE = 1/sqrt(N-3)
#   } else {
#     m = metafor::rma(fz,1/(n-3)) # from wiki page, SE = 1/sqrt(N-3)
#   }
#   mu=m$beta
#   sd = sqrt(m$tau2)
#
#   mean_cor = psych::fisherz2r(mu)
#   lower = upper = sd_r = NA
#   if(!is.null(quantiles)) {
#     if(length(quantiles) != 2) stop('quantiles should have length 2')
#     lower = psych::fisherz2r(qnorm(quantiles[1],mu,sd))
#     upper = psych::fisherz2r(qnorm(quantiles[2],mu,sd))
#   }
#   if(!is.null(n_samp) & !is.na(n_samp) & n_samp > 0 ) {
#     samples = psych::fisherz2r(rnorm(n_samp,mu,sd))
#     sd_r = sd(samples)
#     lower = quantile(samples,quantiles[1])
#     upper = quantile(samples,quantiles[2])
#   }
#
#   result = list(est = mean_cor[1],
#                 sd = sd_r,
#                 lower = lower,
#                 upper = upper)
#   result = na.omit(result)
#   return(result)
# }
#
#
# meta_compare_correlations_old = function(r1,r2,n) {
#   require(metafor)
#   require(psych)
#   # convert r to fz
#   fz1 = psych::fisherz(r1)
#   fz2 = psych::fisherz(r2)
#
#   # fit mean and SD of fz values
#   m = metafor::rma(fz2-fz1,2/(n-3)) # from wiki page, SE = 1/sqrt(N-3)
#   mu=m$beta
#   sd = sqrt(m$tau2)
#
#   # return estimates, p-value, and full model
#
#   list(mu = mu,
#        sd = sd,
#        p = m$pval,
#        model = m)
# }
