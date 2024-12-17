#' Fit reduced and expanded models
#' @title Fit reduced and expanded models
#' @description
#' `run_NI_test()` fits the reduced and expanded models to the data, computes linear combinations of their coefficients, and tests the difference
#' 
#' @param data Data frame containing all the variables in the reduced and expanded models, the cluster variable (if any), and the weight variable (if any).
#' @param reduced Formula object giving the reduced model form. 
#' @param expanded Formula object giving the expanded model form. 
#' @param lincom_var Character containing a string found (only!) in the variables that will be average together. 
#' @param robust Logical indicating whether to use heteroskedasticity-robust standard errors. Default is `TRUE`.
#' @param cluster Character containing the name of the cluster variable in the data OR logical FALSE indicating no clustering. Default is `FALSE`.
#' @param weight Character containing the name of the weight variable in the data OR logical FALSE indicating no weights. Default is `FALSE`.
#' @param null_reduced Logical indicating whether to use the reduced model's residuals in the covariance calculations. Default is `FALSE`; see Appendix B of Bilinski and Hatfield Stat Med (2025) for more detail.
#' @param alpha Numerical value of the type I error rate to use in constructing confidence intervals. Default is 0.05.
#' 
#' @returns 
#' A list containing the treatment effects in the reduced and expanded models, their standard errors, the difference between them, and the uncertainty of the difference (both the standard error and a 1-alpha/2 confidence interval). 
#' \describe{
#'  \item{diff}{Reduced model linear combination minus expanded model linear combination.}
#'  \item{se}{Standard error of this difference.}
#'  \item{CI}{1-alpha/2 confidence interval around the difference.}
#'  \item{tx_r}{Linear combination of estimated coefficients from the reduced model.}
#'  \item{tx_e}{Linear combination of estimated coefficients from the expanded model.}
#'  \item{v_r}{Variance of the linear combination from the reduced model.}
#'  \item{v_e}{Variance of the linear combination from the expanded model.} 
#'  \item{cov_lincom}{Covariance between the two linear combinations.}
#' }
#' 
#' @details
#' Both models must contain at least one variable with `lincom_var` in the name, but there are otherwise no required to have any relationship to the reduced model. A helpful way of adding parameters to your model is to construct them using `update.formula()`, for example, `expanded = update.formula(reduced_model,"~. + treat:year")`
#' Note that this means you should create your own interaction variable in the data if your treatment effects are interaction terms. For instance, if you have variables `post.month` and `treat`, create a variable in the data called `post.by.treat` rather than using a `post.month:treat` interaction term in your model formulas.
#' The linear combination is simply the average of all the coefficients corresponding to `lincom_var` in each model.
#' 
#' @seealso [formula], [update.formula]
#' 
#' @examples
#' data("depcov")
#' 
#' # Reduced model is described in Bilinski and Hatfield
#' # it omits the unemployment control variable from Ankosa-Antwi et al.
#' 
#' red.mod <- as.formula(paste(instype[k],"~",paste(c("fipstate","factor(trend)",
#'   "enact.trt.month","impl.trt.month", # these are treat.by.post
#'   "female","hispanic","white", "asian", "other",
#'   "mar", "student","fpl_ratio","fpl_ratio_2",paste0("age",c(17:25,27:29))),collapse="+")))
#'   
#' # Expanded model adds treatment-group specific linear time trend
#' exp.mod <- update.formula(red.mod,"~ . + fedelig:trend")
#' 
#' # Fit both using robust standard errors, clustering, and weights
#'   grp.lin.yr <- run_NI_test(data=depcov,reduced = red.mod,expanded=exp.mod,lincom_var = 'impl.trt',cluster = 'fipstate',weight = 'weight')
#'
#' # Show the difference and CI
#' grp.lin.yr$diff
#' grp.lin.yr$CI   

#### RUN NON-INFERIORITY TESTS ####
run_NI_test = function(data, reduced, expanded, lincom_var, 
                       robust = TRUE, cluster = FALSE, weight = FALSE, 
                       null_reduced = FALSE, alpha = .05){
  ## Check that both models have the named treatment variable:
  stopifnot("Reduced and expanded models must both contain `lincom_var`."=
              any(grepl(lincom_var,reduced))&any(grepl(lincom_var,expanded)))
  # Check that weight is in the data (or set to FALSE)
  stopifnot("`weight` must either name a variable in the dataset or be set to FALSE."=
              ifelse(is.character(weight),any(grepl(weight,names(data))),TRUE))
  # Check that cluster is in the data (or set to FALSE)
  stopifnot("`cluster` must either name a variable in the dataset or be set to FALSE."=
              ifelse(is.character(cluster),any(grepl(cluster,names(data))),TRUE))
  
  ## If weight is FALSE, set all equal
  if(!is.character(weight)) data$weight = 1
  
  ## Set cluster variable name to cluster
  if(is.character(cluster)) data$cluster = data[,cluster]
  
  # Extract sample size
  n = nrow(data)
  
  # Fit the reduced model
  r = fixest::feols(reduced, data=data, weights=~weight)   # Run reduced model
  Xr = model.matrix(r)                            # Extract model matrix
  ur = r$residuals                                # Extract residuals
  Vr = r$cov.iid/r$sigma2                         # Extract unscaled covariance
  coeff_r = r$coefficients                        # Extract coefficients
  
  # Fit the expanded model  
  e = fixest::feols(expanded, data=data, weights=~weight)  # Run expanded model
  Xe = model.matrix(e)                            # Extract model matrix 
  ue = e$residuals                                # Extract residuals
  Ve = e$cov.iid/e$sigma2                         # Extract unscaled covariance
  coeff_e = e$coefficients                        # Extract coefficients
  
  # extract parameter size
  p_r = ncol(Xr)
  p_e = ncol(Xe)
  
  # Fit variance covariance matrices
  
  # Select residuals to use for reduced model
  # If the null hypothesis is that the reduced model is true: ur
  # Otherwise: ue (see Supplement Lemma 3, Bilinski and Hatfield)
  if(null_reduced){
    ur_for_vcov_r = ur
  }else{
    ur_for_vcov_r = ue
    p_r = ncol(Xe)
  }
  
  # Estimate variance covariance matrices
  if(!robust & !is.character(cluster)){
    # Homoskedastic i.i.d. errors
    vcov_r = Vr*mean(ur_for_vcov_r^2*data$weight)/(n-p_r)*n
    vcov_e = Ve*mean(ue^2*data$weight)/(n-p_e)*n
    cov = Vr*mean(ur_for_vcov_r*ue*data$weight)/(n-p_e)*n
  } else if(is.null(cluster)){
    # Heteroskedasticity robust SEs, without clustering
    
    # add weights
    ur_for_vcov_r = ur_for_vcov_r*data$weight
    ue = ue*data$weight
    
    # calculate variance covariance matrices
    vcov_r = n/(n-1)*eigenMapMatMult(eigenMapMatMult(Vr, eigenMapMatMult(t(Xr*ur_for_vcov_r^2), Xr)), Vr)
    vcov_e = n/(n-1)*eigenMapMatMult(eigenMapMatMult(Ve, eigenMapMatMult(t(Xe*ue^2), Xe)), Ve)
    cov =  n/(n-1)*eigenMapMatMult(eigenMapMatMult(Vr, eigenMapMatMult(t(Xr*ur_for_vcov_r*ue), Xe)), Ve)
  } else if(is.character(cluster)){
    # Heterosketasticity robust SEs, with clustering
    
    # add weights
    ur_for_vcov_r = ur_for_vcov_r*data$weight
    ue = ue*data$weight
    
    # set up for calculating over clusters
    g = as.character(unlist(unique(data$cluster)))
    vcov_r_meat = list()
    vcov_e_meat = list() 
    cov_meat = list() 
    
    # calculate over clusters
    for(i in 1:length(g)){
      # subset observations
      these = which(data$cluster==g[i])
      Xr_sub = collapse::fsubset(Xr, these); ur_for_vcov_r_sub = matrix(collapse::fsubset(ur_for_vcov_r, these), ncol = 1)
      Xe_sub = collapse::fsubset(Xe, these); ue_sub = matrix(collapse::fsubset(ue, these), ncol = 1)
      
      # estimate inner matrices
      vcov_r_meat[[i]] = eigenMapMatMult(eigenMapMatMult(eigenMapMatMult(t(Xr_sub), ur_for_vcov_r_sub), t(ur_for_vcov_r_sub)), Xr_sub)
      vcov_e_meat[[i]] = eigenMapMatMult(eigenMapMatMult(eigenMapMatMult(t(Xe_sub), ue_sub), t(ue_sub)), Xe_sub)
      cov_meat[[i]] = eigenMapMatMult(eigenMapMatMult(eigenMapMatMult(t(Xr_sub), ur_for_vcov_r_sub), t(ue_sub)), Xe_sub)
    }
    
    # calculate over clusters
    vcov_r = eigenMapMatMult(eigenMapMatMult(Vr, Reduce("+", vcov_r_meat)), Vr)*length(g)/(length(g)-1)*(n-1)/(n-p_r)
    vcov_e = eigenMapMatMult(eigenMapMatMult(Ve, Reduce("+", vcov_e_meat)), Ve)*length(g)/(length(g)-1)*(n-1)/(n-p_e)
    cov = eigenMapMatMult(eigenMapMatMult(Vr, Reduce("+", cov_meat)), Ve)*length(g)/(length(g)-1)
    # save space:
    rm(vcov_e_meat,vcov_r_meat); gc()
  }
  
  # Estimate linear combinations
  vec_lincom_r = matrix(0, length(coeff_r), 1)
  comvars_r <- grep(lincom_var,names(coeff_r))
  vec_lincom_r[comvars_r,] = 1/length(comvars_r)
  tx_r = sum(coeff_r*vec_lincom_r)
  
  vec_lincom_e = matrix(0, length(coeff_e), 1)
  comvars_e <- grep(lincom_var,names(coeff_e))
  vec_lincom_e[comvars_e,] = 1/length(comvars_e)
  tx_e = sum(coeff_e*vec_lincom_e)
  
  # calculate difference between reduced and expanded model linear combinations
  diff = tx_r - tx_e
  
  # calculate standard errors
  v_r = eigenMapMatMult(eigenMapMatMult(t(vec_lincom_r), vcov_r),vec_lincom_r)
  v_e = eigenMapMatMult(eigenMapMatMult(t(vec_lincom_e),vcov_e),vec_lincom_e)
  cov_lincom = eigenMapMatMult(eigenMapMatMult(t(vec_lincom_r),cov), vec_lincom_e)
  
  se = sqrt(v_r - 2*cov_lincom + v_e)
  
  # statistical significance
  CI = c(diff - se*qnorm(1-alpha/2), diff + se*qnorm(1-alpha/2))
  
  return(list('diff'=diff, 
              'se'=se, 
              'CI'=CI, 
              'tx_r'=tx_r, 
              'tx_e'=tx_e, 
              'v_r'=v_r, 
              'v_e'=v_e, 
              'cov_lincom'=cov_lincom))
}
