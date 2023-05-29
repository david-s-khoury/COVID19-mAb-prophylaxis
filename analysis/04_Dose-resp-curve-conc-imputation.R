# -------------------------------------------------------------------------
#' 
#' Dose-response curve for prophylactic mAb treatment to prevent sympto-
#' matic SARS-CoV-2 infection
#' Use multiple imputation of the concentration data to include uncertainty
#' in the mAb concentration and in vitro IC50 from the meta-analysis
#' 
# -------------------------------------------------------------------------


# parameters as for fitting to the original data


# impute a data set and analyze it ----------------------------------------

param_best_imputed <- matrix(NA,n_imputations,length(param_best))
hessian_best_imputed <- vector(n_imputations,mode='list')

for(j in 1:n_imputations){
  # sample a random day from the days range (uniformly):
  data_imputed <- data_all
  data_imputed$day <- sample.int.vec(data_imputed$day.min,data_imputed$day.max)
  
  # use the mAb concentration from this day (using data_conc_daily):
  data_imputed <- dplyr::mutate(data_imputed,conc=fun.conc(data_imputed,data_conc_daily)) 
  
  # sample a random IC50 using a normal distribution of log-transformed IC50s using parameters from the meta-analysis:
  data_imputed <- dplyr::mutate(data_imputed,ic50.sampled=sample.ic50s(Means_table,data_imputed))
  
  # make an imputed data set: make fold-IC50 concentration
  data_imputed <- dplyr::mutate(data_imputed,conc.ic50=conc/(ic50.sampled*1e-3))
  
  # analyze the imputed data set in the same way as the original data set: fit a logistic efficacy function only
  data_imputed_fit <- data_imputed[data_imputed$earliest.time==0,] # exclude earliest time interval
  baseline_risks <- data_imputed_fit$e.control/data_imputed_fit$n.control # initialize baseline risk:
  
  RandomInitialEstimate <- NULL
  all_hessians <- vector(n_initial,mode = 'list')
  
  for(i in 1:n_initial){
    # logistic efficacy function:
    par_initial <- c(baseline_risks,LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]]+
                       (UpperB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]]-LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]])*
                       runif(length(param_best),0,1))
    tempfit <- nlm(function(par){nllh(par,data_imputed_fit,model = use_efficacy_model)},par_initial,iterlim=300,hessian = TRUE)
    temp_par <- tail(tempfit$estimate,length(param_best))
    
    output <- data.frame("exitflag"=tempfit$code,
                         "n.iter"=tempfit$iterations,
                         "nllh"=tempfit$minimum)
    output <- cbind(output,t(temp_par))
    
    RandomInitialEstimate <- rbind(RandomInitialEstimate,output)
    all_hessians[[i]] <- tempfit$hessian
  }
  
  # best fit for the current imputed data set:
  param_best_imputed[j,] <- as.numeric(RandomInitialEstimate[which.min(RandomInitialEstimate$nllh),c(4:ncol(RandomInitialEstimate))])
  hessian_best_imputed[[j]] <- all_hessians[[which.min(RandomInitialEstimate$nllh)]]
  
  # print progress:
  print(".")
  if(j%%50==0){
    print(j)
  }
}


# combine analysis of imputed data sets -----------------------------------

# use Rubin's rules to estimate the efficacy parameters and their uncertainty:
param_imputed_est <- colMeans(param_best_imputed) # parameter estimate is the mean of the parameter estimates for the imputed data

# variance-covariance matrix:
var_within <- fun.var.within(hessian_best_imputed) # estimated within imputation variance
var_between <- fun.var.between(param_best_imputed) # estimated between imputation variance
var_total <- var_within+(1+1/nrow(param_best_imputed))*var_between # total variance, taking into account within and between imputation variance

# CIs:
param_imputed_est_ci <- cbind(param_imputed_est-qnorm((1+alpha_CI)/2)*sqrt(diag(var_total)),param_imputed_est+qnorm((1+alpha_CI)/2)*sqrt(diag(var_total)))

# transformed parameters and CIs:
param_imputed_est_transf <- param_imputed_est
param_imputed_est_transf_ci <- param_imputed_est_ci
par_names <- rep("",length(param_best))
if(use_efficacy_model=="logistic"){
  param_imputed_est_transf <- exp(param_imputed_est_transf)
  param_imputed_est_transf[1] <- ifelse(max_transform=="1-exp",1-exp(param_imputed_est[1]),ifelse(max_transform=="exp",exp(param__imputed_est[1]),ifelse(max_transform=="1-abs",1-abs(param_imputed_est[1]),NA)))
  param_imputed_est_transf_ci <- exp(param_imputed_est_transf_ci)
  param_imputed_est_transf_ci[1,] <- c(min(ifelse(max_transform=="1-exp",1-exp(param_imputed_est_ci[1,]),ifelse(max_transform=="exp",exp(param_imputed_est_ci[1,]),ifelse(max_transform=="1-abs",1-abs(param_imputed_est_ci[1,]),NA)))),
                                max(ifelse(max_transform=="1-exp",1-exp(param_imputed_est_ci[1,]),ifelse(max_transform=="exp",exp(param_imputed_est_ci[1,]),ifelse(max_transform=="1-abs",1-abs(param_imputed_est_ci[1,]),NA)))))
  par_names <- c("maximal efficacy","slope","IC50")
}else if(use_efficacy_model=="logistic with max 1"){
  param_imputed_est_transf <- exp(param_imputed_est_transf)
  param_imputed_est_transf_ci <- exp(param_imputed_est_transf_ci)
  par_names <- c("slope","IC50")
}

# parameter estimates table:
table_par <- data.frame(parameter=par_names,
                        estimate=param_imputed_est_transf,
                        CI_lower=param_imputed_est_transf_ci[,1],
                        CI_upper=param_imputed_est_transf_ci[,2])


# save results ------------------------------------------------------------

# save all best fit parameter for imputed data:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(param_best_imputed,param_imputed_est,table_par,param_imputed_est_ci,
     hessian_best_imputed,var_within,var_between,var_total,
     file=glue::glue("output/Best-fit-parameters-imputed-data_{Sys.Date()}_{cur_time}.RData"))
write_xlsx(table_par,glue::glue("output/Imputed-data-parameters_{Sys.Date()}_{cur_time}.xlsx"))


# cleanup -----------------------------------------------------------------

rm(param_best_imputed,hessian_best_imputed,j,data_imputed,baseline_risks,RandomInitialEstimate,
   all_hessians,i,par_initial,tempfit,temp_par,output,param_imputed_est,var_within,var_between,
   var_total,param_imputed_est_ci,param_imputed_est_transf,param_imputed_est_transf_ci,
   par_names,table_par,cur_time)
