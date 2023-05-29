# -------------------------------------------------------------------------
#' 
#' Dose-response curve for prophylactic mAb treatment to prevent sympto-
#' matic SARS-CoV-2 infection
#' 
# -------------------------------------------------------------------------


# parameters for plots and analysis ---------------------------------------

plot_data <- FALSE # plot the efficacy data without the fit


# include concentration fold-IC50 and fold-convalescent -------------------

meta.conv <- Means_table$IC50[Means_table$Antibodies=="convalescent"]
meta.ic50s <- ifelse(data_efficacy$trial=="Levin",Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory=="Wild Type"],
                     ifelse(data_efficacy$trial=="Schmidt (delta)",Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory=="Delta"],
                            ifelse(data_efficacy$trial=="Schmidt (omicron)",Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory=="Omicron/BA.1"],
                                   Means_table$IC50[Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory=="Wild Type"])))


# add variables normalised to ic50 and convalescent to full data table -----------------------------

data_all$meta.ic50=meta.ic50s
data_all$geomean.conc.ic50=data_all$geomean.conc/(meta.ic50s*1e-3)
data_all$conc.min.ic50=data_all$conc.min/(meta.ic50s*1e-3)
data_all$conc.max.ic50=data_all$conc.max/(meta.ic50s*1e-3)
data_all$geomean.conc.conv=data_all$geomean.conc/(meta.conv*meta.ic50s*1e-3)
data_all$conc.min.conv=data_all$conc.min/(meta.conv*meta.ic50s*1e-3)
data_all$conc.max.conv=data_all$conc.max/(meta.conv*meta.ic50s*1e-3)

# remove rows with NA efficacy or NA geometric mean concentration:
data_all <- data_all[!is.na(data_all$efficacy) & !is.na(data_all$geomean.conc),]


# visualize the efficacy data ---------------------------------------------

if(plot_data){
  my_colours <- c("mediumpurple","dodgerblue","darkorange") 
  my_shapes <- c(0,1,2,5,6,11)
  
  p_eff <- ggplot(data_all,aes(x=geomean.conc.ic50,y=efficacy,colour=drug,shape=trial,alpha=as.factor(earliest.time))) +
    # efficacy estimates and 95% CIs:
    geom_errorbar(aes(ymin=efficacy.lower,ymax=efficacy.upper),width=0.05,size=0.25) +
    geom_errorbarh(aes(xmin=conc.min.ic50,xmax=conc.max.ic50),height=2,size=0.25) +
    geom_point(size=2) +
    scale_shape_manual(values=my_shapes, name="Trial:") +
    scale_alpha_manual(values=c(1,0.35), labels=c("no","yes"), name="Earliest time interval:") +
    scale_color_manual(values=my_colours, name="mAb:") +
    # axis, theme, etc.:
    scale_x_log10() +
    scale_y_continuous(minor_breaks = seq(0, 100, 20),breaks=seq(0, 100, 20),labels=seq(0, 100, 20),expand = c(2e-2,1e-1)) +
    coord_cartesian(ylim = c(0,100), xlim = c(20,6e4)) +
    labs(x="Antibody concentration [fold in vitro IC50]", y="Efficacy [%]", title = "Prophylactic mAb treatment ") +
    theme_bw() +
    theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2))
  
  p_eff
  
  # save plot:
  cur_time <- format(Sys.time(),"h%H_m%M_s%S")
  ggsave(glue::glue("output/Protection-by-dose_{Sys.Date()}_{cur_time}.pdf"),height=plot_height,width=plot_width,units="in",plot=p_eff)
  
  # clean up:
  rm(my_colours,my_shapes,p_eff,cur_time)
}


# test for a slope leaving out multiple studies ---------------------------

# reshape data for the statistical model:
# exclude first time interval and select columns
data_stats <- data_all[data_all$earliest.time==0,names(data_all)%in%c("trial","geomean.conc.ic50","n.treatment","e.treatment",
                                             "n.control","e.control")]
data_stats <- cbind(data_stats,trial.time=c(1:nrow(data_stats)))
data_stats.tmp.n <- melt(data_stats[,names(data_stats)%in%c("trial","geomean.conc.ic50","trial.time","n.treatment","n.control")],id=c("trial","geomean.conc.ic50","trial.time"))
data_stats.tmp.n <- dplyr::mutate(data_stats.tmp.n,treatment=ifelse(variable%in%c("n.treatment","e.treatment"),1,0))
data_stats.tmp.n <- data_stats.tmp.n[,names(data_stats.tmp.n)!="variable"]
names(data_stats.tmp.n)[names(data_stats.tmp.n)=="value"] <- "n"
data_stats.tmp.e <- melt(data_stats[,names(data_stats)%in%c("trial","geomean.conc.ic50","trial.time","e.treatment","e.control")],id=c("trial","geomean.conc.ic50","trial.time"))
data_stats.tmp.e <- dplyr::mutate(data_stats.tmp.e,treatment=ifelse(variable%in%c("n.treatment","e.treatment"),1,0))
data_stats.tmp.e <- data_stats.tmp.e[,names(data_stats.tmp.e)!="variable"]
names(data_stats.tmp.e)[names(data_stats.tmp.e)=="value"] <- "events"
data_stats <- dplyr::left_join(data_stats.tmp.n,data_stats.tmp.e,by=c("trial","geomean.conc.ic50","treatment","trial.time"))
data_stats <- dplyr::mutate(data_stats,no_events=n-events)

n_trials <- length(unique(data_stats$trial))
# all subsets of at least two studies:
all_subsets <- unlist(lapply(2:n_trials, combn,x = unique(data_stats$trial), simplify = FALSE),recursive = FALSE)
# fixed_effects <- vector(length(all_subsets),mode='list')
slope_leave_out <- data.frame(matrix(0,length(all_subsets),n_trials+1))
names(slope_leave_out) <- c(unique(data_stats$trial),"pval")

for(i in 1:nrow(slope_leave_out)){
  slope_leave_out[i,which(unique(data_stats$trial)%in%all_subsets[[i]])] <- 1
  
  # compute the p-value for the slope:
  data_stats.tmp <- data_stats[which(data_stats$trial%in%all_subsets[[i]]),]
  tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment + treatment:log10(geomean.conc.ic50) + (1 | trial),
                     control = glmerControl(optimizer="bobyqa"),family=binomial("log"), data=data_stats.tmp)
  
  # est <- exp(fixef(tmp.glmer))
  # tmp.ci <- exp(confint(tmp.glmer))
  # names(fixed_effects)[i] <- paste(all_subsets[[i]],collapse = ", ")
  # fixed_effects[[i]] <- data.frame(parameter=names(est),estimate=est,lower=tmp.ci[2:4,1],upper=tmp.ci[2:4,2])
  
  test.sign <- drop1(tmp.glmer, scope=c("treatment:log10(geomean.conc.ic50)"), test="Chisq")
  slope_leave_out$pval[i] <- test.sign$`Pr(Chi)`[2]
}

# save all p-values:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(slope_leave_out,glue::glue("output/Statistical-models-slope-leave-x-out_{Sys.Date()}_{cur_time}.xlsx"))
# write_xlsx(fixed_effects,glue::glue("output/Statistical-models-slope-leave-x-out-fixed-effects_{Sys.Date()}_{cur_time}.xlsx"))


# test for a slope using the measured concentration -----------------------

# reshape data for the statistical model:
# exclude first time interval and select columns (geometric mean concentration 
# not concentration fold IC50)
data_stats <- data_all[data_all$earliest.time==0,names(data_all)%in%c("trial","geomean.conc","n.treatment","e.treatment",
                                                                      "n.control","e.control")]
data_stats <- cbind(data_stats,trial.time=c(1:nrow(data_stats)))
data_stats.tmp.n <- melt(data_stats[,names(data_stats)%in%c("trial","geomean.conc","trial.time","n.treatment","n.control")],id=c("trial","geomean.conc","trial.time"))
data_stats.tmp.n <- dplyr::mutate(data_stats.tmp.n,treatment=ifelse(variable%in%c("n.treatment","e.treatment"),1,0))
data_stats.tmp.n <- data_stats.tmp.n[,names(data_stats.tmp.n)!="variable"]
names(data_stats.tmp.n)[names(data_stats.tmp.n)=="value"] <- "n"
data_stats.tmp.e <- melt(data_stats[,names(data_stats)%in%c("trial","geomean.conc","trial.time","e.treatment","e.control")],id=c("trial","geomean.conc","trial.time"))
data_stats.tmp.e <- dplyr::mutate(data_stats.tmp.e,treatment=ifelse(variable%in%c("n.treatment","e.treatment"),1,0))
data_stats.tmp.e <- data_stats.tmp.e[,names(data_stats.tmp.e)!="variable"]
names(data_stats.tmp.e)[names(data_stats.tmp.e)=="value"] <- "events"
data_stats <- dplyr::left_join(data_stats.tmp.n,data_stats.tmp.e,by=c("trial","geomean.conc","treatment","trial.time"))
data_stats <- dplyr::mutate(data_stats,no_events=n-events)

n_trials <- length(unique(data_stats$trial))
# all subsets of at least two studies:
all_subsets <- unlist(lapply(2:n_trials, combn,x = unique(data_stats$trial), simplify = FALSE),recursive = FALSE)
# fixed_effects <- vector(length(all_subsets),mode='list')
slope_leave_out <- data.frame(matrix(0,length(all_subsets),n_trials+1))
names(slope_leave_out) <- c(unique(data_stats$trial),"pval")

for(i in 12:nrow(slope_leave_out)){
  slope_leave_out[i,which(unique(data_stats$trial)%in%all_subsets[[i]])] <- 1
  
  # compute the p-value for the slope:
  data_stats.tmp <- data_stats[which(data_stats$trial%in%all_subsets[[i]]),]
  tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment + treatment:log10(geomean.conc) + (1 | trial),
                     control = glmerControl(optimizer="bobyqa"),family=binomial("log"), data=data_stats.tmp)
  
  # est <- exp(fixef(tmp.glmer))
  # tmp.ci <- exp(confint(tmp.glmer))
  # names(fixed_effects)[i] <- paste(all_subsets[[i]],collapse = ", ")
  # fixed_effects[[i]] <- data.frame(parameter=names(est),estimate=est,lower=tmp.ci[2:4,1],upper=tmp.ci[2:4,2])
  
  test.sign <- drop1(tmp.glmer, scope=c("treatment:log10(geomean.conc)"), test="Chisq")
  slope_leave_out$pval[i] <- test.sign$`Pr(Chi)`[2]
}

# save all p-values:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(slope_leave_out,glue::glue("output/Statistical-models-slope-leave-x-out(geomean-conc)_{Sys.Date()}_{cur_time}.xlsx"))
# write_xlsx(fixed_effects,glue::glue("output/Statistical-models-slope-leave-x-out(geomean-conc)-fixed-effects_{Sys.Date()}_{cur_time}.xlsx"))


# maximum-likelihood fit of a dose-response curve -------------------------

# fit to data excluding the earliest time interval:
data_fit <- data_all[data_all$earliest.time==0,]

# initialize baseline risk:
baseline_risks <- data_fit$e.control/data_fit$n.control

# models to use: the efficacy model to use for the analysis and the models to test
if(use_efficacy_model%in%test_models){
  models <- test_models
}else{
  models <- c(use_efficacy_model,test_models)
}

# save output from model fitting for all efficacy functions:
RandomInitialEstimate <- vector(length(models),mode='list')
names(RandomInitialEstimate) <- models
RandomInitialEstimate_full <- vector(length(models),mode='list')
names(RandomInitialEstimate_full) <- models
all_hessians <- vector(length(models),mode = 'list')
names(all_hessians) <- models

for(i in 1:n_initial){ # for each random initial condition
  for(j in 1:length(models)){ # for each model to fit
    # random initial conditions sampled uniformly between the lower and upper bound
    par_initial <- c(baseline_risks,LowerB_all_models[[which(names(LowerB_all_models)==models[j])]]+
                       (UpperB_all_models[[which(names(LowerB_all_models)==models[j])]]-LowerB_all_models[[which(names(LowerB_all_models)==models[j])]])*
                       runif(length(LowerB_all_models[[which(names(LowerB_all_models)==models[j])]]),0,1))
    tempfit <- nlm(function(par){nllh(par,data_fit,model = models[j])},par_initial,iterlim=300,hessian = TRUE)
    par <- tail(tempfit$estimate,length(LowerB_all_models[[which(names(LowerB_all_models)==models[j])]]))
    par_full <- tempfit$estimate
    
    output <- data.frame("exitflag"=tempfit$code,
                         "n.iter"=tempfit$iterations,
                         "nllh"=tempfit$minimum)
    output2 <- cbind(output,t(par_full))
    output <- cbind(output,t(par))
    
    
    RandomInitialEstimate[[j]] <- rbind(RandomInitialEstimate[[j]],output)
    RandomInitialEstimate_full[[j]] <- rbind(RandomInitialEstimate_full[[j]],output2)
    all_hessians[[j]][[i]] <- tempfit$hessian
  }
}

# find and save the best parameter values for all models:
param_best_all <- vector(length(models),mode='list')
names(param_best_all) <- models
param_best_all_full <- vector(length(models),mode='list')
names(param_best_all_full) <- models
hessian_best_all  <- vector(length(models),mode='list')
names(hessian_best_all) <- models
AIC_all  <- vector(length(models),mode='list')
names(AIC_all) <- models

for(j in 1:length(models)){
  param_best_all[[j]] <- as.numeric(RandomInitialEstimate[[j]][which.min(RandomInitialEstimate[[j]]$nllh),c(4:ncol(RandomInitialEstimate[[j]]))])
  param_best_all_full[[j]] <- as.numeric(RandomInitialEstimate_full[[j]][which.min(RandomInitialEstimate_full[[j]]$nllh),c(4:ncol(RandomInitialEstimate_full[[j]]))])
  hessian_best_all[[j]] <- all_hessians[[j]][[which.min(RandomInitialEstimate[[j]]$nllh)]]
  AIC_all[[j]] <- 2*(length(baseline_risks)+length(LowerB_all_models[[which(names(LowerB_all_models)==models[j])]]))+min(RandomInitialEstimate[[j]]$nllh)
}

# best fit parameters for the model to use for the analysis:
param_best <- param_best_all[[which(names(param_best_all)==use_efficacy_model)]]
param_best_full <- param_best_all_full[[which(names(param_best_all_full)==use_efficacy_model)]]
hessian_best <- hessian_best_all[[which(names(hessian_best_all)==use_efficacy_model)]]

# make a table with models and AICs:
table_model_AICs <- data.frame(model=models,
                               AIC=unlist(AIC_all))

# save AIC table, all model best fits, and best fit for the model to use for the analysis:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(table_model_AICs,glue::glue("output/Efficacy-models-AICs_{Sys.Date()}_{cur_time}.xlsx"))
save(param_best_all,hessian_best_all,AIC_all,
     file=glue::glue("output/Best-fit-all-models_{Sys.Date()}_{cur_time}.RData"))
save(param_best,hessian_best,use_efficacy_model,
     file=glue::glue("output/Best-fit-analysis-model_{Sys.Date()}_{cur_time}.RData"))


# compute CIs -------------------------------------------------------------

# compute CIs for specified efficacy function (use_efficacy_model) only:
n_par <- length(param_best) # number of parameters of the model for the analysis

if(CIs_parametric){ # compute CIs using parametric bootstrapping
  covar <- solve(hessian_best)[(nrow(hessian_best)-c((n_par-1):0)),(nrow(hessian_best)-c((n_par-1):0))]
  param_best_ci <- cbind(param_best-qnorm((1+alpha_CI)/2)*sqrt(diag(covar)),
                         param_best+qnorm((1+alpha_CI)/2)*sqrt(diag(covar)))
  
  # make a table with transformed parameter values:
  param_best_transf <- param_best
  param_best_transf_ci <- param_best_ci
  par_names <- rep("",length(param_best))
  if(use_efficacy_model=="logistic"){
    param_best_transf <- exp(param_best_transf)
    param_best_transf[1] <- ifelse(max_transform=="1-exp",1-exp(param_best[1]),ifelse(max_transform=="exp",exp(param_best[1]),ifelse(max_transform=="1-abs",1-abs(param_best[1]),NA)))
    param_best_transf_ci <- exp(param_best_transf_ci)
    param_best_transf_ci[1,] <- c(min(ifelse(max_transform=="1-exp",1-exp(param_best_ci[1,]),ifelse(max_transform=="exp",exp(param_best_ci[1,]),ifelse(max_transform=="1-abs",1-abs(param_best_ci[1,]),NA)))),
                                  max(ifelse(max_transform=="1-exp",1-exp(param_best_ci[1,]),ifelse(max_transform=="exp",exp(param_best_ci[1,]),ifelse(max_transform=="1-abs",1-abs(param_best_ci[1,]),NA)))))
    par_names <- c("maximal efficacy","slope","IC50")
  }else if(use_efficacy_model=="logistic with max 1"){
    param_best_transf <- exp(param_best_transf)
    param_best_transf_ci <- exp(param_best_transf_ci)
    par_names <- c("slope","IC50")
  }
  
  ### goodness of fit test
  degree_of_freedom=nrow(data_fit)*2-length(param_best_full)
  goodness_of_fit_stat<-chisquared_gof_statistic(param_best_full,data_fit,use_efficacy_model)
  chisquared_test_pvalue=1-pchisq(goodness_of_fit_stat,degree_of_freedom)

  goodness_of_fit_table<-data.frame(goodness_of_fit_stat,degree_of_freedom,chisquared_test_pvalue)
  
  table_par_best <- data.frame(parameter=par_names,
                               estimate=param_best_transf,
                               lower=param_best_transf_ci[,1],
                               upper=param_best_transf_ci[,2])
  
  
  
  # save best parameter values:
  cur_time <- format(Sys.time(),"h%H_m%M_s%S")
  write_xlsx(table_par_best,glue::glue("output/Best-fit-parameters-table_{Sys.Date()}_{cur_time}.xlsx"))
  write_xlsx(goodness_of_fit_table,glue::glue("output/Best-fit-goodness-of-fit-test_{Sys.Date()}_{cur_time}.xlsx"))
  
  # for visualization:
  par_bootstr <- rmvnorm(n=n_bootstrap,mean=as.matrix(tail(param_best,n_par)),sigma=covar)
  
  # compute efficacy CIs for a given dose:
  eff_ci_from_dose <- function(d){quantile(sapply(c(1:nrow(par_bootstr)),function(x){efficacy(d,par_bootstr[x,],model=use_efficacy_model)}),
                                           probs = c((1-alpha_CI)/2,(1+alpha_CI)/2))}
  
  # cleanup:
  rm(covar,par_names,param_best_transf,param_best_transf_ci)
  
}else{ # compute CIs using resampling of the data
  
  # fit to data excluding the earliest time interval:
  data_fit <- data_all[data_all$earliest.time==0,]
  
  # for each resampled data set, save the best fit parameter values and exitflags:
  param_resamp <- matrix(NA,n_bootstrap,n_par)
  exitflag_resamp <- matrix(NA,n_bootstrap,1)
  
  # resample data and fit the efficacy function:
  for(j in 1:n_bootstrap){
    # resample data:
    data_resamp <- data_fit[sample.int(nrow(data_fit),nrow(data_fit),replace=TRUE),]
    
    # initialize baseline risk:
    baseline_risks_resamp <- data_resamp$e.control/data_fit$n.control
    
    # save output from model fitting:
    RandomInitialEstimate <- NULL
    
    for(i in 1:n_initial){
      par_initial <- c(baseline_risks,LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]]+
                         (UpperB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]]-LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]])*
                         runif(n_par,0,1))
      tempfit <- nlm(function(par){nllh(par,data_fit,model = use_efficacy_model)},par_initial,iterlim=300,hessian = TRUE)
      par <- tail(tempfit$estimate,n_par)
      
      output <- data.frame("exitflag"=tempfit$code,
                           "n.iter"=tempfit$iterations,
                           "nllh"=tempfit$minimum)
      output <- cbind(output,t(par))
      
      RandomInitialEstimate <- rbind(RandomInitialEstimate,output)
    }
    
    # find best parameter values:
    param_resamp[j,] <- as.numeric(RandomInitialEstimate[which.min(RandomInitialEstimate$nllh),c(4:ncol(RandomInitialEstimate))])
    exitflag_resamp[j] <- RandomInitialEstimate$exitflag[which.min(RandomInitialEstimate$nllh)]
    
    print(glue::glue("j=",j))
  }
  
  # save best parameter values for all resampled data sets:
  cur_time <- format(Sys.time(),"h%H_m%M_s%S")
  save(param_resamp,exitflag_resamp,
       file=glue::glue("output/Best-fit-parameters-resampled-data_{Sys.Date()}_{cur_time}.RData"))
  
  # compute efficacy CIs for a given dose:
  eff_ci_from_dose <- function(d){quantile(sapply(c(1:nrow(param_resamp)),function(x){efficacy(d,param_resamp[x,],model = use_efficacy_model)}),
                                           probs = c((1-alpha_CI)/2,(1+alpha_CI)/2))}
  
  # cleanup:
  rm(exitflag_resamp,data_resamp,baseline_risks_resamp)
}


# Visualize the dose-response curve with CIs ------------------------------

# Visualize the efficacy function fit & CI for the model to use for the analysis

# prepare data:
doses <- 10^(seq(0,5,length.out=1e3)) # doses in fold-IC50
model_fit <- as.numeric(efficacy(doses,param_best,model = use_efficacy_model))
model_ci <- sapply(doses,function(x){eff_ci_from_dose(x)})
model_data <- data.frame(dose=doses,fit=model_fit,lower=model_ci[1,],upper=model_ci[2,])

# save model data:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(model_data,use_efficacy_model,file=glue::glue("output/Model-data-analysis-model_{Sys.Date()}_{cur_time}.RData"))

# visualize:
my_colours <- c("mediumpurple","dodgerblue","darkorange") 
my_shapes <- c(0,1,2,5,6,11)

p_fit <- ggplot(data_all,aes(x=geomean.conc.ic50,y=efficacy,colour=drug,shape=trial,alpha=as.factor(earliest.time))) +
  # model fit:
  geom_ribbon(data=model_data,inherit.aes=FALSE,aes(x=dose,ymin=100*lower, ymax=100*upper),fill="black", alpha = 0.1)+
  geom_line(data=model_data,inherit.aes=FALSE,aes(x=dose,y=100*fit),color="black") +
  # efficacy estimates and 95% CIs:
  geom_errorbar(aes(ymin=efficacy.lower,ymax=efficacy.upper),width=0.05,size=0.25) +
  geom_errorbarh(aes(xmin=conc.min.ic50,xmax=conc.max.ic50),height=2,size=0.25) +
  geom_point(size=2) +
  scale_shape_manual(values=my_shapes, name="Trial:") +
  scale_alpha_manual(values=c(1,0.35), labels=c("no","yes"), name="Earliest time interval:") +
  scale_color_manual(values=my_colours, name="mAb:") +
  # axis, theme, etc.:
  scale_x_log10() +
  scale_y_continuous(minor_breaks = seq(0, 100, 20),breaks=seq(0, 100, 20),labels=seq(0, 100, 20),expand = c(2e-2,1e-1)) +
  coord_cartesian(ylim = c(0,100), xlim = c(20,6e4)) +
  labs(x="Antibody concentration [fold in vitro IC50]", y="Efficacy [%]", title = "Prophylactic mAb treatment") +
  theme_bw() +
  theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2))

p_fit

# save plot:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
ggsave(glue::glue("output/Fig2_dose-response_{Sys.Date()}_{cur_time}.pdf"),height=plot_height,width=plot_width,units="in",plot=p_fit)


# visualize all tested efficacy functions without CIs:

for(j in c(1:length(models))){
  # prepare data:
  model_various <- as.numeric(efficacy(doses,param_best_all[[j]],model = models[j]))
  model_data_various <- data.frame(dose=doses,fit=model_various)
  
  # visualize:
  my_colours <- c("mediumpurple","dodgerblue","darkorange") 
  my_shapes <- c(0,1,2,5,6,11)
  
  p_fit_various <- ggplot(data_all,aes(x=geomean.conc.ic50,y=efficacy,colour=drug,shape=trial,alpha=as.factor(earliest.time))) +
    # model fit:
    geom_line(data=model_data_various,inherit.aes=FALSE,aes(x=dose,y=100*fit),color="black") +
    # efficacy estimates and 95% CIs:
    geom_errorbar(aes(ymin=efficacy.lower,ymax=efficacy.upper),width=0.05,size=0.25) +
    geom_errorbarh(aes(xmin=conc.min.ic50,xmax=conc.max.ic50),height=2,size=0.25) +
    geom_point(size=2) +
    scale_shape_manual(values=my_shapes, name="Trial:") +
    scale_alpha_manual(values=c(1,0.35), labels=c("no","yes"), name="Earliest time interval:") +
    scale_color_manual(values=my_colours, name="mAb:") +
    # axis, theme, etc.:
    scale_x_log10() +
    scale_y_continuous(minor_breaks = seq(0, 100, 20),breaks=seq(0, 100, 20),labels=seq(0, 100, 20),expand = c(2e-2,1e-1)) +
    coord_cartesian(ylim = c(0,100), xlim = c(20,6e4)) +
    labs(x="Antibody concentration [fold in vitro IC50]", y="Efficacy [%]", 
         title = glue::glue("Prophylactic mAb treatment (efficacy model: ",models[j],")")) +
    theme_bw() +
    theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2))
  
  p_fit_various
  
  # save plot:
  cur_time <- format(Sys.time(),"h%H_m%M_s%S")
  ggsave(glue::glue("output/Fig_dose-response-",models[j],"_{Sys.Date()}_{cur_time}.pdf"),height=plot_height,width=plot_width,units="in",plot=p_fit_various)
}


# cleanup -----------------------------------------------------------------

rm(plot_data,data_stats,data_stats.tmp.n,data_stats.tmp.e,n_trials,all_subsets,slope_leave_out,i,data_stats.tmp,
   tmp.glmer,test.sign,cur_time,data_fit,baseline_risks,models,RandomInitialEstimate,all_hessians,j,par_initial,
   tempfit,par,output,param_best_all,hessian_best_all,AIC_all,table_model_AICs,n_par,doses,model_fit,model_ci,
   my_colours,my_shapes,p_fit,model_various,model_data_various,p_fit_various,meta.conv,meta.ic50)

