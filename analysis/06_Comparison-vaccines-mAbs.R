# -------------------------------------------------------------------------
#' 
#' Comparison of the efficacy of prophylactic mAb treatment and
#' vaccinations
#' 
# -------------------------------------------------------------------------


# total efficacy of mAbs and vaccines -------------------------------------

# load data:
data_mAbs_vaccines <- read_excel("raw-data/Data-mAbs-vaccines_2023-05-23.xlsx")
data_mAbs_vaccines$trial <- factor(data_mAbs_vaccines$trial,levels=c("Isa","O'Brien","Herman","Levin","Schmidt (delta)","Pfizer","Moderna"))


# Statistical comparison of mAbs and vaccines:
# reshape data for the model:
data_mab_vac_stats <- melt(data_mAbs_vaccines[,c(1:3,5)],id=c("trial","type"))
data_mab_vac_stats <- dplyr::mutate(data_mab_vac_stats,treatment=ifelse(variable%in%c("n.treatment","e.treatment"),1,0))
data_mab_vac_stats <- data_mab_vac_stats[,-3]
names(data_mab_vac_stats)[3] <- "n"
tmp.events <- cbind(melt(data_mAbs_vaccines[,c(1:2,4,6)],id=c("trial","type")),treatment=c(rep(1,nrow(data_mAbs_vaccines)),rep(0,nrow(data_mAbs_vaccines))))
tmp.events <- tmp.events[,-3]
names(tmp.events)[3] <- "events"
data_mab_vac_stats <- dplyr::left_join(data_mab_vac_stats,tmp.events,by=c("trial","type","treatment"))
data_mab_vac_stats <- dplyr::mutate(data_mab_vac_stats,no_events=n-events)
# add a treatment variable with 3 categories (0=placebo, 1=mAb, 2=vaccine):
data_mab_vac_stats <- dplyr::mutate(data_mab_vac_stats,treatment3=factor(ifelse(treatment==0,0,ifelse(type=="mAb",1,2))))

# statistical model:
model_mAbs_vac <- glmer(cbind(events, no_events) ~ 1 + treatment3 + (1 | trial), control = glmerControl(optimizer="bobyqa"),
                   family=binomial("log"), data=data_mab_vac_stats)
# test if there is a significant difference between the coefficients for mAb and vaccination treatments (treatment3==1 or 2):
comparison.matrix <- matrix(c(0,1,-1),ncol=nrow(summary(model_mAbs_vac)$coefficients),byrow=TRUE)
pval_mAbs_vs_vac <- summary(multcomp::glht(model_mAbs_vac,linfct=comparison.matrix))$test$pvalues[1]


# Average efficacy for mAbs and vaccines:
# mAbs:
tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment + (1 | trial),family=binomial("log"), data=data_mab_vac_stats[data_mab_vac_stats$type=="mAb",])
av_eff_Abs <- 100*(1-exp(fixef(tmp.glmer)[2]))
av_eff_Abs_ci <- as.numeric(sort(100*(1-exp(confint(tmp.glmer))[3,])))

# vaccines:
tmp.glmer <- glmer(cbind(events, no_events) ~ 1 + treatment + (1 | trial),family=binomial("log"), data=data_mab_vac_stats[data_mab_vac_stats$type=="vaccine",])
av_eff_vac <- 100*(1-exp(fixef(tmp.glmer)[2]))
av_eff_vac_ci <- as.numeric(sort(100*(1-exp(confint(tmp.glmer))[3,])))


# save analysis results:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(model_mAbs_vac,pval_mAbs_vs_vac,av_eff_Abs,av_eff_Abs_ci,av_eff_vac,av_eff_vac_ci,
     file=glue::glue("output/MAb-vaccine-comparison_{Sys.Date()}_{cur_time}.RData"))


# Visualize total efficacy of mAbs and vaccines:
my.colours <- c("royalblue3","orangered")
my.shapes <- c(16,17)

p.vaccine.mabs.total <- ggplot(data_mAbs_vaccines,aes(x=trial,y=efficacy,color=type,shape=type)) +
  geom_point() +
  geom_errorbar(aes(ymin=efficacy.lower, ymax=efficacy.upper),width=0.05,size=0.25) +
  scale_color_manual(values=my.colours, name="Treatment:") + 
  scale_shape_manual(values=my.shapes, name="Treatment:") + 
  annotate("segment",x=c(1,6,3,6.5,3),y=c(102,102,102,102,104),xend=c(5,7,3,6.5,6.5),yend=c(102,102,104,104,104)) +
  annotate("text",x=4.75,y=106.3,label="p = 0.002",size = 3) + # see statistical model below for p-value 
  # axis, theme, etc.:
  scale_y_continuous(minor_breaks = seq(0, 100, 5),breaks=seq(0, 100, 10),labels=seq(0, 100, 10),expand = c(2e-2,1e-1)) +
  # coord_cartesian(ylim = c(50,100)) +
  labs(x="Study", y="Efficacy [%]", title = "Overall efficacy of mAbs and vaccines") +
  theme_bw() +
  theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2))

p.vaccine.mabs.total

# save figure:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
ggsave(glue::glue("output/Fig4a_vaccine-mAbs-overall-efficacy_{Sys.Date()}_{cur_time}.pdf"),width=7,height=4,units="in",plot=p.vaccine.mabs.total)


# Visualize mAb and vaccine data with best fits ---------------------------

# load vaccine data:
data_vaccine <- read.csv("raw-data/SummaryTable_Efficacy_NeutRatio_SD_SEM.csv")

# make a combined data frame for vaccine and mAb data:
combined_data_vac_mAbs <- data.frame(type=c(rep("vaccine",nrow(data_vaccine)),rep("mAb",nrow(data_all))),
                                    trial=c(data_vaccine$Study,data_all$trial),
                                    fold.conv=c(data_vaccine$NeutMean/data_vaccine$NeutConv,data_all$geomean.conc.conv),
                                    fold.conv.lower=c(10^(log10(data_vaccine$NeutMean/data_vaccine$NeutConv)-1.96*data_vaccine$SEM),data_all$conc.min.conv),
                                    fold.conv.upper=c(10^(log10(data_vaccine$NeutMean/data_vaccine$NeutConv)+1.96*data_vaccine$SEM),data_all$conc.max.conv),
                                    n.treatment=c(data_vaccine$NumVac,data_all$n.treatment),
                                    e.treatment=c(data_vaccine$InfVac,data_all$e.treatment),
                                    n.control=c(data_vaccine$NumCont,data_all$n.control),
                                    e.control=c(data_vaccine$InfCont,data_all$e.control),
                                    efficacy=c(100*data_vaccine$Efficacy,data_all$efficacy),
                                    efficacy.lower=c(data_vaccine$Lower,data_all$efficacy.lower),
                                    efficacy.upper=c(data_vaccine$Upper,data_all$efficacy.upper),
                                    earliest.time=c(rep(0,nrow(data_vaccine)),data_all$earliest.time))

# exclude the early time points for the fitting:
combined_data_vac_mAbs_fit <- combined_data_vac_mAbs[combined_data_vac_mAbs$earliest.time==0,]

# load best fit model for vaccine data ("Khoury curve"):
load("raw-data/Khoury-curve-data 2022-08-09.RData") # data frame called "LogisticModel_withPoolSD"

# best fit to the mAb data: saved in model_data
# add dose [fold convalescent] to the model data:
model_data <- dplyr::mutate(model_data,dose.conv=dose/Means_table$IC50[Means_table$Antibodies=="convalescent"])

# visualization:
my.colours <- c("royalblue3","orangered") #c("goldenrod1","darkorange","red","red3","lightpink","hotpink","mediumpurple","dodgerblue","turquoise1","palegreen2","springgreen4","darkviolet","purple4","black")
my.shapes <- c(16,17)
dose.range <- c(min(combined_data_vac_mAbs$fold.conv.lower),max(combined_data_vac_mAbs$fold.conv.upper))

p.vaccine.mabs.best.fit <- ggplot(combined_data_vac_mAbs,aes(x=fold.conv,y=efficacy,colour=type,shape=type,alpha=as.factor(earliest.time))) +
  # show the Khoury curve and 95% CI:
  geom_line(data=LogisticModel_withPoolSD[LogisticModel_withPoolSD$NeutRatio_Reported>=log10(min(dose.range)) & LogisticModel_withPoolSD$NeutRatio_Reported<=log10(max(dose.range)),],
            inherit.aes=FALSE,aes(y=100*Efficacy,x=(10^NeutRatio_Reported)),color="orangered") +
  geom_ribbon(data=LogisticModel_withPoolSD[LogisticModel_withPoolSD$NeutRatio_Reported>=log10(min(dose.range)) & LogisticModel_withPoolSD$NeutRatio_Reported<=log10(max(dose.range)),],
              inherit.aes=FALSE,aes(y=100*Efficacy,x=(10^NeutRatio_Reported),ymin=Lower, ymax=Upper),fill="orangered", alpha = 0.3) +
  # model fit:
  geom_ribbon(data=model_data[model_data$dose.conv>=min(dose.range) & model_data$dose.conv<=max(dose.range),],
              inherit.aes=FALSE,aes(x=dose.conv,ymin=100*lower, ymax=100*upper),fill="royalblue3", alpha = 0.1)+
  geom_line(data=model_data[model_data$dose.conv>=min(dose.range) & model_data$dose.conv<=max(dose.range),],
            inherit.aes=FALSE,aes(x=dose.conv,y=100*fit),color="royalblue3") +
  # efficacy estimates and 95% CIs:
  geom_errorbar(aes(ymin=efficacy.lower, ymax=efficacy.upper),width=0.05,size=0.25) +
  geom_errorbarh(aes(xmin=fold.conv.lower, xmax=fold.conv.upper),height=2,size=0.25) +
  geom_point(size=2) +
  scale_shape_manual(values=my.shapes, name="Treatment:") + 
  scale_alpha_manual(values=c(1,0.35), labels=c("no","yes"), name="Earliest time interval:") + 
  scale_color_manual(values=my.colours, name="Treatment:") + 
  # axis, theme, etc.:
  scale_x_log10() +
  scale_y_continuous(minor_breaks = seq(0, 100, 20),breaks=seq(0, 100, 20),labels=seq(0, 100, 20),expand = c(2e-2,1e-1)) +
  coord_cartesian(ylim = c(0,100), xlim = dose.range) + # c(1e-2,100)
  labs(x="Neutralization [fold convalescent]", y="Efficacy [%]", title = "Efficacy of vaccination and prophylactic mAbs") +
  theme_bw() +
  theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2))

p.vaccine.mabs.best.fit

# save figure:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
ggsave(glue::glue("output/Fig_vaccine-mAbs-with-best-fit_{Sys.Date()}_{cur_time}.pdf"),width=7,height=4,units="in",plot=p.vaccine.mabs.best.fit)


# mAb and vaccine data fitting --------------------------------------------

# Compare the vaccine and mAb data by fitting the efficacy model to both data sets simultaneously
# with same/different parameters and comparing different fits with likelihood ratio tests

# specify the parameters that are the same for vaccine and mAb data:
if(use_efficacy_model=="logistic"){
  same.par <- c("m.k.ic50","m.k","m.ic50","k.ic50","m","k","ic50","")
}else if(use_efficacy_model=="logistic with max 1"){
  same.par <- c("k.ic50","k","ic50","")
}

# fit in the same way as the dose-response curve for the mAb data:

# initialize baseline risk:
baseline_risks <- combined_data_vac_mAbs_fit$e.control/combined_data_vac_mAbs_fit$n.control

# save output from model fitting for all models:
RandomInitialEstimate_vac <- vector(length(same.par),mode='list')
names(RandomInitialEstimate_vac) <- paste("same.",same.par,sep="")
RandomInitialEstimate_mAb <- vector(length(same.par),mode='list')
names(RandomInitialEstimate_mAb) <- paste("same.",same.par,sep="")
all_hessians <- vector(length(same.par),mode = 'list')
names(all_hessians) <- paste("same.",same.par,sep="")

for(i in 1:n_initial){ # for each random initial condition
  for(j in 1:length(same.par)){ # for each model to fit
    # index of parameters that are the same for vaccine and mAb data:
    same.par.ind <- get.same.par.index(same.par[j])
    different.par.ind <- get.different.par.index(same.par[j])
    
    # random initial parameters for the vaccine data:
    tempfit <- NA
    while(class(tempfit)!="list"){
      try({
        rand.vac <- LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]]+
          (UpperB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]]-LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]])*
          runif(length(LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]]),0,1)
        rand.vac[get.same.par.index("ic50")] <- log(exp(rand.vac[get.same.par.index("ic50")])/Means_table$IC50[Means_table$Antibodies=="convalescent"]) # convert IC50 from fold-IC50 to fold-convalescent
        # random initial parameters for the mAb data: (the difference to the vaccine parameters)
        rand.mAb <- LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]]+
          (UpperB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]]-LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]])*
          runif(length(LowerB_all_models[[which(names(LowerB_all_models)==use_efficacy_model)]]),0,1)
        rand.mAb[get.same.par.index("ic50")] <- log(exp(rand.mAb[get.same.par.index("ic50")])/Means_table$IC50[Means_table$Antibodies=="convalescent"]) # convert IC50 from fold-IC50 to fold-convalescent
        rand.diff.mAb <- rand.mAb[different.par.ind]-rand.vac[different.par.ind]
        # all initial parameters:
        par_initial <- c(baseline_risks,rand.vac,rand.diff.mAb)
        
        # model fit:
        tempfit <- nlm(function(par){nllh.vac.mAb(par,same.par[j])},par_initial,iterlim=300,hessian = TRUE)},silent = TRUE)
      
    }
    par_est_vac <- tempfit$estimate[nrow(combined_data_vac_mAbs_fit)+c(1:length(rand.vac))]
    par_est_diff_mAb <- rep(0,length(rand.mAb))
    par_est_diff_mAb[different.par.ind] <- tail(tempfit$estimate,length(different.par.ind))
    par_est_mAb <- par_est_vac+par_est_diff_mAb
    
    output <- data.frame("exitflag"=tempfit$code,
                         "n.iter"=tempfit$iterations,
                         "nllh"=tempfit$minimum)
    output_vac <- cbind(output,t(par_est_vac))
    output_mAb <- cbind(output,t(par_est_mAb))
    
    RandomInitialEstimate_vac[[j]] <- rbind(RandomInitialEstimate_vac[[j]],output_vac)
    RandomInitialEstimate_mAb[[j]] <- rbind(RandomInitialEstimate_mAb[[j]],output_mAb)
    all_hessians[[j]][[i]] <- tempfit$hessian
  }
}

# find and save the best parameter values for all models:
param_best_fit_vac <- vector(length(same.par),mode='list')
names(param_best_fit_vac) <- paste("same.",same.par,sep="")
param_best_fit_mAb <- vector(length(same.par),mode='list')
names(param_best_fit_mAb) <- paste("same.",same.par,sep="")
hessian_best_fit_vac_mAb <- vector(length(same.par),mode='list')
names(hessian_best_fit_vac_mAb) <- paste("same.",same.par,sep="") 
  
for(j in 1:length(same.par)){
  param_best_fit_vac[[j]] <- as.numeric(RandomInitialEstimate_vac[[j]][which.min(RandomInitialEstimate_vac[[j]]$nllh),c(4:ncol(RandomInitialEstimate_vac[[j]]))])
  param_best_fit_mAb[[j]] <- as.numeric(RandomInitialEstimate_mAb[[j]][which.min(RandomInitialEstimate_mAb[[j]]$nllh),c(4:ncol(RandomInitialEstimate_mAb[[j]]))])
  hessian_best_fit_vac_mAb[[j]] <- all_hessians[[j]][[which.min(RandomInitialEstimate_vac[[j]]$nllh)]]
}

# save AIC table, all model best fits, and best fit for the model to use for the analysis:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
save(RandomInitialEstimate_vac,RandomInitialEstimate_mAb,all_hessians,
     param_best_fit_vac,param_best_fit_mAb,hessian_best_fit_vac_mAb,
     file=glue::glue("output/Fits-vaccines-and-mAbs_{Sys.Date()}_{cur_time}.RData"))


# comparison of fits with likelihood ratio tests --------------------------

# pairwise comparison of models
lr_comparison_pvals <- data.frame(matrix(NA,length(same.par),length(same.par)))
rownames(lr_comparison_pvals) <- colnames(lr_comparison_pvals) <- paste("same.",same.par,sep="")

for(i in c(1:length(same.par))){
  for(j in c(1:length(same.par))){
    lr_comparison_pvals[i,j] <- lr.pval(same.par1 = same.par[i], same.par2 = same.par[j])
  }
}

# save as a table:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(lr_comparison_pvals,glue::glue("output/Table-likelihood-ratio-test-pvals-vaccine-mAbs_{Sys.Date()}_{cur_time}.xlsx"))


# confidence region for visualization -------------------------------------

# compute the confidence region using parametric bootstrapping

# save the model fit and CIs for each model for vaccine and mAb data:
models_vac_mAb <- vector(length(same.par),mode='list')
names(models_vac_mAb) <- paste("same.",same.par,sep="")
doses <- 10^seq(-3,3,length.out=1e3)

# save the best fit parameters with CIs for each model:
param_best_vac_mAb <- vector(length(same.par),mode='list')
names(param_best_vac_mAb) <- paste("same.",same.par,sep="")

# save the AIC for each fit:
AICs_vac_mAb <- vector(length(same.par),mode='list')
names(AICs_vac_mAb) <- paste("same.",same.par,sep="")

for(j in c(1:length(models_vac_mAb))){
  # covariance matrix:
  covar <- solve(hessian_best_fit_vac_mAb[[j]])[c((length(baseline_risks)+1):nrow(hessian_best_fit_vac_mAb[[j]])),c((length(baseline_risks)+1):nrow(hessian_best_fit_vac_mAb[[j]]))]
  
  # best fit parameters with CIs:
  param_best_vac_mAb[[j]] <- c(param_best_fit_vac[[j]],param_best_fit_mAb[[j]]-param_best_fit_vac[[j]])[which(c(param_best_fit_vac[[j]],param_best_fit_mAb[[j]]-param_best_fit_vac[[j]])!=0)]
  param_ci <- cbind(param_best_vac_mAb[[j]]-qnorm((1+alpha_CI)/2)*sqrt(diag(covar)),
                    param_best_vac_mAb[[j]]+qnorm((1+alpha_CI)/2)*sqrt(diag(covar)))
  param_best_vac_mAb[[j]] <- cbind(param_best_vac_mAb[[j]],param_ci)
  
  # compute AICs:
  AICs_vac_mAb[[j]] <- 2*(nrow(hessian_best_fit_vac_mAb[[j]]))+min(RandomInitialEstimate_mAb[[j]]$nllh)
  
  # model fit for vaccine and mAb data:
  eff_vac <- efficacy(doses,param_best_fit_vac[[j]])
  eff_mAb <- efficacy(doses,param_best_fit_mAb[[j]])
  
  # bootstrap parameters for visualization:
  param_tmp <- c(param_best_fit_vac[[j]],param_best_fit_mAb[[j]]-param_best_fit_vac[[j]])[which(c(param_best_fit_vac[[j]],param_best_fit_mAb[[j]]-param_best_fit_vac[[j]])!=0)]
  par_bootstr <- rmvnorm(n=n_bootstrap,mean=param_tmp,sigma=covar)
  
  # vaccine and mAb specific parameters from the bootstrap:
  par_bootstr_vac <- par_bootstr[,c(1:length(param_best_fit_vac[[j]]))]
  par_bootstr_diff <- matrix(0,nrow(par_bootstr_vac),ncol(par_bootstr_vac))
  if(length(get.different.par.index(same.par[j]))>0){
    par_bootstr_diff[,get.different.par.index(same.par[j])] <- par_bootstr[,c((ncol(par_bootstr_vac)+1):ncol(par_bootstr))]
  }
  par_bootstr_mAb <- par_bootstr_vac+par_bootstr_diff
  
  # compute efficacy CIs for a given dose for vaccines:
  eff_ci_from_dose_vac <- function(d){quantile(sapply(c(1:nrow(par_bootstr_vac)),function(x){efficacy(d,par_bootstr_vac[x,])}),
                                           probs = c((1-alpha_CI)/2,(1+alpha_CI)/2))}
  vac_ci <- sapply(doses,function(x){eff_ci_from_dose_vac(x)})
  
  # compute efficacy CIs for a given dose for vaccines:
  eff_ci_from_dose_mAb <- function(d){quantile(sapply(c(1:nrow(par_bootstr_mAb)),function(x){efficacy(d,par_bootstr_mAb[x,])}),
                                               probs = c((1-alpha_CI)/2,(1+alpha_CI)/2))}
  mAb_ci <- sapply(doses,function(x){eff_ci_from_dose_mAb(x)})
  
  # add model fit and CIs for vaccine and mAb data:
  models_vac_mAb[[j]] <- data.frame(dose=doses,
                                    fit_vac=eff_vac,
                                    fit_vac_lower=vac_ci[1,],
                                    fit_vac_upper=vac_ci[2,],
                                    fit_mAb=eff_mAb,
                                    fit_mAb_lower=mAb_ci[1,],
                                    fit_mAb_upper=mAb_ci[2,])
}

# save the transformed best fit parameters with CIs for each model:
if(use_efficacy_model=="logistic with max 1"){
  param_best_transf_vac_mAb <- param_best_vac_mAb
  for(j in c(1:length(param_best_vac_mAb))){
    param_best_transf_vac_mAb[[j]] <- data.frame(exp(param_best_transf_vac_mAb[[j]]))
    names(param_best_transf_vac_mAb[[j]]) <- c("estimate","CI.lower","CI.upper")
  }
}

# save data:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(param_best_transf_vac_mAb,glue::glue("output/Best-fit-vaccine-mAb-parameters-with-CIs_{Sys.Date()}_{cur_time}.xlsx"))
save(models_vac_mAb,param_best_vac_mAb,AICs_vac_mAb,param_best_transf_vac_mAb,
     file=glue::glue("output/MAb-vaccine-fitting-results_{Sys.Date()}_{cur_time}.RData"))


# visualisation of best fits ----------------------------------------------

# visualize all model fits with CIs:
p.vac.mAbs.all.best.fits <- (vis.vac.mAb.fit("k.ic50") |  vis.vac.mAb.fit("k")) / 
  (vis.vac.mAb.fit("ic50") |  vis.vac.mAb.fit(""))
p.vac.mAbs.all.best.fits

# visualize the best model fit:
p.vac.mAbs.best.model <- vis.vac.mAb.fit("ic50")
p.vac.mAbs.best.model

# vaccine and mAbs comparison with best model fit:
p.vac.mAbs <- p.vaccine.mabs.total | p.vac.mAbs.best.model
p.vac.mAbs

# save figures:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
ggsave(glue::glue("output/FigS_vaccine-mAbs-all-models_{Sys.Date()}_{cur_time}.pdf"),height=17,width=27,units="cm",plot=p.vac.mAbs.all.best.fits)
ggsave(glue::glue("output/Fig4b_vaccine-mAbs-best-model_{Sys.Date()}_{cur_time}.pdf"),width=7,height=4,units="in",plot=p.vac.mAbs.best.model)
ggsave(glue::glue("output/Fig4_vaccine-mAbs-comparison_{Sys.Date()}_{cur_time}.pdf"),height=10,width=30,units="cm",plot=p.vac.mAbs)


# cleanup -----------------------------------------------------------------

rm(data_mAbs_vaccines,data_mab_vac_stats,tmp.events,model_mAbs_vac,comparison.matrix,pval_mAbs_vs_vac,tmp.glmer,
   av_eff_Abs,av_eff_Abs_ci,av_eff_vac,av_eff_vac_ci,cur_time,my.colours,my.shapes,p.vaccine.mabs.total,
   data_vaccine,combined_data_vac_mAbs,combined_data_vac_mAbs_fit,LogisticModel_withPoolSD,doses,model_fit,
   model_ci,model_data,p.vaccine.mabs.best.fit,same.par,baseline_risks,RandomInitialEstimate_vac,RandomInitialEstimate_mAb,
   all_hessians,i,j,same.par.ind,different.par.ind,tempfit,rand.vac,rand.mAb,rand.diff.mAb,par_initial,par_est_vac,
   par_est_diff_mAb,par_est_mAb,output,output_vac,output_mAb,param_best_fit_vac,param_best_fit_mAb,hessian_best_fit_vac_mAb,
   lr_comparison_pvals,doses,AICs_vac_mAb,covar,param_ci,eff_vac,eff_mAb,param_tmp,par_bootstr,par_bootstr_vac,par_bootstr_diff,
   par_bootstr_mAb,eff_ci_from_dose_vac,vac_ci,eff_ci_from_dose_mAb,mAb_ci,p.vac.mAbs.all.best.fits,p.vac.mAbs.best.model,
   p.vac.mAbs)
