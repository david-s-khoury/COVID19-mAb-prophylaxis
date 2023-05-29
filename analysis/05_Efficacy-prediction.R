# -------------------------------------------------------------------------
#' 
#' Prediction of efficacy and duration of protection of prophylactic mAb
#' treatment
#' 
# -------------------------------------------------------------------------


# parameters for analysis -------------------------------------------------

# choose variants for which to predict efficacy:
variants <- c("Wild Type","Delta","Omicron/BA.1","Omicron/BA.2","Omicron/BQ.1")
max.fold.change <- 300 # consider variants with fold-IC50 changes between 1 and max.fold.change


# estimate the half-life of mAbs ------------------------------------------

# peak concentration and time (in mg/L):
peak_evu <- max(data_conc$conc[data_conc$drug=="cilgavimab + tixagevimab"])
peak_time_evu <- data_conc$day[data_conc$drug=="cilgavimab + tixagevimab" & data_conc$conc==peak_evu]
peak_ron <- max(data_conc$conc[data_conc$drug=="casirivimab + imdevimab" & data_conc$trial=="O'Brien & Herman"])
peak_time_ron <- data_conc$day[data_conc$drug=="casirivimab + imdevimab" & data_conc$conc==peak_ron]
peak_adi <- max(data_conc$conc[data_conc$drug=="adintrevimab"])
peak_time_adi <- min(data_conc$day[data_conc$drug=="adintrevimab" & data_conc$conc==peak_adi])

# convert peak concentration to fold-ic50:
peak_evu_ic50 <- peak_evu/(1e-3*Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==unique(data_conc$variant[data_conc$drug=="cilgavimab + tixagevimab"])])
peak_ron_ic50 <- peak_ron/(1e-3*Means_table$IC50[Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==unique(data_conc$variant[data_conc$drug=="casirivimab + imdevimab"])])
peak_adi_ic50 <- peak_adi/(1e-3*Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory=="Wild Type"])
peak_adi_delta_ic50 <- peak_adi/(1e-3*Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory=="Delta"])
peak_adi_omicron_ic50 <- peak_adi/(1e-3*Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory=="Omicron/BA.1"])

# Determine the half-life by fitting to the concentration data from peak:
data_conc <- dplyr::mutate(data_conc, conc.log=log(conc))

# fit a linear model to log-transformed concentration data from peak
lm_ron <- lm(conc.log~day,data=data_conc[data_conc$trial=="O'Brien & Herman" & data_conc$day>=data_conc$day[data_conc$conc==max(data_conc$conc[data_conc$trial=="O'Brien & Herman"])],])
lm_evu <- lm(conc.log~day,data=data_conc[data_conc$trial=="Levin" & is.finite(data_conc$conc.log) & data_conc$day>=data_conc$day[data_conc$conc==max(data_conc$conc[data_conc$trial=="Levin"])],])
lm_adi <- lm(conc.log~day,data=data_conc[data_conc$trial=="Schmidt" & is.finite(data_conc$conc.log) & data_conc$day>=min(data_conc$day[data_conc$conc==max(data_conc$conc[data_conc$trial=="Schmidt"])]),])
half_life_ron_est <- as.numeric(-log(2)/lm_ron$coefficients[2])
half_life_ron_est_ci <- as.numeric(-log(2)/confint(lm_ron)[2,])
half_life_evu_est <- as.numeric(-log(2)/lm_evu$coefficients[2])
half_life_evu_est_ci <- as.numeric(-log(2)/confint(lm_evu)[2,])
half_life_adi_est <- as.numeric(-log(2)/lm_adi$coefficients[2])
half_life_adi_est_ci <- as.numeric(-log(2)/confint(lm_adi)[2,])

# compare the fit to the data with the half-lives from the literature:
half_life_evu <- mean(c(84,89)) # mean of values of tixagevimab and cilgavimab from the Australian Product Information for evusheld 
half_life_ron <- mean(c(31.8,26.9)) # mean of values of casirivimab and imdevimab from the Regeneron fact sheet
half_life_adi <- 141

# make a figure with the concentration data, estimated half-life, and half-life from the literature:
p.conc.half.life.combined <- fun.ab.conc.half.life(ab_name="casirivimab + imdevimab",lin_mod=lm_ron,peak_conc=peak_ron,peak_time=peak_time_ron,half_life=half_life_ron,data_conc) | 
  fun.ab.conc.half.life(ab_name="cilgavimab + tixagevimab",lin_mod=lm_evu,peak_conc=peak_evu,peak_time=peak_time_evu,half_life=half_life_evu,data_conc) | 
  fun.ab.conc.half.life(ab_name="adintrevimab",lin_mod=lm_adi,peak_conc=peak_adi,peak_time=peak_time_adi,half_life=half_life_adi,data_conc)
p.conc.half.life.combined

# save figure:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
ggsave(glue::glue("output/Fig_conc-and-half-life_{Sys.Date()}_{cur_time}.pdf"),width=13.5,height=4,units="in",plot=p.conc.half.life.combined)

# make a table of peak concentration, peak time, estimated half-life, estimated half-life CI, and half-life from the literature:
table_conc_kinetics <- data.frame(trial=c("Schmidt","O'Brien & Herman","Levin"),
                                  mAb=c("adintrevimab","casirivimab + imdevimab","cilgavimab + tixagevimab"),
                                  peak_conc=c(peak_adi,peak_ron,peak_evu),
                                  peak_day=c(peak_time_adi,peak_time_ron,peak_time_evu),
                                  half_life_est=c(half_life_adi_est,half_life_ron_est,half_life_evu_est),
                                  half_life_est_lower=c(half_life_adi_est_ci[1],half_life_ron_est_ci[1],half_life_evu_est_ci[1]),
                                  half_life_est_upper=c(half_life_adi_est_ci[2],half_life_ron_est_ci[2],half_life_evu_est_ci[2]),
                                  half_life_literature=c(half_life_adi,half_life_ron,half_life_evu))

# save peak efficacy table:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(table_conc_kinetics,glue::glue("output/Table-concentration-kinetics_{Sys.Date()}_{cur_time}.xlsx"))


# peak concentration efficacy ---------------------------------------------

# make a table with the peak concentration efficacy against different variants:
peak_efficacies <- data.frame(variant=variants,
                              evu=100*sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==variants[x]),efficacy(peak_evu/(1e-3*Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==variants[x]]),param_best),NA)}),
                              evu.lower=100*sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==variants[x]),eff_ci_from_dose(peak_evu/(1e-3*Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==variants[x]]))[1],NA)}),
                              evu.upper=100*sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==variants[x]),eff_ci_from_dose(peak_evu/(1e-3*Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==variants[x]]))[2],NA)}),
                              ron=100*sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==variants[x]),efficacy(peak_ron/(1e-3*Means_table$IC50[Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==variants[x]]),param_best),NA)}),
                              ron.lower=100*sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==variants[x]),eff_ci_from_dose(peak_ron/(1e-3*Means_table$IC50[Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==variants[x]]))[1],NA)}),
                              ron.upper=100*sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==variants[x]),eff_ci_from_dose(peak_ron/(1e-3*Means_table$IC50[Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==variants[x]]))[2],NA)}),
                              adi=100*sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory==variants[x]),efficacy(peak_adi/(1e-3*Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory==variants[x]]),param_best),NA)}),
                              adi.lower=100*sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory==variants[x]),eff_ci_from_dose(peak_adi/(1e-3*Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory==variants[x]]))[1],NA)}),
                              adi.upper=100*sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory==variants[x]),eff_ci_from_dose(peak_adi/(1e-3*Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory==variants[x]]))[2],NA)}))

# save peak efficacy table:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(peak_efficacies,glue::glue("output/Peak-concentration-efficacy_{Sys.Date()}_{cur_time}.xlsx"))


# in vivo IC50 concentration for mAbs -------------------------------------

# compute the in vivo mAb concentration (in mg/L) that gives 50% protective efficacy
# use the same variants as for the peak concentration efficacy

ic50 <- table_par_best$estimate[table_par_best$parameter=="IC50"] # best fit parameter: efficacy IC50
ic50_ci <- c(table_par_best$lower[table_par_best$parameter=="IC50"],table_par_best$upper[table_par_best$parameter=="IC50"])

# concentration for 50% protection (in mg/L) with CIs and peak in vivo concentration for comparison:
conc_50_protection <- data.frame(variant=variants,
                                 evu=sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==variants[x]),ic50*1e-3*Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==variants[x]],NA)}),
                                 evu.lower=sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==variants[x]),ic50_ci[1]*1e-3*Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==variants[x]],NA)}),
                                 evu.upper=sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==variants[x]),ic50_ci[2]*1e-3*Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory==variants[x]],NA)}),
                                 evu.peak=rep(peak_evu,length(variants)),
                                 ron=sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==variants[x]),ic50*1e-3*Means_table$IC50[Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==variants[x]],NA)}),
                                 ron.lower=sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==variants[x]),ic50_ci[1]*1e-3*Means_table$IC50[Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==variants[x]],NA)}),
                                 ron.upper=sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==variants[x]),ic50_ci[2]*1e-3*Means_table$IC50[Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory==variants[x]],NA)}),
                                 ron.peak=rep(peak_ron,length(variants)),
                                 adi=sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory==variants[x]),ic50*1e-3*Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory==variants[x]],NA)}),
                                 adi.lower=sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory==variants[x]),ic50_ci[1]*1e-3*Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory==variants[x]],NA)}),
                                 adi.upper=sapply(c(1:length(variants)),function(x){ifelse(any(Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory==variants[x]),ic50_ci[2]*1e-3*Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory==variants[x]],NA)}),
                                 adi.peak=rep(peak_adi,length(variants)))

# save data frame:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(conc_50_protection,glue::glue("output/Concentration-50-perc-protection_{Sys.Date()}_{cur_time}.xlsx"))


# time above an efficacy threshold --------------------------------------- 

# compute and visualize the time (since administration) that the protective efficacy is above a
# protective efficacy threshold for different IC50-fold-changes

# compute the time above 50% protection with CIs:
table_time_above_50_eff <- data.frame(variant=variants,
                                      evu=sapply(c(1:length(variants)),function(x){time.above.limit(ab_name="cilgavimab + tixagevimab",variant=variants[x],0.5,peak_evu,lm_evu,param_best,Means_table,ic50.fold.change=1)}),
                                      evu.lower=sapply(c(1:length(variants)),function(x){time.above.limit(ab_name="cilgavimab + tixagevimab",variant=variants[x],0.5,peak_evu,lm_evu,param_best,Means_table,ic50.fold.change=1,CI=TRUE,eff_IC50=max(ic50_ci))}),
                                      evu.upper=sapply(c(1:length(variants)),function(x){time.above.limit(ab_name="cilgavimab + tixagevimab",variant=variants[x],0.5,peak_evu,lm_evu,param_best,Means_table,ic50.fold.change=1,CI=TRUE,eff_IC50=min(ic50_ci))}),
                                      ron=sapply(c(1:length(variants)),function(x){time.above.limit(ab_name="casirivimab + imdevimab",variant=variants[x],0.5,peak_evu,lm_evu,param_best,Means_table,ic50.fold.change=1)}),
                                      ron.lower=sapply(c(1:length(variants)),function(x){time.above.limit(ab_name="casirivimab + imdevimab",variant=variants[x],0.5,peak_evu,lm_evu,param_best,Means_table,ic50.fold.change=1,CI=TRUE,eff_IC50=max(ic50_ci))}),
                                      ron.upper=sapply(c(1:length(variants)),function(x){time.above.limit(ab_name="casirivimab + imdevimab",variant=variants[x],0.5,peak_evu,lm_evu,param_best,Means_table,ic50.fold.change=1,CI=TRUE,eff_IC50=min(ic50_ci))}),
                                      adi=sapply(c(1:length(variants)),function(x){time.above.limit(ab_name="adintrevimab",variant=variants[x],0.5,peak_evu,lm_evu,param_best,Means_table,ic50.fold.change=1)}),
                                      adi.lower=sapply(c(1:length(variants)),function(x){time.above.limit(ab_name="adintrevimab",variant=variants[x],0.5,peak_evu,lm_evu,param_best,Means_table,ic50.fold.change=1,CI=TRUE,eff_IC50=max(ic50_ci))}),
                                      adi.upper=sapply(c(1:length(variants)),function(x){time.above.limit(ab_name="adintrevimab",variant=variants[x],0.5,peak_evu,lm_evu,param_best,Means_table,ic50.fold.change=1,CI=TRUE,eff_IC50=min(ic50_ci))}))

# save data frame:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(table_time_above_50_eff,glue::glue("output/Time-above-50-perc-protection_{Sys.Date()}_{cur_time}.xlsx"))


# Visualization:
fold.changes <- seq(1,max.fold.change,length.out=1e3)

data.vis <- data.frame(mAb=c(rep("cilgavimab + tixagevimab",length(fold.changes)),rep("casirivimab + imdevimab",length(fold.changes)),
                             rep("adintrevimab",length(fold.changes)),rep("adintrevimab (delta)",length(fold.changes)),rep("adintrevimab (omicron)",length(fold.changes))),
                       fold=rep(fold.changes,5),
                       time=c(time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Wild Type",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = fold.changes),
                              time.above.limit(ab_name = "casirivimab + imdevimab",variant="Wild Type",eff_thresh = 0.5,peak_conc = peak_ron,lin_mod = lm_ron,par_eff = param_best,Means_table,ic50.fold.change = fold.changes),
                              time.above.limit(ab_name = "adintrevimab",variant="Wild Type",eff_thresh = 0.5,peak_conc = peak_adi,lin_mod = lm_adi,par_eff = param_best,Means_table,ic50.fold.change = fold.changes),
                              time.above.limit(ab_name = "adintrevimab",variant="Delta",eff_thresh = 0.5,peak_conc = peak_adi,lin_mod = lm_adi,par_eff = param_best,Means_table,ic50.fold.change = fold.changes),
                              time.above.limit(ab_name = "adintrevimab",variant="Omicron/BA.1",eff_thresh = 0.5,peak_conc = peak_adi,lin_mod = lm_adi,par_eff = param_best,Means_table,ic50.fold.change = fold.changes)),
                       time.lower=c(time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Wild Type",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = fold.changes,CI=TRUE,eff_IC50 = max(ic50_ci)),
                                    time.above.limit(ab_name = "casirivimab + imdevimab",variant="Wild Type",eff_thresh = 0.5,peak_conc = peak_ron,lin_mod = lm_ron,par_eff = param_best,Means_table,ic50.fold.change = fold.changes,CI=TRUE,eff_IC50 = max(ic50_ci)),
                                    time.above.limit(ab_name = "adintrevimab",variant="Wild Type",eff_thresh = 0.5,peak_conc = peak_adi,lin_mod = lm_adi,par_eff = param_best,Means_table,ic50.fold.change = fold.changes,CI=TRUE,eff_IC50 = max(ic50_ci)),
                                    time.above.limit(ab_name = "adintrevimab",variant="Delta",eff_thresh = 0.5,peak_conc = peak_adi,lin_mod = lm_adi,par_eff = param_best,Means_table,ic50.fold.change = fold.changes,CI=TRUE,eff_IC50 = max(ic50_ci)),
                                    time.above.limit(ab_name = "adintrevimab",variant="Omicron/BA.1",eff_thresh = 0.5,peak_conc = peak_adi,lin_mod = lm_adi,par_eff = param_best,Means_table,ic50.fold.change = fold.changes,CI=TRUE,eff_IC50 = max(ic50_ci))),
                       time.upper=c(time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Wild Type",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = fold.changes,CI=TRUE,eff_IC50 = min(ic50_ci)),
                                    time.above.limit(ab_name = "casirivimab + imdevimab",variant="Wild Type",eff_thresh = 0.5,peak_conc = peak_ron,lin_mod = lm_ron,par_eff = param_best,Means_table,ic50.fold.change = fold.changes,CI=TRUE,eff_IC50 = min(ic50_ci)),
                                    time.above.limit(ab_name = "adintrevimab",variant="Wild Type",eff_thresh = 0.5,peak_conc = peak_adi,lin_mod = lm_adi,par_eff = param_best,Means_table,ic50.fold.change = fold.changes,CI=TRUE,eff_IC50 = min(ic50_ci)),
                                    time.above.limit(ab_name = "adintrevimab",variant="Delta",eff_thresh = 0.5,peak_conc = peak_adi,lin_mod = lm_adi,par_eff = param_best,Means_table,ic50.fold.change = fold.changes,CI=TRUE,eff_IC50 = min(ic50_ci)),
                                    time.above.limit(ab_name = "adintrevimab",variant="Omicron/BA.1",eff_thresh = 0.5,peak_conc = peak_adi,lin_mod = lm_adi,par_eff = param_best,Means_table,ic50.fold.change = fold.changes,CI=TRUE,eff_IC50 = min(ic50_ci))))

change.evu.BA1 <- Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory=="Omicron/BA.1"]/Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory=="Wild Type"]
change.evu.BA2 <- Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory=="Omicron/BA.2"]/Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory=="Wild Type"]

data.fold.drops <- data.frame(fold.name=c("1","5","10","50","Omicron/BA.1 (cil.+tix.)","Omicron/BA.2 (cil.+tix.)"),
                              fold=c(1,5,10,50,change.evu.BA1,change.evu.BA2),
                              time.evu=c(time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Wild Type",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = c(1,5,10,50,change.evu.BA1,change.evu.BA2))),
                              time.ron=c(time.above.limit(ab_name = "casirivimab + imdevimab",variant="Wild Type",eff_thresh = 0.5,peak_conc = peak_ron,lin_mod = lm_ron,par_eff = param_best,Means_table,ic50.fold.change = c(1,5,10,50)),NA,NA),
                              time.adi=c(time.above.limit(ab_name = "adintrevimab",variant="Wild Type",eff_thresh = 0.5,peak_conc = peak_adi,lin_mod = lm_adi,par_eff = param_best,Means_table,ic50.fold.change = c(1,5,10,50)),NA,NA),
                              time.adi.delta=c(time.above.limit(ab_name = "adintrevimab",variant="Delta",eff_thresh = 0.5,peak_conc = peak_adi,lin_mod = lm_adi,par_eff = param_best,Means_table,ic50.fold.change = c(1,5,10,50)),NA,NA),
                              time.adi.omicron=c(time.above.limit(ab_name = "adintrevimab",variant="Omicron/BA.1",eff_thresh = 0.5,peak_conc = peak_adi,lin_mod = lm_adi,par_eff = param_best,Means_table,ic50.fold.change = c(1,5,10,50)),NA,NA))

my.colours.fold.escape <- c("black","olivedrab4","olivedrab3","olivedrab1","red4","red")
my.colours <- c("dodgerblue","darkorange","mediumpurple")
vis.x.max <- 300

# plot the days above 50% protection by level of escape [fold-increase in the IC50]
p.variants <- ggplot(data.vis[data.vis$mAb%in%c("cilgavimab + tixagevimab","casirivimab + imdevimab","adintrevimab"),],aes(x=fold,y=time,color=mAb)) +
   geom_ribbon(inherit.aes=FALSE,aes(x=fold,y=time,ymin=time.lower,ymax=time.upper,fill=mAb),alpha = 0.1) +
   geom_line() +
   # add points showing different fold-drops:
   geom_point(data=data.fold.drops,inherit.aes = FALSE, aes(x=fold,y=time.evu,color=factor(fold)), shape=16) +
   geom_point(data=data.fold.drops,inherit.aes = FALSE, aes(x=fold,y=time.ron,color=factor(fold)), shape=16) +
   geom_point(data=data.fold.drops,inherit.aes = FALSE, aes(x=fold,y=time.adi,color=factor(fold)), shape=16) +
   # geom_point(data=data.fold.drops,inherit.aes = FALSE, aes(x=fold,y=time.adi.omicron,color=factor(fold)), shape=16) +
   annotate("rect",xmin=1,xmax=vis.x.max,ymin=-10,ymax=30,alpha=1,fill="gray95") +
   # annotate("segment",x=c(0,change.evu.BA1),y=c(time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Omicron/BA.1",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = 1),0),xend=c(change.evu.BA1,change.evu.BA1),
   #          yend=c(time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Omicron/BA.1",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = 1),time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Omicron/BA.1",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = 1)),
   #          color="red4",linetype="dashed") +
   annotate("segment",x=c(0,change.evu.BA2),y=c(time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Omicron/BA.2",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = 1),0),xend=c(change.evu.BA2,change.evu.BA2),
            yend=c(time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Omicron/BA.2",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = 1),time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Omicron/BA.2",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = 1)),
            color="red",linetype="dashed") +
   # annotate("text",x=1.5,y=100,label=paste(round(time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Omicron/BA.1",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = 1),1)," days",sep=""),size = 3,color="red4") +
   annotate("text",x=1.5,y=10+ceiling(time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Omicron/BA.2",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = 1)),
            label=paste(round(time.above.limit(ab_name = "cilgavimab + tixagevimab",variant="Omicron/BA.2",eff_thresh = 0.5,peak_conc = peak_evu,lin_mod = lm_evu,par_eff = param_best,Means_table,ic50.fold.change = 1),1)," days",sep=""),size = 3,color="red") +
   # axis, theme, etc.:
   scale_colour_manual(values = c('casirivimab + imdevimab'=my.colours[1],'cilgavimab + tixagevimab'=my.colours[2],
                                  'adintrevimab (delta)'=my.colours[3],#'adintrevimab (omicron)'=my.colours[4],
                                  'Ancestral'=my.colours.fold.escape[1],'5-fold'=my.colours.fold.escape[2],'10-fold'=my.colours.fold.escape[3],
                                  '50-fold'=my.colours.fold.escape[4],'BA.1 (for cil.+tix.)'=my.colours.fold.escape[5],
                                  'BA.2 (for cil.+tix.)'=my.colours.fold.escape[6]), name="",
                       guide = guide_legend(byrow=TRUE,ncol=2,title="mAb:",title.vjust = 1,
                                            override.aes = list(linetype = c(rep("solid",3),rep("blank",6)),
                                                                fill = c(my.colours[c(1:3)],rep("white",6)),
                                                                shape = c(rep(NA,3),rep(16,6))))) +
   scale_fill_manual(values=my.colours[c(3,1,2)]) +
   scale_y_continuous(name = "Days above 50% protection",breaks=c(30,seq(0,1100,by=100)),minor_breaks = seq(50,1200,by=50),expand = c(0,0)) +
   scale_x_log10(breaks=c(seq(1,10,by=1),seq(20,100,by=10),150,200,250,300),minor_breaks=c(seq(1,10,by=1),15,seq(20,100,by=5),150,200,250,300),expand = c(1e-2,1e-2)) +
   coord_cartesian(xlim = c(1,vis.x.max),ylim=c(-5,1010)) +
   labs(x="Level of escape [fold-increase in the IC50]", title = "Time from treatment to 50% protection") +
   theme_bw() +
   guides(fill="none") +
   theme(legend.key.size = unit(0.4, 'cm'),panel.grid = element_line(colour="gray95",size = 0.2),legend.position = 'bottom')

p.variants

# save plot:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
ggsave(glue::glue("output/Fig_time-above-50-perc-efficacy-by_fold-ic50-change_{Sys.Date()}_{cur_time}.pdf"),height=12,width=15,units="cm",plot=p.variants)


# visualize protective efficacy over time ---------------------------------

# parameters or visualization:
day_max_ron <- 1e3 #500
day_max_evu <- 1e3
day_max_adi <- 1e3 #2.5e3

p.efficacy.by.time.combined <- ( vis.eff.by.time(day_max = day_max_ron, ab_name = "casirivimab + imdevimab", lin_mod = lm_ron) | 
                                    vis.eff.by.time(day_max = day_max_evu, ab_name = "cilgavimab + tixagevimab", lin_mod = lm_evu)) / 
   ( vis.eff.by.time(day_max = day_max_adi, ab_name = "adintrevimab", lin_mod = lm_adi) | p.variants )

p.efficacy.by.time.combined

# save figure:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
ggsave(glue::glue("output/Fig3_efficacy-by-time-with-time-above-50-percent-protection_{Sys.Date()}_{cur_time}.pdf"),height=17,width=20,units="cm",plot=p.efficacy.by.time.combined)


# fold-increase in IC50 that leads to 30 days of 50% protection -----------

# which fold-change in the mAb IC50 leads to 30 days of 50% protective efficacy and fold-changes for different variants: 
table_ic50_change <- data.frame(mAb=c("adintrevimab","casirivimab + imdevimab","cilgavimab + tixagevimab"),
                                fold.change=c(fun.fold.change(ab_name = "adintrevimab", lin_mod = lm_adi),
                                              fun.fold.change(ab_name = "casirivimab + imdevimab", lin_mod = lm_ron),
                                              fun.fold.change(ab_name = "cilgavimab + tixagevimab", lin_mod = lm_evu)),
                                fold.change.lower=c(fun.fold.change(ab_name = "adintrevimab", lin_mod = lm_adi, CI=TRUE, eff_IC50 = max(ic50_ci)),
                                                    fun.fold.change(ab_name = "casirivimab + imdevimab", lin_mod = lm_ron, CI=TRUE, eff_IC50 = max(ic50_ci)),
                                                    fun.fold.change(ab_name = "cilgavimab + tixagevimab", lin_mod = lm_evu, CI=TRUE, eff_IC50 = max(ic50_ci))),
                                fold.change.upper=c(fun.fold.change(ab_name = "adintrevimab", lin_mod = lm_adi, CI=TRUE, eff_IC50 = min(ic50_ci)),
                                                    fun.fold.change(ab_name = "casirivimab + imdevimab", lin_mod = lm_ron, CI=TRUE, eff_IC50 = min(ic50_ci)),
                                                    fun.fold.change(ab_name = "cilgavimab + tixagevimab", lin_mod = lm_evu, CI=TRUE, eff_IC50 = min(ic50_ci))),
                                change.Omicron.BA.1=c(Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory=="Omicron/BA.1"]/Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory=="Wild Type"],
                                                      Means_table$IC50[Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory=="Omicron/BA.1"]/Means_table$IC50[Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory=="Wild Type"],
                                                      Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory=="Omicron/BA.1"]/Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory=="Wild Type"]),
                                change.Omicron.BA.2=c(Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory=="Omicron/BA.2"]/Means_table$IC50[Means_table$Antibodies=="adintrevimab" & Means_table$VariantCategory=="Wild Type"],
                                                      Means_table$IC50[Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory=="Omicron/BA.2"]/Means_table$IC50[Means_table$Antibodies=="casirivimab + imdevimab" & Means_table$VariantCategory=="Wild Type"],
                                                      Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory=="Omicron/BA.2"]/Means_table$IC50[Means_table$Antibodies=="cilgavimab + tixagevimab" & Means_table$VariantCategory=="Wild Type"]))

# save table:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
write_xlsx(table_ic50_change,glue::glue("output/Table-fold-IC50-change-30days-50perc-protection_{Sys.Date()}_{cur_time}.xlsx"))


# cleanup -----------------------------------------------------------------

rm(variants,max.fold.change,peak_evu,peak_time_evu,peak_ron,peak_time_ron,peak_adi,peak_time_adi,
   peak_evu_ic50,peak_ron_ic50,peak_adi_delta_ic50,peak_adi_omicron_ic50,lm_ron,lm_evu,
   lm_adi,half_life_ron_est,half_life_ron_est_ci,half_life_evu_est,half_life_evu_est_ci,half_life_adi_est,
   half_life_adi_est_ci,half_life_evu,half_life_ron,half_life_adi,p.conc.half.life.combined,cur_time,
   table_conc_kinetics,peak_efficacies,ic50,ic50_ci,conc_50_protection,table_time_above_50_eff,data.vis,
   change.evu.BA1,change.evu.BA2,data.fold.drops,my.colours.fold.escape,my.colours,vis.x.max,p.variants,
   day_max_ron,day_max_evu,day_max_adi,p.efficacy.by.time.combined,table_ic50_change)
