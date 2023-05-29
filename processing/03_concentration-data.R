# -------------------------------------------------------------------------
#' 
#' Process data of in vivo mAb concentration of prophylactic mAb treatment
#' of SARS-CoV-2
#' 
# -------------------------------------------------------------------------


# load data ---------------------------------------------------------------

# load concentration data and compute the geometric mean concentration:
obrien_conc <- read_excel("raw-data/mAbs in vivo concentration_FINAL_2023-05-16.xlsx",sheet = "O'Brien")
herman_conc <- read_excel("raw-data/mAbs in vivo concentration_FINAL_2023-05-16.xlsx",sheet = "Herman")
isa_conc <- read_excel("raw-data/mAbs in vivo concentration_FINAL_2023-05-16.xlsx",sheet = "Isa")
levin_conc <- read_excel("raw-data/mAbs in vivo concentration_FINAL_2023-05-16.xlsx",sheet = "Levin")
schmidt_conc <- read_excel("raw-data/mAbs in vivo concentration_FINAL_2023-05-16.xlsx",sheet = "Schmidt")


# compute the geometric mean concentration --------------------------------

# O'Brien study:
geomean_conc_obrien <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="O'Brien",])),
                              function(x){fun_geomean_conc(obrien_conc$Day_corrected[obrien_conc$Drug=="casirivimab"],
                                                           obrien_conc$Conc_geomean[obrien_conc$Drug=="casirivimab"]+obrien_conc$Conc_geomean[obrien_conc$Drug=="imdevimab"],
                                                           data_efficacy$Day_min[which(data_efficacy$trial=="O'Brien")[x]],data_efficacy$Day_max[which(data_efficacy$trial=="O'Brien")[x]])})
obrien_conc_daily <- approx(obrien_conc$Day_corrected[obrien_conc$Drug=="casirivimab"],
                            obrien_conc$Conc_geomean[obrien_conc$Drug=="casirivimab"]+obrien_conc$Conc_geomean[obrien_conc$Drug=="imdevimab"],
                            xout=c(0:max(obrien_conc$Day_corrected[obrien_conc$Drug=="casirivimab"])))$y
conc_min_obrien <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="O'Brien",])),function(x){min(obrien_conc_daily[c(data_efficacy$Day_min[which(data_efficacy$trial=="O'Brien")[x]]:data_efficacy$Day_max[which(data_efficacy$trial=="O'Brien")[x]])])})
conc_max_obrien <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="O'Brien",])),function(x){max(obrien_conc_daily[c(data_efficacy$Day_min[which(data_efficacy$trial=="O'Brien")[x]]:data_efficacy$Day_max[which(data_efficacy$trial=="O'Brien")[x]])])})

# Herman study:
# use O'Brien concentration data where possible (as in vivo measurement), then use Herman concentration data
geomean_conc_herman <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Herman",])),
                              function(x){fun_geomean_conc(c(obrien_conc$Day_corrected[obrien_conc$Drug=="casirivimab"],herman_conc$Day[herman_conc$Day>max(obrien_conc$Day_corrected)]),
                                                           c(obrien_conc$Conc_geomean[obrien_conc$Drug=="casirivimab"]+obrien_conc$Conc_geomean[obrien_conc$Drug=="imdevimab"],herman_conc$Conc_geomean[herman_conc$Day>max(obrien_conc$Day_corrected)]),
                                                           data_efficacy$Day_min[which(data_efficacy$trial=="Herman")[x]],data_efficacy$Day_max[which(data_efficacy$trial=="Herman")[x]])})
herman_conc_daily <- approx(c(obrien_conc$Day_corrected[obrien_conc$Drug=="casirivimab"],herman_conc$Day[herman_conc$Day>max(obrien_conc$Day_corrected)]),
                            c(obrien_conc$Conc_geomean[obrien_conc$Drug=="casirivimab"]+obrien_conc$Conc_geomean[obrien_conc$Drug=="imdevimab"],herman_conc$Conc_geomean[herman_conc$Day>max(obrien_conc$Day_corrected)]),
                            xout=c(0:max(herman_conc$Day_corrected)))$y
conc_min_herman <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Herman",])),function(x){min(herman_conc_daily[c(data_efficacy$Day_min[which(data_efficacy$trial=="Herman")[x]]:data_efficacy$Day_max[which(data_efficacy$trial=="Herman")[x]])])})
conc_max_herman <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Herman",])),function(x){max(herman_conc_daily[c(data_efficacy$Day_min[which(data_efficacy$trial=="Herman")[x]]:data_efficacy$Day_max[which(data_efficacy$trial=="Herman")[x]])])})

# Isa study:
geomean_conc_isa <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Isa",])),
                           function(x){fun_geomean_conc(isa_conc$Day[isa_conc$Drug=="casirivimab + imdevimab"],isa_conc$Conc_geomean[isa_conc$Drug=="casirivimab + imdevimab"],
                                                        data_efficacy$Day_min[which(data_efficacy$trial=="Isa")[x]],data_efficacy$Day_max[which(data_efficacy$trial=="Isa")[x]])})
isa_conc_daily <- approx(isa_conc$Day[isa_conc$Drug=="casirivimab + imdevimab"],isa_conc$Conc_geomean[isa_conc$Drug=="casirivimab + imdevimab"],
                         xout=c(0:max(isa_conc$Day_corrected)))$y
conc_min_isa <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Isa",])),function(x){min(isa_conc_daily[c(data_efficacy$Day_min[which(data_efficacy$trial=="Isa")[x]]:data_efficacy$Day_max[which(data_efficacy$trial=="Isa")[x]])])})
conc_max_isa <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Isa",])),function(x){max(isa_conc_daily[c(data_efficacy$Day_min[which(data_efficacy$trial=="Isa")[x]]:data_efficacy$Day_max[which(data_efficacy$trial=="Isa")[x]])])})

# Levin study:
geomean_conc_levin <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Levin",])),
                             function(x){fun_geomean_conc(levin_conc$Day[levin_conc$Drug=="cilgavimab + tixagevimab"],levin_conc$Conc_geomean[levin_conc$Drug=="cilgavimab + tixagevimab"],
                                                          data_efficacy$Day_min[which(data_efficacy$trial=="Levin")[x]],data_efficacy$Day_max[which(data_efficacy$trial=="Levin")[x]])})
levin_conc_daily <- approx(levin_conc$Day[levin_conc$Drug=="cilgavimab + tixagevimab"],levin_conc$Conc_geomean[levin_conc$Drug=="cilgavimab + tixagevimab"],
                             xout=c(0:max(levin_conc$Day_corrected)))$y
conc_min_levin <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Levin",])),function(x){min(levin_conc_daily[c(data_efficacy$Day_min[which(data_efficacy$trial=="Levin")[x]]:data_efficacy$Day_max[which(data_efficacy$trial=="Levin")[x]])])})
conc_max_levin <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Levin",])),function(x){max(levin_conc_daily[c(data_efficacy$Day_min[which(data_efficacy$trial=="Levin")[x]]:data_efficacy$Day_max[which(data_efficacy$trial=="Levin")[x]])])})

# Schmidt study (delta & omicron):
geomean_conc_schmidt_delta <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Schmidt (delta)",])),
                                     function(x){fun_geomean_conc(schmidt_conc$Day_corrected,schmidt_conc$Conc_geomean,
                                                                  data_efficacy$Day_min[which(data_efficacy$trial=="Schmidt (delta)")[x]],data_efficacy$Day_max[which(data_efficacy$trial=="Schmidt (delta)")[x]])})
schmidt_conc_daily <- approx(schmidt_conc$Day_corrected,schmidt_conc$Conc_geomean,xout=c(1:max(schmidt_conc$Day_corrected)))$y
conc_min_schmidt_delta <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Schmidt (delta)",])),function(x){min(schmidt_conc_daily[c(data_efficacy$Day_min[which(data_efficacy$trial=="Schmidt (delta)")[x]]:data_efficacy$Day_max[which(data_efficacy$trial=="Schmidt (delta)")[x]])])})
conc_max_schmidt_delta <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Schmidt (delta)",])),function(x){max(schmidt_conc_daily[c(data_efficacy$Day_min[which(data_efficacy$trial=="Schmidt (delta)")[x]]:data_efficacy$Day_max[which(data_efficacy$trial=="Schmidt (delta)")[x]])])})

geomean_conc_schmidt_omicron <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Schmidt (omicron)",])),
                                       function(x){fun_geomean_conc(schmidt_conc$Day_corrected,schmidt_conc$Conc_geomean,
                                                                    data_efficacy$Day_min[which(data_efficacy$trial=="Schmidt (omicron)")[x]],data_efficacy$Day_max[which(data_efficacy$trial=="Schmidt (omicron)")[x]])})
conc_min_schmidt_omicron <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Schmidt (omicron)",])),function(x){min(schmidt_conc_daily[c(data_efficacy$Day_min[which(data_efficacy$trial=="Schmidt (omicron)")[x]]:data_efficacy$Day_max[which(data_efficacy$trial=="Schmidt (omicron)")[x]])])})
conc_max_schmidt_omicron <- sapply(c(1:nrow(data_efficacy[data_efficacy$trial=="Schmidt (omicron)",])),function(x){max(schmidt_conc_daily[c(data_efficacy$Day_min[which(data_efficacy$trial=="Schmidt (omicron)")[x]]:data_efficacy$Day_max[which(data_efficacy$trial=="Schmidt (omicron)")[x]])])})



# combine raw concentration data ------------------------------------------

data_conc <- data.frame(trial=c(rep("O'Brien & Herman",length(obrien_conc$Day_corrected[obrien_conc$Drug=="casirivimab"])+length(herman_conc$Day[herman_conc$Day>max(obrien_conc$Day_corrected)])),
                                rep("Isa",length(isa_conc$Day[isa_conc$Drug=="casirivimab + imdevimab"])),rep("Levin",length(levin_conc$Day[levin_conc$Drug=="cilgavimab + tixagevimab"])),rep("Schmidt",nrow(schmidt_conc))),
                        day=c(obrien_conc$Day_corrected[obrien_conc$Drug=="casirivimab"],herman_conc$Day[herman_conc$Day>max(obrien_conc$Day_corrected)],
                              isa_conc$Day[isa_conc$Drug=="casirivimab + imdevimab"],levin_conc$Day[levin_conc$Drug=="cilgavimab + tixagevimab"],schmidt_conc$Day_corrected),
                        conc=c(obrien_conc$Conc_geomean[obrien_conc$Drug=="casirivimab"]+obrien_conc$Conc_geomean[obrien_conc$Drug=="imdevimab"],herman_conc$Conc_geomean[herman_conc$Day>max(obrien_conc$Day_corrected)],
                               isa_conc$Conc_geomean[isa_conc$Drug=="casirivimab + imdevimab"],levin_conc$Conc_geomean[levin_conc$Drug=="cilgavimab + tixagevimab"],schmidt_conc$Conc_geomean),
                        drug=c(rep("casirivimab + imdevimab",length(obrien_conc$Day_corrected[obrien_conc$Drug=="casirivimab"])+length(herman_conc$Day[herman_conc$Day>max(obrien_conc$Day_corrected)])+length(isa_conc$Day[isa_conc$Drug=="casirivimab + imdevimab"])),
                               rep("cilgavimab + tixagevimab",length(levin_conc$Day[levin_conc$Drug=="cilgavimab + tixagevimab"])),rep("adintrevimab",nrow(schmidt_conc))),
                        dose.admin=c(rep(1200,length(obrien_conc$Day_corrected[obrien_conc$Drug=="casirivimab"])+length(herman_conc$Day[herman_conc$Day>max(obrien_conc$Day_corrected)])+length(isa_conc$Day[isa_conc$Drug=="casirivimab + imdevimab"])),
                                     rep(300,length(levin_conc$Day[levin_conc$Drug=="cilgavimab + tixagevimab"])),rep(300,nrow(schmidt_conc))),
                        administration=c(rep("SC",length(obrien_conc$Day_corrected[obrien_conc$Drug=="casirivimab"])+length(herman_conc$Day[herman_conc$Day>max(obrien_conc$Day_corrected)])+length(isa_conc$Day[isa_conc$Drug=="casirivimab + imdevimab"])),
                                         rep("IM",length(levin_conc$Day[levin_conc$Drug=="cilgavimab + tixagevimab"])),rep("IM",nrow(schmidt_conc))),
                        variant=c(rep(unique(c(obrien_conc$dominant_variant,herman_conc$dominant_variant)),length(obrien_conc$Day_corrected[obrien_conc$Drug=="casirivimab"])+length(herman_conc$Day[herman_conc$Day>max(obrien_conc$Day_corrected)])),
                                  rep(unique(isa_conc$dominant_variant),length(isa_conc$Day[isa_conc$Drug=="casirivimab + imdevimab"])),rep(unique(levin_conc$dominant_variant),length(levin_conc$Day[levin_conc$Drug=="cilgavimab + tixagevimab"])),
                                  rep(glue::glue(unique(schmidt_conc$dominant_variant_delta)," and ",unique(schmidt_conc$dominant_variant_omicron)),nrow(schmidt_conc))))


# make daily concentration data frame for all studies ---------------------

data_conc_daily <- data.frame(trial=c(rep("O'Brien & Herman",length(herman_conc_daily)),rep("Isa",length(isa_conc_daily)),
                                      rep("Levin",length(levin_conc_daily)),rep("Schmidt",length(schmidt_conc_daily))),
                              day=c(c(0:max(herman_conc$Day_corrected)),c(0:max(isa_conc$Day_corrected)),
                                    c(0:max(levin_conc$Day_corrected)),c(1:max(schmidt_conc$Day_corrected))),
                              conc=c(herman_conc_daily,isa_conc_daily,levin_conc_daily,schmidt_conc_daily),
                              drug=c(rep("casirivimab + imdevimab",length(herman_conc_daily)),rep("casirivimab",length(isa_conc_daily)),
                                     rep("cilgavimab + tixagevimab",length(levin_conc_daily)),rep("adintrevimab",length(schmidt_conc_daily))))


# combine all concentration and efficacy data -----------------------------

data_all <- data.frame(trial=data_efficacy$trial,
                       day=(data_efficacy$Day_max+data_efficacy$Day_min)/2,
                       day.min=data_efficacy$Day_min,
                       day.max=data_efficacy$Day_max,
                       geomean.conc=c(geomean_conc_obrien,geomean_conc_herman,geomean_conc_isa,geomean_conc_levin,geomean_conc_schmidt_delta,geomean_conc_schmidt_omicron),
                       conc.min=c(conc_min_obrien,conc_min_herman,conc_min_isa,conc_min_levin,conc_min_schmidt_delta,conc_min_schmidt_omicron),
                       conc.max=c(conc_max_obrien,conc_max_herman,conc_max_isa,conc_max_levin,conc_max_schmidt_delta,conc_max_schmidt_omicron),
                       n.treatment=data_efficacy$N_treatment,
                       e.treatment=data_efficacy$Events_treatment,
                       n.control=data_efficacy$N_control,
                       e.control=data_efficacy$Events_control,
                       efficacy=data_efficacy$Efficacy,
                       efficacy.lower=data_efficacy$Efficacy_lower,
                       efficacy.upper=data_efficacy$Efficacy_upper,
                       earliest.time=as.numeric(data_efficacy$Day_min%in%c(0,1)),
                       drug=ifelse(data_efficacy$trial=="Levin","cilgavimab + tixagevimab",ifelse(data_efficacy$trial%in%c("Schmidt (delta)","Schmidt (omicron)"),"adintrevimab","casirivimab + imdevimab")),
                       dose.admin=ifelse(data_efficacy$trial%in%c("Levin","Schmidt (delta)","Schmidt (omicron)"),300,1200),
                       administration=ifelse(data_efficacy$trial%in%c("Levin","Schmidt (delta)","Schmidt (omicron)"),"IM","SC")
                       )



# cleanup -----------------------------------------------------------------

rm(obrien_conc,herman_conc,isa_conc,levin_conc,schmidt_conc,
   geomean_conc_obrien,obrien_conc_daily,conc_min_obrien,conc_max_obrien,
   geomean_conc_herman,herman_conc_daily,conc_min_herman,conc_max_herman,
   geomean_conc_isa,isa_conc_daily,conc_min_isa,conc_max_isa,
   geomean_conc_levin,levin_conc_daily,conc_min_levin,conc_max_levin,
   geomean_conc_schmidt_delta,schmidt_conc_daily,conc_min_schmidt_delta,conc_max_schmidt_delta,
   geomean_conc_schmidt_omicron,conc_min_schmidt_omicron,conc_max_schmidt_omicron)

