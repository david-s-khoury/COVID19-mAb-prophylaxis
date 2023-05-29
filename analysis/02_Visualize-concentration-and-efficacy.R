# -------------------------------------------------------------------------
#' 
#' Visualization of the concentration and efficacy data of prophylactic
#' mAb treatment to prevent symptomatic SARS-CoV-2 infection
#' 
# -------------------------------------------------------------------------


# use the data "data_efficacy"


# Visualization with standard parameters ----------------------------------

# combine all plots
p.conc.eff.combined <- (conc.eff.plot("Levin") | conc.eff.plot("O'Brien & Herman")) / (conc.eff.plot("Isa") | conc.eff.plot("Schmidt"))
p.conc.eff.combined 

# save figure:
cur_time <- format(Sys.time(),"h%H_m%M_s%S")
ggsave(glue::glue("output/Fig1_concentration-and-efficacy_{Sys.Date()}_{cur_time}.pdf"),height=15,width=16,units="cm",plot=p.conc.eff.combined)


# cleanup -----------------------------------------------------------------

rm(p.conc.eff.combined,cur_time)
