# -------------------------------------------------------------------------
#' 
#' Process data on the efficacy of prophylactic mAb treatment for the pre-
#' vention of symptomatic SARS-CoV-2 infections.
#' 
# -------------------------------------------------------------------------


# load data ---------------------------------------------------------------

obrien_data <- read_excel("raw-data/Efficacy over time_FINAL_2023-05-16.xlsx",sheet="O'Brien")
herman_data <- read_excel("raw-data/Efficacy over time_FINAL_2023-05-16.xlsx",sheet="Herman")
isa_data <- read_excel("raw-data/Efficacy over time_FINAL_2023-05-16.xlsx",sheet="Isa")
levin_data <- read_excel("raw-data/Efficacy over time_FINAL_2023-05-16.xlsx",sheet="Levin")
schmidt_delta_data <- read_excel("raw-data/Efficacy over time_FINAL_2023-05-16.xlsx",sheet="Schmidt (delta)")
schmidt_omicron_data <- read_excel("raw-data/Efficacy over time_FINAL_2023-05-16.xlsx",sheet="Schmidt (omicron)")


# combine data from different studies -------------------------------------

obrien_data$Week <- NA # change week to month in the O'Brien data in order to be able to combine the data
names(obrien_data) <- names(herman_data)
names(schmidt_delta_data) <- names(herman_data)
names(schmidt_omicron_data) <- names(herman_data)
data_efficacy <- rbind(data.frame(obrien_data,trial=rep("O'Brien",nrow(obrien_data))),
                       data.frame(herman_data,trial=rep("Herman",nrow(herman_data))),
                       data.frame(isa_data,trial=rep("Isa",nrow(isa_data))),
                       data.frame(levin_data,trial=rep("Levin",nrow(levin_data))),
                       data.frame(schmidt_delta_data,trial=rep("Schmidt (delta)",nrow(schmidt_delta_data))),
                       data.frame(schmidt_omicron_data,trial=rep("Schmidt (omicron)",nrow(schmidt_omicron_data))))


# cleanup -----------------------------------------------------------------

rm(obrien_data,herman_data,isa_data,levin_data,schmidt_delta_data,schmidt_omicron_data)

