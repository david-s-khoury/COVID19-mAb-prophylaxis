# -------------------------------------------------------------------------
#' 
#' COVID-19 mAb prophylaxis
#' 
#' Analysis of the efficacy of prophylactic mAb treatment for preventing
#' symptomatic SARS-CoV-2 infection
#' 
# -------------------------------------------------------------------------


# choose analysis ---------------------------------------------------------

# choose analysis to do and how to do it:
run_stanford_db_meta_regression <- FALSE # change to TRUE if you want to re-run regression
CIs_parametric <- TRUE # compute confidence region using parametric bootstrapping (TRUE) or by resampling the data (FALSE)
do_imputation <- FALSE # change to TRUE if you want to fit the model with imputing the mAb concentration and meta-analysis IC50
max_transform <- "1-exp" # choose transform of maximum efficacy in logistic efficacy function: "1-exp","exp","1-abs" (describes how the max is obtained from the fitted parameter)
test_models <- c("logistic","single hit","threshold","double logistic","logistic with slope 1","logistic with max 1") # specify efficacy functions to test/compare: "logistic","single hit","powerlaw","threshold","slope threshold","double logistic","logistic with slope 1","logistic with max 1"
use_efficacy_model <- "logistic with max 1" # efficacy function to use for the full analysis; values: "logistic" or one of the possible values of other_models


# setup -------------------------------------------------------------------

source("setup.R")


# processing --------------------------------------------------------------

# load data and clean data

source("processing/01_ic50-stanford-data.R")
source("processing/02_efficacy-data.R")
source("processing/03_concentration-data.R")



# analysis ----------------------------------------------------------------

source("analysis/01_IC50-meta-analysis.R") 
source("analysis/02_Visualize-concentration-and-efficacy.R")
source("analysis/03_Dose-response-curve.R")
if(do_imputation){
  source("analysis/04_Dose-resp-curve-conc-imputation.R")
}
source("analysis/05_Efficacy-prediction.R")
source("analysis/06_Comparison-vaccines-mAbs.R")

