# -------------------------------------------------------------------------
#' 
#' This file is used to set up everything that is needed across the project,
#' it loads libraries, creates functions, sets themes and defaults.
#' 
# -------------------------------------------------------------------------


# load packages -----------------------------------------------------------

library(readxl)
library(writexl)
library(ggplot2)
library(ggtext)
library(patchwork)
library(reshape2)
library(Matrix)
library(lme4)
library(mvtnorm)
library(stringr)
library(lmec)
library(nlme)


# load local functions ----------------------------------------------------

source("helper-functions/general-helper-functions.R")


# set defaults ------------------------------------------------------------


alpha_CI <- 0.95 # confidence level for confidence intervals

plot_width <- 7
plot_height <- 4

n_bootstrap <- 1e5 # number of bootstraps

# parameters for model fitting:
n_initial <- 100 # number of random initial parameter values for the model fitting

# number of data sets imputed with concentration data
n_imputations <- 1e2

# upper and lower bounds for all parameters in all other efficacy functions:
all_models <- c("logistic","single hit","powerlaw","threshold","slope threshold","double logistic","logistic with slope 1","logistic with max 1")

LowerB_all_models <- vector(length(all_models),mode='list')
names(LowerB_all_models) <- all_models
LowerB_all_models$logistic <- c(ifelse(max_transform=="exp",log(0.6),ifelse(max_transform=="1-exp",log(1-0.6),ifelse(max_transform=="1-abs",1-0.6,par[1]))),log(0.1),log(10)) 
LowerB_all_models$`single hit` <- c(0)
LowerB_all_models$powerlaw <- c(0,0)
LowerB_all_models$threshold <- c(50,0,0)
LowerB_all_models$`slope threshold` <- c(50,0,0,0)
LowerB_all_models$`double logistic` <- c(50,ifelse(max_transform=="exp",log(0.6),ifelse(max_transform=="1-exp",log(1-0.6),ifelse(max_transform=="1-abs",1-0.6,par[1]))),
                                         log(0.1),log(0.1))
LowerB_all_models$`logistic with slope 1` <- c(ifelse(max_transform=="exp",log(0.6),ifelse(max_transform=="1-exp",log(1-0.6),ifelse(max_transform=="1-abs",1-0.6,par[1]))),log(10)) 
LowerB_all_models$`logistic with max 1` <- c(log(0.1),log(10)) 

UpperB_all_models <- vector(length(all_models),mode='list')
names(UpperB_all_models) <- all_models
UpperB_all_models$logistic <- c(ifelse(max_transform=="exp",log(0.99),ifelse(max_transform=="1-exp",log(1-0.99),ifelse(max_transform=="1-abs",1-0.99,par[1]))),log(100),log(5000)) # maximum, slope, IC50
UpperB_all_models$`single hit` <- c(1e3)
UpperB_all_models$powerlaw<- c(1e3,10)
UpperB_all_models$threshold  <- c(3.5e4,1,1)
UpperB_all_models$`slope threshold` <- c(3.5e4,1,1,1e2)
UpperB_all_models$`double logistic` <- c(3.5e4,ifelse(max_transform=="exp",log(0.99),ifelse(max_transform=="1-exp",log(1-0.99),ifelse(max_transform=="1-abs",1-0.99,par[1]))),
                                         log(100),log(100))
UpperB_all_models$`logistic with slope 1` <- c(ifelse(max_transform=="exp",log(0.99),ifelse(max_transform=="1-exp",log(1-0.99),ifelse(max_transform=="1-abs",1-0.99,par[1]))),log(5000)) 
UpperB_all_models$`logistic with max 1` <- c(log(100),log(5000)) 
