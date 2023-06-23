############## IDA METHOD 1A: Individual Participants Cluster Uncorrected Linear Regression ##############

setwd("/projects/p31385/crystallized")

library(methods)
library(brms)        # bayesian models
library(estimatr)    # robust standard error regression
library(lme4)        # Frequentist MLM
library(broom.mixed) # summaries of models
library(bootpredictlme4)    # for calculating prediction intervals
library(rstan)       # bayes underpinnings
library(tidybayes)   # pretty bayes draws and plots
library(plyr)        # data wrangling
library(tidyverse)   # data wrangling

sessionInfo()

jobid = as.integer(Sys.getenv("PBS_ARRAYID"))
print(jobid)
args <- read.table("scripts/cluster/args/bayesian/ipd1a_bayesian.txt", header = F, stringsAsFactors = F)[jobid,]
print(args)
