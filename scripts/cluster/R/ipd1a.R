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

jobid = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))+1
print(jobid)
args <- read.table("scripts/cluster/args/bayesian/ipd1a_bayesian.txt", header = F, stringsAsFactors = F)[jobid,]
print(args)

### MODEL FUNCTION
ipd1a_mod_fun <- function(trait, outcome, type, mod, cov){
  ## load the data
  load(sprintf("data/one_stage/%s_%s.RData", trait, outcome))
  
  ## model formula 
  if (cov == "all") cv <- c("age", "gender", "education")
  if (!cov %in% c("all", "none")) cv <- cov
  rhs <- "p_value"
  rhs <- if(cov != "none") c(rhs, cv) else rhs
  if(mod != "none"){rhs <- c(rhs, paste("p_value", mod, sep = "*"))}
  rhs <- paste(rhs, collapse = " + ")
  f <- paste("o_value ~ ", rhs, collapse = "")
  
  ## compiled Bayesian model to speed up processing and avoid crashing
  if(type == "Bayesian") load("results/1a_ipd_reg/bayes_sample_mod.RData")
  
  ## run the models & save
  m <- if(type == "Frequentist"){do.call("lm", list(formula = f, data = quote(d)))} else {update(m, formula = f, newdata = d)}
  save(m, file = sprintf("results/1a_ipd_reg/%s/models/%s_%s_%s_%s.RData"
                         , type, outcome, trait, mod, cov))
  
  ## extract model terms and confidence intervals & save
  fx <- tidy(m, conf.int = T) %>%
    select(term, estimate, conf.low, conf.high)
  save(fx, file = sprintf("results/1a_ipd_reg/%s/summary/%s_%s_%s_%s.RData"
                          , type, outcome, trait, mod, cov))
  
  ## get simple effects for moderator tests
  if(mod != "none"){
    pred.fx <- ipd1a_simpeff_fun(m, mod, type)
    save(pred.fx, file = sprintf("results/1a_ipd_reg/%s/predicted/%s_%s_%s_%s.RData"
                                 , type, outcome, trait, mod, cov))
  }
  
  ## clean up the local function environment
  rm(list = c("d", "f", "rhs", "m", "fx", "cv"))
  gc()
}

##### SIMPLE EFFECTS FUNCTION
ipd1a_simpeff_fun <- function(m, moder, type){
  d <- if(type == "Bayesian") m$data else m$model
  d <- d %>% select(-o_value, -p_value)
  cols <- colnames(d)
  md_cl <- class(d[,moder])
  if(any(sapply(d, class) == "numeric")){
    msd <- d %>%
      select_if(is.numeric) %>%
      pivot_longer(everything()
                   , names_to = "item"
                   , values_to = "value") %>%
      group_by(item) %>%
      summarize_at(vars(value), lst(mean, sd)) %>%
      ungroup()
  }
  if(any(sapply(d, class) == "factor")){
    fct_lev <- d %>% 
      select_if(is.factor) %>%
      summarize_all(~list(unique(.)))
  }
  d <- d %>% select(-one_of(moder))
  md_levs <- if(md_cl == "numeric") with(msd, c(mean[item == moder] - sd[item == moder], mean[item == moder], mean[item == moder] + sd[item == moder])) else unique(fct_lev[,moder][[1]])
  
  mod_frame <- expand.grid(
    p_value = seq(0,10,.5)
    , modvalue = md_levs
    , stringsAsFactors = F
  ) %>% setNames(c("p_value", moder))
  
  if(ncol(d) > 0){
    if(any(sapply(d, class) == "numeric")){
      mod_frame <- tibble(mod_frame, d %>% select_if(is.numeric) %>% summarize_all(mean))
    } 
    if(any(sapply(d, class) == "factor")){
      mod_frame <- tibble(mod_frame, d %>% select_if(is.factor) %>% summarize_all(~levels(.)[1]))
    }
  }
  
  pred.fx <- if(type == "Bayesian"){
    bind_cols(
      mod_frame, 
      fitted(m, newdata = mod_frame) %>% data.frame
    ) %>%
      select(one_of(colnames(m$data)), pred = Estimate, lower = Q2.5, upper = Q97.5)
  } else {
    bind_cols(
      mod_frame, 
      predict(m, newdata = mod_frame, interval = "confidence") %>% data.frame 
    ) %>%
      select(one_of(colnames(m$model)), pred = fit, lower = lwr, upper = upr)
  }
  
  rm(list = c("m", "mod_frame", "d", "md_levs"))
  gc()
  return(pred.fx)
}

ipd1a_mod_fun(args[,1], args[,2], args[,3], args[,4], args[,5])