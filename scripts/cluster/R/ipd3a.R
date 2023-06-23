############## IDA METHOD 3A: Two-Stage Individual Participant Linear Regression ##############

setwd("/projects/p31385/crystallized")

library(methods)
library(brms)        # bayesian models
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
args <- read.table("scripts/cluster/args/bayesian/ipd3_bayesian.txt"
                   , header = T, stringsAsFactors = F)[jobid,]
print(args)

ipd3_study_mod_fun <- function(trait, outcome, type, mod, study, cov){
  ## load the data
  load(sprintf("data/two_stage/%s_%s_%s.RData", trait, outcome, study))
  
  ## compiled Bayesian model to speed up processing and avoid crashing
  if(type == "Bayesian") load("results/3_ipd_meta/bayes_sample_mod.RData")
  
  ## model formula 
  if (cov == "all") cv <- c("age", "gender", "education")
  if (!cov %in% c("all", "none")) cv <- cov
  rhs <- "p_value"
  rhs <- if(cov != "none") c(rhs, cv) else rhs
  if(mod != "none"){rhs <- c(rhs, paste("p_value", mod, sep = "*"))}
  rhs <- paste(rhs, collapse = " + ")
  f <- paste("o_value ~ ", rhs, collapse = "")
  
  ## run the models & save
  m <- if(type == "Frequentist"){do.call("lm", list(formula = f, data = quote(d)))} else {update(m, formula = f, newdata = d
                                                                                                 , refresh = F, iter = 2000, warmup = 1000)}
  summary(m)
  save(m, file = sprintf("results/3_ipd_meta/%s/studyModels/%s_%s_%s_%s_%s.RData"
                         , type, outcome, trait, mod, cov, study))
  
  ## extract model terms and confidence intervals & save
  fx <- tidy(m, conf.int = T) %>%
    select(term, estimate, conf.low, conf.high); fx
  save(fx, file = sprintf("results/3_ipd_meta/%s/studySummary/%s_%s_%s_%s_%s.RData"
                          , type, outcome, trait, mod, cov, study))
  
  ## calculate effect sizes for random effects meta analysis
  es <- ipd3_es_fun(m, type, mod); es
  save(es, file = sprintf("results/3_ipd_meta/%s/studyEffects/%s_%s_%s_%s_%s.RData"
                          , type, outcome, trait, mod, cov, study))
  
  ## calculate simple effects  
  if(mod != "none"){
    pred.rx <- ipd3_study_simpeff_fun(m, mod, type)
    save(pred.rx, file = sprintf("results/3_ipd_meta/%s/studyPredicted/%s_%s_%s_%s_%s.RData"
                                 , type, outcome, trait, mod, cov, study))
  }
  
  ## clean up the local function environment
  rm(list = c("d", "f", "m", "fx", "es", "rhs"))
  gc()
}

ipd3_es_fun <- function(m, type, mod){
  ## extract model features needed for meta-analysis
  ts <- insight::clean_parameters(m)$Cleaned_Parameter
  ts <- if(mod == "none") "p_value" else ts[grepl("p_value.", ts)]
  ts <- str_replace(ts, "[.]", ":")
  ## standardize the model 
  # ms <- standardize(m)
  ## get standardized model coefficients and standard errors
  ## for bayesian models this is the sd of the posterior estimates
  es <- if(type == "Frequentist"){ 
    summary(m)$coef[ts, c("Estimate", "Std. Error")] %>% setNames(c("Estimate", "SEI"))
  } else {
    fixef(m)[ts, c("Estimate", "Est.Error")] %>% setNames(c("Estimate", "SEI"))
  }
  ## format to standardized format
  es <- es %>%
    t() %>% 
    as.data.frame() %>% 
    mutate(ni = if(type == "Frequentist") {nrow(m$model)} else {nrow(m$data)})
  return(es)
}

ipd3_study_simpeff_fun <- function(m, moder, type){
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
  
  mod_frame <- crossing(
    p_value = seq(0,10,.5)
    , modvalue = md_levs
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
ipd3_study_mod_fun(args[,1], args[,2], args[,4], args[,5], args[,3], args[,6])
