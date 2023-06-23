
############## IDA METHOD 1B: Individual Participants Cluster Corrected Linear Regression ##############

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

jobid = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")) + 1
print(jobid)
args <- read.table("scripts/cluster/args/bayesian/ipd2c_bayesian.txt"
                   , header = T, stringsAsFactors = F)[jobid,]
print(args)

ipd2c_mod_fun <- function(trait, outcome, type, mod, cov){
  ## load the data
  load(sprintf("data/one_stage/%s_%s.RData", trait, outcome))
  
  ## compiled Bayesian model to speed up processing and avoid crashing
  load("results/2c_ipd_melsm/bayes_sample_mod.RData")
  
  stdyModers <- c("continent", "country", "scale", "baseAge", "baseYear", "predInt")
  
  ## formula 
  if (cov == "all") cv <- c("age", "gender", "education")
  if (!cov %in% c("all", "none")) cv <- cov
  rhs <- "p_value"
  rhs <- if(cov != "none") c(rhs, cv) else rhs
  if(mod != "none"){rhs <- c(rhs, paste("p_value", mod, sep = "*"))}
  re <- if(mod == "none" | mod %in% stdyModers) "(p_value | study)" else paste(paste("(p_value", mod, sep = " * "), "| study)")
  rhs <- paste(c(rhs, re), collapse = " + ")
  f1 <- paste("o_value ~ ", rhs, collapse = "")
  f2 <- paste("sigma ~ ", rhs, collapse = "")
  f <- bf(f1, f2)
  
  ## run the models & save
  m <- update(m, formula = f, newdata = d, iter = 2000, warmup = 1000, cores = 12)
  save(m, file = sprintf("results/2c_ipd_melsm/%s/models/%s_%s_%s_%s.RData"
                         , type, outcome, trait, mod, cov))
  
  ## extract model terms and confidence intervals & save
  fx <- fixef(m, probs = c(0.025, 0.975)) %>% data.frame %>% 
    rownames_to_column("term") %>%
    mutate(group = ifelse(grepl("sigma", term), "sigma", "b")
           , term = str_remove_all(term, "sigma_")) %>%
    as_tibble  %>%
    select(term, group, estimate = Estimate, conf.low = Q2.5, conf.high = Q97.5)
  rx <- std_eff_fun(m, type)
  save(fx, rx, file = sprintf("results/2c_ipd_melsm/%s/summary/%s_%s_%s_%s.RData"
                              , type, outcome, trait, mod, cov))
  
  ## extract heterogeneity estimates
  het <- ipd2c_hetero_fun(m, type)
  save(het, file = sprintf("results/2c_ipd_melsm/%s/heterogeneity/%s_%s_%s_%s.RData"
                           , type, outcome, trait, mod, cov))
  
  ## simple effects for moderators
  if(mod != "none"){
    pred.fx <- ipd2c_fx_pred_fun(m, mod, type)
    pred.rx <- if(mod %in% stdyModers) NULL else ipd2c_rx_pred_fun(m, mod, type)
    save(pred.fx, pred.rx, file = sprintf("results/2c_ipd_melsm/%s/predicted/%s_%s_%s_%s.RData"
                                          , type, outcome, trait, mod, cov))
  }
  
  ## clean up the local function environment
  rm(list = c("d", "f", "rhs", "m", "fx", "rx", "het"))
  gc()
}

std_eff_fun <- function(m, type){
  stdyModers <- c("continent", "country", "scale", "baseAge", "baseYear", "predInt")
  if(type == "Frequentist"){
    # coef function gives fixed + random effect estimates but not SE's
    # so we'll take those estimates and get their SE's from the parameters package
    coef(m)$study %>%
      data.frame() %>%
      rownames_to_column("study") %>%
      mutate(term = "estimate") %>%
      full_join(
        parameters::standard_error(m, effects = "random")$study %>%
          data.frame() %>%
          rownames_to_column("study") %>%
          mutate(term = "SE")) %>%
      select(study, term, Intercept = X.Intercept., p_value) %>%
      pivot_longer(c(-study, -term), names_to = "names", values_to = "estimate") %>%
      pivot_wider(names_from = "term", values_from = "estimate") %>%
      mutate(conf.low = estimate - 2*SE, conf.high = estimate + 2*SE)
  } else {
    coef(m, probs = c(0.025, 0.975))[[1]] %>% array_tree(3) %>% 
      tibble(term = names(.), data = .) %>% 
      mutate(data = map(data, ~(.) %>% data.frame %>% 
                          rownames_to_column("study"))) %>% 
      unnest(data) %>% 
      mutate(group = ifelse(grepl("sigma", term), "sigma", "b")
             , term = str_remove_all(term, "sigma_")) %>%
      as_tibble  %>%
      select(term, study, group, estimate = Estimate, conf.low = Q2.5, conf.high = Q97.5)
  }
}

ipd2c_hetero_fun <- function(m, type){
  stdyModers <- c("continent", "country", "scale", "baseAge", "baseYear", "predInt")
  if(type == "Frequentist"){
    tidy(m, effects = "ran_pars", conf.int = T, nsim = 100, conf.method = "bot")
  } else {
    VarCorr(m)$study$sd %>%
      data.frame() %>%
      rownames_to_column("term") %>%
      mutate(group = ifelse(grepl("sigma", term), "sigma", "b")
             , term = str_remove_all(term, "sigma_")) %>%
      mutate_at(vars(Estimate, Q2.5, Q97.5), ~.^2) %>%
      select(term, group, estimate = Estimate, conf.low = Q2.5, conf.high = Q97.5)
  }
  args$conf.method = "boot"; args$nsim <- 10
  do.call(tidy, args) %>%
    select(group, term, estimate, conf.low, conf.high) %>%
    separate(term, c("est", "term"), sep = "__") %>%
    mutate_at(vars(estimate:conf.high), ~ifelse(est == "sd", .^2, .)) %>%
    mutate(est = ifelse(est == "sd", "var", est))
}

predIntlme4 <- function(m, mod_frame, ref){
  b <- bootMer(m, FUN = function(x) 
    lme4:::predict.merMod(x, newdata = mod_frame , re.form = ref)
    , nsim = 100)
  ci <- apply(b$t, 2, quantile, probs = c(.05/2, 1 - .05/2)) %>% t()
  data.frame(pred = predict(m, newdata = mod_frame, re.form = ref), ci) %>% 
    setNames(c("pred", "lower", "upper")) %>% as_tibble()
}

ipd2c_fx_pred_fun <- function(m, moder, type){
  stdyModers <- c("continent", "country", "scale", "baseAge", "baseYear", "predInt")
  d <- if(type == "Bayesian") m$data else m@frame
  d <- d %>% select(-o_value, -p_value, -study)
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
      fitted(m
             , newdata = mod_frame
             , re_formula = NA) %>% data.frame
    ) %>%
      select(one_of(colnames(m$data)), pred = Estimate, lower = Q2.5, upper = Q97.5)
  } else {
    bind_cols(
      mod_frame, 
      predIntlme4(m, mod_frame, NA)
    ) #%>%
    # select(one_of(colnames(m@frame)), pred, lower, upper)
  }
  
  rm(list = c("m", "mod_frame", "d", "md_levs"))
  gc()
  return(pred.fx)
}

ipd2c_rx_pred_fun <- function(m, moder, type){
  stdyModers <- c("continent", "country", "scale", "baseAge", "baseYear", "predInt")
  d <- if(type == "Bayesian") m$data else m@frame
  d <- d %>% select(-o_value, -p_value)
  cols <- colnames(d)
  md_cl <- class(d[,moder])
  if(any(sapply(d, class) == "numeric")){
    msd <- d %>%
      group_by(study) %>%
      select_if(is.numeric) %>%
      pivot_longer(-study
                   , names_to = "item"
                   , values_to = "value") %>%
      group_by(study, item) %>%
      summarize_at(vars(value), lst(mean, sd), na.rm = T) %>%
      ungroup() 
  }
  if(any(sapply(d, class) == "factor")){
    fct_lev <- d %>% 
      group_by(study) %>%
      select_if(is.factor) %>%
      summarize_all(unique) %>%
      ungroup()
  }
  d <- d %>% select(-one_of(moder))
  
  mod_frame <- if(md_cl == "numeric") {
    crossing(
      p_value = seq(0,10,.5),
      study = unique(d$study)
    ) %>% full_join(
      msd %>% 
        filter(item == moder) %>%
        mutate(lower = mean - sd, upper = mean + sd) %>%
        select(-sd) %>% 
        pivot_longer(cols = c(mean, lower, upper)
                     , names_to = "meas"
                     , values_to = "modvalue") %>%
        pivot_wider(names_from = "item", values_from = "modvalue") %>%
        select(study, one_of(moder))
    )
  } else {
    crossing(
      p_value = seq(0,10,.5)
      , mod_value = unique(fct_lev[,moder])
      , study = unique(d$study)
    ) %>%
      setNames(c("p_value", moder, "study"))
  }
  
  if(ncol(d) > 0){
    if(any(sapply(d, class) == "numeric")){
      mod_frame <- d %>% 
        group_by(study) %>% 
        select_if(is.numeric) %>% 
        summarize_all(mean, na.rm = T) %>%
        ungroup() %>%
        full_join(mod_frame)
    } 
    if(any(sapply(d, class) == "factor")){
      mod_frame <- d %>% 
        group_by(study) %>%
        select_if(is.factor) %>% 
        summarize_all(levels) %>%
        ungroup() %>%
        full_join(mod_frame) 
    }
  }
  
  pred.rx <- if(type == "Bayesian"){
    bind_cols(
      mod_frame, 
      fitted(m
             , newdata = mod_frame) %>% data.frame
    ) %>%
      select(one_of(colnames(m$data)), pred = Estimate, lower = Q2.5, upper = Q97.5)
  } else {
    bind_cols(
      mod_frame, 
      predIntlme4(m, mod_frame, NULL)
    ) 
  }
  
  rm(list = c("m", "mod_frame", "d"))
  gc()
  return(pred.rx)
}

ipd2c_mod_fun(args[,1], args[,2], args[,3], args[,4], args[,5])