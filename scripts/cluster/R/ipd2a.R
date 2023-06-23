############## IDA METHOD 1B: Individual Participants Cluster Corrected Linear Regression ##############

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
args <- read.table("scripts/cluster/args/bayesian/ipd2a_bayesian.txt"
                   , header = T, stringsAsFactors = F)[jobid,]
print(args)

ipd2a_mod_fun <- function(trait, outcome, type, mod, cov){
  ## load the data
  load(sprintf("data/one_stage/%s_%s.RData", trait, outcome))
  
  stdyModers <- c("continent", "country", "scale", "baseAge", "baseYear", "predInt")
  
  ## Applt effects codes  
  d <- contr_fun(d)
  
  ## compiled Bayesian model to speed up processing and avoid crashing
  if(type == "Bayesian") load("results/2a_ipd_dc/bayes_sample_mod.RData")
  
  # get the model formula
  ## model formula 
  if (cov == "all") cv <- c("age", "gender", "education")
  if (!cov %in% c("all", "none")) cv <- cov
  rhs <- c("p_value", "study", "p_value:study")
  rhs <- if(cov != "none") c(rhs, cv) else rhs
  if(!mod %in% c("none", stdyModers)){rhs <- c(rhs, paste("p_value", mod, "study", sep = "*"))}
  if(mod %in% stdyModers){rhs <- c(rhs, mod, paste("p_value", mod, sep = ":"))}
  rhs <- paste(rhs, collapse = " + ")
  f <- paste("o_value ~ ", rhs, collapse = "")
  
  ## run the models & save
  m <- if(type == "Frequentist"){do.call("lm", list(formula = f, data = quote(d)))} else {update(m, formula = f, newdata = d, iter = 2000, warmup = 1000, cores = 4)}
  closeAllConnections()
  save(m, file = sprintf("results/2a_ipd_dc/%s/models/%s_%s_%s_%s.RData"
                         , type, outcome, trait, mod, cov))
  
  ## extract model terms and confidence intervals & save
  fx <- tidy(m, conf.int = T) %>%
    select(term, estimate, conf.low, conf.high)
  rx <- if(!mod %in% stdyModers) std_eff_fun(m, type, mod) else NA
  save(fx, rx, file = sprintf("results/2a_ipd_dc/%s/summary/%s_%s_%s_%s.RData"
                              , type, outcome, trait, mod, cov))
  
  ## clean up the local function environment
  rm(list = c("d", "f", "rhs", "m", "fx", "rx"))
  gc()
}

contr_fun <- function(d){
  stdyModers <- c("continent", "country", "scale", "baseAge", "baseYear", "predInt")
  ## effect code the data  
  d <- d %>%
    mutate(study = str_remove_all(study, "-")
           , study = factor(study))
  std <- rownames(contrasts(d$study))
  cntr <- contr.sum(length(std)); rownames(cntr) <- std; colnames(cntr) <- std[1:(length(std)-1)]
  contrasts(d$study) <- cntr
  return(d)
}

std_eff_fun <- function(m, type, mod){
  if(type == "Bayesian"){
    std <- as.character(unique(m$data$study)) # vector of studies
    nc <- nrow(fixef(m)); nr <- length(std) # number of rows and columns
    trms <- rownames(fixef(m))
  } else{
    std <- as.character(unique(m$model$study)) # vector of studies
    nc <- length(coef(m)); nr <- length(std) # number of rows and columns
    trms <- names(coef(m))
  }
  
  modtrm <- if(mod %in% c("none", stdyModers)) "p_value" else trms[grepl(paste0("p_value:", mod), trms)]
  stdtrms <- trms[grepl("p_value:study", trms)]
  stdtrms <- if (mod %in% c("none", stdyModers)) stdtrms else stdtrms[grepl(mod, stdtrms)]
  
  # create character contrasts 
  cntrm <- paste(stdtrms, "+", modtrm, " = 0")
  cntrm <- c(cntrm, paste0(" - ", stdtrms, collapse = "") %>% 
               paste(modtrm, ., collapse = "") %>% 
               paste(., "= 0", collapse = ""))
  names(cntrm) <- std
  
  # Run the contrasts and rename to match overall model coefficient tables 
  h <- if(type == "Bayesian") {
    hypothesis(m, cntrm)$hypothesis %>% # brms hypothesis function
      select(study = Hypothesis, estimate = Estimate, 
             conf.low = CI.Lower, conf.high = CI.Upper) %>%
      mutate(term = modtrm) %>%
      as_tibble()
  } else {
    (multcomp::glht(m, cntrm) %>% # multcomp hypothesis function
       confint(., calpha = multcomp::univariate_calpha()))$confint %>%
      data.frame() %>% 
      mutate(study = std, 
             term = modtrm) %>%
      rename(estimate = Estimate, conf.low = lwr, conf.high = upr) %>%
      as_tibble()
  }
  return(h)
}

fx_lm_pred_fun <- function(m, newdata){
  tt <- terms(m)
  Terms <- delete.response(tt)
  mf <- model.frame(Terms, newdata, xlev = m$xlevels)
  X <- model.matrix(Terms, mf, contrasts.arg = m$contrasts)
  X[,grepl("study", colnames(X))]  <- 0
  p <- m$rank
  p1 <- seq_len(p)
  piv <- m$qr$pivot[p1] 
  beta <- m$coefficients
  predictor <- drop(X[, piv, drop = FALSE] %*% beta[piv])
  
  w <- m$weights
  res.var <- {
    r <- m$residuals
    rss <- sum(r^2)
    df <- m$df.residual
    rss/df
  }
  
  XRinv <- X[, piv] %*% qr.solve(qr.R(m$qr)[p1, p1])
  ip <- drop(XRinv^2 %*% rep(res.var, p))
  
  tfrac <- qt((.05)/2, df)
  hwid <- tfrac * sqrt(ip)
  predictor <- cbind(predictor, predictor + hwid %o% c(1, -1))
  colnames(predictor) <- c("fit", "lwr", "upr")
  return(predictor)
}

fx_bayes_pred_fun <- function(m, newdata){
  ps <- posterior_samples(m, "^b", as.matrix = T) 
  X <- prepare_predictions(
    m, newdata %>% mutate(o_value = 1)
  )$dpars$mu$fe$X[,1:ncol(ps), drop = FALSE]
  X[,grepl("study", colnames(X))]  <- 0
  predictor <- X %*% t(ps)
  predictor <- posterior_summary(t(predictor))
  return(predictor)
}

ipd2a_simpeff_fun <- function(m, moder, type){
  d <- if(type == "Bayesian") m$data else m$model
  std_lev <- rev(levels(d$study))[1]
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
  md_levs <- if(md_cl == "numeric"){
    if(moder %in% c("age", "baseAge", "baseYear")) {
      c(-10, 0, 10)
    } else if (moder %in% c("predInt", "education")) {
      c(-5, 0, 5) 
    } else {
      with(msd, c(mean[item == moder] - sd[item == moder], mean[item == moder], mean[item == moder] + sd[item == moder]))
    }
  } else { 
    unique(fct_lev[,moder][[1]])
  }
  
  mod_frame <- crossing(
    p_value = seq(0,10,.5)
    , modvalue = md_levs
  ) %>% setNames(c("p_value", moder)) %>%
    mutate(study = std_lev)
  
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
      mod_frame %>% select(-study), 
      fx_bayes_pred_fun(m, newdata = mod_frame) %>% data.frame
    ) %>%
      select(one_of(colnames(m$data)), pred = Estimate, lower = Q2.5, upper = Q97.5)
  } else {
    bind_cols(
      mod_frame %>% select(-study), 
      fx_lm_pred_fun(m, newdata = mod_frame) %>% data.frame 
    ) %>%
      select(one_of(colnames(m$model)), pred = fit, lower = lwr, upper = upr)
  }
  
  rm(list = c("m", "mod_frame", "d", "md_levs"))
  gc()
  return(pred.fx)
}

crossing_fun <- function(df, mod, mod_lev){
  stdyModers <- c("continent", "country", "scale", "baseAge", "baseYear", "predInt")
  pred.rx <- crossing(
    p_value = seq(0,10,.5),
    mod_value = mod_lev
  ) %>%
    setNames(c("p_value", mod))
  return(pred.rx)
}

ipd2a_rx_pred_fun <- function(m, moder, type){
  d <- if(type == "Bayesian") m$data else m$model
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
      summarize_all(~list(unique(.))) %>%
      ungroup() %>%
      data.frame()
  }
  
  d <- d %>% select(-one_of(moder))
  
  md_levs <- if(md_cl == "numeric"){
    if(moder %in% c("age", "baseAge", "baseYear")) {
      c(-10, 0, 10)
    } else if (moder %in% c("predInt", "education")) {
      c(-5, 0, 5) 
    } else {
      with(msd, c(mean[item == moder] - sd[item == moder], mean[item == moder], mean[item == moder] + sd[item == moder]))
    }
  } else { 
    unique(fct_lev[,moder][[1]])
  }
  
  mod_frame <- if(md_cl == "numeric") {
    crossing(
      p_value = seq(0,10,.5),
      modvalue = md_levs,
      study = unique(d$study)
    ) %>% setNames(c("p_value", moder, "study"))#%>% full_join(
    #   msd %>% 
    #     filter(item == moder) %>%
    #     mutate(lower = mean - sd, upper = mean + sd) %>%
    #     select(-sd) %>% 
    #     pivot_longer(cols = c(mean, lower, upper)
    #                  , names_to = "meas"
    #                  , values_to = "modvalue") %>%
    #     pivot_wider(names_from = "item", values_from = "modvalue") %>%
    #     select(study, one_of(moder))
    # )
  } else {
    crossing(
      p_value = seq(0,10,.5)
      , mod_value = md_levs#unique(fct_lev[,moder][[1]])
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
        summarize_all(~levels(.)) %>%
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
      predict(m, newdata = mod_frame, interval = "confidence") %>% data.frame 
    ) %>%
      select(one_of(colnames(m$model)), pred = fit, lower = lwr, upper = upr)
  }
  
  rm(list = c("m", "mod_frame", "d", "md_levs"))
  gc()
  return(pred.rx)
}

stdyModers <- c("continent", "country", "scale", "baseAge", "baseYear", "predInt")

ipd2a_mod_fun(args[,1], args[,2], args[,3], args[,4], args[,5])