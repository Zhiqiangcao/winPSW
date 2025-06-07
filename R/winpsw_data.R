#' Generate data with ordinal outcomes and corresponding covariates in a randomized clinical trial 
#' 
#' This function can generate ordinal outcomes with three levels and 6 covariates with quadratic terms of continuous covariates in a randomized clinical trial 
#' 
#' @param n sample size of the data set
#' @param betatrt a vector of regression coefficients of linear predictors for ordinal outcome model in treatment group
#' @param betactrl a vector of regression coefficients of linear predictors for ordinal outcomes model in control group
#' @param b01 intercept for ordinal outcome model in treatment group
#' @param b02 intercept for ordinal outcome model in control group
#' @param sds a vector of standard errors for three continuous covariates
#' @param bernparam a vector of probabilities for generating three binary covariates
#' @param pe the probability for generating treatment variable in a randomized clinical trial
#' 
#' @return A data.frame including ordinal outcome, treatment, covariates x1-x6, and quadratic
#'         terms of 3 continuous covariates
#' 
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}, Scott Zuo \email{scott.zuo@@northwestern.edu} and Fan Li \email{fan.f.li@@yale.edu} 
#' 
#' @importFrom stats rnorm 
#' @importFrom stats rbinom
#' 
#' @examples
#' set.seed(123456)
#' n = 200
#' betatrt = c(1,-1,1,-1,1,-1)
#' betactrl = c(0.5,-0.5,0.5,-0.5,0.5,-0.5)
#' b01 = 1
#' b02 = 0.05
#' sds = c(0.3,0.4,0.5)
#' bernparam = c(0.75,0.50,0.25)
#' pe = 0.5
#' mydata = winpsw_data(n,betatrt,betactrl,b01,b02,sds,bernparam,pe)
#' print(mydata)
#' 

winpsw_data = function(n,betatrt,betactrl,b01,b02,sds,bernparam,pe){
  #six covariates: 3 continuous and 3 binary
  x1 = rnorm(n, mean = 1, sd = sds[1])
  x2 = rnorm(n, mean = 0.9, sd = sds[2])
  x3 = rnorm(n, mean = 0.8, sd = sds[3])
  
  x4 = rbinom(n,1,bernparam[1])
  x5 = rbinom(n,1,bernparam[2])
  x6 = rbinom(n,1,bernparam[3])
  
  df_cov = data.frame(x1, x2, x3, x4, x5, x6)
  treatment_assignment = rbinom(n, 1, pe) #pe=0.5 for balanced design; otherwise, it is unbalanced
  
  trt_cov = df_cov[treatment_assignment == 1, ]
  ctrl_cov = df_cov[treatment_assignment == 0, ]
  
  #construct quadratic terms of 3 continuous covariates
  trt_cov_quad = trt_cov[,1:3]^2   
  ctrl_cov_quad = ctrl_cov[,1:3]^2
  
  betatrt_quad = 2*betatrt[1:3]     #regression coefficients of quadratic terms
  betactrl_quad = 2*betactrl[1:3]
  
  combined_trt_quad = cbind(trt_cov, trt_cov_quad)
  combined_ctrl_quad = cbind(ctrl_cov, ctrl_cov_quad)
  
  logodds1_trt = b01 + as.matrix(combined_trt_quad) %*% c(betatrt, betatrt_quad) + 1
  logodds2_trt = b02 + as.matrix(combined_trt_quad) %*% c(betatrt, betatrt_quad) + 1
  
  logodds1_ctrl = b01 + as.matrix(combined_ctrl_quad) %*% c(betactrl, betactrl_quad)
  logodds2_ctrl = b02 + as.matrix(combined_ctrl_quad) %*% c(betactrl, betactrl_quad)
  
  #order: first < second < third, the larger is more favorable
  inv_logit = function(logit) exp(logit)/(1 + exp(logit))
  
  ## Probability Trt
  prob_2to3_trt = inv_logit(logodds1_trt)
  prob_3_trt = inv_logit(logodds2_trt)
  prob_1_trt = 1 - prob_2to3_trt
  prob_2_trt = prob_2to3_trt - prob_3_trt
  
  ## Probability Ctrl
  prob_2to3_ctrl = inv_logit(logodds1_ctrl)
  prob_3_ctrl = inv_logit(logodds2_ctrl)
  prob_1_ctrl = 1 - prob_2to3_ctrl
  prob_2_ctrl = prob_2to3_ctrl - prob_3_ctrl
  
  #vector of outcomes as a 3 factor with ordered levels
  outcomes3lvl = factor(c("first", "second", "third"), 
                         levels = c("first", "second", "third"), 
                         ordered = TRUE)
  
  #generate random outcomes
  n1 = nrow(trt_cov)
  outcomes_trt = numeric(n1)
  for (i in 1:n1) {
    outcomes_trt[i] = sample(
      outcomes3lvl, 
      size = 1,
      prob = c(prob_1_trt[i], prob_2_trt[i], prob_3_trt[i])
    )
  }
  
  n0 = nrow(ctrl_cov)
  outcomes_ctrl = numeric(n0)
  for (i in 1:n0) {
    outcomes_ctrl[i] = sample(
      outcomes3lvl, 
      size = 1,
      prob = c(prob_1_ctrl[i], prob_2_ctrl[i], prob_3_ctrl[i])
    )
  }
  
  ######################combined data#######################################
  df_trt = data.frame(outcomes_trt)
  df_ctrl = data.frame(outcomes_ctrl)
  colnames(df_trt) =  colnames(df_ctrl) = "outcome"
  
  df_comb = rbind(df_trt, df_ctrl)
  df_comb$treatment = c(rep(1,n1), rep(0,n0))
  
  df_comb$x1 = c(x1[treatment_assignment == 1], x1[treatment_assignment == 0])
  df_comb$x2 = c(x2[treatment_assignment == 1], x2[treatment_assignment == 0])
  df_comb$x3 = c(x3[treatment_assignment == 1], x3[treatment_assignment == 0])
  df_comb$x4 = c(x4[treatment_assignment == 1], x4[treatment_assignment == 0])
  df_comb$x5 = c(x5[treatment_assignment == 1], x5[treatment_assignment == 0])
  df_comb$x6 = c(x6[treatment_assignment == 1], x6[treatment_assignment == 0])
  df_comb$x1s = (df_comb$x1)^2  #quadratic terms
  df_comb$x2s = (df_comb$x2)^2
  df_comb$x3s = (df_comb$x3)^2
  return (df_comb)
}
