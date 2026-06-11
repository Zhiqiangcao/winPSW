###simulation of unadjusted, ipw and ow
rm(list=ls())
library(MASS)
library(nnet)
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)
library(mvtnorm)


####data generation
#vector of outcomes as a factor with ordered levels
outcomes_3lvl = factor(c("first", "second", "third"), 
                       levels = c("first", "second", "third"), 
                       ordered = TRUE)

#Order: first < second < third, the larger the better
inv_logit = function(logit) exp(logit)/(1 + exp(logit))

#generate simulation data (quadratic model)
gen_data = function(n,outcomes_3lvl,bi_trt,bi_ctrl,b01,b02,sds,bern_param,pe){
  # covariates
  x1 = rnorm(n, mean = 1, sd = sds[1])
  x2 = rnorm(n, mean = 0.9, sd = sds[2])
  x3 = rnorm(n, mean = 0.8, sd = sds[3])
  
  x4 = rbinom(n,1,bern_param[1])
  x5 = rbinom(n,1,bern_param[2])
  x6 = rbinom(n,1,bern_param[3])
  
  df_cov = data.frame(x1, x2, x3, x4, x5, x6)
  treatment_assignment = rbinom(n, 1, pe) #pe=0.5 for balanced; pe=0.7 for unbalanced
  
  trt_cov = df_cov[treatment_assignment == 1, ]
  ctrl_cov = df_cov[treatment_assignment == 0, ]
  
  trt_cov_quad = trt_cov^2   #quadratic covariates of continuous covariates
  ctrl_cov_quad = ctrl_cov^2
  
  bi_trt_quad = 2*bi_trt     #regression coefficients of quadratic terms
  bi_ctrl_quad = 2*bi_ctrl
  
  combined_trt_quad = cbind(trt_cov, trt_cov_quad)
  combined_ctrl_quad = cbind(ctrl_cov, ctrl_cov_quad)
  
  logodds1_trt = b01 + as.matrix(combined_trt_quad) %*% c(bi_trt, bi_trt_quad) + trt_eff1
  logodds2_trt = b02 + as.matrix(combined_trt_quad) %*% c(bi_trt, bi_trt_quad) + trt_eff1
  
  logodds1_ctrl = b01 + as.matrix(combined_ctrl_quad) %*% c(bi_ctrl, bi_ctrl_quad)
  logodds2_ctrl = b02 + as.matrix(combined_ctrl_quad) %*% c(bi_ctrl, bi_ctrl_quad)
  
  ## Probability trt
  prob_2to3_trt = inv_logit(logodds1_trt)
  prob_3_trt = inv_logit(logodds2_trt)
  prob_1_trt = 1 - prob_2to3_trt
  prob_2_trt = prob_2to3_trt - prob_3_trt
  
  ## Probability ctrl
  prob_2to3_ctrl = inv_logit(logodds1_ctrl)
  prob_3_ctrl = inv_logit(logodds2_ctrl)
  prob_1_ctrl = 1 - prob_2to3_ctrl
  prob_2_ctrl = prob_2to3_ctrl - prob_3_ctrl
  
  #generate random outcomes
  n1 = nrow(trt_cov)
  outcomes_trt <- numeric(n1)
  for (i in 1:n1) {
    outcomes_trt[i] = sample(
      outcomes_3lvl, 
      size = 1,
      prob = c(prob_1_trt[i], prob_2_trt[i], prob_3_trt[i])
    )
  }
  
  n0 = nrow(ctrl_cov)
  outcomes_ctrl = numeric(n0)
  for (i in 1:n0) {
    outcomes_ctrl[i] = sample(
      outcomes_3lvl, 
      size = 1,
      prob = c(prob_1_ctrl[i], prob_2_ctrl[i], prob_3_ctrl[i])
    )
  }
  
  ####################################################
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
  df_comb$x1s = (df_comb$x1)^2  #quadratic covariates
  df_comb$x2s = (df_comb$x2)^2
  df_comb$x3s = (df_comb$x3)^2
  df_comb$x4s = (df_comb$x4)^2  #quadratic covariates
  df_comb$x5s = (df_comb$x5)^2
  df_comb$x6s = (df_comb$x6)^2
  return (df_comb)
}

gen_data_int = function(n,outcomes_3lvl,bi_trt,bi_ctrl,b01,b02,sds,bern_param,pe){
  # covariates
  x1 = rnorm(n, mean = 1, sd = sds[1])
  x2 = rnorm(n, mean = 0.9, sd = sds[2])
  x3 = rnorm(n, mean = 0.8, sd = sds[3])
  
  x4 = rbinom(n,1,bern_param[1])
  x5 = rbinom(n,1,bern_param[2])
  x6 = rbinom(n,1,bern_param[3])
  
  df_cov = data.frame(x1, x2, x3, x4, x5, x6)
  treatment_assignment = rbinom(n, 1, pe) #pe=0.5 for balanced; pe=0.7 for unbalanced
  
  trt_cov = df_cov[treatment_assignment == 1, ]
  ctrl_cov = df_cov[treatment_assignment == 0, ]
  
  #construct interaction terms
  trt_cov_interactions <- data.frame(matrix(nrow = nrow(trt_cov), ncol = 0))
  ctrl_cov_interactions <- data.frame(matrix(nrow = nrow(ctrl_cov), ncol = 0))
  
  for (i in 1:(ncol(trt_cov)-1)) {
    for (j in (i+1):ncol(trt_cov)) {
      trt_cov_interactions[, paste0(names(trt_cov)[i], "_", names(trt_cov)[j])] <- trt_cov[,i] * trt_cov[,j]
      ctrl_cov_interactions[, paste0(names(ctrl_cov)[i], "_x_", names(ctrl_cov)[j])] <- ctrl_cov[,i] * ctrl_cov[,j]
    }
  }
  
  #coefficients of interaction terms
  bi_trt_interactions <- numeric()
  bi_ctrl_interactions <- numeric()
  
  for (i in 1:(length(bi_trt)-1)) {
    for (j in (i+1):length(bi_trt)) {
      if(sign(bi_trt[i]) == sign(bi_trt[j])){
        bi_trt_interactions <- c(bi_trt_interactions, bi_trt[i] + bi_trt[j])
        bi_ctrl_interactions <- c(bi_ctrl_interactions, bi_ctrl[i] + bi_ctrl[j])
      }
      else {
        bi_trt_interactions <- c(bi_trt_interactions, 0.25*bi_trt[i] * bi_trt[j])
        bi_ctrl_interactions <- c(bi_ctrl_interactions, 0.25*bi_ctrl[i] * bi_ctrl[j])
      }
    }
  }
  
  combined_trt_interact = cbind(trt_cov, trt_cov_interactions)
  combined_ctrl_interact = cbind(ctrl_cov, ctrl_cov_interactions)
  
  logodds1_trt = b01 + as.matrix(combined_trt_interact) %*% c(bi_trt, bi_trt_interactions) + trt_eff1
  logodds2_trt = b02 + as.matrix(combined_trt_interact) %*% c(bi_trt, bi_trt_interactions) + trt_eff1
  
  logodds1_ctrl = b01 + as.matrix(combined_ctrl_interact) %*% c(bi_ctrl, bi_ctrl_interactions)
  logodds2_ctrl = b02 + as.matrix(combined_ctrl_interact) %*% c(bi_ctrl, bi_ctrl_interactions)
  
  
  ## Probability trt
  prob_2to3_trt = inv_logit(logodds1_trt)
  prob_3_trt = inv_logit(logodds2_trt)
  prob_1_trt = 1 - prob_2to3_trt
  prob_2_trt = prob_2to3_trt - prob_3_trt
  
  ## Probability ctrl
  prob_2to3_ctrl = inv_logit(logodds1_ctrl)
  prob_3_ctrl = inv_logit(logodds2_ctrl)
  prob_1_ctrl = 1 - prob_2to3_ctrl
  prob_2_ctrl = prob_2to3_ctrl - prob_3_ctrl
  
  #generate random outcomes
  n1 = nrow(trt_cov)
  outcomes_trt <- numeric(n1)
  for (i in 1:n1) {
    outcomes_trt[i] = sample(
      outcomes_3lvl, 
      size = 1,
      prob = c(prob_1_trt[i], prob_2_trt[i], prob_3_trt[i])
    )
  }
  
  n0 = nrow(ctrl_cov)
  outcomes_ctrl = numeric(n0)
  for (i in 1:n0) {
    outcomes_ctrl[i] = sample(
      outcomes_3lvl, 
      size = 1,
      prob = c(prob_1_ctrl[i], prob_2_ctrl[i], prob_3_ctrl[i])
    )
  }
  
  ##################################################################
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
  #interaction terms
  df_comb$x12 = df_comb$x1*df_comb$x2
  df_comb$x13 = df_comb$x1*df_comb$x3
  df_comb$x14 = df_comb$x1*df_comb$x4
  df_comb$x15 = df_comb$x1*df_comb$x5
  df_comb$x16 = df_comb$x1*df_comb$x6
  df_comb$x23 = df_comb$x2*df_comb$x3
  df_comb$x24 = df_comb$x2*df_comb$x4
  df_comb$x25 = df_comb$x2*df_comb$x5
  df_comb$x26 = df_comb$x2*df_comb$x6
  df_comb$x34 = df_comb$x3*df_comb$x4
  df_comb$x35 = df_comb$x3*df_comb$x5
  df_comb$x36 = df_comb$x3*df_comb$x6
  df_comb$x45 = df_comb$x4*df_comb$x5
  df_comb$x46 = df_comb$x4*df_comb$x6
  df_comb$x56 = df_comb$x5*df_comb$x6
  return (df_comb)
}

#cal P value of log(WR) and WD based on wald type z normal test
cal_p_value = function(ws,ws_se){
  #ws=(WR, WD); ws_se are the corresponding standard error of ws
  wr = ws[1]; wd = ws[2]
  se_wr = ws_se[1] 
  se_wd = ws_se[2]
  log_wr = log(wr)
  se_log_wr = (1/abs(wr))*se_wr
  #p-value
  z_wr = log_wr/se_log_wr
  p_wr = 2*pnorm(-abs(z_wr))
  z_wd = wd/se_wd
  p_wd = 2*pnorm(-abs(z_wd))
  res = c(p_wr,p_wd)
  return(res)
}

#simulation parameter settings
setwd("C:/Users/82655/Dropbox/research/covariate adjustment win ratio/simulation_code")
scenarios = read.table("WR_Params.txt", header=TRUE, sep="")
scenarios
sds = c(scenarios$sd_x1, scenarios$sd_x2, scenarios$sd_x3)
bern_param = c(scenarios$param_x4, scenarios$param_x5, scenarios$param_x6)

trt_eff1 = 0 #no treatment effect. scenarios$trt_eff1
bi_trt = c(1,-1,1,-1,1,-1)*scenarios$bi_trt
bi_ctrl = bi_trt  #c(1,-1,1,-1,1,-1)*scenarios$bi_ctrl

b01 = 1
b02 = 0.05
n = 300  #sample size
N = 1000 #simulation times
#keep estimation results
wr = wd = matrix(0,N,7)
#for balanced or unbalanced design
quad = 1
choose = 1
if(quad==1){ #quadratic model
  if(choose == 1){
    pe = 0.5
  }else{
    pe = 0.7
  }
}else{     #interaction model
  if(choose == 1){
    pe = 0.5
  }else{
    pe = 0.7
  }
}

#input source code
source("function_setup_new1.R")
kk = k = 1
#main simulation
while(k <= N){
  set.seed(123456+kk)
  cat("iter=",k,"kk=",kk,"\n")
  if(quad==1){
    mydata = gen_data(n,outcomes_3lvl,bi_trt,bi_ctrl,b01,b02,sds,bern_param,pe)
  }else{
    mydata = gen_data_int(n,outcomes_3lvl,bi_trt,bi_ctrl,b01,b02,sds,bern_param,pe)
  }
  mydata1 = mydata[mydata$treatment==1,]
  mydata0 = mydata[mydata$treatment==0,]
  ltrt = length(unique(mydata1$outcome))
  lctrl = length(unique(mydata0$outcome))
  if(ltrt !=3 | lctrl != 3){
    kk = kk+1
    k = k
  }else{
    #unadjusted estimation
    unadj_res = unadj_new_estimation(data=mydata,outcomevar="outcome",treatment="treatment")
    est1 = unadj_res$EST; se1 = unadj_res$SE
    #p-value for log(WR) and WD
    p_value_res1 = cal_p_value(est1[3:4],se1[3:4])
    if(p_value_res1[1]<0.05) wr[k,1] = 1
    if(p_value_res1[2]<0.05) wd[k,1] = 1
    #ipw estimation
    ipw_res = ipw_estimation(data=mydata,outcomevar="outcome",treatment="treatment",
                             covariate=c("x1","x2","x3","x4","x5","x6"))
    est2 = ipw_res$EST; se2 = ipw_res$SE
    #p-value for log(WR) and WD
    p_value_res2 = cal_p_value(est2[3:4],se2[3:4])
    if(p_value_res2[1]<0.05) wr[k,2] = 1
    if(p_value_res2[2]<0.05) wd[k,2] = 1
    #ow estimation
    ow_res = ow_estimation(data=mydata,outcomevar="outcome",treatment="treatment",
                           covariate=c("x1","x2","x3","x4","x5","x6"))
    est3 = ow_res$EST; se3 = ow_res$SE
    #p-value for log(WR) and WD
    p_value_res3 = cal_p_value(est3[3:4],se3[3:4])
    if(p_value_res3[1]<0.05) wr[k,3] = 1
    if(p_value_res3[2]<0.05) wd[k,3] = 1
    #aipw
    aipw_res = aipw_estimation(data=mydata,outcomevar="outcome",treatment="treatment",
                               covariate=c("x1","x2","x3","x4","x5","x6"),
                               covariate_reg=c("x1","x2","x3","x4","x5","x6",
                                               "x1s","x2s","x3s","x4s","x5s","x6s"))
    est4 = aipw_res$EST; se4 = aipw_res$SE
    #p-value for log(WR) and WD
    p_value_res4 = cal_p_value(est4[3:4],se4[3:4])
    if(p_value_res4[1]<0.05) wr[k,4] = 1
    if(p_value_res4[2]<0.05) wd[k,4] = 1
    #aow
    aow_res = aow_estimation(data=mydata,outcomevar="outcome",treatment="treatment",
                               covariate=c("x1","x2","x3","x4","x5","x6"),
                               covariate_reg=c("x1","x2","x3","x4","x5","x6",
                                               "x1s","x2s","x3s","x4s","x5s","x6s"))
    est5 = aow_res$EST; se5 = aow_res$SE
    #p-value for log(WR) and WD
    p_value_res5 = cal_p_value(est5[3:4],se5[3:4])
    if(p_value_res5[1]<0.05) wr[k,5] = 1
    if(p_value_res5[2]<0.05) wd[k,5] = 1
    
    #aipw_mis
    aipw_mis_res = aipw_estimation(data=mydata,outcomevar="outcome",treatment="treatment",
                               covariate=c("x1","x2","x3","x4","x5","x6"),
                               covariate_reg=c("x1","x2","x3","x4","x5","x6"))
    est6 = aipw_mis_res$EST; se6 = aipw_mis_res$SE
    #p-value for log(WR) and WD
    p_value_res6 = cal_p_value(est6[3:4],se6[3:4])
    if(p_value_res6[1]<0.05) wr[k,6] = 1
    if(p_value_res6[2]<0.05) wd[k,6] = 1
    
    #aow_mis
    aow_mis_res = aow_estimation(data=mydata,outcomevar="outcome",treatment="treatment",
                               covariate=c("x1","x2","x3","x4","x5","x6"),
                               covariate_reg=c("x1","x2","x3","x4","x5","x6"))
    est7 = aow_mis_res$EST; se7 = aow_mis_res$SE
    #p-value for log(WR) and WD
    p_value_res7 = cal_p_value(est7[3:4],se7[3:4])
    if(p_value_res7[1]<0.05) wr[k,7] = 1
    if(p_value_res7[2]<0.05) wd[k,7] = 1
    k = k + 1
    kk = kk +1
  }
}

#summary estimation results
#for unadjusted
em_rr_wr = apply(wr,2,mean)
em_rr_wd = apply(wd,2,mean)
res = rbind(em_rr_wr,em_rr_wd)
res = as.matrix(res)
colnames(res) = c("UNADJ","IPW","OW","AIPW","AOW","AIPW-mis","AOW-mis")
row.names(res) = c("log(WR)","WD")

write.csv(res,file=paste0("emp_type_one_error_quad_",quad,"_choose_",choose, "_n_", n, ".csv"))  

