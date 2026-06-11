###simulation of unadjusted, ipw and ow
rm(list=ls())
library(MASS)
library(nnet)
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

gen_data = function(n,outcomes_3lvl,mean_set,covariance_set,bi_trt,bi_ctrl,b01,b02,pe){
  # covariates
  x = rmvnorm(n,mean=mean_set,sigma=covariance_set)
  x1 <- x[,1]
  x2 <- x[,2]
  x3 <- x[,3]
  x4 <- x[,4]
  x5 <- x[,5]
  x6 <- x[,6]
  
  df_cov <- data.frame(x1, x2, x3, x4, x5, x6)
  treatment_assignment <- rbinom(n, 1, pe) #0.5 is for balanced design; 0.7 is for unbalanced
  
  trt_cov <- df_cov[treatment_assignment == 1, ]
  ctrl_cov <- df_cov[treatment_assignment == 0, ]
  
  logodds1_trt <- b01 + as.matrix(trt_cov) %*% bi_trt + trt_eff1
  logodds2_trt <- b02 + as.matrix(trt_cov) %*% bi_trt + trt_eff1
  
  logodds1_ctrl <- b01 + as.matrix(ctrl_cov) %*% bi_ctrl
  logodds2_ctrl <- b02 + as.matrix(ctrl_cov) %*% bi_ctrl
  
  ## Probability Trt
  prob_2to3_trt <- inv_logit(logodds1_trt)
  prob_3_trt <- inv_logit(logodds2_trt)
  prob_1_trt <- 1 - prob_2to3_trt
  prob_2_trt <- prob_2to3_trt - prob_3_trt
  
  ## Probability Ctrl
  prob_2to3_ctrl <- inv_logit(logodds1_ctrl)
  prob_3_ctrl <- inv_logit(logodds2_ctrl)
  prob_1_ctrl <- 1 - prob_2to3_ctrl
  prob_2_ctrl <- prob_2to3_ctrl - prob_3_ctrl
  
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
covariance_set = matrix(c(1,0.8,0.8,0.2,0.2,0.2,0.8,1,0.8,0.2,0.2,0.2,
                          0.8,0.8,1,0.2,0.2,0.2,0.2,0.2,0.2,1,0.8,0.8,
                          0.2,0.2,0.2,0.8,1,0.8,0.2,0.2,0.2,0.8,0.8,1),byrow = T,6,6)
#mean vector for x
mean_set = seq(0,1,by=0.2)
trt_eff1 = 0 #no treatment effect. scenarios$trt_eff1
bi_trt = c(1,-1,1,-1,1,-1)*scenarios$bi_trt
bi_ctrl = bi_trt  #c(1,-1,1,-1,1,-1)*scenarios$bi_ctrl

b01 = 1
b02 = 0.05
n = 400  #sample size
N = 1000 #simulation times
#keep estimation results
wr = wd = matrix(0,N,7)
#for balanced or unbalanced design
choose = 2
if(choose == 1){
  pe = 0.5
}else{
  pe = 0.7
}


#input source code
source("function_setup_new1.R")


kk = k = 1
#main simulation
while(k <= N){
  set.seed(100000+kk)
  cat("iter=",k,"kk=",kk,"\n")
  mydata = gen_data(n,outcomes_3lvl,mean_set,covariance_set,bi_trt,bi_ctrl,b01,b02,pe)
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
                               covariate_reg=c("x1","x2","x3","x4","x5","x6"))
    est4 = aipw_res$EST; se4 = aipw_res$SE
    #p-value for log(WR) and WD
    p_value_res4 = cal_p_value(est4[3:4],se4[3:4])
    if(p_value_res4[1]<0.05) wr[k,4] = 1
    if(p_value_res4[2]<0.05) wd[k,4] = 1
    #aow
    aow_res = aow_estimation(data=mydata,outcomevar="outcome",treatment="treatment",
                             covariate=c("x1","x2","x3","x4","x5","x6"),
                             covariate_reg=c("x1","x2","x3","x4","x5","x6"))
    est5 = aow_res$EST; se5 = aow_res$SE
    #p-value for log(WR) and WD
    p_value_res5 = cal_p_value(est5[3:4],se5[3:4])
    if(p_value_res5[1]<0.05) wr[k,5] = 1
    if(p_value_res5[2]<0.05) wd[k,5] = 1
    
    #aipw_mis
    aipw_mis_res = aipw_estimation(data=mydata,outcomevar="outcome",treatment="treatment",
                                   covariate=c("x1","x2","x3","x4","x5","x6"),
                                   covariate_reg=c("x1","x4"))
    est6 = aipw_mis_res$EST; se6 = aipw_mis_res$SE
    #p-value for log(WR) and WD
    p_value_res6 = cal_p_value(est6[3:4],se6[3:4])
    if(p_value_res6[1]<0.05) wr[k,6] = 1
    if(p_value_res6[2]<0.05) wd[k,6] = 1
    
    #aow_mis
    aow_mis_res = aow_estimation(data=mydata,outcomevar="outcome",treatment="treatment",
                                 covariate=c("x1","x2","x3","x4","x5","x6"),
                                 covariate_reg=c("x1","x4"))
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

write.csv(res,file=paste0("emp_type_one_error_choose_",choose, "_n_", n, ".csv"))  

