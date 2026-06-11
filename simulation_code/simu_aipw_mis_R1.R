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
  treatment_assignment <- rbinom(n, 1, 0.5) #0.5 is for balanced design; 0.7 is for unbalanced
  
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


#simulation parameter settings
scenarios = read.table("C:/Users/82655/Dropbox/research/covariate adjustment win ratio/WR_Params.txt", header=TRUE, sep="")
scenarios
sds = c(scenarios$sd_x1, scenarios$sd_x2, scenarios$sd_x3)
bern_param = c(scenarios$param_x4, scenarios$param_x5, scenarios$param_x6)

trt_eff1 = scenarios$trt_eff1 
bi_trt = c(1,-1,1,-1,1,-1)*scenarios$bi_trt
bi_ctrl = c(1,-1,1,-1,1,-1)*scenarios$bi_ctrl
#variance-covariance matrix for x
covariance_set = matrix(c(1,0.8,0.8,0.2,0.2,0.2,0.8,1,0.8,0.2,0.2,0.2,
                          0.8,0.8,1,0.2,0.2,0.2,0.2,0.2,0.2,1,0.8,0.8,
                          0.2,0.2,0.2,0.8,1,0.8,0.2,0.2,0.2,0.8,0.8,1),byrow = T,6,6)
#mean vector for x
mean_set = seq(0,1,by=0.2)
b01 = 1
b02 = 0.05
n = 400  #sample size
N = 1000 #simulation times
#keep estimation results
est1 = est2 = est3 = matrix(0,N,4)
se1 = se2 = se3 = matrix(0,N,4)
cova1 = cova2 = cova3 = matrix(0,N,4)

#for balanced or unbalanced design
choose = 1
if(choose == 1){
  WP_trt_sim_true = 0.3736102
  WP_ctrl_sim_true = 0.2429874
  WR_sim_true = 1.546113
  WD_sim_true = 0.1306228
  pe = 0.5
}else{
  WP_trt_sim_true = 0.3726256 
  WP_ctrl_sim_true = 0.2437491
  WR_sim_true = 1.537628  
  WD_sim_true = 0.1288766
  pe = 0.7
}

true_val = c(WP_trt_sim_true,WP_ctrl_sim_true,WR_sim_true,WD_sim_true)
#input source code
setwd("C:/Users/82655/Dropbox/research/covariate adjustment win ratio/simulation_code")
source("function_setup_new1.R")

kk = k = 1
#main simulation
while(k <= N){
  set.seed(123456+kk)
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
    aipw_res = aipw_estimation(data=mydata,outcomevar="outcome",treatment="treatment",
                               covariate=c("x1","x2","x3","x4","x5","x6"),
                               covariate_reg=c("x1","x4"))
    est1[k,] = aipw_res$EST
    se1[k,] = aipw_res$SE
    low1 = aipw_res$conf_low
    high1 = aipw_res$conf_high
    cova1[k,] = as.numeric(low1<=true_val&true_val<=high1)
    k = k + 1
    kk = kk +1
  }
  
}

#summary estimation results
wp_est = apply(est1,2,mean)
wp_sd = apply(est1,2,sd)
wp_se = apply(se1,2,mean)
wp_cr = apply(cova1,2,mean)
res1 = data.frame(true_val,wp_est,wp_sd,wp_se,wp_cr)
row.names(res1) = c("WP","LP","WR","WD")
res1
write.csv(res1,file=paste0("aipw_mis_choose_", choose, "_n_", n, ".csv"))  