###simulation of aipw-mis in Tables 1-2 of main text and Tables S2-S3 of Web Appendix
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
gen_data_quad = function(n,outcomes_3lvl,bi_trt,bi_ctrl,b01,b02,sds,bern_param,pe){
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

#interaction model
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

#generate data in Web Appendix F  
gen_data = function(n,outcomes_3lvl,mean_set,covariance_set,bi_trt,bi_ctrl,b01,b02,pe){
  # covariates
  x = rmvnorm(n,mean=mean_set,sigma=covariance_set)
  x1 <- x[,1]; x2 <- x[,2]
  x3 <- x[,3]; x4 <- x[,4]
  x5 <- x[,5]; x6 <- x[,6]
  
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

#simulation parameter settings
scenarios = read.table("C:/Users/82655/Dropbox/research/covariate adjustment win ratio/R_package_update/WR_Params.txt", header=TRUE, sep="")
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
n = 200  #sample size
N = 1000 #simulation times

#keep estimation results
est1 = se1 = cova1 = matrix(0,N,4)

#input source code
source("C:/Users/82655/Dropbox/research/covariate adjustment win ratio/R_package_update/winstat_aipw_new.R")

for(scen in 1:6){
  if(scen == 1){ #quadratic model, balanced
    WP_trt_sim_true = 0.2956395
    WP_ctrl_sim_true = 0.2044269 
    WR_sim_true = 1.45289
    WD_sim_true = 0.09121263
    pe = 0.5
    trt_eff1 = 1
    covariate_reg=c("x1","x2","x3","x4","x5","x6")
  }else if(scen == 2){  #quadratic model, unbalanced
    WP_trt_sim_true = 0.295409 
    WP_ctrl_sim_true = 0.2046325
    WR_sim_true = 1.449744 
    WD_sim_true = 0.09077645
    pe = 0.7
    covariate_reg=c("x1","x2","x3","x4","x5","x6")
    trt_eff1 = 1
  }else if(scen == 3){ #interaction model, balanced
    WP_trt_sim_true = 0.3174456 
    WP_ctrl_sim_true = 0.1917862
    WR_sim_true = 1.662504
    WD_sim_true = 0.1256594
    pe = 0.5
    trt_eff1 = 2
    covariate_reg=c("x1","x2","x3","x4","x5","x6")
  }else if(scen == 4){ #interaction model, unbalanced
    WP_trt_sim_true = 0.3167853
    WP_ctrl_sim_true = 0.1920353 
    WR_sim_true = 1.657127
    WD_sim_true = 0.12475
    pe = 0.7
    trt_eff1 = 2
    covariate_reg=c("x1","x2","x3","x4","x5","x6")
  }else if(scen==5){ #correlation model, balanced
    WP_trt_sim_true = 0.3736102
    WP_ctrl_sim_true = 0.2429874
    WR_sim_true = 1.546113
    WD_sim_true = 0.1306228
    pe = 0.5
    covariate_reg=c("x1","x4")
    trt_eff1 = 1
  }else if(scen==6){  #correlation model, unbalanced
    WP_trt_sim_true = 0.3726256 
    WP_ctrl_sim_true = 0.2437491
    WR_sim_true = 1.537628  
    WD_sim_true = 0.1288766
    pe = 0.7
    covariate_reg=c("x1","x4")
    trt_eff1 = 1
  }
  
  true_val = c(WP_trt_sim_true,WP_ctrl_sim_true,WR_sim_true,WD_sim_true)
  
  #main simulation
  kk = k = 1
  #main simulation
  while(k <= N){
    cat("scen",scen,"iter=",k,"\n")
    set.seed(123456+kk)
    if(scen<=2){
      mydata = gen_data_quad(n,outcomes_3lvl,bi_trt,bi_ctrl,b01,b02,sds,bern_param,pe)
    }else if(scen==3|scen==4){
      mydata = gen_data_int(n,outcomes_3lvl,bi_trt,bi_ctrl,b01,b02,sds,bern_param,pe)
    }else{
      mydata = gen_data(n,outcomes_3lvl,mean_set,covariance_set,bi_trt,bi_ctrl,b01,b02,pe)
    }
    
    mydata1 = mydata[mydata$treatment==1,]
    mydata0 = mydata[mydata$treatment==0,]
    ltrt = length(unique(mydata1$outcome))
    lctrl = length(unique(mydata0$outcome))
    if(ltrt !=3 | lctrl != 3){
      kk = kk+1
      k = k
    }else{
      #AIPW estimation
      aipw_res = winstat_aipw_new(data=mydata,outcomevar="outcome",treatment="treatment",
                                  covariate=c("x1","x2","x3","x4","x5","x6"),
                                  covariate_reg=covariate_reg)
      est1[k,] = aipw_res$EST
      se1[k,] = aipw_res$SE
      low1 = aipw_res$conf_low
      high1 = aipw_res$conf_high
      cova1[k,] = as.numeric(low1<=true_val&true_val<=high1)
      
      
      k = k + 1
      kk = kk +1
    }
  }
  wp_est = apply(est1,2,mean)
  wp_sd = apply(est1,2,sd)
  wp_se = apply(se1,2,mean)
  wp_cr = apply(cova1,2,mean)
  res1 = data.frame(true_val,wp_est,wp_sd=wp_sd,wp_se=wp_se,wp_cr=wp_cr)
  row.names(res1) = c("WP","LP","WR","WD")
  
  final_res = res1
  write.csv(final_res,paste0('scen_',scen,'_aipw_mis_','n=',n, ".csv"))  
}

