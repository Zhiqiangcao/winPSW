#' Estimate win statistics by unadjusted method 
#' 
#' The unadjusted estimator (i.e., without covariate adjustment) is realized by three methods introduced in paper (Cao et al., 2026). 
#' i.e., Bebu and Lachin method (BLM), Bias-reduced method (BRM) and the influence function method (IFM). 
#' In addition to estimating point estimates of four win statistics, that is, win probability (WP), loss probability (LP), 
#' win ratio (WR) and win difference (WD), their standard errors as well as 95% confidence intervals will also be returned
#' 
#' @param data A dataset of randomized controlled trials (RCT), which should include ordinal outcome and treatment.
#' @param outcomevar the name of the outcome in RCT data 
#' @param treatment the name of the treatment in RCT data
#' @param method the method used to estimate variance of WP and LP. Three methods are provided, one is the default influence function method (IFM), the second is the Bebu and Lachin method (BLM),
#' and the third is Bias-reduced method (BRM) 
#'  
#' @return A data.frame containing point and interval estimates of WP, LP, WR and WD
#' @export
#' 
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}, Scott Zuo \email{scott.zuo@@northwestern.edu} and Fan Li \email{fan.f.li@@yale.edu} 
#' 
#' @importFrom stats qnorm
#' 
#' @references Cao Z., Zuo S., Ryan M.M., Davis-Plourde K., Heagerty P., Tong G. and Li F. Covariate-adjusted win statistics in 
#' randomized clinical trials. Canadian Journal of Statistics, in press. 2026;0(0):1-27.
#'
#' Bebu I. and Lachin J. Large sample inference for a win ratio analysis of a composite outcome based on prioritized components.
#' Biostatistics 2016; 17(1): 178-187.
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
#' # example 1: using the default IFM method
#' res1 = winstat_unadj_new(data=mydata, outcomevar="outcome", treatment="treatment")
#' print(res1)
#' # example 2: Bebu and Lachin method
#' res2 = winstat_unadj_new(data=mydata, outcomevar="outcome", treatment="treatment",method="BLM")
#' print(res2)
#' # example 3: Bias-reduced method
#' res3 = winstat_unadj_new(data=mydata, outcomevar="outcome", treatment="treatment",method="BRM")
#' print(res3)
#' 


winstat_unadj_new = function(data, outcomevar, treatment, method="IFM"){
  if(method=="BLM"){
    #input treatment and outcome names in data
    treatments = data[,treatment]
    outcomes = data[,outcomevar]
    #classify treatment and control groups
    outcomes_trt = outcomes[treatments==1]
    outcomes_ctrl = outcomes[treatments==0]
    
    df_trt = data.frame(outcomes_trt)
    df_ctrl = data.frame(outcomes_ctrl)
    n1 = nrow(df_trt); n0 = nrow(df_ctrl)
    
    # Initialize matrices
    trtwin = ctrlwin = matrix(0, nrow = n1, ncol = n0)
    
    # Define the pairwise comparison function
    compare_rows_unadj = function(i) {
      trt_row = df_trt[i,1]
      trt_vec = ctrl_vec = numeric(n0)
      for (j in 1:n0) {
        if (trt_row > df_ctrl[j,1]) {
          trt_vec[j] = 1
          ctrl_vec[j] = 0
        } else if (trt_row < df_ctrl[j,1]) {
          trt_vec[j] = 0
          ctrl_vec[j] = 1
        }
      }
      return(list(trt = trt_vec, ctrl = ctrl_vec))
    }
    results = apply(matrix(1:n1,nrow=1),2,compare_rows_unadj)
    
    # Aggregate results
    for (i in 1:length(results)) {
      trtwin[i, ] = results[[i]]$trt
      ctrlwin[i, ] = results[[i]]$ctrl
    }
    # Calculate win proportions
    trt_winpr = sum(as.vector(trtwin), na.rm = TRUE) / (n1*n0)
    ctrl_winpr = sum(as.vector(ctrlwin), na.rm = TRUE) / (n1*n0)
    
    #point estimates of WP, LP, WR and WD
    WP_trt_unadj = trt_winpr
    WP_ctrl_unadj = ctrl_winpr
    WR_unadj = trt_winpr / ctrl_winpr
    WD_unadj = trt_winpr-ctrl_winpr
    
    ###########Theoretical Variance of BL ###########
    colnames(df_trt) = colnames(df_ctrl) = outcomevar
    
    #probability of ijik 
    p_ijik_parallel = function(i) {
      temp_inner1 = temp_inner2 = temp_inner3 = numeric(n0)
      for (j in 1:n0) {
        wij_unadjusted = 1
        temp_inner1[j] = wij_unadjusted * as.numeric(df_trt[i,outcomevar] > df_ctrl[j,outcomevar])
        temp_inner2[j] = wij_unadjusted * as.numeric(df_trt[i,outcomevar] < df_ctrl[j,outcomevar])
      }
      return(list(unadj1 = (sum(temp_inner1))^2, unadj2 = (sum(temp_inner2))^2,
                  unadj3 = sum(temp_inner1)*sum(temp_inner2)))
    }
    
    # Parallelize the outer loop
    results_ijik = apply(matrix(1:n1,nrow=1),2,p_ijik_parallel)
    # Aggregate results
    ijik_unadj1 = do.call(rbind, lapply(results_ijik, function(x) x$unadj1))
    ijik_unadj2 = do.call(rbind, lapply(results_ijik, function(x) x$unadj2))
    term1_trt1 = sum(as.vector(ijik_unadj1))/(n1 * n0 * (n0 - 1))
    term2_trt1 = (1/(n0-1)) * WP_trt_unadj
    p_ijik_trt_results = term1_trt1 - term2_trt1
    term1_ctrl1 = sum(as.vector(ijik_unadj2))/(n1 * n0 * (n0 - 1))
    term2_ctrl1 = (1/(n0-1)) * WP_ctrl_unadj
    p_ijik_ctrl_results = term1_ctrl1 - term2_ctrl1
    #for covariance
    ijik_unadj3 = do.call(rbind, lapply(results_ijik, function(x) x$unadj3))
    p_ijik_cov_results = sum(as.vector(ijik_unadj3))/(n1 * n0 * (n0 - 1))
    
    
    #probability of ijkj 
    p_ijkj_parallel = function(j) {
      temp_inner1 = temp_inner2 = temp_inner3 = numeric(n1)
      for (i in 1:n1) {
        wij_unadjusted = 1
        temp_inner1[i] = wij_unadjusted * as.numeric(df_trt[i,outcomevar] > df_ctrl[j,outcomevar])
        temp_inner2[i] = wij_unadjusted * as.numeric(df_trt[i,outcomevar] < df_ctrl[j,outcomevar])
      }
      return(list(unadj1 = (sum(temp_inner1))^2, unadj2 = (sum(temp_inner2))^2,
                  unadj3 = sum(temp_inner1)*sum(temp_inner2)))
    }
    
    # Parallelize the outer loop
    results_ijkj = apply(matrix(1:n0,nrow=1),2,p_ijkj_parallel)
    # Aggregate results
    ijkj_unadj1 = do.call(rbind, lapply(results_ijkj, function(x) x$unadj1))
    ijkj_unadj2 = do.call(rbind, lapply(results_ijkj, function(x) x$unadj2))
    term1_trt2 = sum(as.vector(ijkj_unadj1))/(n0 * n1 * (n1 - 1))
    term2_trt2 = (1/(n1-1)) * WP_trt_unadj
    p_ijkj_trt_results = term1_trt2 - term2_trt2
    term1_ctrl2 = sum(as.vector(ijkj_unadj2))/(n0 * n1 * (n1 - 1))
    term2_ctrl2 = (1/(n1-1)) * WP_ctrl_unadj
    p_ijkj_ctrl_results = term1_ctrl2 - term2_ctrl2
    #for covariance
    ijkj_unadj3 = do.call(rbind, lapply(results_ijkj, function(x) x$unadj3))
    p_ijkj_cov_results = sum(as.vector(ijkj_unadj3))/(n0 * n1 * (n1 - 1))
    
    #variance and covariance of WP and LP
    var_trt_unadj = (1/n1)*p_ijik_trt_results + (1/n0)*p_ijkj_trt_results - (n1+n0)/(n1*n0)*(WP_trt_unadj)^2
    var_ctrl_unadj = (1/n1)*p_ijik_ctrl_results + (1/n0)*p_ijkj_ctrl_results - (n1+n0)/(n0*n1)*(WP_ctrl_unadj)^2 
    cov_unadj = (1/n1)*p_ijik_cov_results + (1/n0)*p_ijkj_cov_results-(n1+n0)/(n1*n0)*(WP_trt_unadj*WP_ctrl_unadj)
  }else if(method == "BRM"){ #bias reduced method
    treatments = data[,treatment]
    outcomes = data[,outcomevar]
    outcomes_trt = outcomes[treatments==1]
    outcomes_ctrl = outcomes[treatments==0]
    
    df_trt = data.frame(outcomes_trt)
    df_ctrl = data.frame(outcomes_ctrl)
    n1 = nrow(df_trt); n0 = nrow(df_ctrl)
    
    # Initialize matrices
    trtwin = matrix(0, nrow = n1, ncol = n0)
    ctrlwin = matrix(0, nrow = n1, ncol = n0)
    
    # Define the comparison function
    compare_rows_unadj = function(i) {
      trt_row = df_trt[i,1]
      trt_vec = ctrl_vec = numeric(n0)
      
      for (j in 1:n0) {
        if (trt_row > df_ctrl[j,1]) {
          trt_vec[j] = 1
          ctrl_vec[j] = 0
        } else if (trt_row < df_ctrl[j,1]) {
          trt_vec[j] = 0
          ctrl_vec[j] = 1
        }
      }
      return(list(trt = trt_vec, ctrl = ctrl_vec))
    }
    
    results = apply(matrix(1:n1,nrow=1),2,compare_rows_unadj)
    
    # Aggregate results
    for (i in 1:length(results)) {
      trtwin[i, ] = results[[i]]$trt
      ctrlwin[i, ] = results[[i]]$ctrl
    }
    # Calculate win proportions
    trt_winpr = sum(as.vector(trtwin), na.rm = TRUE) / (n1*n0)
    ctrl_winpr = sum(as.vector(ctrlwin), na.rm = TRUE) / (n1*n0)
    
    #win probability (WP_trt), loss probability (WP_ctrl), win ratio and win difference
    WP_trt_unadj = trt_winpr
    WP_ctrl_unadj = ctrl_winpr
    WR_unadj = trt_winpr / ctrl_winpr
    WD_unadj = trt_winpr-ctrl_winpr
    
    
    ###########Theoretical Variance###########
    colnames(df_trt) = colnames(df_ctrl) = outcomevar
    
    #probability of ijik 
    p_ijik_parallel = function(i) {
      temp_inner1 = temp_inner2 = temp_inner3 = numeric(n0)
      for (j in 1:n0) {
        wij_unadjusted = 1
        temp_inner1[j] = wij_unadjusted * as.numeric(df_trt[i,outcomevar] > df_ctrl[j,outcomevar])
        temp_inner2[j] = wij_unadjusted * as.numeric(df_trt[i,outcomevar] < df_ctrl[j,outcomevar])
      }
      return(list(unadj1 = (sum(temp_inner1))^2, unadj2 = (sum(temp_inner2))^2,
                  unadj3 = sum(temp_inner1)*sum(temp_inner2)))
    }
    
    # Parallelize the outer loop
    results_ijik = apply(matrix(1:n1,nrow=1),2,p_ijik_parallel)
    # Aggregate results
    ijik_unadj1 = do.call(rbind, lapply(results_ijik, function(x) x$unadj1))
    ijik_unadj2 = do.call(rbind, lapply(results_ijik, function(x) x$unadj2))
    term1_trt1 = sum(as.vector(ijik_unadj1))/(n1 * n0^2)
    term2_trt1 = 0
    p_ijik_trt_results = term1_trt1 - term2_trt1
    term1_ctrl1 = sum(as.vector(ijik_unadj2))/(n1 * n0^2)
    term2_ctrl1 = 0
    p_ijik_ctrl_results = term1_ctrl1 - term2_ctrl1
    #for covariance
    ijik_unadj3 = do.call(rbind, lapply(results_ijik, function(x) x$unadj3))
    p_ijik_cov_results = sum(as.vector(ijik_unadj3))/(n1 * n0^2)
    
    
    #probability of ijkj 
    p_ijkj_parallel = function(j) {
      temp_inner1 = temp_inner2 = temp_inner3 = numeric(n1)
      for (i in 1:n1) {
        wij_unadjusted = 1
        temp_inner1[i] = wij_unadjusted * as.numeric(df_trt[i,outcomevar] > df_ctrl[j,outcomevar])
        temp_inner2[i] = wij_unadjusted * as.numeric(df_trt[i,outcomevar] < df_ctrl[j,outcomevar])
      }
      return(list(unadj1 = (sum(temp_inner1))^2, unadj2 = (sum(temp_inner2))^2,
                  unadj3 = sum(temp_inner1)*sum(temp_inner2)))
    }
    
    # Parallelize the outer loop
    results_ijkj = apply(matrix(1:n0,nrow=1),2,p_ijkj_parallel)
    # Aggregate results
    ijkj_unadj1 = do.call(rbind, lapply(results_ijkj, function(x) x$unadj1))
    ijkj_unadj2 = do.call(rbind, lapply(results_ijkj, function(x) x$unadj2))
    term1_trt2 = sum(as.vector(ijkj_unadj1))/(n0 * n1^2)
    term2_trt2 = 0
    p_ijkj_trt_results = term1_trt2 - term2_trt2
    term1_ctrl2 = sum(as.vector(ijkj_unadj2))/(n0 * n1^2)
    term2_ctrl2 = 0
    p_ijkj_ctrl_results = term1_ctrl2 - term2_ctrl2
    #for covariance
    ijkj_unadj3 = do.call(rbind, lapply(results_ijkj, function(x) x$unadj3))
    p_ijkj_cov_results = sum(as.vector(ijkj_unadj3))/(n0 * n1^2)
    
    #use no weight 
    var_trt_unadj = (1/n1)*p_ijik_trt_results + (1/n0)*p_ijkj_trt_results - (n1+n0)/(n1*n0)*(WP_trt_unadj)^2
    var_ctrl_unadj = (1/n1)*p_ijik_ctrl_results + (1/n0)*p_ijkj_ctrl_results - (n1+n0)/(n0*n1)*(WP_ctrl_unadj)^2 
    cov_unadj = (1/n1)*p_ijik_cov_results + (1/n0)*p_ijkj_cov_results-(n1+n0)/(n1*n0)*(WP_trt_unadj*WP_ctrl_unadj)
  }else if(method == "IFM"){  #IFM
    treatments = data[,treatment]
    outcomes = data[,outcomevar]
    n = length(outcomes)
    # Define a_w^{UNADJ}
    #pairwise comparison
    compare_aw_unadj = function(i) {
      aw_unadj1_vec = aw_unadj0_vec = numeric(n)
      for (j in 1:n) {
        if (outcomes[i] > outcomes[j]) {
          aw_unadj1_vec[j] = treatments[i] * (1-treatments[j])
          aw_unadj0_vec[j] = 0
        } else if (outcomes[i] < outcomes[j]) {
          aw_unadj1_vec[j] = 0
          aw_unadj0_vec[j] = treatments[i] * (1-treatments[j])
        }
      }
      return(list(unadj1 = aw_unadj1_vec, unadj0 = aw_unadj0_vec))
    }
    # Parallelize the outer loop
    results_aw_unadj = apply(matrix(1:n,nrow=1),2,compare_aw_unadj)
    # Aggregate results
    g_aw_unadj1 = do.call(rbind, lapply(results_aw_unadj, function(x) x$unadj1))
    g_aw_unadj0 = do.call(rbind, lapply(results_aw_unadj, function(x) x$unadj0))
    
    #A^{IPW}
    A_unadj1 = sum(as.vector(g_aw_unadj1), na.rm = TRUE)/(n*(n-1))
    A_unadj0 = sum(as.vector(g_aw_unadj0), na.rm = TRUE)/(n*(n-1))
    
    #B^{UNADJ}
    B_unadj = (1/(n*(n-1))) * sum(treatments)*sum(1-treatments)
    
    # Calculate tau values, win and loss probabilities
    tau1_unadj = A_unadj1/B_unadj
    tau0_unadj = A_unadj0/B_unadj
    
    #point estimates of WP, LP, WR and WD
    WP_trt_unadj = tau1_unadj
    WP_ctrl_unadj = tau0_unadj
    WR_unadj = tau1_unadj/tau0_unadj
    WD_unadj = tau1_unadj-tau0_unadj
    
    ###########unadj Theoretical Variance###########
    #calculate a1^{unadj0} and b1^{unadj}
    aw1_b1_unadj_fun = function(i) {
      #treatment
      hij = treatments[i] * (1 - treatments)
      hji = (1-treatments[i]) * treatments
      aw_unadj1 = (1/2)*(hij*as.numeric(outcomes[i]>outcomes)+hji*as.numeric(outcomes[i] < outcomes))
      sum_aw1_unadj_trt = sum(aw_unadj1)
      b_unadj1 = (1/2) * (hij + hji)
      sum_b1_unadj_trt = sum(b_unadj1)
      
      #control
      aw_unadj0 = (1/2)*(hij*as.numeric(outcomes[i]<outcomes)+hji*as.numeric(outcomes[i] > outcomes))
      sum_aw1_unadj_ctrl = sum(aw_unadj0)
      b_unadj0 = b_unadj1 
      sum_b1_unadj_ctrl = sum(b_unadj0)
      return(list(aw1_trt = (1/(n-1))*sum_aw1_unadj_trt, 
                  aw1_ctrl = (1/(n-1))*sum_aw1_unadj_ctrl,
                  b1_trt = (1/(n-1))*sum_b1_unadj_trt,
                  b1_ctrl = (1/(n-1))*sum_b1_unadj_ctrl))
    }
    
    results_aw1_b1 = apply(matrix(1:n,nrow=1),2,aw1_b1_unadj_fun)
    # Aggregate results
    aw1_unadj_trt = as.numeric(do.call(rbind, lapply(results_aw1_b1, function(x) x$aw1_trt)))
    aw1_unadj_ctrl = as.numeric(do.call(rbind, lapply(results_aw1_b1, function(x) x$aw1_ctrl)))
    b1_unadj_trt = as.numeric(do.call(rbind, lapply(results_aw1_b1, function(x) x$b1_trt)))
    b1_unadj_ctrl = as.numeric(do.call(rbind, lapply(results_aw1_b1, function(x) x$b1_ctrl)))
    
    #influence function
    inf_unadj_trt = inf_unadj_ctrl = numeric(n)
    for(i in 1:n){
      #use (A.8) of web appendix to compute influence function of unadjustment
      q_unadj_i_trt = WP_trt_unadj + (aw1_unadj_trt[i] - WP_trt_unadj*b1_unadj_trt[i])/B_unadj
      q_unadj_i_ctrl = WP_ctrl_unadj + (aw1_unadj_ctrl[i] - WP_ctrl_unadj*b1_unadj_ctrl[i])/B_unadj
      
      #If we use (A.9) of web appendix to compute influence function of unadjustment, then we should use the following
      #q_unadj_i_trt = WP_trt_unadj + (aw1_unadj_trt[i]/B_unadj - WP_trt_unadj)
      #q_unadj_i_ctrl = WP_ctrl_unadj + (aw1_unadj_ctrl[i]/B_unadj - WP_ctrl_unadj)
      #note: results from (A.8) and (A.9) are quite similar.
      #in Table S1, we use (A.9) to compute variance of unadj estimator
      #and in other tables of paper, we use (A.8) to compute variance of unadj estimator
      
      temp_trt = 2*(q_unadj_i_trt - WP_trt_unadj) 
      temp_ctrl = 2*(q_unadj_i_ctrl - WP_ctrl_unadj) 
      
      inf_unadj_trt[i] = temp_trt
      inf_unadj_ctrl[i] = temp_ctrl
    }
    m_inf_trt = mean(inf_unadj_trt)
    m_inf_ctrl = mean(inf_unadj_ctrl)
    #variance and covariance of WP and LP
    unadj_Var_trt = sum((inf_unadj_trt-m_inf_trt)^2)/n
    unadj_Var_ctrl = sum((inf_unadj_ctrl-m_inf_ctrl)^2)/n
    unadj_cov = sum((inf_unadj_trt-m_inf_trt)*(inf_unadj_ctrl-m_inf_ctrl))/n
    
    var_trt_unadj = unadj_Var_trt/n
    var_ctrl_unadj = unadj_Var_ctrl/n
    cov_unadj = unadj_cov/n
  }
  
  #approximation variance of WR
  var_ratio = function(uW,uL,varW,varL,covWL){
    WR_appro_var = (uW^2/uL^2)*(varW/uW^2-2*covWL/(uW*uL)+varL/uL^2)
    return(WR_appro_var)
  }
  
  ###Variance of WR and WD
  var_WR_unadj = var_ratio(WP_trt_unadj,WP_ctrl_unadj,var_trt_unadj,
                           var_ctrl_unadj,cov_unadj)
  var_WD_unadj = var_trt_unadj + var_ctrl_unadj - 2*cov_unadj
  
  #################Unadjusted Estimator Output#########################
  Estimation = c(WP_trt_unadj,WP_ctrl_unadj,WR_unadj,WD_unadj)
  SE = c(sqrt(var_trt_unadj),sqrt(var_ctrl_unadj),sqrt(var_WR_unadj),sqrt(var_WD_unadj))
  conf_low = Estimation-qnorm(0.975)*SE
  conf_high = Estimation+qnorm(0.975)*SE
  res = data.frame(EST=Estimation,SE=SE,conf_low=conf_low,conf_high=conf_high)
  row.names(res) = c("WP","LP","WR","WD")
  return(res)
  
}