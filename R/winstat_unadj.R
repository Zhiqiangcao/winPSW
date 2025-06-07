#' Estimate win statistics by unadjusted method 
#' 
#' The unadjusted estimator (i.e., without covariate adjustment) is realized by two methods introduced in paper (Cao et al., 2025). When using the Bebu and Lachin method (BLM), 
#' the point estimate is based on formula (4) of paper, and variance estimate is based on formula (5) of paper. When using the influence function method (IFM), 
#' the point estimate is based on formula (9) of paper with weights =1, and variance estimate is based on formula (18) of paper. 
#' In addition to estimating point estimates of four win statistics, that is, win probability (WP), loss probability (LP), 
#' win ratio (WR) and win difference (WD), their standard errors as well as 95% confidence intervals will also be returned
#' 
#' @param data A dataset of randomized controlled trials (RCT), which should include ordinal outcome and treatment.
#' @param outcomevar the name of the outcome in RCT data 
#' @param treatment the name of the treatment in RCT data
#' @param method the method used to estimate variance of WP and LP. Two methods are provided, one is the default influence function method (IFM), the other is the Bebu and Lachin method (BLM) 
#'  
#' @return A data.frame containing point and interval estimates of WP, LP, WR and WD
#' @export
#' 
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}, Scott Zuo \email{scott.zuo@@northwestern.edu} and Fan Li \email{fan.f.li@@yale.edu} 
#' 
#' @importFrom stats qnorm
#' 
#' @references Cao Z., Zuo S., Ryan M.M., Davis-Plourde K., Heagerty P., Tong G. and Li F. Covariate-adjusted win statistics in 
#' randomized clinical trials. under review. 2025;0(0):1-20.
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
#' res1 = winstat_unadj(data=mydata, outcomevar="outcome", treatment="treatment")
#' print(res1)
#' # example 2: Bebu and Lachin method
#' res2 = winstat_unadj(data=mydata, outcomevar="outcome", treatment="treatment",method="BLM")
#' print(res2)
#' 


winstat_unadj = function(data, outcomevar, treatment, method="IFM"){
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
  }else{
    #treatments and outcomes
    treatments = data[,treatment]
    outcomes = data[,outcomevar]
    n = length(outcomes)
    #pairwise comparison
    compare_rows_unadj = function(i) {
      g_hfunc_unadj1_vec = g_hfunc_unadj2_vec = numeric(n)
      for (j in 1:n) {
        if (outcomes[i] > outcomes[j]) {
          g_hfunc_unadj1_vec[j] = treatments[i] * (1-treatments[j])
          g_hfunc_unadj2_vec[j] = 0
        } else if (outcomes[i] < outcomes[j]) {
          g_hfunc_unadj1_vec[j] = 0
          g_hfunc_unadj2_vec[j] = treatments[i] * (1-treatments[j])
        }
      }
      return(list(unadj1 = g_hfunc_unadj1_vec, unadj2 = g_hfunc_unadj2_vec))
    }
    # Parallelize the outer loop
    results_unadj = apply(matrix(1:n,nrow=1),2,compare_rows_unadj)
    
    # Aggregate results
    g_hfunc_unadj1 = do.call(rbind, lapply(results_unadj, function(x) x$unadj1))
    g_hfunc_unadj2 = do.call(rbind, lapply(results_unadj, function(x) x$unadj2))
    
    adjusted_sum_i_neq_j = sum(treatments)*sum(1-treatments)
    g1_constant_denominator_unadj = (1/(n*(n-1))) * adjusted_sum_i_neq_j
    
    # Calculate tau values
    tau1_unadj = sum(as.vector(g_hfunc_unadj1), na.rm = TRUE)/adjusted_sum_i_neq_j
    tau0_unadj = sum(as.vector(g_hfunc_unadj2), na.rm = TRUE)/adjusted_sum_i_neq_j
    
    #point estimates of WP, LP, WR and WD
    WP_trt_unadj = tau1_unadj
    WP_ctrl_unadj = tau0_unadj
    WR_unadj = tau1_unadj/tau0_unadj
    WD_unadj = tau1_unadj-tau0_unadj
    
    ##combine treatment and control
    g1_unadj_fun = function(i) {
      pairwise_comparisons1 = (1/2) * (
        (treatments[i] * (1 - treatments) * (outcomes[i] > outcomes)) + 
          ((1-treatments[i]) * treatments * (outcomes[i] < outcomes))
      )
      sum_g1_unadj_trt = sum(pairwise_comparisons1) / g1_constant_denominator_unadj
      pairwise_comparisons0 = (1/2) * (
        (treatments[i] * (1 - treatments) * (outcomes[i] < outcomes))+ 
          ((1-treatments[i]) * treatments * (outcomes[i] > outcomes))
      )
      sum_g1_unadj_ctrl = sum(pairwise_comparisons0) / g1_constant_denominator_unadj
      return(list(trt_g1=(1/(n-1))*sum_g1_unadj_trt, ctrl_g1=(1/(n-1))*sum_g1_unadj_ctrl))
    }
    # Parallelize the outer loop
    results_g1 = apply(matrix(1:n,nrow=1),2,g1_unadj_fun)
    
    # Aggregate results
    g1_unadj_trt = as.numeric(do.call(rbind, lapply(results_g1, function(x) x$trt_g1)))
    g1_unadj_ctrl = as.numeric(do.call(rbind, lapply(results_g1, function(x) x$ctrl_g1)))
    
    #variance and covariance of WP and LP
    var_trt_sum = var_ctrl_sum = cov_sum = 0
    for(i in 1:n){
      temp_trt = 2 * (g1_unadj_trt[i] - WP_trt_unadj)
      temp_ctrl = 2 * (g1_unadj_ctrl[i] - WP_ctrl_unadj)
      var_trt_sum = var_trt_sum + temp_trt^2
      var_ctrl_sum = var_ctrl_sum + temp_ctrl^2
      cov_sum = cov_sum + temp_trt*temp_ctrl
    }
    unadj_Var_trt = var_trt_sum/n
    unadj_Var_ctrl = var_ctrl_sum/n
    unadj_cov = cov_sum/n
    
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
