#' Estimate win statistics by OW method 
#' 
#' The OW (i.e., overlap weight) method is a Hajek-type estimator by leveraging U-statistic theory, see section 3 of paper (Cao et al., 2025). 
#' The point estimate is based on formula (9) of paper with OW weights, and variance estimate is based on formula (15) of paper with OW weights. 
#' In addition to estimating point estimates of four win statistics, that is, win probability (WP), loss probability (LP),  
#' win ratio (WR) and win difference (WD), their standard errors as well as 95% confidence intervals will also be returned
#' 
#' @param data A dataset of randomized controlled trials (RCT), which should include ordinal outcome, treatment and covariates
#' @param outcomevar the name of the outcome in RCT data 
#' @param treatment the name of the treatment in RCT data
#' @param covariate the name of vector covariates used in the propensity score model
#' 
#' @return A data.frame containing point and interval estimates of WP, LP, WR and WD
#' 
#' @export
#' @author Zhiqiang Cao \email{zcaoae@@connect.ust.hk}, Scott Zuo \email{scott.zuo@@northwestern.edu} and Fan Li \email{fan.f.li@@yale.edu} 
#' 
#' @importFrom stats glm
#' @importFrom stats as.formula  
#' @importFrom stats fitted
#' @importFrom stats model.matrix
#' @importFrom stats qnorm
#' @importFrom stats binomial 
#' @importFrom MASS ginv
#' 
#' @references Cao Z., Zuo S., Ryan M.M., Davis-Plourde K., Heagerty P., Tong G. and Li F. Covariate-adjusted win statistics in 
#' randomized clinical trials. under review. 2025;0(0):1-20.
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
#' 
#' # example 1 
#' res1 = winstat_ow(data=mydata, outcomevar="outcome", treatment="treatment",
#'        covariate=paste0("x",1:6,seq=""))
#' print(res1)
#' 



winstat_ow = function(data,outcomevar,treatment,covariate){
  n = dim(data)[1]  #sample size
  ps.formula = as.formula(paste(treatment,"~",paste(covariate,collapse="+"),sep=""))
  # calculate propensity scores
  PropScore = glm(ps.formula,data = data,family=binomial(link="logit"))
  pi_func = fitted(PropScore)
  
  #simplify notations--outcomes, treatments, covariates
  outcomes = data[,outcomevar]
  treatments = data[,treatment]
  covariates = model.matrix(ps.formula,data=data) #including intercept as one covariate
  p = dim(covariates)[2]
  
  compare_rows_OW = function(i) {
    g_hfunc_OW1_vec = g_hfunc_OW0_vec = numeric(n)
    for (j in 1:(n)) {
      hij = (treatments[i] * (1-treatments[j]) * (1-pi_func[i]) * pi_func[j])
      if (outcomes[i] > outcomes[j]) {
        g_hfunc_OW1_vec[j] = hij
        g_hfunc_OW0_vec[j] = 0
      } else if (outcomes[i] < outcomes[j]) {
        g_hfunc_OW1_vec[j] = 0
        g_hfunc_OW0_vec[j] = hij
      }
    }
    return(list(OW1 = g_hfunc_OW1_vec, OW0 = g_hfunc_OW0_vec))
  }
  # Parallelize the outer loop
  results_OW = apply(matrix(1:n,nrow=1),2,compare_rows_OW)
  
  # Aggregate results
  g_hfunc_OW1 = do.call(rbind, lapply(results_OW, function(x) x$OW1))
  g_hfunc_OW0 = do.call(rbind, lapply(results_OW, function(x) x$OW0))
  
  adjusted_sum_i_neq_j = sum(treatments*(1-pi_func),na.rm = TRUE)*sum((1-treatments)*pi_func,na.rm = TRUE)
  # Calculate tau values
  tau1_OW = sum(as.vector(g_hfunc_OW1), na.rm = TRUE)/adjusted_sum_i_neq_j
  tau0_OW = sum(as.vector(g_hfunc_OW0), na.rm = TRUE)/adjusted_sum_i_neq_j
  
  #point estimates of WP, LP, WR and WD
  WP_trt_ow = tau1_OW
  WP_ctrl_ow = tau0_OW
  WR_ow = tau1_OW / tau0_OW
  WD_ow = tau1_OW - tau0_OW
  
  ###########OW Theoretical Variance###########
  #g1_OW
  g1_constant_denominator_OW = (1/(n*(n-1))) * adjusted_sum_i_neq_j
  g1_OW_fun = function(i) {
    #for treatment
    pairwise_comparisons_trt = (1/2) * (
      (((1 - pi_func[i]) * pi_func) * treatments[i] * (1 - treatments) * (outcomes[i] > outcomes)) +
        (((1 - pi_func) * pi_func[i]) * treatments * (1 - treatments[i]) * (outcomes > outcomes[i]))
    )
    sum_g1_OW_trt = sum(pairwise_comparisons_trt) / g1_constant_denominator_OW
    #for control
    pairwise_comparisons_ctrl = (1/2) * (
      (((1 - pi_func[i]) * pi_func) * treatments[i] * (1 - treatments) * (outcomes[i] < outcomes)) +
        (((1 - pi_func) * pi_func[i]) * treatments * (1 - treatments[i]) * (outcomes < outcomes[i]))
    )
    sum_g1_OW_ctrl = sum(pairwise_comparisons_ctrl) / g1_constant_denominator_OW
    return(list(trt_g1 = (1/(n-1))*sum_g1_OW_trt, 
                ctrl_g1 = (1/(n-1))*sum_g1_OW_ctrl))
  }
  
  results_g1 = apply(matrix(1:n,nrow=1),2,g1_OW_fun)
  # Aggregate results
  g1_OW_trt = as.numeric(do.call(rbind, lapply(results_g1, function(x) x$trt_g1)))
  g1_OW_ctrl = as.numeric(do.call(rbind, lapply(results_g1, function(x) x$ctrl_g1)))
  
  ### B^{ow}, 1*p vector
  dBeta_B_OW = function() {
    X = as.matrix(covariates)
    p = dim(X)[2]
    B = matrix(0, nrow = 1, ncol = p) # Initialize B as a zero matrix with 1 row and 6 columns
    # Pre-compute the values that are not dependent on i
    ps = pi_func
    Z = treatments
    # Compute sum_j1 and sum_j2 outside the loop
    sum_j1 = sum((1 - Z) * ps)
    sum_j2 = colSums((1 - Z) * ps * (1 - ps) * X)
    for (i in 1:n) {
      ps_i = ps[i]
      Z_i = Z[i]
      X_i = X[i, ]
      term1 = (-Z_i * ps_i * (1-ps_i)) * X_i * sum_j1
      term2 = (Z_i * (1 - ps_i)) * sum_j2
      B = B + term1 + term2
    }
    Bn = B / (n * (n - 1))
    return(Bn)
  }
  dBeta_B_OW_val = dBeta_B_OW()
  
  ### C^{OW}, numeric
  dBeta_C_OW = function() {
    ps =  pi_func
    Z = treatments
    sum_i = sum(Z * (1 - ps))
    sum_j = sum((1 - Z) * ps)
    Cn = (1 / (n * (n - 1))) * sum_i * sum_j
    return(Cn)
  }
  dBeta_C_OW_val = dBeta_C_OW()
  
  calculate_An_OW = function(outcomes,treatments,covariates,pi_func,adjusted_sum_i_neq_j,dBeta_B_OW_val,dBeta_C_OW_val){
    An_first_f = function(i){
      Xi = as.matrix(covariates[i,])
      p = dim(Xi)[2]
      numerator_vec1 = numerator_vec0 = numeric(p)
      for(j in 1:n){
        Xj = as.matrix(covariates[j,])
        e_Xi = pi_func[i]
        e_Xj = pi_func[j]
        hij = ((1-e_Xj) * Xj - e_Xi * Xi) * (treatments[i] * (1 - treatments[j])) * 
          (e_Xj * (1 - e_Xi))
        #term1 = ((1-e_Xj) * Xj - e_Xi * Xi) * (treatments[i] * (1 - treatments[j])) * 
        #  (e_Xj * (1 - e_Xi)) * (outcomes[i] > outcomes[j]) 
        term1 = hij * (outcomes[i] > outcomes[j])
        #term0 = ((1-e_Xj) * Xj - e_Xi * Xi) * (treatments[i] * (1 - treatments[j])) * 
        #  (e_Xj * (1 - e_Xi)) * (outcomes[i] < outcomes[j]) 
        term0 = hij * (outcomes[i] < outcomes[j]) 
        numerator_vec1 = numerator_vec1 + term1
        numerator_vec0 = numerator_vec0 + term0
      }
      return(list(numerator1 = numerator_vec1, numerator0 = numerator_vec0))
    }
    An_first_temp = apply(matrix(1:n,nrow=1),2,An_first_f)
    total_numerator1 = Reduce("+", lapply(An_first_temp, `[[`, "numerator1"))
    total_numerator0 = Reduce("+", lapply(An_first_temp, `[[`, "numerator0"))
    An_first_trt = total_numerator1/adjusted_sum_i_neq_j
    An_first_ctrl = total_numerator0/adjusted_sum_i_neq_j
    An_trt = An_first_trt-t(dBeta_B_OW_val)*(WP_trt_ow/dBeta_C_OW_val)
    An_ctrl = An_first_ctrl-t(dBeta_B_OW_val)*(WP_ctrl_ow/dBeta_C_OW_val)
    return (list(An_trt=An_trt,An_ctrl=An_ctrl))
  }
  An_OW = calculate_An_OW(outcomes,treatments,covariates,pi_func,adjusted_sum_i_neq_j,dBeta_B_OW_val,dBeta_C_OW_val)
  
  An_trt_OW = matrix(An_OW$An_trt,ncol=1)
  An_ctrl_OW = matrix(An_OW$An_ctrl,ncol=1)
  
  ### E_beta0beta0
  X = as.matrix(covariates)
  ps = pi_func
  Ebeta0 = crossprod(sqrt(ps*(1-ps)) * X) / n
  #generalized inverse
  inv_Ebeta0 = ginv(Ebeta0)
  
  ###l_beta 
  l_beta_fun = function(i) {
    X_i = as.numeric(covariates[i,])
    Z_i = treatments[i]
    ps_i = pi_func[i]
    l_beta = inv_Ebeta0 %*% X_i * (Z_i - ps_i)
    return(l_beta)
  }
  
  #variance and covariance of WP and LP
  var_trt_sum = var_ctrl_sum = cov_sum = 0
  An_t_trt = t(An_trt_OW)
  An_t_ctrl = t(An_ctrl_OW)
  for(i in 1:n){
    l_beta_fun_i = l_beta_fun(i)
    temp_trt = 2*(g1_OW_trt[i] - WP_trt_ow) + An_t_trt %*% l_beta_fun_i
    temp_ctrl = 2*(g1_OW_ctrl[i] - WP_ctrl_ow) + An_t_ctrl %*% l_beta_fun_i
    var_trt_sum = var_trt_sum + temp_trt^2
    var_ctrl_sum = var_ctrl_sum + temp_ctrl^2
    cov_sum = cov_sum + temp_trt*temp_ctrl
  }
  OW_Var_trt = var_trt_sum/n
  OW_Var_ctrl = var_ctrl_sum/n
  OW_cov = cov_sum/n
  
  var_trt_ow = OW_Var_trt/n
  var_ctrl_ow = OW_Var_ctrl/n
  cov_ow = OW_cov/n
  
  #approximation variance of win ratio
  var_ratio = function(uW,uL,varW,varL,covWL){
    WR_appro_var = (uW^2/uL^2)*(varW/uW^2-2*covWL/(uW*uL)+varL/uL^2)
    return(WR_appro_var)
  }
  
  ###Variance of WR and WD
  var_WR_ow = var_ratio(WP_trt_ow,WP_ctrl_ow,var_trt_ow,var_ctrl_ow,cov_ow)
  var_WD_ow = var_trt_ow + var_ctrl_ow - 2*cov_ow
  #################OW Estimator Output#########################
  Estimation = c(WP_trt_ow,WP_ctrl_ow,WR_ow,WD_ow)
  SE = c(sqrt(var_trt_ow),sqrt(var_ctrl_ow),sqrt(var_WR_ow),sqrt(var_WD_ow))
  conf_low = Estimation-qnorm(0.975)*SE
  conf_high = Estimation+qnorm(0.975)*SE
  res = data.frame(EST=Estimation,SE=SE,conf_low=conf_low,conf_high=conf_high)
  row.names(res) = c("WP","LP","WR","WD")
  return(res)
}

