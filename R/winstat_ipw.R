#' Estimate win statistics by IPW method 
#' 
#' 
#' The IPW (i.e., the inverse probability weight) method is a Hajek-type estimator by leveraging U-statistic theory, see section 3 of paper (Cao et al., 2025). 
#' The point estimate is based on formula (9) of paper with IPW weights, and variance estimate is based on formula (15) of paper with IPW weights. 
#' In addition to estimating point estimates of four win statistics, that is, win probability (WP), loss probability (LP),  
#' win ratio (WR) and win difference (WD), their standard errors as well as 95% confidence intervals will also be returned
#' 
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
#' res1 = winstat_ipw(data=mydata, outcomevar="outcome", treatment="treatment",
#'        covariate=paste0("x",1:6,seq=""))
#' print(res1)
#' 


winstat_ipw = function(data,outcomevar,treatment,covariate){
  n = dim(data)[1]  #sample size
  ps.formula = as.formula(paste(treatment,"~",paste(covariate,collapse="+"),sep=""))
  # calculate propensity scores
  PropScore = glm(ps.formula,data = data,family=binomial(link="logit"))
  pi_func = fitted(PropScore)
  
  #simplify notations--outcomes, treatments, covariates
  outcomes = data[,outcomevar]
  treatments = data[,treatment]
  covariates = model.matrix(ps.formula,data=data) #including intercept as one covariate
  p = dim(covariates)[2] #number of covariates
  
  # Define the IPW comparison function
  compare_rows_IPW = function(i) {
    g_hfunc_IPW1_vec = g_hfunc_IPW0_vec = numeric(n)
    for (j in 1:n) {
      hij = (treatments[i] * (1-treatments[j]) * 1 / (pi_func[i] * (1-pi_func[j])))
      if (outcomes[i] > outcomes[j]) {
        g_hfunc_IPW1_vec[j] = hij
        g_hfunc_IPW0_vec[j] = 0
      } else if (outcomes[i] < outcomes[j]) {
        g_hfunc_IPW1_vec[j] = 0
        g_hfunc_IPW0_vec[j] = hij
      }
    }
    return(list(IPW1 = g_hfunc_IPW1_vec, IPW0 = g_hfunc_IPW0_vec))
  }
  # Parallelize the outer loop
  results_IPW = apply(matrix(1:n,nrow=1),2,compare_rows_IPW)
  
  # Aggregate results
  g_hfunc_IPW1 = do.call(rbind, lapply(results_IPW, function(x) x$IPW1))
  g_hfunc_IPW0 = do.call(rbind, lapply(results_IPW, function(x) x$IPW0))
  
  # Pre-calculate constants
  adjusted_sum_i_neq_j = sum(treatments/pi_func,na.rm = TRUE)*sum((1-treatments)/(1-pi_func),na.rm = TRUE)
  
  # Calculate tau values
  tau1_IPW = sum(as.vector(g_hfunc_IPW1), na.rm = TRUE)/adjusted_sum_i_neq_j
  tau0_IPW = sum(as.vector(g_hfunc_IPW0), na.rm = TRUE)/adjusted_sum_i_neq_j
  
  #point estimates of WP, LP, WR and WD
  WP_trt_ipw = tau1_IPW
  WP_ctrl_ipw = tau0_IPW
  WR_ipw = tau1_IPW / tau0_IPW
  WD_ipw = tau1_IPW - tau0_IPW
  
  ###########IPW Theoretical Variance###########
  ###g1_IPW_trt###
  g1_constant_denominator_IPW = (1/(n*(n-1))) * adjusted_sum_i_neq_j
  
  g1_IPW_fun = function(i) {
    #treatment
    pairwise_comparisons1 = (1/2) * (
      ((treatments[i] * (1 - treatments) * as.numeric(outcomes[i] > outcomes)) / 
         (pi_func[i] * (1 - pi_func))) + 
        (((1-treatments[i]) * treatments * as.numeric(outcomes[i] < outcomes)) / 
           (pi_func * (1 - pi_func[i]))))
    # Sum and adjust according to the formula
    sum_g1_IPW_trt = sum(pairwise_comparisons1) / g1_constant_denominator_IPW
    #control
    pairwise_comparisons0 = (1/2) * (
      ((treatments[i] * (1 - treatments) * as.numeric(outcomes[i] < outcomes)) / 
         (pi_func[i] * (1 - pi_func))) + 
        (((1-treatments[i]) * treatments * as.numeric(outcomes[i] > outcomes)) / 
           (pi_func * (1 - pi_func[i]))))
    # Sum and adjust according to the formula
    sum_g1_IPW_ctrl = sum(pairwise_comparisons0) / g1_constant_denominator_IPW
    return(list(trt_g1 = (1/(n-1))*sum_g1_IPW_trt, 
                ctrl_g1 = (1/(n-1))*sum_g1_IPW_ctrl))
  }
  
  results_g1 = apply(matrix(1:n,nrow=1),2,g1_IPW_fun)
  # Aggregate results
  g1_IPW_trt = as.numeric(do.call(rbind, lapply(results_g1, function(x) x$trt_g1)))
  g1_IPW_ctrl = as.numeric(do.call(rbind, lapply(results_g1, function(x) x$ctrl_g1)))
  
  ### B^{IPW}, 1*p vector
  dBeta_B = function() {
    X = as.matrix(covariates)
    # Initialize B as a zero matrix with 1 row and p columns
    B = matrix(0, nrow = 1, ncol = p) 
    # Pre-compute the values that are not dependent on i
    ps = pi_func
    Z = treatments
    # Compute sum_j1 and sum_j2 outside the loop
    sum_j1 = sum((1 - Z) / (1 - ps))
    sum_j2 = colSums((1 - Z) * ps * X / (1 - ps))
    for (i in 1:n) {
      ps_i = ps[i]
      Z_i = Z[i]
      X_i = X[i, ]
      term1 = (-Z_i * (1 - ps_i) / ps_i) * X_i * sum_j1
      term2 = (Z_i / ps_i) * sum_j2
      B = B + term1 + term2
    }
    Bn = B / (n * (n - 1))
    return (Bn)
  }
  dBeta_B_val = dBeta_B()
  
  ### C^{IPW}, numeric 
  dBeta_C = function() {
    ps = pi_func
    Z = treatments
    sum_i = sum(Z / ps)
    sum_j = sum((1 - Z) / (1 - ps))
    Cn = (1 / (n * (n - 1))) * sum_i * sum_j
    return(Cn)
  }
  dBeta_C_val = dBeta_C()
  
  ##calculate An
  calculate_An_IPW = function(outcomes,treatments,covariates,pi_func,adjusted_sum_i_neq_j,dBeta_B_val,dBeta_C_val){
    An_first_f = function(i){
      Xi = as.matrix(covariates[i,])
      numerator_vec1 = numerator_vec0 = numeric(p)
      for(j in 1:n){
        Xj = as.matrix(covariates[j,])
        e_Xi = pi_func[i]
        e_Xj = pi_func[j]
        hij = (e_Xj * Xj - (1 - e_Xi) * Xi) * (treatments[i] * (1 - treatments[j]))/ (e_Xi * (1 - e_Xj))
        #term1 = (e_Xj * Xj - (1 - e_Xi) * Xi) * (treatments[i] * (1 - treatments[j])) *
        #  (outcomes[i] > outcomes[j]) / (e_Xi * (1 - e_Xj))
        term1 = hij * (outcomes[i] > outcomes[j])
        #term2 = (e_Xj * Xj - (1 - e_Xi) * Xi) * (treatments[i] * (1 - treatments[j])) *
        #  (outcomes[i] < outcomes[j]) / (e_Xi * (1 - e_Xj))
        term0 = hij * (outcomes[i] < outcomes[j])
        numerator_vec1 = numerator_vec1 + term1
        numerator_vec0 = numerator_vec0 + term0
      }
      return(list(numerator1 = numerator_vec1,numerator0 = numerator_vec0))
    }
    An_first_temp = apply(matrix(1:n,nrow=1),2,An_first_f)
    total_numerator1 = Reduce("+", lapply(An_first_temp, `[[`, "numerator1"))
    total_numerator0 = Reduce("+", lapply(An_first_temp, `[[`, "numerator0"))
    An_first_trt = total_numerator1/adjusted_sum_i_neq_j
    An_first_ctrl = total_numerator0/adjusted_sum_i_neq_j
    An_trt = An_first_trt-t(dBeta_B_val)*(WP_trt_ipw/dBeta_C_val)
    An_ctrl = An_first_ctrl-t(dBeta_B_val)*(WP_ctrl_ipw/dBeta_C_val)
    return (list(An_trt=An_trt,An_ctrl=An_ctrl))
  }
  An_IPW = calculate_An_IPW(outcomes,treatments,covariates,pi_func,adjusted_sum_i_neq_j,dBeta_B_val,dBeta_C_val)
  
  An_trt_IPW = matrix(An_IPW$An_trt,ncol=1)
  An_ctrl_IPW = matrix(An_IPW$An_ctrl,ncol=1)
  
  ### E_beta0beta0
  X = as.matrix(covariates)
  ps = pi_func
  Ebeta0 = crossprod(sqrt(ps*(1-ps)) * X) / n
  #generalized inverse
  inv_Ebeta0 = ginv(Ebeta0)
  ###influence function of \beta, i.e., l_beta 
  l_beta_fun = function(i) {
    X_i = as.numeric(covariates[i,])
    Z_i = treatments[i]
    ps_i = pi_func[i]
    l_beta = inv_Ebeta0 %*% X_i * (Z_i - ps_i)
    return(l_beta)
  }
  
  #variance and covariance of WP and LP
  var_trt_sum = var_ctrl_sum = cov_sum = 0
  An_t_trt = t(An_trt_IPW)
  An_t_ctrl = t(An_ctrl_IPW)
  for(i in 1:n){
    l_beta_fun_i = l_beta_fun(i)
    temp_trt = 2*(g1_IPW_trt[i] - WP_trt_ipw) + An_t_trt %*% l_beta_fun_i
    temp_ctrl = 2*(g1_IPW_ctrl[i] - WP_ctrl_ipw) + An_t_ctrl %*% l_beta_fun_i
    var_trt_sum = var_trt_sum + temp_trt^2
    var_ctrl_sum = var_ctrl_sum + temp_ctrl^2
    cov_sum = cov_sum + temp_trt*temp_ctrl
  }
  IPW_Var_trt = var_trt_sum/n
  IPW_Var_ctrl = var_ctrl_sum/n
  IPW_cov = cov_sum/n
  
  var_trt_ipw = IPW_Var_trt/n
  var_ctrl_ipw = IPW_Var_ctrl/n
  cov_ipw = IPW_cov/n
  
  #approximation variance of WR 
  var_ratio = function(uW,uL,varW,varL,covWL){
    WR_appro_var = (uW^2/uL^2)*(varW/uW^2-2*covWL/(uW*uL)+varL/uL^2)
    return(WR_appro_var)
  }
  
  #Variance of WR and WD
  var_WR_ipw = var_ratio(WP_trt_ipw,WP_ctrl_ipw,var_trt_ipw,var_ctrl_ipw,cov_ipw)
  var_WD_ipw = var_trt_ipw + var_ctrl_ipw - 2*cov_ipw
  #################IPW Estimator Output#########################
  Estimation = c(WP_trt_ipw,WP_ctrl_ipw,WR_ipw,WD_ipw)
  SE = c(sqrt(var_trt_ipw),sqrt(var_ctrl_ipw),sqrt(var_WR_ipw),sqrt(var_WD_ipw))
  conf_low = Estimation-qnorm(0.975)*SE
  conf_high = Estimation+qnorm(0.975)*SE
  res = data.frame(EST=Estimation,SE=SE,conf_low=conf_low,conf_high=conf_high)
  row.names(res) = c("WP","LP","WR","WD")
  return(res)
}
