#' Estimate win statistics by IPW method 
#' 
#' 
#' The IPW (i.e., the inverse probability weight) method is a Hajek-type estimator by leveraging U-statistic theory, see section 3 of paper (Cao et al., 2026). 
#' The point estimate is based on formula (9) of paper with IPW weights, and variance estimate is based on formula (A.11) of Web Appendix with IPW weights. 
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
#' randomized clinical trials. Canadian Journal of Statistics, in press. 2026;0(0):1-27.
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
#' res1 = winstat_ipw_new(data=mydata, outcomevar="outcome", treatment="treatment",
#'        covariate=paste0("x",1:6,seq=""))
#' print(res1)
#' 


winstat_ipw_new = function(data,outcomevar,treatment,covariate){
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
  
  # Define a_w^{IPW}
  compare_aw_IPW = function(i) {
    aw_IPW1_vec = aw_IPW0_vec = numeric(n)
    for (j in 1:n) {
      hij = (treatments[i] * (1-treatments[j]) * 1 / (pi_func[i] * (1-pi_func[j])))
      if (outcomes[i] > outcomes[j]) {
        aw_IPW1_vec[j] = hij
        aw_IPW0_vec[j] = 0
      } else if (outcomes[i] < outcomes[j]) {
        aw_IPW1_vec[j] = 0
        aw_IPW0_vec[j] = hij
      }
    }
    return(list(IPW1 = aw_IPW1_vec, IPW0 = aw_IPW0_vec))
  }
  # Parallelize the outer loop
  results_aw_IPW = apply(matrix(1:n,nrow=1),2,compare_aw_IPW)
  
  # Aggregate results
  g_aw_IPW1 = do.call(rbind, lapply(results_aw_IPW, function(x) x$IPW1))
  g_aw_IPW0 = do.call(rbind, lapply(results_aw_IPW, function(x) x$IPW0))
  #A^{IPW}
  A_IPW1 = sum(as.vector(g_aw_IPW1), na.rm = TRUE)/(n*(n-1))
  A_IPW0 = sum(as.vector(g_aw_IPW0), na.rm = TRUE)/(n*(n-1))
  
  #B^{IPW}
  B_IPW = (1/(n*(n-1)))*sum(treatments/pi_func,na.rm = TRUE)*sum((1-treatments)/(1-pi_func),na.rm = TRUE)
  
  # Calculate tau values, win and loss probabilities
  tau1_IPW = A_IPW1/B_IPW
  tau0_IPW = A_IPW0/B_IPW
  
  #point estimates of WP, LP, WR and WD
  WP_trt_ipw = tau1_IPW
  WP_ctrl_ipw = tau0_IPW
  WR_ipw = tau1_IPW / tau0_IPW
  WD_ipw = tau1_IPW - tau0_IPW
  
  ###########IPW Theoretical Variance###########
  #calculate a1^{IPW} and b1^{IPW}
  aw1_b1_IPW_fun = function(i) {
    #treatment
    hij = ((treatments[i]*(1-treatments)/(pi_func[i] * (1 - pi_func))))
    hji = ((1-treatments[i]) * treatments/(pi_func * (1 - pi_func[i])))
    aw_ipw1 = (1/2)*(hij*as.numeric(outcomes[i]>outcomes)+hji*as.numeric(outcomes[i] < outcomes))
    sum_aw1_IPW_trt = sum(aw_ipw1)
    b_ipw1 = (1/2) * (hij + hji)
    sum_b1_IPW_trt = sum(b_ipw1)
    
    #control
    aw_ipw0 = (1/2)*(hij*as.numeric(outcomes[i]<outcomes)+hji*as.numeric(outcomes[i] > outcomes))
    sum_aw1_IPW_ctrl = sum(aw_ipw0)
    b_ipw0 = b_ipw1 
    sum_b1_IPW_ctrl = sum(b_ipw0)
    return(list(aw1_trt = (1/(n-1))*sum_aw1_IPW_trt, 
                aw1_ctrl = (1/(n-1))*sum_aw1_IPW_ctrl,
                b1_trt = (1/(n-1))*sum_b1_IPW_trt,
                b1_ctrl = (1/(n-1))*sum_b1_IPW_ctrl))
  }
  
  results_aw1_b1 = apply(matrix(1:n,nrow=1),2,aw1_b1_IPW_fun)
  # Aggregate results
  aw1_IPW_trt = as.numeric(do.call(rbind, lapply(results_aw1_b1, function(x) x$aw1_trt)))
  aw1_IPW_ctrl = as.numeric(do.call(rbind, lapply(results_aw1_b1, function(x) x$aw1_ctrl)))
  b1_IPW_trt = as.numeric(do.call(rbind, lapply(results_aw1_b1, function(x) x$b1_trt)))
  b1_IPW_ctrl = as.numeric(do.call(rbind, lapply(results_aw1_b1, function(x) x$b1_ctrl)))
  
  ###U_n\nabla_{\beta}a_w^{\text{IPW}}(\widehat{\beta})
  ###U_n\nabla_{\beta}b^{\text{IPW}}(\widehat{\beta})
  d_aw_b_IPW = function(outcomes,treatments,covariates,pi_func){
    cal_aw_b_derivative = function(i){
      Xi = as.matrix(covariates[i,])
      aw_vec1 = aw_vec0 = b_vec1 = b_vec0 = numeric(p)
      for(j in 1:n){
        if(j != i){
          Xj = as.matrix(covariates[j,])
          e_Xi = pi_func[i]
          e_Xj = pi_func[j]
          hij = (e_Xj*Xj-(1-e_Xi)*Xi)*(treatments[i]*(1-treatments[j]))/(e_Xi*(1-e_Xj))
          term1 = hij * as.numeric(outcomes[i] > outcomes[j])
          term0 = hij * as.numeric(outcomes[i] < outcomes[j])
          aw_vec1 = aw_vec1 + term1
          aw_vec0 = aw_vec0 + term0
          b_vec1 = b_vec1 + hij
          b_vec0 = b_vec0 + hij
        }
      }
      return(list(aw1 = aw_vec1,aw0 = aw_vec0,
                  b1 = b_vec1, b0 = b_vec0))
    }
    aw_b_derivative_temp = apply(matrix(1:n,nrow=1),2,cal_aw_b_derivative)
    total_aw_der1 = Reduce("+", lapply(aw_b_derivative_temp , `[[`, "aw1"))
    total_aw_der0 = Reduce("+", lapply(aw_b_derivative_temp , `[[`, "aw0"))
    total_b_der1 = Reduce("+", lapply(aw_b_derivative_temp , `[[`, "b1"))
    total_b_der0 = Reduce("+", lapply(aw_b_derivative_temp , `[[`, "b0"))
    
    der_aw_IPW_trt = total_aw_der1/(n*(n-1))
    der_aw_IPW_ctrl = total_aw_der0/(n*(n-1))
    der_b_IPW_trt = total_b_der1/(n*(n-1))
    der_b_IPW_ctrl = total_b_der0/(n*(n-1))
    
    return (list(der_aw_IPW_trt = der_aw_IPW_trt,der_aw_IPW_ctrl = der_aw_IPW_ctrl,
                 der_b_IPW_trt = der_b_IPW_trt, der_b_IPW_ctrl = der_b_IPW_ctrl))
  }
  
  der_aw_b_IPW_res = d_aw_b_IPW(outcomes,treatments,covariates,pi_func)
  
  der_aw_IPW_trt = matrix(der_aw_b_IPW_res$der_aw_IPW_trt, ncol=1)
  der_aw_IPW_ctrl = matrix(der_aw_b_IPW_res$der_aw_IPW_ctrl, ncol=1)
  der_b_IPW_trt = matrix(der_aw_b_IPW_res$der_b_IPW_trt, ncol=1)
  der_b_IPW_ctrl = matrix(der_aw_b_IPW_res$der_b_IPW_ctrl, ncol=1)
  
  G_beta_IPW_trt = (der_aw_IPW_trt - WP_trt_ipw*der_b_IPW_trt)/B_IPW
  G_beta_IPW_ctrl = (der_aw_IPW_ctrl - WP_ctrl_ipw*der_b_IPW_ctrl)/B_IPW
  
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
  
  #influence function
  G_beta_IPW_trt = t(G_beta_IPW_trt)
  G_beta_IPW_ctrl = t(G_beta_IPW_ctrl)
  inf_ipw_trt = inf_ipw_ctrl = numeric(n)
  for(i in 1:n){
    l_beta_fun_i = l_beta_fun(i)
    q_IPW_i_trt = WP_trt_ipw + (aw1_IPW_trt[i] - WP_trt_ipw*b1_IPW_trt[i])/B_IPW
    q_IPW_i_ctrl = WP_ctrl_ipw + (aw1_IPW_ctrl[i] - WP_ctrl_ipw*b1_IPW_ctrl[i])/B_IPW
    
    temp_trt = 2*(q_IPW_i_trt - WP_trt_ipw) + G_beta_IPW_trt %*% l_beta_fun_i
    temp_ctrl = 2*(q_IPW_i_ctrl - WP_ctrl_ipw) + G_beta_IPW_ctrl %*% l_beta_fun_i
    
    inf_ipw_trt[i] = temp_trt
    inf_ipw_ctrl[i] = temp_ctrl
  }
  m_inf_trt = mean(inf_ipw_trt)
  m_inf_ctrl = mean(inf_ipw_ctrl)
  #variance and covariance of WP and LP
  IPW_Var_trt = sum((inf_ipw_trt-m_inf_trt)^2)/n
  IPW_Var_ctrl = sum((inf_ipw_ctrl-m_inf_ctrl)^2)/n
  IPW_cov = sum((inf_ipw_trt-m_inf_trt)*(inf_ipw_ctrl-m_inf_ctrl))/n
  
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
