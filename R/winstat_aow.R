#' Estimate win statistics by AOW method
#' 
#' The AOW (i.e., augmented overlap weight) method is a Hajek-type estimator by leveraging U-statistic theory, see section 4 of paper (Cao et al., 2025). 
#' The point estimate is based on formula (19) of paper with OW weights, and variance estimate is based on formula (20) of paper with OW weights.
#' In addition to estimating point estimates of four win statistics, that is, win probability (WP), loss probability (LP), 
#' win ratio (WR) and win difference (WD), their standard errors as well as 95% confidence intervals will also be returned
#' 
#' @param data A dataset of randomized controlled trials (RCT), which should include ordinal outcome, treatment and covariates
#' @param outcomevar the name of the outcome in RCT data
#' @param treatment the name of the treatment in RCT data
#' @param covariate the name of vector covariates used in the propensity score model
#' @param covariate_reg the name of vector covariates used in the outcome model (i.e., ordinal logistic regression model)
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
#' @importFrom stats predict
#' @importFrom stats binomial 
#' @importFrom MASS ginv
#' @importFrom MASS polr
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
#' res1 = winstat_aow(data=mydata, outcomevar="outcome", treatment="treatment",
#'        covariate=paste0("x",1:6,seq=""),
#'        covariate_reg=c(paste0("x",1:6,seq=""),paste0("x",1:3,"s",seq="")))
#' print(res1)
#' 




winstat_aow = function(data,outcomevar,treatment,covariate,covariate_reg){
  #data is a data.frame including outcome, treatment, covariate in PS model
  #and covariate_reg in outcome model
  n = dim(data)[1]  #sample size
  ps.formula = as.formula(paste(treatment,"~",paste(covariate,collapse="+"),sep=""))
  # calculate propensity scores
  PropScore = glm(ps.formula,data = data,family=binomial(link="logit"))
  pi_func = fitted(PropScore)
  
  #simplify notations--outcomes, treatments, covariates
  outcomes = data[,outcomevar]
  treatments = data[,treatment]
  covariates = model.matrix(ps.formula,data=data) #including intercept as covariate
  covariates_reg = as.matrix(data[,covariate_reg])
  p = dim(covariates)[2]
  po = dim(covariates_reg)[2]
  
  ### mu-function
  lp_start = rep(0, po)
  th_start = c(-1, 1)
  start_values = c(lp_start, th_start)
  model.formual = as.formula(paste("factor(",paste(outcomevar),")","~",paste(covariate_reg,collapse="+"),sep=""))
  dr_trt = polr(model.formual,data = data[treatments==1,],start = start_values)
  cond_prob_trt = predict(dr_trt, newdata = data, type = "probs")
  
  dr_ctrl = polr(model.formual, data = data[treatments == 0,],start = start_values)
  cond_prob_ctrl = predict(dr_ctrl, newdata = data, type = "probs")
  
  
  #number of categories in ordinal outcome
  L = length(unique(outcomes))
  
  compare_rows_AOW = function(i) {
    AOW_num1_vec = numeric(n)
    AOW_denom1_vec = numeric(n)
    AOW_num0_vec = numeric(n)
    AOW_denom0_vec = numeric(n)
    AOW_mu1_vec = numeric(n)
    AOW_mu0_vec = numeric(n)
    AOW_mu_denom1_vec = numeric(n)
    AOW_mu_denom0_vec = numeric(n)
    
    for(j in 1:n){
      if(i != j){
        mu_ij_trt = 0
        for(l in 2:L){
          temp1 = cond_prob_trt[i, l] * (sum(cond_prob_ctrl[j,1:(l-1)]))
          mu_ij_trt = mu_ij_trt+temp1
        }
        
        AOW_num1_vec[j] = (treatments[i] * (1 - treatments[j])) * ((1-pi_func[i]) * pi_func[j]) *
          ((outcomes[i] > outcomes[j]) - mu_ij_trt)
        AOW_denom1_vec[j] = (treatments[i] * (1 - treatments[j])) * ((1-pi_func[i]) * pi_func[j])
        
        #mu_j_i
        mu_ij_ctrl = 0
        for(l in 1:(L-1)){
          temp2 = cond_prob_trt[i, l] * (sum(cond_prob_ctrl[j,(l+1):L]))
          mu_ij_ctrl = mu_ij_ctrl+temp2
        }
        AOW_num0_vec[j] = (treatments[i] * (1 - treatments[j])) * ((1-pi_func[i]) * pi_func[j]) *
          ((outcomes[i] < outcomes[j]) - mu_ij_ctrl)
        #AOW_denom0_vec[j] = (treatments[i] * (1 - treatments[j])) * ((1-pi_func[i]) * pi_func[j])
        AOW_denom0_vec[j] = AOW_denom1_vec[j]
        
        hij = pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])
        AOW_mu1_vec[j] = hij * mu_ij_trt 
        AOW_mu0_vec[j] = hij * mu_ij_ctrl
        
        AOW_mu_denom1_vec[j] = hij
        AOW_mu_denom0_vec[j] = hij
      }
    }
    
    return(list(AOW_num1 = AOW_num1_vec, AOW_denom1 = AOW_denom1_vec, 
                AOW_num0 = AOW_num0_vec, AOW_denom0 = AOW_denom0_vec, 
                AOW_mu1 = AOW_mu1_vec, AOW_mu_denom1 = AOW_mu_denom1_vec,
                AOW_mu0 = AOW_mu0_vec, AOW_mu_denom0 = AOW_mu_denom0_vec))
  }
  
  results_AOW = apply(matrix(1:n,nrow=1),2,compare_rows_AOW)
  
  
  # Combine the results
  tau_AOW_num1 = do.call(rbind, lapply(results_AOW, function(x) x$AOW_num1))
  tau_AOW_denom1 = do.call(rbind, lapply(results_AOW, function(x) x$AOW_denom1))
  tau_AOW_num0 = do.call(rbind, lapply(results_AOW, function(x) x$AOW_num0))
  tau_AOW_denom0 = do.call(rbind, lapply(results_AOW, function(x) x$AOW_denom0))
  tau_AOW_mu1 = do.call(rbind, lapply(results_AOW, function(x) x$AOW_mu1))
  tau_AOW_mu_denom1 = do.call(rbind, lapply(results_AOW, function(x) x$AOW_mu_denom1))
  tau_AOW_mu0 = do.call(rbind, lapply(results_AOW, function(x) x$AOW_mu0))
  tau_AOW_mu_denom0 = do.call(rbind, lapply(results_AOW, function(x) x$AOW_mu_denom0))
  
  tau1_AOW = (sum(as.vector(tau_AOW_num1), na.rm = T) / sum(as.vector(tau_AOW_denom1), na.rm = T)) +
    (sum(as.vector(tau_AOW_mu1), na.rm = T)/sum(as.vector(tau_AOW_mu_denom1), na.rm = T))
  
  tau0_AOW = (sum(as.vector(tau_AOW_num0), na.rm = T) / sum(as.vector(tau_AOW_denom0), na.rm = T)) +
    (sum(as.vector(tau_AOW_mu0), na.rm = T)/sum(as.vector(tau_AOW_mu_denom0), na.rm = T))
  
  #point estimates of WP, LP, WR and WD
  WP_trt_aow = tau1_AOW
  WP_ctrl_aow = tau0_AOW
  WR_aow = tau1_AOW / tau0_AOW
  WD_aow = tau1_AOW - tau0_AOW
  
  ########### AOW Theoretical Variance###########
  mu_func_trt = function(a, b, cond_prob_trt, cond_prob_ctrl) {
    mu_a_b = 0
    for(l in 2:L){
      temp1 = cond_prob_trt[a, l] * (sum(cond_prob_ctrl[b,1:(l-1)]))
      mu_a_b = mu_a_b+temp1
    }
    return(mu_a_b)
  }
  
  mu_func_ctrl = function(a, b, cond_prob_trt, cond_prob_ctrl) {
    mu_a_b = 0
    for(l in 1:(L-1)){
      temp2 = cond_prob_trt[a, l] * (sum(cond_prob_ctrl[b,(l+1):L]))
      mu_a_b = mu_a_b+temp2
    }
    return(mu_a_b)
  }
  
  
  ###g1_AOW_trt###
  adjusted_denom1_term = sum(treatments*(1-pi_func))*sum((1-treatments)*pi_func)
  g1_denominator1_AOW = (1/(n*(n-1))) * adjusted_denom1_term
  adjusted_denom2_term = sum(sapply(1:n, function(j) {sum(pi_func[-j]*(1-pi_func[-j])) * 
      (pi_func[j]) * (1-pi_func[j])}))
  g1_denominator2_AOW = (1/(n*(n-1))) * adjusted_denom2_term
  
  g1_AOW_fun = function(i){
    #for treatment
    pairwise_comparisons_AOW1_term1 = (1/2) * (
      ((treatments[i] * (1 - treatments) * (1 - pi_func[i]) * pi_func) * 
         ((outcomes[i] > outcomes) - sapply(1:n, function(b) mu_func_trt(i, b, cond_prob_trt, cond_prob_ctrl)))) +
        (((1 - treatments[i]) * treatments * pi_func[i] * (1 - pi_func)) * 
           ((outcomes > outcomes[i]) - sapply(1:n, function(a) mu_func_trt(a, i, cond_prob_trt, cond_prob_ctrl))))
    )
    pairwise_comparisons_AOW1_term2 = (1/2) * pi_func[i] * (1 - pi_func[i]) *
      pi_func * (1 - pi_func) * (sapply(1:n, function(b) mu_func_trt(i, b, cond_prob_trt, cond_prob_ctrl)) +
                                   sapply(1:n, function(a) mu_func_trt(a, i, cond_prob_trt, cond_prob_ctrl)))
    # Sum and adjust according to the formula
    sum_g1_AOW_trt <- (sum(pairwise_comparisons_AOW1_term1) / g1_denominator1_AOW) +
      (sum(pairwise_comparisons_AOW1_term2) / g1_denominator2_AOW)
    #for control 
    pairwise_comparisons_AOW0_term1 = (1/2) * (
      ((treatments[i] * (1 - treatments) * (1 - pi_func[i]) * pi_func) * 
         ((outcomes[i] < outcomes) - sapply(1:n, function(b) mu_func_ctrl(i, b, cond_prob_trt, cond_prob_ctrl)))) +
        (((1 - treatments[i]) * treatments * pi_func[i] * (1 - pi_func)) * 
           ((outcomes < outcomes[i]) - sapply(1:n, function(a) mu_func_ctrl(a, i, cond_prob_trt, cond_prob_ctrl))))
    )
    pairwise_comparisons_AOW0_term2 = (1/2) * pi_func[i] * (1 - pi_func[i]) *
      pi_func * (1 - pi_func) * (sapply(1:n, function(b) mu_func_ctrl(i, b, cond_prob_trt, cond_prob_ctrl)) +
                                   sapply(1:n, function(a) mu_func_ctrl(a, i, cond_prob_trt, cond_prob_ctrl)))
    
    sum_g1_AOW_ctrl = (sum(pairwise_comparisons_AOW0_term1) / g1_denominator1_AOW) +
      (sum(pairwise_comparisons_AOW0_term2) / g1_denominator2_AOW)
    
    return(list(trt_g1 = (1/(n-1))*sum_g1_AOW_trt, 
                ctrl_g1 = (1/(n-1))*sum_g1_AOW_ctrl))
  }
  
  results_g1 = apply(matrix(1:n,nrow=1),2,g1_AOW_fun)
  # Aggregate results
  g1_AOW_trt = as.numeric(do.call(rbind, lapply(results_g1, function(x) x$trt_g1)))
  g1_AOW_ctrl = as.numeric(do.call(rbind, lapply(results_g1, function(x) x$ctrl_g1)))
  
  ### B ###
  ### B1_trt & B1_ctrl is a p*1 matrix
  B1_AOW = function(outcomes,treatments,covariates,pi_func, cond_prob_trt, cond_prob_ctrl) {
    B1_value_f = function(j){
      Xj = as.matrix(covariates[j,])
      numerator_vec1 = numerator_vec0 = numeric(p)
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        Xi = as.matrix(covariates[i,])
        e_Xi = pi_func[i]
        e_Xj = pi_func[j]
        mu_ij_trt = mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
        hij_xixj = (treatments[i] * (1 - treatments[j])) * (1 - e_Xi) * e_Xj * 
          ((1-e_Xj) * Xj - e_Xi * Xi)
        #term1 = (treatments[i] * (1 - treatments[j])) * (1 - e_Xi) * e_Xj * 
        #  ((1-e_Xj) * Xj - e_Xi * Xi) * ((outcomes[i] > outcomes[j]) - mu_ij_trt)
        term1 = hij_xixj * ((outcomes[i] > outcomes[j]) - mu_ij_trt)
        numerator_vec1 = numerator_vec1 + term1
        mu_ij_ctrl = mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
        #term2 = (treatments[i] * (1 - treatments[j])) * (1 - e_Xi) * e_Xj * 
        #  ((1-e_Xj) * Xj - e_Xi * Xi) * ((outcomes[i] < outcomes[j]) - mu_ij_ctrl)
        term0 = hij_xixj * ((outcomes[i] < outcomes[j]) - mu_ij_ctrl)
        numerator_vec0 = numerator_vec0 + term0
      }
      return(list(numerator1 = numerator_vec1, numerator0 = numerator_vec0))
    }
    B1_values = apply(matrix(1:n,nrow=1),2,B1_value_f)
    
    total_numerator1 = Reduce("+", lapply(B1_values, `[[`, "numerator1"))
    total_numerator0 = Reduce("+", lapply(B1_values, `[[`, "numerator0"))
    denominator = adjusted_denom1_term 
    B1_trt = total_numerator1 / denominator
    B1_ctrl = total_numerator0 / denominator
    return(list(B1_trt = B1_trt,B1_ctrl = B1_ctrl))
  }
  B1 = B1_AOW(outcomes,treatments,covariates,pi_func, cond_prob_trt, cond_prob_ctrl) 
  
  ### B2_trt & B2_ctrl is a p*1 matrix
  B2_AOW = function(outcomes,treatments,pi_func, cond_prob_trt, cond_prob_ctrl) {
    B2_value_f = function(j){
      numerator_vec1 = numerator_vec0 = 0
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        e_Xi = pi_func[i]
        e_Xj = pi_func[j]
        mu_ij_trt = mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
        hij_xixj = (treatments[i] * (1 - treatments[j])) * (1-e_Xi) * e_Xj
        #term1 = (treatments[i] * (1 - treatments[j])) * (1-e_Xi) * e_Xj *
        #  ((outcomes[i] > outcomes[j]) - mu_ij_trt)
        term1 = hij_xixj * ((outcomes[i] > outcomes[j]) - mu_ij_trt)
        numerator_vec1 = numerator_vec1 + term1
        mu_ij_ctrl = mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
        #term2 = (treatments[i] * (1 - treatments[j])) * (1-e_Xi) * e_Xj *
        #  ((outcomes[i] < outcomes[j]) - mu_ij_ctrl)
        term0 = hij_xixj * ((outcomes[i] < outcomes[j]) - mu_ij_ctrl)
        numerator_vec0 = numerator_vec0 + term0
      }
      return(list(numerator1 = numerator_vec1, numerator0 = numerator_vec0))
    }
    
    B2_values = apply(matrix(1:n,nrow=1),2,B2_value_f)
    total_numerator1 = Reduce("+", lapply(B2_values, `[[`, "numerator1"))
    total_numerator0 = Reduce("+", lapply(B2_values, `[[`, "numerator0"))
    denominator = adjusted_denom1_term 
    B2_trt = total_numerator1 / denominator
    B2_ctrl = total_numerator0 / denominator
    return(list(B2_trt = B2_trt,B2_ctrl = B2_ctrl))
  }
  B2 = B2_AOW(outcomes,treatments,pi_func, cond_prob_trt, cond_prob_ctrl) 
  
  ### B3 ###
  B3_AOW_fun = function(treatments,covariates,pi_func) {
    B3_value_f = function(j){
      Xj = as.matrix(covariates[j,])
      numerator_vec = numeric(p)
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        Xi = as.matrix(covariates[i,])
        e_Xi = pi_func[i]
        e_Xj = pi_func[j]
        term = (treatments[i] * (1 - treatments[j])) * (1 - e_Xi) * e_Xj * 
          ((1-e_Xj) * Xj - e_Xi * Xi)
        numerator_vec = numerator_vec + term
      }
      return(list(numerator = numerator_vec))
    }
    B3_values = apply(matrix(1:n,nrow=1),2,B3_value_f)
    
    total_numerator = Reduce("+", lapply(B3_values, `[[`, "numerator"))
    return(total_numerator)
  }
  B3 = B3_AOW_fun(treatments,covariates,pi_func)
  
  ### B4 is a numeric
  B4_AOW_fun = function(treatments,pi_func) {
    B4_term_f = function(j){
      term_val = 0
      for (i in 1:n) {
        if (i != j) {
          term_val = term_val + treatments[i] * (1-pi_func[i])
        }
      }
      term_val = term_val * (1 - treatments[j]) * pi_func[j]
      return(list(B4_term_val = term_val))
    }
    B4_term = apply(matrix(1:n,nrow=1),2,B4_term_f)
    B4_AOW = Reduce("+", lapply(B4_term, `[[`, "B4_term_val"))
    return(B4_AOW)  
  }
  ##note: B4=adjusted_denom1_term
  B4 = adjusted_denom1_term
  
  ### B8 is a numeric, and B8 is the denominator of B5 and B6
  B8_AOW_fun = function(pi_func) {
    B8_term_f = function(j){
      term_val = 0
      for (i in 1:n) {
        if (i != j) {
          term_val = term_val + pi_func[i] * (1-pi_func[i])
        }
      }
      term_val = term_val * pi_func[j] * (1-pi_func[j])
      return(list(B8_term_val = term_val))
    }
    B8_term = apply(matrix(1:n,nrow=1),2,B8_term_f)
    B8_AOW = Reduce("+", lapply(B8_term, `[[`, "B8_term_val"))
    return(B8_AOW)  
  }
  #note: B8 = adjusted_denom2_term
  B8 = adjusted_denom2_term
  
  ### B5, p*1 vector
  B5_AOW = function(covariates, pi_func, cond_prob_trt, cond_prob_ctrl) {
    B5_value_f = function(j){
      numerator_vec1 = numerator_vec0 = numeric(p)
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        if (i != j) {
          Xi = as.matrix(covariates[i,])
          e_Xi = pi_func[i]
          e_Xj = pi_func[j]
          h_ij = e_Xi * (1-e_Xi) * e_Xj * (1-e_Xj)
          mu_ij_trt = mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
          term1 = h_ij * ((1 - 2 * e_Xi) * Xi) * mu_ij_trt
          numerator_vec1 = numerator_vec1 + term1
          mu_ij_ctrl = mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
          term0 = h_ij * ((1 - 2 * e_Xi) * Xi) * mu_ij_ctrl
          numerator_vec0 = numerator_vec0 + term0
        }
      }
      return(list(numerator1 = numerator_vec1, numerator0 = numerator_vec0))
    }
    B5_values = apply(matrix(1:n,nrow=1),2,B5_value_f)
    
    total_numerator1 = Reduce("+", lapply(B5_values, `[[`, "numerator1"))
    total_numerator0 = Reduce("+", lapply(B5_values, `[[`, "numerator0"))
    denominator = adjusted_denom2_term
    B5_trt = total_numerator1 / denominator
    B5_ctrl = total_numerator0 / denominator
    return(list(B5_trt=B5_trt,B5_ctrl=B5_ctrl))
  }
  B5 = B5_AOW(covariates, pi_func, cond_prob_trt, cond_prob_ctrl)
  
  ### B6 is a numeric
  B6_AOW = function(pi_func, cond_prob_trt, cond_prob_ctrl) {
    B6_value_f = function(j){
      numerator_val1 = numerator_val0 = 0
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        if (i != j) {
          e_Xi = pi_func[i]
          e_Xj = pi_func[j]
          h_ij = e_Xi * (1-e_Xi) * e_Xj * (1-e_Xj)
          mu_ij_trt = mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
          term1 = h_ij * mu_ij_trt
          numerator_val1 = numerator_val1 + term1
          mu_ij_ctrl = mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
          term0 = h_ij * mu_ij_ctrl
          numerator_val0 = numerator_val0 + term0
        }
      }
      return(list(numerator1 = numerator_val1, numerator0 = numerator_val0))
    }
    B6_values = apply(matrix(1:n,nrow=1),2,B6_value_f)
    
    total_numerator1 = Reduce("+", lapply(B6_values, `[[`, "numerator1"))
    total_numerator0 = Reduce("+", lapply(B6_values, `[[`, "numerator0"))
    denominator = adjusted_denom2_term
    B6_trt = total_numerator1 / denominator
    B6_ctrl = total_numerator0 / denominator
    return(list(B6_trt=B6_trt,B6_ctrl=B6_ctrl))
  }
  B6 = B6_AOW(pi_func, cond_prob_trt, cond_prob_ctrl)
  
  ### B7 is p*1 vector
  B7_AOW_fun = function(covariates, pi_func) {
    B7_term_f = function(j){
      term1_vec =  term2_vec = numeric(p)
      Xj = as.matrix(covariates[j,])
      for (i in 1:n) {
        if (i != j) {
          Xi = as.matrix(covariates[i,])
          term1_vec = term1_vec + pi_func[i] * (1-pi_func[i]) * ((1-2*pi_func[i]) * Xi)
          term2_vec = term2_vec + (pi_func[i] * (1-pi_func[i]))
        }
      }
      term1_vec = term1_vec * pi_func[j] * (1 - pi_func[j])
      term2_vec = term2_vec * pi_func[j] * (1 - pi_func[j]) * (1 - 2*pi_func[j]) * Xj
      return(list(term1 = term1_vec,term2 = term2_vec))
    }
    B7_term = apply(matrix(1:n,nrow=1),2,B7_term_f)
    total_term1 = Reduce("+", lapply(B7_term, `[[`, "term1"))
    total_term2 = Reduce("+", lapply(B7_term, `[[`, "term2"))
    return(total_term1 + total_term2)  
  }
  B7 = B7_AOW_fun(covariates,pi_func)
  
  #Bn^{trt}
  B1_trt = B1$B1_trt; B2_trt = B2$B2_trt
  B5_trt = B5$B5_trt; B6_trt = B6$B6_trt
  Bn_trt = B1_trt-B3*(B2_trt/B4)+B5_trt-B7*(B6_trt/B8)
  #Bn^{ctrl}
  B1_ctrl = B1$B1_ctrl; B2_ctrl = B2$B2_ctrl
  B5_ctrl = B5$B5_ctrl; B6_ctrl = B6$B6_ctrl
  Bn_ctrl = B1_ctrl-B3*(B2_ctrl/B4)+B5_ctrl-B7*(B6_ctrl/B8)
  
  ### E_beta0beta0 p*p
  X = as.matrix(covariates)
  ps = pi_func
  Ebeta0 = crossprod(sqrt(ps*(1-ps)) * X) / n
  
  ###influence function of \beta, l_beta 
  inv_Ebeta0 = ginv(Ebeta0)
  l_beta_fun = function(i) {
    X_i = as.matrix(covariates[i,])
    Z_i = treatments[i]
    ps_i = pi_func[i]
    l_beta = inv_Ebeta0 %*% X_i * (Z_i - ps_i)
    return(l_beta)
  }
  
  ### E_theta_theta
  thresholds_trt = dr_trt$zeta
  thresholds_ctrl = dr_ctrl$zeta
  
  E_thetatheta_fun = function(outcomes, treatments, covariates_reg, cond_prob_trt, cond_prob_ctrl, thresholds_trt, thresholds_ctrl) {
    F_matrices = lapply(1:n, function(i) {
      Xi = as.matrix(covariates_reg[i,])
      thresholds = if(treatments[i] == 1) thresholds_trt else thresholds_ctrl
      cond_prob = if(treatments[i] == 1) cond_prob_trt else cond_prob_ctrl
      # Calculate delta values based on the outcome for individual i
      delta = as.numeric(outcomes[i] == 1:L) #a 1*L vector
      if(L==3){
        # Calculate the c values for individual i
        c1 = cond_prob[i, 1]
        c2 = cond_prob[i, 1] + cond_prob[i,2]
        # Calculate the F matrix components for individual i
        F11 = -(delta[1] + delta[2]) * (c1 - c1^2) - delta[2] * exp(thresholds[2] - thresholds[1]) *
          ((exp(thresholds[2] - thresholds[1]) - 1)^(-2))
        F12 = delta[2] * exp(thresholds[2] - thresholds[1]) * ((exp(thresholds[2] - thresholds[1]) - 1)^(-2))
        F21 = F12
        F13 = (delta[1] + delta[2]) * (c1 - c1^2) * t(Xi)
        F31 = t(F13)
        F22 = -(delta[2] + delta[3]) * (c2 - c2^2) - delta[2] * exp(thresholds[1] - thresholds[2]) * 
          (1 - exp(thresholds[1] - thresholds[2]))^(-2)
        F23 = (delta[2] + delta[3]) * (c2 - c2^2) * t(Xi)
        F32 = t(F23)
        F33 = -((delta[1] + delta[2]) * (c1 - c1^2) + (delta[2] + delta[3]) * (c2 - c2^2)) * (Xi %*% t(Xi))
        
        # Construct the F matrix for individual i
        top_left = rbind(cbind(F11, F12), cbind(F21, F22))
        top_2row = cbind(top_left, rbind(F13, F23))
        bottom_row = cbind(cbind(F31, F32), F33)
        F_matrix = as.matrix(rbind(top_2row, bottom_row))
      }else{ #in this case L>3
        k = L-1+po
        FM = matrix(0,k,k)
        #the first row
        c1 = cond_prob[i,1]
        f11 = -sum(delta[1:2])*(c1-c1^2)-delta[2]*exp(thresholds[2]-thresholds[1])*((exp(thresholds[2]-thresholds[1])-1)^(-2))
        f12 = delta[2]*exp(thresholds[2]-thresholds[1])*((exp(thresholds[2]-thresholds[1])-1)^(-2))
        #3:(L-1) columns
        f1lm1 = rep(0,L-3)
        f1l = sum(delta[1:2])*(c1 - c1^2)*t(Xi) #1*po
        FM[1,] = c(f11,f12,f1lm1,f1l)
        #2:(L-2) row
        for(l in 2:(L-2)){
          for(j in 1:k){
            if (j < l) FM[l,j] = FM[j,l]
            else if (j == l){
              cl = sum(cond_prob[i,1:l])
              temp1 = sum(delta[l:(l+1)])*(cl-cl^2)
              temp2 = delta[l]*exp(thresholds[l-1]-thresholds[l])*(1-exp(thresholds[l-1]-thresholds[l]))^(-2)
              temp3 = delta[l+1]*exp(thresholds[l+1]-thresholds[l])*((exp(thresholds[l+1]-thresholds[l])-1)^(-2))
              FM[l,j] = -temp1-temp2-temp3
            }
            else if (j==(l+1)) FM[l,j] = temp3
            else if (j<=(L-1)) FM[l,j] = 0
            else FM[l,L:k] = temp1*t(Xi)
          }
        }
        #the L-1 row
        for(j in 1:(L-2)){
          FM[L-1,j] = FM[j,L-1]
        }
        cLm1 = sum(cond_prob[i,1:(L-1)])
        temp1 = sum(delta[(L-1):L])*(cLm1-cLm1^2)
        temp2 = delta[L-1]*exp(thresholds[L-2]-thresholds[L-1])*(1-exp(thresholds[L-2]-thresholds[L-1]))^(-2)
        FM[L-1,L-1] = -temp1-temp2
        FM[L-1,L:k] = temp1*t(Xi)
        #row--L:k, 1:(L-1) columns
        for(l in L:k){
          for(j in 1:(L-1)){
            FM[l,j] = FM[j,l]
          }
        }
        #row--L:k, L:k columns
        factor_val = 0
        for(j in 1:(L-1)){
          cj = sum(cond_prob[i,1:j])
          temp = sum(delta[j:(j+1)])*(cj-cj^2)
          factor_val = factor_val+temp
        }
        FM[L:k,L:k] = -factor_val* (Xi %*% t(Xi))
        F_matrix = FM
      }
      return(F_matrix)
    })
    # Compute the average of the F matrices across all individuals
    E_thetatheta = Reduce("+", F_matrices) / n
    return(E_thetatheta)
  }
  E_thetatheta = E_thetatheta_fun(outcomes, treatments, covariates_reg, cond_prob_trt, cond_prob_ctrl, thresholds_trt, thresholds_ctrl)
  
  ### mu(Xi, Xj) p*1
  mu_tilde_AOW_fun = function(i, j, covariates_reg, cond_prob_trt, cond_prob_ctrl) {
    # Extract Xi and Xj as row vectors
    Xi = as.matrix(covariates_reg[i,])
    Xj = as.matrix(covariates_reg[j,])
    if(L==3){
      # Calculate the c values for individuals i and j
      ci_1 = cond_prob_trt[i, 1]
      ci_2 = cond_prob_trt[i, 1] + cond_prob_trt[i, 2]
      cj_1 = cond_prob_ctrl[j, 1]
      cj_2 = cond_prob_ctrl[j, 1] + cond_prob_ctrl[j, 2]
      
      # Construct the mu_tilde vector according to the formula
      mu_tilde_vector_1 = -ci_1 * (1 - ci_1) * cj_1 + cj_1 * (1 - cj_1) * (ci_2 - ci_1)
      mu_tilde_vector_2 = -ci_2 * (1 - ci_2) * (cj_2 - cj_1) + (1 - ci_2) * cj_2 * (1 - cj_2)
      mu_tilde_vector_3 = (cj_1 * (ci_1 - ci_1^2) + (ci_2 - ci_2^2) * (cj_2 - cj_1)) * Xi - 
        ((ci_2 - ci_1) * (cj_1 - cj_1^2) + (1 - ci_2) * (cj_2 - cj_2^2)) * Xj
      mu_tilde_vector = rbind(mu_tilde_vector_1, mu_tilde_vector_2,
                              mu_tilde_vector_3)
    }else{
      mu_tilde_vector_l = matrix(0,L-1,1)
      #the first element
      ci1 = cond_prob_trt[i,1]; ci2 = sum(cond_prob_trt[i,1:2])
      cj1 = cond_prob_ctrl[j,1];
      mu_tilde_vector_l[1] = -(ci1-ci1^2)*cj1+(ci2-ci1)*(cj1-cj1^2)
      #the 2 to (L-2) element
      for(l in 2:(L-2)){
        cil = sum(cond_prob_trt[i,1:l]); cilp1 = sum(cond_prob_trt[i,1:(l+1)])
        cjlm1 = sum(cond_prob_ctrl[j,1:(l-1)]); cjl = sum(cond_prob_ctrl[j,1:l])
        mu_tilde_vector_l[l] = (cil-cil^2)*(cjlm1-cjl)+(cilp1-cil)*(cjl-cjl^2)
      }
      #the (L-1)th element
      ciLm1 = sum(cond_prob_trt[i,1:(L-1)]); cjLm2 = sum(cond_prob_ctrl[j,1:(L-2)]) 
      cjLm1 = sum(cond_prob_ctrl[j,1:(L-1)])
      mu_tilde_vector_l[L-1] = (ciLm1-ciLm1^2)*(cjLm2-cjLm1)+(1-ciLm1)*(cjLm1-cjLm1^2)
      #calculate the final term
      factor1 = factor2 = 0
      for(l in 1:(L-2)){
        cil = sum(cond_prob_trt[i,1:l]); cilp1 = sum(cond_prob_trt[i,1:(l+1)])
        cjl = sum(cond_prob_ctrl[j,1:l])
        temp1 = ((cil-cil^2)-(cilp1-cilp1^2))*cjl
        factor1 = factor1 + temp1
        temp2 = (cilp1-cil)*(cjl-cjl^2)
        factor2 = factor2 + temp2
      }
      ciLm1 = sum(cond_prob_trt[i,1:(L-1)]) 
      cjLm1 = sum(cond_prob_ctrl[j,1:(L-1)])
      factor1 = factor1 + (ciLm1-ciLm1^2)*cjLm1
      factor2 = factor2 + (1-ciLm1)*(cjLm1-cjLm1^2)
      mu_tilde_vector_final = factor1*Xi - factor2*Xj
      mu_tilde_vector = rbind(mu_tilde_vector_l,mu_tilde_vector_final)
    }
    return(mu_tilde_vector)
  }
  
  ### Cn p*1 vector
  Cn_AOW_fun = function(treatments, covariates_reg, pi_func, cond_prob_trt, cond_prob_ctrl) {
    # Calculate the first term of Cn
    Cn_term_f = function(j){
      numerator_vec1 = numerator_vec2 = 0
      # Loop through all control indices and calculate contributions to numerator and denominator
      for (i in 1:n) {
        if(i != j){
          e_Xi = pi_func[i]
          e_Xj = pi_func[j]
          mu_tilde_ij = mu_tilde_AOW_fun(i, j, covariates_reg, cond_prob_trt, cond_prob_ctrl)
          term1 = (treatments[i] * (1 - treatments[j])) * (1-e_Xi) * e_Xj * (-mu_tilde_ij)
          numerator_vec1 = numerator_vec1 + term1
          
          h_ij = e_Xi * (1-e_Xi) * e_Xj * (1-e_Xj)
          term2 = h_ij * (mu_tilde_ij)
          numerator_vec2 = numerator_vec2 + term2
        }
      }
      return(list(numerator1 = numerator_vec1,numerator2 = numerator_vec2))
    }
    Cn_term = apply(matrix(1:n,nrow=1),2,Cn_term_f)
    term1_numerator = Reduce("+", lapply(Cn_term, `[[`, "numerator1"))
    term2_numerator = Reduce("+", lapply(Cn_term, `[[`, "numerator2"))
    term1_denominator = adjusted_denom1_term
    first_term = term1_numerator / term1_denominator
    term2_denominator = adjusted_denom2_term
    second_term = term2_numerator / term2_denominator
    Cn_value = first_term + second_term
    return(Cn_value)
  }
  # Use the function to calculate Cn
  Cn_value = Cn_AOW_fun(treatments, covariates_reg, pi_func, cond_prob_trt, cond_prob_ctrl)
  
  #influence function of \theta, l_tilde_theta
  E_thetatheta_inv = ginv(E_thetatheta)
  l_theta_AOW_fun = function(i,outcomes,treatments,covariates_reg, cond_prob_trt, cond_prob_ctrl, 
                              thresholds_trt, thresholds_ctrl, E_thetatheta_inv) {
    Xi = as.matrix(covariates_reg[i,])
    cond_prob = if(treatments[i] == 1) cond_prob_trt else cond_prob_ctrl
    delta = as.numeric(outcomes[i] == 1:L)
    thresholds = if (treatments[i] == 1) thresholds_trt else thresholds_ctrl
    if(L==3){
      c1 = cond_prob[i, 1]
      c2 = c1 + cond_prob[i, 2]
      # Construct the l_tilde vector according to the formula
      l_tilde_vector_1 = delta[1] - (delta[1] + delta[2]) * c1 - delta[2] * (exp(thresholds[2] - thresholds[1])-1)^(-1)
      l_tilde_vector_2 = -(delta[2] + delta[3]) * c2 + delta[2] * (1 - exp(thresholds[1] - thresholds[2]))^(-1)
      l_tilde_vector_3 = (-(delta[1] + delta[2]) + (delta[1] + delta[2]) * c1 + (delta[2] + delta[3]) * c2) * Xi
      l_tilde_vector = rbind(l_tilde_vector_1, l_tilde_vector_2, l_tilde_vector_3)
    }else{
      l_tilde_vector_l = matrix(0,L-1,1)
      c1 = cond_prob[i,1]
      l_tilde_vector_l[1] = delta[1]-sum(delta[1:2])*c1-delta[2]*(exp(thresholds[2]-thresholds[1])-1)^(-1)
      for(l in 2:(L-2)){
        cl = sum(cond_prob[i,1:l])
        l_tilde_vector_l[l] = -sum(delta[l:(l+1)])*cl+delta[l]*1/(1-exp(thresholds[l-1]-thresholds[l]))-delta[l+1]*1/(exp(thresholds[l+1]-thresholds[l])-1)
      }
      cLm1 = sum(cond_prob[i,1:(L-1)])
      l_tilde_vector_l[L-1] = -sum(delta[(L-1):L])*cLm1+delta[L-1]*1/(1-exp(thresholds[L-2]-thresholds[L-1]))
      factor_val = 0
      for(l in 1:(L-1)){
        cl = sum(cond_prob[i,1:l])
        temp = sum(delta[l:(l+1)])*cl
        factor_val = factor_val+temp
      }
      factor_val = factor_val-sum(delta[1:(L-1)])
      l_tilde_vector_finall = factor_val*Xi
      l_tilde_vector = rbind(l_tilde_vector_l,l_tilde_vector_finall)
    }
    l_tilde_value = E_thetatheta_inv %*% l_tilde_vector
    return(l_tilde_value)
  }
  
  #variance and covariance of WP and LP
  var_trt_sum = var_ctrl_sum = cov_sum = 0
  Bn_t_trt = t(Bn_trt)
  Bn_t_ctrl = t(Bn_ctrl)
  Cn_t = t(Cn_value)
  for(i in 1:n){
    l_theta = l_theta_AOW_fun(i,outcomes,treatments,covariates_reg,cond_prob_trt, cond_prob_ctrl, 
                              thresholds_trt, thresholds_ctrl, E_thetatheta_inv)
    l_beta = l_beta_fun(i)
    temp_trt = (2*(g1_AOW_trt[i]-WP_trt_aow) + Bn_t_trt %*% l_beta + Cn_t %*% l_theta)
    temp_ctrl = (2*(g1_AOW_ctrl[i]-WP_ctrl_aow) + Bn_t_ctrl %*% l_beta + Cn_t %*% l_theta)
    var_trt_sum = var_trt_sum + temp_trt^2
    var_ctrl_sum = var_ctrl_sum + temp_ctrl^2
    cov_sum = cov_sum + temp_trt*temp_ctrl
  }
  AOW_Var_trt = var_trt_sum/n
  AOW_Var_ctrl = var_ctrl_sum/n
  AOW_cov = cov_sum/n
  
  var_trt_aow = AOW_Var_trt/n
  var_ctrl_aow = AOW_Var_ctrl/n
  cov_aow = AOW_cov/n
  
  #approximation variance of WR
  var_ratio = function(uW,uL,varW,varL,covWL){
    WR_appro_var = (uW^2/uL^2)*(varW/uW^2-2*covWL/(uW*uL)+varL/uL^2)
    return(WR_appro_var)
  }
  
  ###Variance of WR and WD
  var_WR_aow = var_ratio(WP_trt_aow,WP_ctrl_aow,var_trt_aow,var_ctrl_aow,cov_aow)
  var_WD_aow = var_trt_aow + var_ctrl_aow - 2*cov_aow
  #################AOW Estimator Output#########################
  Estimation = c(WP_trt_aow,WP_ctrl_aow,WR_aow,WD_aow)
  SE = c(sqrt(var_trt_aow),sqrt(var_ctrl_aow),sqrt(var_WR_aow),sqrt(var_WD_aow))
  conf_low = Estimation-qnorm(0.975)*SE
  conf_high = Estimation+qnorm(0.975)*SE
  res = data.frame(EST=Estimation,SE=SE,conf_low=conf_low,conf_high=conf_high)
  row.names(res) = c("WP","LP","WR","WD")
  return(res)
}

