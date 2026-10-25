#' Estimate win statistics by AOW method
#'  
#' The AOW (i.e., augmented overlap weight) method is a Hajek-type estimator by leveraging U-statistic theory,  see section 4 of paper (Cao et al., 2026).
#' The point estimate is based on formula (18) of paper with OW weights, and variance estimate is based on formula (A.20) of Web Appendix with OW weights. 
#' In addition to estimating point estimates of four win statistics, that is, win probability (WP), loss probability (LP), 
#' win ratio (WR) and win difference (WD), their standard errors as well as 95% confidence intervals will also be returned
#' 
#' @param data A dataset of randomized controlled trials (RCT), which should include ordinal outcome, treatment and covariates.
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
#' res1 = winstat_aow_new(data=mydata, outcomevar="outcome", treatment="treatment",
#'       covariate=paste0("x",1:6,seq=""),
#'       covariate_reg=c(paste0("x",1:6,seq=""),paste0("x",1:3,"s",seq="")))
#' print(res1)
#' 


winstat_aow_new = function(data,outcomevar,treatment,covariate,covariate_reg){
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
  
  #number of categories in ordinal outcome
  L = length(unique(outcomes))
  
  ### mu-function
  lp_start = rep(0, po)
  th_start = c(-1, 1)
  start_values = c(lp_start, th_start)
  model.formual = as.formula(paste("factor(",paste(outcomevar),")","~",paste(covariate_reg,collapse="+"),sep=""))
  dr_trt = polr(model.formual,data = data[treatments==1,], start=start_values)
  cond_prob_trt = predict(dr_trt, newdata = data, type = "probs")
  
  dr_ctrl = polr(model.formual, data = data[treatments == 0,], start=start_values)
  cond_prob_ctrl = predict(dr_ctrl, newdata = data, type = "probs")
  
  compare_rows_AOW = function(i) {
    aw_aow_trt_vec = aw_aow_ctrl_vec = numeric(n)
    b_aow_trt_vec = b_aow_ctrl_vec = numeric(n)
    cw_aow_trt_vec = cw_aow_ctrl_vec = numeric(n)
    d_aow_trt_vec = d_aow_ctrl_vec = numeric(n)
    
    for(j in 1:n){
      if(i != j){
        #aw
        #mu_ij (i>j)
        mu_ij_trt = 0
        for(l in 2:L){
          temp1 = cond_prob_trt[i,l] * (sum(cond_prob_ctrl[j,1:(l-1)]))
          mu_ij_trt = mu_ij_trt+temp1
        }
        kij = (treatments[i] * (1 - treatments[j])) * ((1-pi_func[i]) * pi_func[j])
        aw_aow_trt_vec[j] = kij*(as.numeric(outcomes[i]>outcomes[j]) - mu_ij_trt)
        
        #mu_ji(j>i)
        mu_ij_ctrl = 0
        for(l in 1:(L-1)){
          temp2 = cond_prob_trt[i, l] * (sum(cond_prob_ctrl[j,(l+1):L]))
          mu_ij_ctrl = mu_ij_ctrl+temp2
        }
        aw_aow_ctrl_vec[j] = kij*(as.numeric(outcomes[i]<outcomes[j]) - mu_ij_ctrl)
        
        ##b
        b_aow_trt_vec[j] = kij
        b_aow_ctrl_vec[j] = b_aow_trt_vec[j]
        
        ##cw
        hij = pi_func[i] * pi_func[j] * (1-pi_func[i]) *(1-pi_func[j])
        cw_aow_trt_vec[j] = hij * mu_ij_trt 
        cw_aow_ctrl_vec[j] = hij * mu_ij_ctrl 
        
        ##d
        d_aow_trt_vec[j] = hij
        d_aow_ctrl_vec[j] = d_aow_trt_vec[j]
      }
    }
    
    return(list(aw_trt = aw_aow_trt_vec, aw_ctrl = aw_aow_ctrl_vec, 
                b_trt = b_aow_trt_vec, b_ctrl = b_aow_ctrl_vec,
                cw_trt = cw_aow_trt_vec, cw_ctrl = cw_aow_ctrl_vec,
                d_trt = d_aow_trt_vec, d_ctrl = d_aow_ctrl_vec))
  }
  
  results_AOW = apply(matrix(1:n,nrow=1),2,compare_rows_AOW)
  
  
  # Combine the results
  g_aw_AOW1 = do.call(rbind, lapply(results_AOW, function(x) x$aw_trt))
  g_aw_AOW0 = do.call(rbind, lapply(results_AOW, function(x) x$aw_ctrl))
  g_b_AOW1 = do.call(rbind, lapply(results_AOW, function(x) x$b_trt))
  g_b_AOW0 = do.call(rbind, lapply(results_AOW, function(x) x$b_ctrl))
  g_cw_AOW1 = do.call(rbind, lapply(results_AOW, function(x) x$cw_trt))
  g_cw_AOW0 = do.call(rbind, lapply(results_AOW, function(x) x$cw_ctrl))
  g_d_AOW1 = do.call(rbind, lapply(results_AOW, function(x) x$d_trt))
  g_d_AOW0 = do.call(rbind, lapply(results_AOW, function(x) x$d_ctrl))
  
  A_AOW1 = sum(as.vector(g_aw_AOW1), na.rm = TRUE)/(n*(n-1))
  A_AOW0 = sum(as.vector(g_aw_AOW0), na.rm = TRUE)/(n*(n-1))
  B_AOW1 = sum(as.vector(g_b_AOW1), na.rm = TRUE)/(n*(n-1))
  B_AOW0 = sum(as.vector(g_b_AOW0), na.rm = TRUE)/(n*(n-1))
  C_AOW1 = sum(as.vector(g_cw_AOW1), na.rm = TRUE)/(n*(n-1))
  C_AOW0 = sum(as.vector(g_cw_AOW0), na.rm = TRUE)/(n*(n-1))
  D_AOW1 = sum(as.vector(g_d_AOW1), na.rm = TRUE)/(n*(n-1))
  D_AOW0 = sum(as.vector(g_d_AOW0), na.rm = TRUE)/(n*(n-1))
  
  rho1_aow_trt = A_AOW1/B_AOW1
  rho1_aow_ctrl = A_AOW0/B_AOW0
  rho2_aow_trt = C_AOW1/D_AOW1
  rho2_aow_ctrl = C_AOW0/D_AOW0
  
  tau1_AOW = rho1_aow_trt + rho2_aow_trt 
  tau0_AOW = rho1_aow_ctrl + rho2_aow_ctrl
  
  #point estimates of WP, LP, WR and WD
  WP_trt_aow = tau1_AOW
  WP_ctrl_aow = tau0_AOW
  WR_aow = tau1_AOW / tau0_AOW
  WD_aow = tau1_AOW - tau0_AOW
  
  ########### AOW Theoretical Variance###########
  mu_func_trt = function(a, b, cond_prob_trt, cond_prob_ctrl) {
    mu_a_b = 0
    for(l in 2:L){
      temp1 = cond_prob_trt[a,l] * (sum(cond_prob_ctrl[b,1:(l-1)]))
      mu_a_b = mu_a_b+temp1
    }
    return(mu_a_b)
  }
  
  mu_func_ctrl = function(a, b, cond_prob_trt, cond_prob_ctrl) {
    mu_a_b = 0
    for(l in 1:(L-1)){
      temp2 = cond_prob_trt[a,l] * (sum(cond_prob_ctrl[b,(l+1):L]))
      mu_a_b = mu_a_b+temp2
    }
    return(mu_a_b)
  }
  
  #calculate a1^{AOW}, b1^{AOW}, c1^{AOW} and d1^{AOW}
  aw1_b1_cw1_d1_AOW_fun = function(i){
    #for treatment
    mu_ij = sapply(1:n, function(b) mu_func_trt(i, b, cond_prob_trt, cond_prob_ctrl)) #i>j
    mu_ji = sapply(1:n, function(b) mu_func_ctrl(i, b, cond_prob_trt, cond_prob_ctrl)) #i<j
    
    kij = (treatments[i] * (1 - treatments) * (1 - pi_func[i]) * pi_func)
    kji = ((1 - treatments[i]) * treatments * pi_func[i] * (1 - pi_func))
    aw_aow1 = (1/2)*(kij*(as.numeric(outcomes[i] > outcomes)-mu_ij)
                     +kji*(as.numeric(outcomes > outcomes[i])-mu_ji))
    b_aow1 = (1/2) * (kij + kji)
    hij = pi_func[i] * (1 - pi_func[i]) * pi_func * (1 - pi_func)  
    cw_aow1 = (1/2) * hij * (mu_ij + mu_ji)
    d_aow1 = (1/2) *(hij + hij)
    sum_aw1_AOW_trt = sum(aw_aow1)
    sum_b1_AOW_trt = sum(b_aow1)
    sum_cw1_AOW_trt = sum(cw_aow1)
    sum_d1_AOW_trt = sum(d_aow1)
    
    #for control 
    aw_aow0 = (1/2)*(kij*(as.numeric(outcomes[i] < outcomes)-mu_ji)
                     +kji*(as.numeric(outcomes < outcomes[i])-mu_ij)) 
    b_aow0 = b_aow1
    cw_aow0 = cw_aow1
    d_aow0 = d_aow1
    sum_aw1_AOW_ctrl = sum(aw_aow0)
    sum_b1_AOW_ctrl = sum(b_aow0)
    sum_cw1_AOW_ctrl = sum(cw_aow0)
    sum_d1_AOW_ctrl = sum(d_aow0)
    
    return(list(aw1_trt = (1/(n-1))*sum_aw1_AOW_trt, 
                aw1_ctrl = (1/(n-1))*sum_aw1_AOW_ctrl,
                b1_trt = (1/(n-1))*sum_b1_AOW_trt,
                b1_ctrl = (1/(n-1))*sum_b1_AOW_ctrl,
                cw1_trt = (1/(n-1))*sum_cw1_AOW_trt, 
                cw1_ctrl = (1/(n-1))*sum_cw1_AOW_ctrl,
                d1_trt = (1/(n-1))*sum_d1_AOW_trt,
                d1_ctrl = (1/(n-1))*sum_d1_AOW_ctrl))
  }
  
  results_aw1_b1_cw1_d1 = apply(matrix(1:n,nrow=1),2,aw1_b1_cw1_d1_AOW_fun)
  # Aggregate results
  aw1_AOW_trt = as.numeric(do.call(rbind, lapply(results_aw1_b1_cw1_d1, function(x) x$aw1_trt)))
  aw1_AOW_ctrl = as.numeric(do.call(rbind, lapply(results_aw1_b1_cw1_d1, function(x) x$aw1_ctrl)))
  b1_AOW_trt = as.numeric(do.call(rbind, lapply(results_aw1_b1_cw1_d1, function(x) x$b1_trt)))
  b1_AOW_ctrl = as.numeric(do.call(rbind, lapply(results_aw1_b1_cw1_d1, function(x) x$b1_ctrl)))
  cw1_AOW_trt = as.numeric(do.call(rbind, lapply(results_aw1_b1_cw1_d1, function(x) x$cw1_trt)))
  cw1_AOW_ctrl = as.numeric(do.call(rbind, lapply(results_aw1_b1_cw1_d1, function(x) x$cw1_ctrl)))
  d1_AOW_trt = as.numeric(do.call(rbind, lapply(results_aw1_b1_cw1_d1, function(x) x$d1_trt)))
  d1_AOW_ctrl = as.numeric(do.call(rbind, lapply(results_aw1_b1_cw1_d1, function(x) x$d1_ctrl)))
  
  ###U_n\nabla_{\beta}a_w^{\text{AIPW}}(\widehat{\eta})
  ###U_n\nabla_{\beta}b^{\text{AIPW}}(\widehat{\beta})
  ###U_n\nabla_{\beta}c_w^{\text{AIPW}}(\widehat{\eta})
  ###U_n\nabla_{\beta}d^{\text{AIPW}}(\widehat{\beta})
  d_aw_b_cw_d_AOW =function(outcomes,treatments,covariates,pi_func, cond_prob_trt, cond_prob_ctrl) {
    cal_aw_b_cw_d_derivative = function(i){
      Xi = matrix(covariates[i,])
      aw_vec1 = aw_vec0 = b_vec1 = b_vec0 = numeric(p)
      cw_vec1 = cw_vec0 = d_vec1 = d_vec0 = numeric(p)
      for (j in 1:n) {
        if(j != i){
          Xj = matrix(covariates[j,])
          e_Xi = pi_func[i]
          e_Xj = pi_func[j]
          mu_ij_trt = mu_func_trt(i, j, cond_prob_trt, cond_prob_ctrl)
          mu_ij_ctrl = mu_func_ctrl(i, j, cond_prob_trt, cond_prob_ctrl)
          kij_xixj = (treatments[i]*(1-treatments[j]))*(1-e_Xi)*e_Xj*((1-e_Xj)*Xj-e_Xi*Xi)
          term1 = kij_xixj*(as.numeric(outcomes[i] > outcomes[j])-mu_ij_trt)
          term0 = kij_xixj*(as.numeric(outcomes[i] < outcomes[j])-mu_ij_ctrl)  
          aw_vec1 = aw_vec1 + term1
          aw_vec0 = aw_vec0 + term0
          b_vec1 = b_vec1 + kij_xixj
          b_vec0 = b_vec0 + kij_xixj
          hij_xixj = e_Xi*e_Xj*(1-e_Xi)*(1-e_Xj)*((1-2*e_Xi)*Xi+(1-2*e_Xj)*Xj)
          term1f = hij_xixj*mu_ij_trt
          term0f = hij_xixj*mu_ij_ctrl
          cw_vec1 = cw_vec1 + term1f
          cw_vec0 = cw_vec0 + term0f
          d_vec1 = d_vec1 + hij_xixj
          d_vec0 = d_vec0 + hij_xixj
        }
      }
      return(list(aw1 = aw_vec1, aw0 = aw_vec0,
                  b1 = b_vec1, b0 = b_vec0,
                  cw1 = cw_vec1, cw0 = cw_vec0,
                  d1 = d_vec1, d0 = d_vec0))
    }
    aw_b_cw_d_derivative_temp = apply(matrix(1:n,nrow=1),2,cal_aw_b_cw_d_derivative)
    total_aw_der1 = Reduce("+", lapply(aw_b_cw_d_derivative_temp, `[[`, "aw1"))
    total_aw_der0 = Reduce("+", lapply(aw_b_cw_d_derivative_temp, `[[`, "aw0"))
    total_b_der1 = Reduce("+", lapply(aw_b_cw_d_derivative_temp, `[[`, "b1"))
    total_b_der0 = Reduce("+", lapply(aw_b_cw_d_derivative_temp, `[[`, "b0"))
    total_cw_der1 = Reduce("+", lapply(aw_b_cw_d_derivative_temp, `[[`, "cw1"))
    total_cw_der0 = Reduce("+", lapply(aw_b_cw_d_derivative_temp, `[[`, "cw0"))
    total_d_der1 = Reduce("+", lapply(aw_b_cw_d_derivative_temp, `[[`, "d1"))
    total_d_der0 = Reduce("+", lapply(aw_b_cw_d_derivative_temp, `[[`, "d0"))
    
    der_aw_AOW_trt = total_aw_der1/(n*(n-1))
    der_aw_AOW_ctrl = total_aw_der0/(n*(n-1))
    der_b_AOW_trt = total_b_der1/(n*(n-1))
    der_b_AOW_ctrl = total_b_der0/(n*(n-1))
    der_cw_AOW_trt = total_cw_der1/(n*(n-1))
    der_cw_AOW_ctrl = total_cw_der0/(n*(n-1))
    der_d_AOW_trt = total_d_der1/(n*(n-1))
    der_d_AOW_ctrl = total_d_der0/(n*(n-1))
    return(list(der_aw_AOW_trt = der_aw_AOW_trt, der_aw_AOW_ctrl = der_aw_AOW_ctrl,
                der_b_AOW_trt = der_b_AOW_trt, der_b_AOW_ctrl = der_b_AOW_ctrl,
                der_cw_AOW_trt = der_cw_AOW_trt, der_cw_AOW_ctrl = der_cw_AOW_ctrl,
                der_d_AOW_trt = der_d_AOW_trt, der_d_AOW_ctrl = der_d_AOW_ctrl))
  }
  der_aw_b_cw_d_AOW_res = d_aw_b_cw_d_AOW(outcomes,treatments,covariates,pi_func,cond_prob_trt,cond_prob_ctrl)
  
  der_aw_AOW_trt = matrix(der_aw_b_cw_d_AOW_res$der_aw_AOW_trt, nrow=1)
  der_aw_AOW_ctrl = matrix(der_aw_b_cw_d_AOW_res$der_aw_AOW_ctrl, nrow=1)
  der_b_AOW_trt = matrix(der_aw_b_cw_d_AOW_res$der_b_AOW_trt, nrow=1)
  der_b_AOW_ctrl = matrix(der_aw_b_cw_d_AOW_res$der_b_AOW_ctrl, nrow=1)
  der_cw_AOW_trt = matrix(der_aw_b_cw_d_AOW_res$der_cw_AOW_trt, nrow=1)
  der_cw_AOW_ctrl = matrix(der_aw_b_cw_d_AOW_res$der_cw_AOW_ctrl, nrow=1)
  der_d_AOW_trt = matrix(der_aw_b_cw_d_AOW_res$der_d_AOW_trt, nrow=1)
  der_d_AOW_ctrl = matrix(der_aw_b_cw_d_AOW_res$der_d_AOW_ctrl, nrow=1)
  
  G_beta_AOW_trt = (der_aw_AOW_trt-rho1_aow_trt*der_b_AOW_trt)/B_AOW1+(der_cw_AOW_trt-rho2_aow_trt*der_d_AOW_trt)/D_AOW1
  G_beta_AOW_ctrl = (der_aw_AOW_ctrl-rho1_aow_ctrl*der_b_AOW_ctrl)/B_AOW0+(der_cw_AOW_ctrl-rho2_aow_ctrl*der_d_AOW_ctrl)/D_AOW0
  
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
  
  ###U_n\nabla_{\alpha}a_w^{\text{AOW}}(\widehat{\eta})
  ###U_n\nabla_{\alpha}c_w^{\text{AOW}}(\widehat{\eta})
  d_aw_cw_AOW = function(treatments, covariates_reg, pi_func, cond_prob_trt, cond_prob_ctrl) {
    # Calculate the first term of Cn
    cal_aw_cw_derivative = function(i){
      aw_vec = cw_vec = 0
      for(j in 1:n){
        if(j != i){
          e_Xi = pi_func[i]
          e_Xj = pi_func[j]
          mu_tilde_ij = mu_tilde_AOW_fun(i, j, covariates_reg, cond_prob_trt, cond_prob_ctrl)
          term1 = (treatments[i]*(1-treatments[j]))*(1-e_Xi)*e_Xj*(-mu_tilde_ij)
          aw_vec = aw_vec + term1
          term2 = e_Xi * (1-e_Xi) * e_Xj * (1-e_Xj)*(mu_tilde_ij)
          cw_vec = cw_vec + term2
        }
      }
      return(list(aw_vec = aw_vec, cw_vec = cw_vec))
    }
    der_aw_cw_aow_res = apply(matrix(1:n,nrow=1),2,cal_aw_cw_derivative)
    der_aw_aow = Reduce("+", lapply(der_aw_cw_aow_res, `[[`, "aw_vec"))
    der_cw_aow = Reduce("+", lapply(der_aw_cw_aow_res, `[[`, "cw_vec"))
    der_aw_alpha_AOW = der_aw_aow/(n*(n-1))
    der_cw_alpha_AOW = der_cw_aow/(n*(n-1))
    return(list(der_aw_alpha_AOW = der_aw_alpha_AOW,
                der_cw_alpha_AOW = der_cw_alpha_AOW))
  }
  
  d_aw_cw_AOW_res = d_aw_cw_AOW(treatments, covariates_reg, pi_func, cond_prob_trt, cond_prob_ctrl)
  der_aw_alpha_AOW = matrix(d_aw_cw_AOW_res$der_aw_alpha_AOW, ncol=1)
  der_cw_alpha_AOW = matrix(d_aw_cw_AOW_res$der_cw_alpha_AOW, ncol=1)
  
  G_alpha_AOW = der_aw_alpha_AOW/B_AOW1+der_cw_alpha_AOW/D_AOW1
  
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
  G_alpha_AOW = t(G_alpha_AOW)
  inf_aow_trt = inf_aow_ctrl = numeric(n)
  for(i in 1:n){
    l_theta = l_theta_AOW_fun(i,outcomes,treatments,covariates_reg,cond_prob_trt, cond_prob_ctrl, 
                              thresholds_trt, thresholds_ctrl, E_thetatheta_inv)
    l_beta = l_beta_fun(i)
    q_AOW_i_trt = WP_trt_aow+(aw1_AOW_trt[i]-rho1_aow_trt*b1_AOW_trt[i])/B_AOW1
                  +(cw1_AOW_trt[i]-rho2_aow_trt*d1_AOW_trt[i])/D_AOW1
    q_AOW_i_ctrl = WP_ctrl_aow+(aw1_AOW_ctrl[i]-rho1_aow_ctrl*b1_AOW_ctrl[i])/B_AOW0
                  +(cw1_AOW_ctrl[i]-rho2_aow_ctrl*d1_AOW_ctrl[i])/D_AOW0
    temp_trt = 2*(q_AOW_i_trt - WP_trt_aow) + G_beta_AOW_trt %*% l_beta + G_alpha_AOW %*% l_theta
    temp_ctrl = 2*(q_AOW_i_ctrl - WP_ctrl_aow) + G_beta_AOW_ctrl %*% l_beta + G_alpha_AOW %*% l_theta
    inf_aow_trt[i] = temp_trt
    inf_aow_ctrl[i] = temp_ctrl
  }
  m_inf_trt = mean(inf_aow_trt)
  m_inf_ctrl = mean(inf_aow_ctrl)
  #variance and covariance of WP and LP
  AOW_Var_trt = sum((inf_aow_trt-m_inf_trt)^2)/n
  AOW_Var_ctrl = sum((inf_aow_ctrl-m_inf_ctrl)^2)/n
  AOW_cov = sum((inf_aow_trt-m_inf_trt)*(inf_aow_ctrl-m_inf_ctrl))/n
  
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

