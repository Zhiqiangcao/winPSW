rm(list=ls())
library(xtable)
data1 = read.csv("C:/Users/82655/Dropbox/research/covariate adjustment win ratio/simulation results/ba_corr_r1.csv",header=T)
data2 = read.csv("C:/Users/82655/Dropbox/research/covariate adjustment win ratio/simulation results/unba_corr_r1.csv",header=T)

####bias calculation
#balanced design
unadj_200 = as.numeric(data1[1:4,4])-as.numeric(data1[1:4,3])
ipw_200 = as.numeric(data1[6:9,4])-as.numeric(data1[6:9,3])
ow_200 = as.numeric(data1[11:14,4])-as.numeric(data1[11:14,3])
aipw_200 = as.numeric(data1[16:19,4])-as.numeric(data1[16:19,3])
aow_200 = as.numeric(data1[21:24,4])-as.numeric(data1[21:24,3])
aipw_m_200 = as.numeric(data1[26:29,4])-as.numeric(data1[26:29,3])
aow_m_200 = as.numeric(data1[31:34,4])-as.numeric(data1[31:34,3])

unadj_300 = as.numeric(data1[36:39,4])-as.numeric(data1[36:39,3])
ipw_300 = as.numeric(data1[41:44,4])-as.numeric(data1[41:44,3])
ow_300 = as.numeric(data1[46:49,4])-as.numeric(data1[46:49,3])
aipw_300 = as.numeric(data1[51:54,4])-as.numeric(data1[51:54,3])
aow_300 = as.numeric(data1[56:59,4])-as.numeric(data1[56:59,3])
aipw_m_300 = as.numeric(data1[61:64,4])-as.numeric(data1[61:64,3])
aow_m_300 = as.numeric(data1[66:69,4])-as.numeric(data1[66:69,3])

unadj_400 = as.numeric(data1[71:74,4])-as.numeric(data1[71:74,3])
ipw_400 = as.numeric(data1[76:79,4])-as.numeric(data1[76:79,3])
ow_400 = as.numeric(data1[81:84,4])-as.numeric(data1[81:84,3])
aipw_400 = as.numeric(data1[86:89,4])-as.numeric(data1[86:89,3])
aow_400 = as.numeric(data1[91:94,4])-as.numeric(data1[91:94,3])
aipw_m_400 = as.numeric(data1[96:99,4])-as.numeric(data1[96:99,3])
aow_m_400 = as.numeric(data1[101:104,4])-as.numeric(data1[101:104,3])

bias_200 = cbind(unadj_200,ipw_200,ow_200,aipw_200,aow_200,aipw_m_200,aow_m_200)/as.numeric(data1[1:4,3])
bias_300 = cbind(unadj_300,ipw_300,ow_300,aipw_300,aow_300,aipw_m_300,aow_m_300)/as.numeric(data1[1:4,3])
bias_400 = cbind(unadj_400,ipw_400,ow_400,aipw_400,aow_400,aipw_m_400,aow_m_400)/as.numeric(data1[1:4,3])
ba_bias = rbind(bias_200,bias_300,bias_400)

#unbalanced design
unadj_200 = as.numeric(data2[1:4,4])-as.numeric(data2[1:4,3])
ipw_200 = as.numeric(data2[6:9,4])-as.numeric(data2[6:9,3])
ow_200 = as.numeric(data2[11:14,4])-as.numeric(data2[11:14,3])
aipw_200 = as.numeric(data2[16:19,4])-as.numeric(data2[16:19,3])
aow_200 = as.numeric(data2[21:24,4])-as.numeric(data2[21:24,3])
aipw_m_200 = as.numeric(data2[26:29,4])-as.numeric(data2[26:29,3])
aow_m_200 = as.numeric(data2[31:34,4])-as.numeric(data2[31:34,3])

unadj_300 = as.numeric(data2[36:39,4])-as.numeric(data2[36:39,3])
ipw_300 = as.numeric(data2[41:44,4])-as.numeric(data2[41:44,3])
ow_300 = as.numeric(data2[46:49,4])-as.numeric(data2[46:49,3])
aipw_300 = as.numeric(data2[51:54,4])-as.numeric(data2[51:54,3])
aow_300 = as.numeric(data2[56:59,4])-as.numeric(data2[56:59,3])
aipw_m_300 = as.numeric(data2[61:64,4])-as.numeric(data2[61:64,3])
aow_m_300 = as.numeric(data2[66:69,4])-as.numeric(data2[66:69,3])

unadj_400 = as.numeric(data2[71:74,4])-as.numeric(data2[71:74,3])
ipw_400 = as.numeric(data2[76:79,4])-as.numeric(data2[76:79,3])
ow_400 = as.numeric(data2[81:84,4])-as.numeric(data2[81:84,3])
aipw_400 = as.numeric(data2[86:89,4])-as.numeric(data2[86:89,3])
aow_400 = as.numeric(data2[91:94,4])-as.numeric(data2[91:94,3])
aipw_m_400 = as.numeric(data2[96:99,4])-as.numeric(data2[96:99,3])
aow_m_400 = as.numeric(data2[101:104,4])-as.numeric(data2[101:104,3])

bias_un_200 = cbind(unadj_200,ipw_200,ow_200,aipw_200,aow_200,aipw_m_200,aow_m_200)/as.numeric(data2[1:4,3])
bias_un_300 = cbind(unadj_300,ipw_300,ow_300,aipw_300,aow_300,aipw_m_300,aow_m_300)/as.numeric(data2[1:4,3])
bias_un_400 = cbind(unadj_400,ipw_400,ow_400,aipw_400,aow_400,aipw_m_400,aow_m_400)/as.numeric(data2[1:4,3])
unba_bias = rbind(bias_un_200,bias_un_300,bias_un_400)
#keep WR and WD
overall_bias = rbind(ba_bias[c(3:4,7:8,11:12),],unba_bias[c(3:4,7:8,11:12),])*100
xtable(overall_bias,digits = 2)

#balanced design
####relative efficiency
#wr
cal_ref = function(data1){
  ipw_eff_wr_200 = as.numeric(data1[3,5]^2)/as.numeric(data1[8,5]^2)
  ow_eff_wr_200 = as.numeric(data1[3,5]^2)/as.numeric(data1[13,5]^2)
  aipw_eff_wr_200 = as.numeric(data1[3,5]^2)/as.numeric(data1[18,5]^2)
  aow_eff_wr_200 = as.numeric(data1[3,5]^2)/as.numeric(data1[23,5]^2)
  aipw_m_eff_wr_200 = as.numeric(data1[3,5]^2)/as.numeric(data1[28,5]^2)
  aow_m_eff_wr_200 = as.numeric(data1[3,5]^2)/as.numeric(data1[33,5]^2)
  wr_eff_200 = cbind(ipw_eff_wr_200,ow_eff_wr_200,aipw_eff_wr_200,aow_eff_wr_200,aipw_m_eff_wr_200,aow_m_eff_wr_200)
  
  ipw_eff_wr_300 = as.numeric(data1[38,5]^2)/as.numeric(data1[43,5]^2)
  ow_eff_wr_300 = as.numeric(data1[38,5]^2)/as.numeric(data1[48,5]^2)
  aipw_eff_wr_300 = as.numeric(data1[38,5]^2)/as.numeric(data1[53,5]^2)
  aow_eff_wr_300 = as.numeric(data1[38,5]^2)/as.numeric(data1[58,5]^2)
  aipw_m_eff_wr_300 = as.numeric(data1[38,5]^2)/as.numeric(data1[63,5]^2)
  aow_m_eff_wr_300 = as.numeric(data1[38,5]^2)/as.numeric(data1[68,5]^2)
  wr_eff_300 = cbind(ipw_eff_wr_300,ow_eff_wr_300,aipw_eff_wr_300,aow_eff_wr_300,aipw_m_eff_wr_300,aow_m_eff_wr_300)
  
  
  ipw_eff_wr_400 = as.numeric(data1[73,5]^2)/as.numeric(data1[78,5]^2)
  ow_eff_wr_400 = as.numeric(data1[73,5]^2)/as.numeric(data1[83,5]^2)
  aipw_eff_wr_400 = as.numeric(data1[73,5]^2)/as.numeric(data1[88,5]^2)
  aow_eff_wr_400 = as.numeric(data1[73,5]^2)/as.numeric(data1[93,5]^2)
  aipw_m_eff_wr_400 = as.numeric(data1[73,5]^2)/as.numeric(data1[98,5]^2)
  aow_m_eff_wr_400 = as.numeric(data1[73,5]^2)/as.numeric(data1[103,5]^2)
  wr_eff_400 = cbind(ipw_eff_wr_400,ow_eff_wr_400,aipw_eff_wr_400,aow_eff_wr_400,aipw_m_eff_wr_400,aow_m_eff_wr_400)
  
  #wd
  ipw_eff_wd_200 = as.numeric(data1[4,5]^2)/as.numeric(data1[9,5]^2)
  ow_eff_wd_200 = as.numeric(data1[4,5]^2)/as.numeric(data1[14,5]^2)
  aipw_eff_wd_200 = as.numeric(data1[4,5]^2)/as.numeric(data1[19,5]^2)
  aow_eff_wd_200 = as.numeric(data1[4,5]^2)/as.numeric(data1[24,5]^2)
  aipw_m_eff_wd_200 = as.numeric(data1[4,5]^2)/as.numeric(data1[29,5]^2)
  aow_m_eff_wd_200 = as.numeric(data1[4,5]^2)/as.numeric(data1[34,5]^2)
  wd_eff_200 = cbind(ipw_eff_wd_200,ow_eff_wd_200,aipw_eff_wd_200,aow_eff_wd_200,aipw_m_eff_wd_200,aow_m_eff_wd_200)
  
  ipw_eff_wd_300 = as.numeric(data1[39,5]^2)/as.numeric(data1[44,5]^2)
  ow_eff_wd_300 = as.numeric(data1[39,5]^2)/as.numeric(data1[49,5]^2)
  aipw_eff_wd_300 = as.numeric(data1[39,5]^2)/as.numeric(data1[54,5]^2)
  aow_eff_wd_300 = as.numeric(data1[39,5]^2)/as.numeric(data1[59,5]^2)
  aipw_m_eff_wd_300 = as.numeric(data1[39,5]^2)/as.numeric(data1[64,5]^2)
  aow_m_eff_wd_300 = as.numeric(data1[39,5]^2)/as.numeric(data1[69,5]^2)
  wd_eff_300 = cbind(ipw_eff_wd_300,ow_eff_wd_300,aipw_eff_wd_300,aow_eff_wd_300,aipw_m_eff_wd_300,aow_m_eff_wd_300)
  
  
  ipw_eff_wd_400 = as.numeric(data1[74,5]^2)/as.numeric(data1[79,5]^2)
  ow_eff_wd_400 = as.numeric(data1[74,5]^2)/as.numeric(data1[84,5]^2)
  aipw_eff_wd_400 = as.numeric(data1[74,5]^2)/as.numeric(data1[89,5]^2)
  aow_eff_wd_400 = as.numeric(data1[74,5]^2)/as.numeric(data1[94,5]^2)
  aipw_m_eff_wd_400 = as.numeric(data1[74,5]^2)/as.numeric(data1[99,5]^2)
  aow_m_eff_wd_400 = as.numeric(data1[74,5]^2)/as.numeric(data1[104,5]^2)
  wd_eff_400 = cbind(ipw_eff_wd_400,ow_eff_wd_400,aipw_eff_wd_400,aow_eff_wd_400,aipw_m_eff_wd_400,aow_m_eff_wd_400)
  
  ipw_wr_eff = c(wr_eff_200[1],wr_eff_300[1],wr_eff_400[1])
  ow_wr_eff = c(wr_eff_200[2],wr_eff_300[2],wr_eff_400[2])
  aipw_wr_eff = c(wr_eff_200[3],wr_eff_300[3],wr_eff_400[3])
  aow_wr_eff = c(wr_eff_200[4],wr_eff_300[4],wr_eff_400[4])
  aipw_mis_wr_eff = c(wr_eff_200[5],wr_eff_300[5],wr_eff_400[5])
  aow_mis_wr_eff = c(wr_eff_200[6],wr_eff_300[6],wr_eff_400[6])
  wr_eff = data.frame(ipw_wr_eff,ow_wr_eff,aipw_wr_eff,aow_wr_eff,aipw_mis_wr_eff,aow_mis_wr_eff)
  
  ipw_wd_eff = c(wd_eff_200[1],wd_eff_300[1],wd_eff_400[1])
  ow_wd_eff = c(wd_eff_200[2],wd_eff_300[2],wd_eff_400[2])
  aipw_wd_eff = c(wd_eff_200[3],wd_eff_300[3],wd_eff_400[3])
  aow_wd_eff = c(wd_eff_200[4],wd_eff_300[4],wd_eff_400[4])
  aipw_mis_wd_eff = c(wd_eff_200[5],wd_eff_300[5],wd_eff_400[5])
  aow_mis_wd_eff = c(wd_eff_200[6],wd_eff_300[6],wd_eff_400[6])
  wd_eff = data.frame(ipw_wd_eff,ow_wd_eff,aipw_wd_eff,aow_wd_eff,aipw_mis_wd_eff,aow_mis_wd_eff)
  final_res = list(wr_eff,wd_eff)
  return(final_res)
}

quad_ba_res = cal_ref(data1)
quad_unba_res = cal_ref(data2)

#####est_var/MC Var
#wr
unadj_emc_wr_200 = as.numeric(data1[3,6]^2)/as.numeric(data1[3,5]^2)
ipw_emc_wr_200 = as.numeric(data1[8,6]^2)/as.numeric(data1[8,5]^2)
ow_emc_wr_200 = as.numeric(data1[13,6]^2)/as.numeric(data1[13,5]^2)
aipw_emc_wr_200 = as.numeric(data1[18,6]^2)/as.numeric(data1[18,5]^2)
aow_emc_wr_200 = as.numeric(data1[23,6]^2)/as.numeric(data1[23,5]^2)
aipw_m_emc_wr_200 = as.numeric(data1[28,6]^2)/as.numeric(data1[28,5]^2)
aow_m_emc_wr_200 = as.numeric(data1[33,6]^2)/as.numeric(data1[33,5]^2)
wr_emc_200 = cbind(unadj_emc_wr_200,ipw_emc_wr_200,ow_emc_wr_200,aipw_emc_wr_200,aow_emc_wr_200,aipw_m_emc_wr_200,aow_m_emc_wr_200)

unadj_emc_wr_300 = as.numeric(data1[38,6]^2)/as.numeric(data1[38,5]^2)
ipw_emc_wr_300 = as.numeric(data1[43,6]^2)/as.numeric(data1[43,5]^2)
ow_emc_wr_300 = as.numeric(data1[48,6]^2)/as.numeric(data1[48,5]^2)
aipw_emc_wr_300 = as.numeric(data1[53,6]^2)/as.numeric(data1[53,5]^2)
aow_emc_wr_300 = as.numeric(data1[58,6]^2)/as.numeric(data1[58,5]^2)
aipw_m_emc_wr_300 = as.numeric(data1[63,6]^2)/as.numeric(data1[63,5]^2)
aow_m_emc_wr_300 = as.numeric(data1[68,6]^2)/as.numeric(data1[68,5]^2)
wr_emc_300 = cbind(unadj_emc_wr_300,ipw_emc_wr_300,ow_emc_wr_300,aipw_emc_wr_300,aow_emc_wr_300,aipw_m_emc_wr_300,aow_m_emc_wr_300)

unadj_emc_wr_400 = as.numeric(data1[73,6]^2)/as.numeric(data1[73,5]^2)
ipw_emc_wr_400 = as.numeric(data1[78,6]^2)/as.numeric(data1[78,5]^2)
ow_emc_wr_400 = as.numeric(data1[83,6]^2)/as.numeric(data1[83,5]^2)
aipw_emc_wr_400 = as.numeric(data1[88,6]^2)/as.numeric(data1[88,5]^2)
aow_emc_wr_400 = as.numeric(data1[93,6]^2)/as.numeric(data1[93,5]^2)
aipw_m_emc_wr_400 = as.numeric(data1[98,6]^2)/as.numeric(data1[98,5]^2)
aow_m_emc_wr_400 = as.numeric(data1[103,6]^2)/as.numeric(data1[103,5]^2)
wr_emc_400 = cbind(unadj_emc_wr_400,ipw_emc_wr_400,ow_emc_wr_400,aipw_emc_wr_400,aow_emc_wr_400,aipw_m_emc_wr_400,aow_m_emc_wr_400)

#wd
unadj_emc_wd_200 = as.numeric(data1[4,6]^2)/as.numeric(data1[4,5]^2)
ipw_emc_wd_200 = as.numeric(data1[9,6]^2)/as.numeric(data1[9,5]^2)
ow_emc_wd_200 = as.numeric(data1[14,6]^2)/as.numeric(data1[14,5]^2)
aipw_emc_wd_200 = as.numeric(data1[19,6]^2)/as.numeric(data1[19,5]^2)
aow_emc_wd_200 = as.numeric(data1[24,6]^2)/as.numeric(data1[24,5]^2)
aipw_m_emc_wd_200 = as.numeric(data1[29,6]^2)/as.numeric(data1[29,5]^2)
aow_m_emc_wd_200 = as.numeric(data1[34,6]^2)/as.numeric(data1[34,5]^2)
wd_emc_200 = cbind(unadj_emc_wd_200,ipw_emc_wd_200,ow_emc_wd_200,aipw_emc_wd_200,aow_emc_wd_200,aipw_m_emc_wd_200,aow_m_emc_wd_200)

unadj_emc_wd_300 = as.numeric(data1[39,6]^2)/as.numeric(data1[39,5]^2)
ipw_emc_wd_300 = as.numeric(data1[44,6]^2)/as.numeric(data1[44,5]^2)
ow_emc_wd_300 = as.numeric(data1[49,6]^2)/as.numeric(data1[49,5]^2)
aipw_emc_wd_300 = as.numeric(data1[54,6]^2)/as.numeric(data1[54,5]^2)
aow_emc_wd_300 = as.numeric(data1[59,6]^2)/as.numeric(data1[59,5]^2)
aipw_m_emc_wd_300 = as.numeric(data1[64,6]^2)/as.numeric(data1[64,5]^2)
aow_m_emc_wd_300 = as.numeric(data1[69,6]^2)/as.numeric(data1[69,5]^2)
wd_emc_300 = cbind(unadj_emc_wd_300,ipw_emc_wd_300,ow_emc_wd_300,aipw_emc_wd_300,aow_emc_wd_300,aipw_m_emc_wd_300,aow_m_emc_wd_300)

unadj_emc_wd_400 = as.numeric(data1[74,6]^2)/as.numeric(data1[74,5]^2)
ipw_emc_wd_400 = as.numeric(data1[79,6]^2)/as.numeric(data1[79,5]^2)
ow_emc_wd_400 = as.numeric(data1[84,6]^2)/as.numeric(data1[84,5]^2)
aipw_emc_wd_400 = as.numeric(data1[89,6]^2)/as.numeric(data1[89,5]^2)
aow_emc_wd_400 = as.numeric(data1[94,6]^2)/as.numeric(data1[94,5]^2)
aipw_m_emc_wd_400 = as.numeric(data1[99,6]^2)/as.numeric(data1[99,5]^2)
aow_m_emc_wd_400 = as.numeric(data1[104,6]^2)/as.numeric(data1[104,5]^2)
wd_emc_400 = cbind(unadj_emc_wd_400,ipw_emc_wd_400,ow_emc_wd_400,aipw_emc_wd_400,aow_emc_wd_400,aipw_m_emc_wd_400,aow_m_emc_wd_400)

######coerage rate
#wr
unadj_cr_wr_200 = as.numeric(data1[3,7])
ipw_cr_wr_200 = as.numeric(data1[8,7])
ow_cr_wr_200 = as.numeric(data1[13,7])
aipw_cr_wr_200 = as.numeric(data1[18,7])
aow_cr_wr_200 = as.numeric(data1[23,7])
aipw_m_cr_wr_200 = as.numeric(data1[28,7])
aow_m_cr_wr_200 = as.numeric(data1[33,7])
wr_cr_200 = cbind(unadj_cr_wr_200,ipw_cr_wr_200,ow_cr_wr_200,aipw_cr_wr_200,aow_cr_wr_200,aipw_m_cr_wr_200,aow_m_cr_wr_200)

unadj_cr_wr_300 = as.numeric(data1[38,7])
ipw_cr_wr_300 = as.numeric(data1[43,7])
ow_cr_wr_300 = as.numeric(data1[48,7])
aipw_cr_wr_300 = as.numeric(data1[53,7])
aow_cr_wr_300 = as.numeric(data1[58,7])
aipw_m_cr_wr_300 = as.numeric(data1[63,7])
aow_m_cr_wr_300 = as.numeric(data1[68,7])
wr_cr_300 = cbind(unadj_cr_wr_300,ipw_cr_wr_300,ow_cr_wr_300,aipw_cr_wr_300,aow_cr_wr_300,aipw_m_cr_wr_300,aow_m_cr_wr_300)

unadj_cr_wr_400 = as.numeric(data1[73,7])
ipw_cr_wr_400 = as.numeric(data1[78,7])
ow_cr_wr_400 = as.numeric(data1[83,7])
aipw_cr_wr_400 = as.numeric(data1[88,7])
aow_cr_wr_400 = as.numeric(data1[93,7])
aipw_m_cr_wr_400 = as.numeric(data1[98,7])
aow_m_cr_wr_400 = as.numeric(data1[103,7])
wr_cr_400 = cbind(unadj_cr_wr_400,ipw_cr_wr_400,ow_cr_wr_400,aipw_cr_wr_400,aow_cr_wr_400,aipw_m_cr_wr_400,aow_m_cr_wr_400)

#wd
unadj_cr_wd_200 = as.numeric(data1[4,7])
ipw_cr_wd_200 = as.numeric(data1[9,7])
ow_cr_wd_200 = as.numeric(data1[14,7])
aipw_cr_wd_200 = as.numeric(data1[19,7])
aow_cr_wd_200 = as.numeric(data1[24,7])
aipw_m_cr_wd_200 = as.numeric(data1[29,7])
aow_m_cr_wd_200 = as.numeric(data1[34,7])
wd_cr_200 = cbind(unadj_cr_wd_200,ipw_cr_wd_200,ow_cr_wd_200,aipw_cr_wd_200,aow_cr_wd_200,aipw_m_cr_wd_200,aow_m_cr_wd_200)

unadj_cr_wd_300 = as.numeric(data1[39,7])
ipw_cr_wd_300 = as.numeric(data1[44,7])
ow_cr_wd_300 = as.numeric(data1[49,7])
aipw_cr_wd_300 = as.numeric(data1[54,7])
aow_cr_wd_300 = as.numeric(data1[59,7])
aipw_m_cr_wd_300 = as.numeric(data1[64,7])
aow_m_cr_wd_300 = as.numeric(data1[69,7])
wd_cr_300 = cbind(unadj_cr_wd_300,ipw_cr_wd_300,ow_cr_wd_300,aipw_cr_wd_300,aow_cr_wd_300,aipw_m_cr_wd_300,aow_m_cr_wd_300)

unadj_cr_wd_400 = as.numeric(data1[74,7])
ipw_cr_wd_400 = as.numeric(data1[79,7])
ow_cr_wd_400 = as.numeric(data1[84,7])
aipw_cr_wd_400 = as.numeric(data1[89,7])
aow_cr_wd_400 = as.numeric(data1[94,7])
aipw_m_cr_wd_400 = as.numeric(data1[99,7])
aow_m_cr_wd_400 = as.numeric(data1[104,7])
wd_cr_400 = cbind(unadj_cr_wd_400,ipw_cr_wd_400,ow_cr_wd_400,aipw_cr_wd_400,aow_cr_wd_400,aipw_m_cr_wd_400,aow_m_cr_wd_400)

#balanced results (do not include effiency)
ba_wr = rbind(cbind(wr_emc_200,wr_cr_200),
              cbind(wr_emc_300,wr_cr_300),
              cbind(wr_emc_400,wr_cr_400))
ba_wd = rbind(cbind(wd_emc_200,wd_cr_200),
              cbind(wd_emc_300,wd_cr_300),
              cbind(wd_emc_400,wd_cr_400))
ba_res = rbind(ba_wr,ba_wd)


############unbalanced design
#####est_var/MC Var
#wr
unadj_emc_wr_200 = as.numeric(data2[3,6]^2)/as.numeric(data2[3,5]^2)
ipw_emc_wr_200 = as.numeric(data2[8,6]^2)/as.numeric(data2[8,5]^2)
ow_emc_wr_200 = as.numeric(data2[13,6]^2)/as.numeric(data2[13,5]^2)
aipw_emc_wr_200 = as.numeric(data2[18,6]^2)/as.numeric(data2[18,5]^2)
aow_emc_wr_200 = as.numeric(data2[23,6]^2)/as.numeric(data2[23,5]^2)
aipw_m_emc_wr_200 = as.numeric(data2[28,6]^2)/as.numeric(data2[28,5]^2)
aow_m_emc_wr_200 = as.numeric(data2[33,6]^2)/as.numeric(data2[33,5]^2)
wr_emc_200 = cbind(unadj_emc_wr_200,ipw_emc_wr_200,ow_emc_wr_200,aipw_emc_wr_200,aow_emc_wr_200,aipw_m_emc_wr_200,aow_m_emc_wr_200)

unadj_emc_wr_300 = as.numeric(data2[38,6]^2)/as.numeric(data2[38,5]^2)
ipw_emc_wr_300 = as.numeric(data2[43,6]^2)/as.numeric(data2[43,5]^2)
ow_emc_wr_300 = as.numeric(data2[48,6]^2)/as.numeric(data2[48,5]^2)
aipw_emc_wr_300 = as.numeric(data2[53,6]^2)/as.numeric(data2[53,5]^2)
aow_emc_wr_300 = as.numeric(data2[58,6]^2)/as.numeric(data2[58,5]^2)
aipw_m_emc_wr_300 = as.numeric(data2[63,6]^2)/as.numeric(data2[63,5]^2)
aow_m_emc_wr_300 = as.numeric(data2[68,6]^2)/as.numeric(data2[68,5]^2)
wr_emc_300 = cbind(unadj_emc_wr_300,ipw_emc_wr_300,ow_emc_wr_300,aipw_emc_wr_300,aow_emc_wr_300,aipw_m_emc_wr_300,aow_m_emc_wr_300)

unadj_emc_wr_400 = as.numeric(data2[73,6]^2)/as.numeric(data2[73,5]^2)
ipw_emc_wr_400 = as.numeric(data2[78,6]^2)/as.numeric(data2[78,5]^2)
ow_emc_wr_400 = as.numeric(data2[83,6]^2)/as.numeric(data2[83,5]^2)
aipw_emc_wr_400 = as.numeric(data2[88,6]^2)/as.numeric(data2[88,5]^2)
aow_emc_wr_400 = as.numeric(data2[93,6]^2)/as.numeric(data2[93,5]^2)
aipw_m_emc_wr_400 = as.numeric(data2[98,6]^2)/as.numeric(data2[98,5]^2)
aow_m_emc_wr_400 = as.numeric(data2[103,6]^2)/as.numeric(data2[103,5]^2)
wr_emc_400 = cbind(unadj_emc_wr_400,ipw_emc_wr_400,ow_emc_wr_400,aipw_emc_wr_400,aow_emc_wr_400,aipw_m_emc_wr_400,aow_m_emc_wr_400)

#wd
unadj_emc_wd_200 = as.numeric(data2[4,6]^2)/as.numeric(data2[4,5]^2)
ipw_emc_wd_200 = as.numeric(data2[9,6]^2)/as.numeric(data2[9,5]^2)
ow_emc_wd_200 = as.numeric(data2[14,6]^2)/as.numeric(data2[14,5]^2)
aipw_emc_wd_200 = as.numeric(data2[19,6]^2)/as.numeric(data2[19,5]^2)
aow_emc_wd_200 = as.numeric(data2[24,6]^2)/as.numeric(data2[24,5]^2)
aipw_m_emc_wd_200 = as.numeric(data2[29,6]^2)/as.numeric(data2[29,5]^2)
aow_m_emc_wd_200 = as.numeric(data2[34,6]^2)/as.numeric(data2[34,5]^2)
wd_emc_200 = cbind(unadj_emc_wd_200,ipw_emc_wd_200,ow_emc_wd_200,aipw_emc_wd_200,aow_emc_wd_200,aipw_m_emc_wd_200,aow_m_emc_wd_200)

unadj_emc_wd_300 = as.numeric(data2[39,6]^2)/as.numeric(data2[39,5]^2)
ipw_emc_wd_300 = as.numeric(data2[44,6]^2)/as.numeric(data2[44,5]^2)
ow_emc_wd_300 = as.numeric(data2[49,6]^2)/as.numeric(data2[49,5]^2)
aipw_emc_wd_300 = as.numeric(data2[54,6]^2)/as.numeric(data2[54,5]^2)
aow_emc_wd_300 = as.numeric(data2[59,6]^2)/as.numeric(data2[59,5]^2)
aipw_m_emc_wd_300 = as.numeric(data2[64,6]^2)/as.numeric(data2[64,5]^2)
aow_m_emc_wd_300 = as.numeric(data2[69,6]^2)/as.numeric(data2[69,5]^2)
wd_emc_300 = cbind(unadj_emc_wd_300,ipw_emc_wd_300,ow_emc_wd_300,aipw_emc_wd_300,aow_emc_wd_300,aipw_m_emc_wd_300,aow_m_emc_wd_300)

unadj_emc_wd_400 = as.numeric(data2[74,6]^2)/as.numeric(data2[74,5]^2)
ipw_emc_wd_400 = as.numeric(data2[79,6]^2)/as.numeric(data2[79,5]^2)
ow_emc_wd_400 = as.numeric(data2[84,6]^2)/as.numeric(data2[84,5]^2)
aipw_emc_wd_400 = as.numeric(data2[89,6]^2)/as.numeric(data2[89,5]^2)
aow_emc_wd_400 = as.numeric(data2[94,6]^2)/as.numeric(data2[94,5]^2)
aipw_m_emc_wd_400 = as.numeric(data2[99,6]^2)/as.numeric(data2[99,5]^2)
aow_m_emc_wd_400 = as.numeric(data2[104,6]^2)/as.numeric(data2[104,5]^2)
wd_emc_400 = cbind(unadj_emc_wd_400,ipw_emc_wd_400,ow_emc_wd_400,aipw_emc_wd_400,aow_emc_wd_400,aipw_m_emc_wd_400,aow_m_emc_wd_400)

######coerage rate
#wr
unadj_cr_wr_200 = as.numeric(data2[3,7])
ipw_cr_wr_200 = as.numeric(data2[8,7])
ow_cr_wr_200 = as.numeric(data2[13,7])
aipw_cr_wr_200 = as.numeric(data2[18,7])
aow_cr_wr_200 = as.numeric(data2[23,7])
aipw_m_cr_wr_200 = as.numeric(data2[28,7])
aow_m_cr_wr_200 = as.numeric(data2[33,7])
wr_cr_200 = cbind(unadj_cr_wr_200,ipw_cr_wr_200,ow_cr_wr_200,aipw_cr_wr_200,aow_cr_wr_200,aipw_m_cr_wr_200,aow_m_cr_wr_200)

unadj_cr_wr_300 = as.numeric(data2[38,7])
ipw_cr_wr_300 = as.numeric(data2[43,7])
ow_cr_wr_300 = as.numeric(data2[48,7])
aipw_cr_wr_300 = as.numeric(data2[53,7])
aow_cr_wr_300 = as.numeric(data2[58,7])
aipw_m_cr_wr_300 = as.numeric(data2[63,7])
aow_m_cr_wr_300 = as.numeric(data2[68,7])
wr_cr_300 = cbind(unadj_cr_wr_300,ipw_cr_wr_300,ow_cr_wr_300,aipw_cr_wr_300,aow_cr_wr_300,aipw_m_cr_wr_300,aow_m_cr_wr_300)

unadj_cr_wr_400 = as.numeric(data2[73,7])
ipw_cr_wr_400 = as.numeric(data2[78,7])
ow_cr_wr_400 = as.numeric(data2[83,7])
aipw_cr_wr_400 = as.numeric(data2[88,7])
aow_cr_wr_400 = as.numeric(data2[93,7])
aipw_m_cr_wr_400 = as.numeric(data2[98,7])
aow_m_cr_wr_400 = as.numeric(data2[103,7])
wr_cr_400 = cbind(unadj_cr_wr_400,ipw_cr_wr_400,ow_cr_wr_400,aipw_cr_wr_400,aow_cr_wr_400,aipw_m_cr_wr_400,aow_m_cr_wr_400)

#wd
unadj_cr_wd_200 = as.numeric(data2[4,7])
ipw_cr_wd_200 = as.numeric(data2[9,7])
ow_cr_wd_200 = as.numeric(data2[14,7])
aipw_cr_wd_200 = as.numeric(data2[19,7])
aow_cr_wd_200 = as.numeric(data2[24,7])
aipw_m_cr_wd_200 = as.numeric(data2[29,7])
aow_m_cr_wd_200 = as.numeric(data2[34,7])
wd_cr_200 = cbind(unadj_cr_wd_200,ipw_cr_wd_200,ow_cr_wd_200,aipw_cr_wd_200,aow_cr_wd_200,aipw_m_cr_wd_200,aow_m_cr_wd_200)

unadj_cr_wd_300 = as.numeric(data2[39,7])
ipw_cr_wd_300 = as.numeric(data2[44,7])
ow_cr_wd_300 = as.numeric(data2[49,7])
aipw_cr_wd_300 = as.numeric(data2[54,7])
aow_cr_wd_300 = as.numeric(data2[59,7])
aipw_m_cr_wd_300 = as.numeric(data2[64,7])
aow_m_cr_wd_300 = as.numeric(data2[69,7])
wd_cr_300 = cbind(unadj_cr_wd_300,ipw_cr_wd_300,ow_cr_wd_300,aipw_cr_wd_300,aow_cr_wd_300,aipw_m_cr_wd_300,aow_m_cr_wd_300)

unadj_cr_wd_400 = as.numeric(data2[74,7])
ipw_cr_wd_400 = as.numeric(data2[79,7])
ow_cr_wd_400 = as.numeric(data2[84,7])
aipw_cr_wd_400 = as.numeric(data2[89,7])
aow_cr_wd_400 = as.numeric(data2[94,7])
aipw_m_cr_wd_400 = as.numeric(data2[99,7])
aow_m_cr_wd_400 = as.numeric(data2[104,7])
wd_cr_400 = cbind(unadj_cr_wd_400,ipw_cr_wd_400,ow_cr_wd_400,aipw_cr_wd_400,aow_cr_wd_400,aipw_m_cr_wd_400,aow_m_cr_wd_400)

#balanced results (do not include effiency)
unba_wr = rbind(cbind(wr_emc_200,wr_cr_200),
              cbind(wr_emc_300,wr_cr_300),
              cbind(wr_emc_400,wr_cr_400))
unba_wd = rbind(cbind(wd_emc_200,wd_cr_200),
              cbind(wd_emc_300,wd_cr_300),
              cbind(wd_emc_400,wd_cr_400))
unba_res = rbind(unba_wr,unba_wd)

overall_res = rbind(ba_res,unba_res)
xtable(overall_res,digits = 3)
