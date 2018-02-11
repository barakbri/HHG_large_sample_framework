#This is a power simulation for the non parametric tests of independence, for small sample sizes
setwd('~')
library(HHG)

MODE_RUN_SMALL_N_SIMULATION = FALSE
MODE_PLOT_SMALL_N_RESULTS = FALSE

generate_linear_data = function(N,slope){
  ret = list()
  ret$X = rnorm(N,0,1)
  ret$Y = rnorm(N,0,1) +  slope * ret$X
  return(ret)
}

N_vec = seq(10,40,10)

NULL_TABLE_SIZE = 2000
POWER_REPS = 2000

SmallN_res_list = list()
SmallN_res_list$Power_result_slope_1 = rep(0,length(N_vec))
SmallN_res_list$Power_result_Spearman_slope_1 = rep(0,length(N_vec))
SmallN_res_list$Power_result_slope_0_5 = rep(0,length(N_vec))
SmallN_res_list$Power_result_Spearman_slope_0_5 = rep(0,length(N_vec))
SmallN_res_list$Power_result_slope_0_25 = rep(0,length(N_vec))
SmallN_res_list$Power_result_Spearman_slope_0_25 = rep(0,length(N_vec))
SmallN_res_list$T1E_result = rep(0,length(N_vec))
VERBOSE = T
alpha = 0.05

set.seed(1)
for(i in 1:length(N_vec)){
  #generate Null table 
  current_N = N_vec[i] 
  if(VERBOSE)
    print(paste0('Doing N: ',current_N))
  nt = Fast.independence.test.nulltable(current_N, nr.perm = NULL_TABLE_SIZE)
  
  
  #with slope 1
  for(b in 1:POWER_REPS){
    if(VERBOSE & b %% (ceiling(POWER_REPS/10)) == 1)
      print(paste0('Power rep: ',b))
    data = generate_linear_data(current_N,slope = 1)
    res = Fast.independence.test(data$X,data$Y,nt)
    SmallN_res_list$Power_result_slope_1[i] = SmallN_res_list$Power_result_slope_1[i] + (1 * (res$MinP.pvalue<=alpha)) / POWER_REPS
    SmallN_res_list$Power_result_Spearman_slope_1[i] = SmallN_res_list$Power_result_Spearman_slope_1[i] + 1 * (cor.test(data$X,data$Y,method = 'spearman')$p.value<=alpha) / POWER_REPS
  }
  
  #with slope 0.5
  for(b in 1:POWER_REPS){
    if(VERBOSE & b %% (ceiling(POWER_REPS/10)) == 1)
      print(paste0('Power rep: ',b))
    data = generate_linear_data(current_N,slope = 0.5)
    res = Fast.independence.test(data$X,data$Y,nt)
    SmallN_res_list$Power_result_slope_0_5[i] = SmallN_res_list$Power_result_slope_0_5[i] + (1 * (res$MinP.pvalue<=alpha)) / POWER_REPS
    SmallN_res_list$Power_result_Spearman_slope_0_5[i] = SmallN_res_list$Power_result_Spearman_slope_0_5[i] + 1 * (cor.test(data$X,data$Y,method = 'spearman')$p.value<=alpha) / POWER_REPS
  }
  
  #with slope 0.25
  for(b in 1:POWER_REPS){
    if(VERBOSE & b %% (ceiling(POWER_REPS/10)) == 1)
      print(paste0('Power rep: ',b))
    data = generate_linear_data(current_N,slope = 0.25)
    res = Fast.independence.test(data$X,data$Y,nt)
    SmallN_res_list$Power_result_slope_0_25[i] = SmallN_res_list$Power_result_slope_0_25[i] + (1 * (res$MinP.pvalue<=alpha)) / POWER_REPS
    SmallN_res_list$Power_result_Spearman_slope_0_25[i] = SmallN_res_list$Power_result_Spearman_slope_0_25[i] + 1 * (cor.test(data$X,data$Y,method = 'spearman')$p.value<=alpha) / POWER_REPS
  }
  
  #slope 0
  for(b in 1:POWER_REPS){
    if(VERBOSE & b %% (ceiling(POWER_REPS/10)) == 1)
      print(paste0('Power rep: ',b))
    data = generate_linear_data(current_N,slope = 0)
    res = Fast.independence.test(data$X,data$Y,nt)
    SmallN_res_list$T1E_result[i] = SmallN_res_list$T1E_result[i] + (1 * (res$MinP.pvalue<=alpha)) / POWER_REPS
  }
  
}


save(SmallN_res_list,file = 'SmallN_res.RData')

pdf('SmallN_Ind_Results.pdf', width = 6, height = 6)
plot(N_vec,SmallN_res_list$Power_result_slope_0_5,ylim = c(0,1),col = 'blue',lty = 1,type = 'l', xlab = 'N (Sample size)',ylab = 'Power',lwd = 2)
lines(N_vec,SmallN_res_list$Power_result_Spearman_slope_0_5,ylim = c(0,1),col = 'blue',lty = 2,type = 'l',lwd = 2)
lines(N_vec,SmallN_res_list$Power_result_slope_0_25,ylim = c(0,1),col = 'red',lty = 1,type = 'l',lwd = 2)
lines(N_vec,SmallN_res_list$Power_result_Spearman_slope_0_25,ylim = c(0,1),col = 'red',lty = 2,type = 'l',lwd = 2)
lines(N_vec,SmallN_res_list$Power_result_slope_1,ylim = c(0,1),col = 'green',lty = 1,type = 'l',lwd = 2)
lines(N_vec,SmallN_res_list$Power_result_Spearman_slope_1,ylim = c(0,1),col = 'green',lty = 2,type = 'l',lwd = 2)
lines(N_vec,SmallN_res_list$T1E_result,ylim = c(0,1),col = 'black', lty = 1, type = 'l',lwd = 2)
abline(h = 0.05,col = 'grey',lty = 2,lwd = 2)
dev.off()

