#This is a power simulation for the non parametric tests of independence, for small sample sizes
# See brill-heller-heller.R for an index used for running all scripts.

#setwd('~') #set in main script, disabled here

library(HHG)

# This is a function used for generating a dataset with X & Y from a linear setting, with a required slope.
generate_linear_data = function(N,slope){
  ret = list()
  ret$X = rnorm(N,0,1)
  ret$Y = rnorm(N,0,1) +  slope * ret$X
  return(ret)
}

# This is the vector of sample sizes considered for the power study: 10,20,30,40
N_vec = seq(10,40,10)

NULL_TABLE_SIZE = 2000 # we produce null tables of size 2000
POWER_REPS = 2000 # power is estimated over 2000 data generations

# this is used to store the results for small N
# vectors are for a test and slope, and store all results by N
SmallN_res_list = list()
SmallN_res_list$Power_result_slope_1 = rep(0,length(N_vec)) #ADP results with slope 1
SmallN_res_list$Power_result_Spearman_slope_1 = rep(0,length(N_vec)) # spearman, slope 1
SmallN_res_list$Power_result_slope_0_5 = rep(0,length(N_vec)) # ADP, slope 0.5
SmallN_res_list$Power_result_Spearman_slope_0_5 = rep(0,length(N_vec)) #spearman, slope 0.5
SmallN_res_list$Power_result_slope_0_25 = rep(0,length(N_vec)) # ADP, slope 0.25
SmallN_res_list$Power_result_Spearman_slope_0_25 = rep(0,length(N_vec)) #spearma, slope 0.25
SmallN_res_list$T1E_result = rep(0,length(N_vec)) # Type I error for ADP, over all sample sizes

VERBOSE = T #should messages be printed
alpha = 0.05 #Required Type I error rate

set.seed(1)
#iterate over sample sizes
for(i in 1:length(N_vec)){
  #generate Null table 
  current_N = N_vec[i] 
  if(VERBOSE)
    print(paste0('Doing N: ',current_N))
  nt = Fast.independence.test.nulltable(current_N, nr.perm = NULL_TABLE_SIZE) # null table will be generated with the default parameters, m.max = 10, N.atoms = N (as N<=40)
  
  
  #with slope 1 - generate data, test, store results
  for(b in 1:POWER_REPS){
    if(VERBOSE & b %% (ceiling(POWER_REPS/10)) == 1)
      print(paste0('Power rep: ',b))
    data = generate_linear_data(current_N,slope = 1)
    res = Fast.independence.test(data$X,data$Y,nt)
    SmallN_res_list$Power_result_slope_1[i] = SmallN_res_list$Power_result_slope_1[i] + (1 * (res$MinP.pvalue<=alpha)) / POWER_REPS
    SmallN_res_list$Power_result_Spearman_slope_1[i] = SmallN_res_list$Power_result_Spearman_slope_1[i] + 1 * (cor.test(data$X,data$Y,method = 'spearman')$p.value<=alpha) / POWER_REPS
  }
  
  #with slope 0.5 - generate data, test, store results
  for(b in 1:POWER_REPS){
    if(VERBOSE & b %% (ceiling(POWER_REPS/10)) == 1)
      print(paste0('Power rep: ',b))
    data = generate_linear_data(current_N,slope = 0.5)
    res = Fast.independence.test(data$X,data$Y,nt)
    SmallN_res_list$Power_result_slope_0_5[i] = SmallN_res_list$Power_result_slope_0_5[i] + (1 * (res$MinP.pvalue<=alpha)) / POWER_REPS
    SmallN_res_list$Power_result_Spearman_slope_0_5[i] = SmallN_res_list$Power_result_Spearman_slope_0_5[i] + 1 * (cor.test(data$X,data$Y,method = 'spearman')$p.value<=alpha) / POWER_REPS
  }
  
  #with slope 0.25 - generate data, test, store results
  for(b in 1:POWER_REPS){
    if(VERBOSE & b %% (ceiling(POWER_REPS/10)) == 1)
      print(paste0('Power rep: ',b))
    data = generate_linear_data(current_N,slope = 0.25)
    res = Fast.independence.test(data$X,data$Y,nt)
    SmallN_res_list$Power_result_slope_0_25[i] = SmallN_res_list$Power_result_slope_0_25[i] + (1 * (res$MinP.pvalue<=alpha)) / POWER_REPS
    SmallN_res_list$Power_result_Spearman_slope_0_25[i] = SmallN_res_list$Power_result_Spearman_slope_0_25[i] + 1 * (cor.test(data$X,data$Y,method = 'spearman')$p.value<=alpha) / POWER_REPS
  }
  
  #slope 0 - generate data, test, store results
  for(b in 1:POWER_REPS){
    if(VERBOSE & b %% (ceiling(POWER_REPS/10)) == 1)
      print(paste0('Power rep: ',b))
    data = generate_linear_data(current_N,slope = 0)
    res = Fast.independence.test(data$X,data$Y,nt)
    SmallN_res_list$T1E_result[i] = SmallN_res_list$T1E_result[i] + (1 * (res$MinP.pvalue<=alpha)) / POWER_REPS
  }
  
}

#save the results to file
save(SmallN_res_list,file = 'SmallN_res.RData')

#plot to file
pdf('SmallN_Ind_Results.pdf', width = 6, height = 6)
#plot slope 0.5
plot(N_vec,SmallN_res_list$Power_result_slope_0_5,ylim = c(0,1),col = 'blue',lty = 1,type = 'l', xlab = 'N (Sample size)',ylab = 'Power',lwd = 2)
lines(N_vec,SmallN_res_list$Power_result_Spearman_slope_0_5,ylim = c(0,1),col = 'blue',lty = 2,type = 'l',lwd = 2)
#plot slope 0.25
lines(N_vec,SmallN_res_list$Power_result_slope_0_25,ylim = c(0,1),col = 'red',lty = 1,type = 'l',lwd = 2)
lines(N_vec,SmallN_res_list$Power_result_Spearman_slope_0_25,ylim = c(0,1),col = 'red',lty = 2,type = 'l',lwd = 2)
#plot slope 1
lines(N_vec,SmallN_res_list$Power_result_slope_1,ylim = c(0,1),col = 'green',lty = 1,type = 'l',lwd = 2)
lines(N_vec,SmallN_res_list$Power_result_Spearman_slope_1,ylim = c(0,1),col = 'green',lty = 2,type = 'l',lwd = 2)
#plot type I error
lines(N_vec,SmallN_res_list$T1E_result,ylim = c(0,1),col = 'black', lty = 1, type = 'l',lwd = 2)
#plot wanted Type I error for comparison
abline(h = alpha,col = 'grey',lty = 2,lwd = 2)
dev.off()

