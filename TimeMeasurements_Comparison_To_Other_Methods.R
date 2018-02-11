#This is a script for comparing time measurements to other methods
setwd('~')
library(microbenchmark)
library(HHG)
library(dHSIC)
library(minerva)
library(energy)

MODE_TIME_RESULTS_BY_N = TRUE
MODE_TIME_RESULTS_BY_NR_ATOMS = TRUE

if(MODE_TIME_RESULTS_BY_N){
  
  MICROBENCHMARK_REPETITION = 100
  N_VEC = seq(500,5000,500)
  NULL_TABLE_SIZE = 1000
  PERMUTATIONS_FOR_TEST = 1000
  time_results = data.frame(Test = NA, N = NA,Time = NA)
  
  row = 1
  
  microbenchmark_wrapper_Fast_NA_45_ML = function(){
    Fast.independence.test(x,y,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 45)
  }
  
  microbenchmark_wrapper_dcov = function(){
    energy::dcov.test(x,y,R = PERMUTATIONS_FOR_TEST)
  }
  
  microbenchmark_wrapper_dHSIC = function(){
    dHSIC::dhsic.test(x, y, method="permutation", B = PERMUTATIONS_FOR_TEST)
  }
  
  microbenchmark_wrapper_MIC = function(){
    mic_res = minerva::mine(x,y)
    #MIC if very time expensive for full 1000 permutation test and time repetitions. We will compute a single statistic and multiply time by 1001, to imitate a full null table.
    #mic_res = minerva::mine(x,y)
    #mine_permutations = rep(NA,1000)
    #for(b in 1:1000){print(b);mine_permutations[b] = minerva::mine(x,sample(y))$MIC}
    #mic_pvalue = (1+sum(mine_permutations>=mic_res))/1001
    #mic_pvalue 
  }
  
  
  for(n_i in 1:length(N_VEC)){
    
    
    set.seed(1)
    current_N = N_VEC[n_i]
    x = rnorm(current_N)
    y = x + rnorm(current_N)
    
    print(paste0('Measuring Time for sample ',current_N))    
    
    print('MXL')
    mcb_F45_ML = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_45_ML(),replications = MICROBENCHMARK_REPETITION)$elapsed/MICROBENCHMARK_REPETITION
    print('dCOV')
    mcb_dcov = rbenchmark::benchmark(microbenchmark_wrapper_dcov(),replications = MICROBENCHMARK_REPETITION)$elapsed/MICROBENCHMARK_REPETITION
    print('dHSIC')
    mcb_dHSIC = rbenchmark::benchmark(microbenchmark_wrapper_dHSIC(),replications = MICROBENCHMARK_REPETITION)$elapsed/MICROBENCHMARK_REPETITION
    print('MIC')
    mcb_MIC = rbenchmark::benchmark(microbenchmark_wrapper_MIC(),replications = MICROBENCHMARK_REPETITION)$elapsed/MICROBENCHMARK_REPETITION * (NULL_TABLE_SIZE + 1)
    
    time_results[row,] = c('MXL , 45 Atoms',current_N,mcb_F45_ML) ; row = row + 1
    time_results[row,] = c('dCOV',current_N,mcb_dcov) ; row = row + 1
    time_results[row,] = c('dHSIC',current_N,mcb_dHSIC) ; row = row + 1
    time_results[row,] = c('MIC',current_N,mcb_MIC) ; row = row + 1
  }
  
  time_results$N = as.numeric(time_results$N)
  time_results$Time = as.numeric(time_results$Time)
  save(time_results,file = 'Time_Results_Comparison.RData')
  ggplot(time_results,aes(x = log(N),y = log(Time),color = Test))+geom_point()+geom_line()
  
  
}

if(MODE_TIME_RESULTS_BY_NR_ATOMS){
  
  #we measure times for ADP in this section
  set.seed(1)
  #nr atoms, in a lattice for simulation
  atoms_ind= c(10:30)*10
  #measured times
  time_ind_m_10 = rep(NA,length(atoms_ind))
  time_ind_m_15 = rep(NA,length(atoms_ind))
  
  
  microbenchmark_wrapper_MMAX_10 = function(atoms_i){
    res = HHG::hhg.univariate.ind.stat(1:300,sample(1:300),'ADP-EQP-ML',mmax = 10,nr.atoms = atoms_ind[atoms_i])
  }
  
  microbenchmark_wrapper_MMAX_15 = function(atoms_i){
    res = HHG::hhg.univariate.ind.stat(1:300,sample(1:300),'ADP-EQP-ML',mmax = 15,nr.atoms = atoms_ind[atoms_i])
  }
  
  
  
  #number of repetitions per setting
  nr.reps= 50
  for(atoms_i in 1:length(atoms_ind)){
    print(paste0('doing atoms_i: ',atoms_i ))
    time = rbenchmark::benchmark(microbenchmark_wrapper_MMAX_10(atoms_i),replications = nr.reps)$elapsed/nr.reps
    time_ind_m_10[atoms_i] = time
    
    time = rbenchmark::benchmark(microbenchmark_wrapper_MMAX_15(atoms_i),replications = nr.reps)$elapsed/nr.reps
    time_ind_m_15[atoms_i] = time
  }
  
  #save results and plot to file
  time_results = list(time_ind_m_10 = time_ind_m_10,time_ind_m_15 = time_ind_m_15)
  save(time_results,file = 'time_results.RData')
  PLOT_TO_FILE = T
  if(PLOT_TO_FILE){
    pdf('RunningTimeComparison.pdf',
        width     = 3.25,
        height    = 3.25,
        pointsize = 4
        #units     = "in",
        #res       = 1200,
        #pointsize = 4)
    )
  }
  plot(log(atoms_ind),log(time_ind_m_10),col='black',xlab='ln(Nr.Atoms)',ylab='ln(Time[Sec]/ 1[Sec])',pch='X',cex=0.8)
  points(log(atoms_ind),log(time_ind_m_15),col='red',cex=1)
  model_ind = lm(log(time_ind_m_10)~log(atoms_ind))
  abline(model_ind$coefficients[1],model_ind$coefficients[2],col='black',lty=2,lwd=1)
  model_ind_ml = lm(log(time_ind_m_15)~log(atoms_ind))
  abline(model_ind_ml$coefficients[1],model_ind_ml$coefficients[2],col='red',lty=3,lwd=1)
  ind_string = paste0('Linear fit, m.max = 10: ',round(model_ind$coefficients[1],2)," + x*",round(model_ind$coefficients[2],2),' ,R^2: ',round(summary(model_ind)$r.squared,6))
  ind_string_ml = paste0('Linear fit, m.max=15: ',round(model_ind_ml$coefficients[1],2)," + x*",round(model_ind_ml$coefficients[2],2),' ,R^2: ',round(summary(model_ind_ml)$r.squared,6))
  text(ind_string,x=4.9,y=4,col='black',cex=0.7)
  text(ind_string_ml,x=4.9,y=3.5,col='red',cex=0.7)
  if(PLOT_TO_FILE){
    dev.off() 
  }
}


if(F){
  
  load('time_results.RData') #=> time_results
  time_results_by_atoms = time_results
  load('Time_Results_Comparison.RData')  #=> time_results
  time_results$logN = log(time_results$N)  
  time_results$logT = log(time_results$Time)
  ind_MXL = which(time_results$Test == "MXL , 45 Atoms")
  ind_dCOV = which(time_results$Test == "dCOV")
  ind_dHSIC = which(time_results$Test == "dHSIC")
  ind_MIC = which(time_results$Test == "MIC")
  
  logtime_MXL = time_results$logT[ind_MXL]
  logn_MXL = time_results$logN[ind_MXL]
  
  logtime_dCOV = time_results$logT[ind_dCOV]
  logn_dCOV = time_results$logN[ind_dCOV]
  
  logtime_dHSIC = time_results$logT[ind_dHSIC]
  logn_dHSIC = time_results$logN[ind_dHSIC]
  
  logtime_MIC = time_results$logT[ind_MIC]
  logn_MIC = time_results$logN[ind_MIC]
  
  x_lim = c( 0.97*min(time_results$logN), 1.03 * max(time_results$logN))
  y_lim = c( 0.97*min(time_results$logT), 1.03 * max(time_results$logT))
  
  
  pdf(file = './Combined_Time_Chart.pdf',width = 8,height = 4,pointsize = 8)
  
  par(mfrow = c(1,2))
  lwd_param = 1
  plot(logn_MXL,logtime_MXL, xlim = x_lim, ylim = y_lim,type = 'b', xlab = 'ln(N)', ylab = 'ln(Time[Sec]/1[Sec])', lwd = lwd_param)
  lines(logn_dCOV,logtime_dCOV,col = 'red',type = 'b', lwd = lwd_param)
  lines(logn_dHSIC,logtime_dHSIC, col = 'blue',type = 'b', lwd = lwd_param)
  lines(logn_MIC,logtime_MIC, col = 'green',type = 'b', lwd = lwd_param)
  legend(x = 7.7, y = 2.2,col = c('black','red','blue','green'),legend = c("MXL , 45 Atoms","dCOV","dHSIC","MIC"),lty = 1, lwd = lwd_param ,pch = 1)
  
  atoms_ind= c(10:30)*10
  attach(time_results_by_atoms)
  plot(log(atoms_ind),log(time_ind_m_10),col='black',xlab='ln(Nr.Atoms)',ylab='ln(Time[Sec]/ 1[Sec])', pch='X', cex=0.8)
  points(log(atoms_ind),log(time_ind_m_15),col='red',cex=1)
  model_ind = lm(log(time_ind_m_10)~log(atoms_ind))
  abline(model_ind$coefficients[1],model_ind$coefficients[2],col='black',lty=2,lwd = 1)
  model_ind_ml = lm(log(time_ind_m_15)~log(atoms_ind))
  abline(model_ind_ml$coefficients[1],model_ind_ml$coefficients[2],col='red',lty=3,lwd = 1)
  ind_string = paste0('Linear fit, m.max = 10: ',round(model_ind$coefficients[1],2)," + x*",round(model_ind$coefficients[2],2),' ,R^2: ',round(summary(model_ind)$r.squared,6))
  ind_string_ml = paste0('Linear fit, m.max=15: ',round(model_ind_ml$coefficients[1],2)," + x*",round(model_ind_ml$coefficients[2],2),' ,R^2: ',round(summary(model_ind_ml)$r.squared,6))
  text(ind_string,x=5.1,y=4,col='black',cex=0.7)
  text(ind_string_ml,x=5.1,y=3.5,col='red',cex=0.7)
  par(mfrow = c(1,1))
  
  dev.off()
}