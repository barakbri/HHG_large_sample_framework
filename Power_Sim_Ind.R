setwd('~')
NR_CORES = detectCores() - 1 # currently set to 7, for reproducability
#Power Simulation for independence settings

library(ggplot2)
library(grid)
library(rbenchmark)
library(HHG)
library(dHSIC)
library(minerva)
library(energy)

library(parallel)
library(doParallel)
library(foreach)
library(doRNG)

#Subsection A: Scenarios



# Linear
datagenLine = function(n, m = n) {
  if (m == 0) {
    sig = 0
  } else if (m == 30) {
    sig = 0.6
  } else if (m == 50) {
    sig = 1
  } else if (m == 100) {
    sig = 1.5
  } else if (m == 300) {
    sig = 2
  } else if(m == 500){
    sig = 2.5
  } else if(m == 700){
    sig = 3.0
  } else if (m == 2500) {
    sig = 5
  } 
  
  x = runif(n, 0, 1)
  y = x 
  y = y + rnorm(n) * sig
  return (rbind(x, y))
}

# Exponential [2 ^ x]
datagenExp2x = function(n, m = n) {
  if (m == 0) {
    sig = 0
  } else if (m == 30) {
    sig = 10 ^ 2.7
  } else if (m == 50) {
    sig = 10 ^ 2.7
  } else if (m == 100) {
    sig = 10 ^ 2.8
  } else if (m == 300) {
    sig = 10 ^ 2.9
  } else if (m == 500) {
    sig = 10 ^ 3.0
  } else if (m == 700) {
    sig = 10 ^ 3.1
  } else if (m == 2500) {
    sig = 10 ^ 3.6
  }
  
  x = runif(n, 0, 10)
  y = 2^x 
  y = y + rnorm(n) * sig
  return (rbind(x, y))
}


# Sine
datagenSine = function(n, m = n) {
  if (m == 0) {
    sig = 0
  } else if (m == 30) {
    sig = 0.3
  } else if (m == 50) {
    sig = 0.3
  } else if (m == 100) {
    sig = 0.5
  } else if (m == 300) {
    sig = 1.5
  } else if (m == 500) {
    sig = 1.75
  } else if (m == 700) {
    sig = 2.0
  } else if (m == 2500) {
    sig = 5
  }
  
  x = runif(n, 0, 1)
  y = sin(8 * pi * x)
  y = y + rnorm(n) * sig    
  return (rbind(x, y))
}

# Circles
datagenCircles = function(n, m = n) {
  if (m == 0) {
    sig = 0
  } else if (m == 30) {
    sig = 0.01
  } else if (m == 50) {
    sig = 0.03
  } else if (m == 100) {
    sig = 0.05
  } else if (m == 300) {
    sig = 0.2
  } else if (m == 2500) {
    sig = 0.55
  }
  
  cx = sample(c(-1, 1), n, replace = T)
  cy = sample(c(-1, 1), n, replace = T)
  r = sqrt(2)
  phi = runif(n, 0, 2*pi)
  x = cx + r * cos(phi)
  y = cy + r * sin(phi)
  y = y + rnorm(n) * sig
  x = x + rnorm(n) * sig
  return (rbind(x, y))
}  


Scenario_list = list(datagenLine,datagenExp2x,datagenCircles,datagenSine)
Scenario_names = c('Line','Exp2x','Circles','Sine')

#Subsection B: Declarations

NR_SCENARIOS = length(Scenario_list)
N_vec = c(500,1000,1500,2000,2500)
NULL_TABLE_SIZE = 1000
POWER_REPETITIONS = 2000 #Number of realizations for power evaluation
PERMUTATIONS_FOR_TEST = NULL_TABLE_SIZE #Number of permutations for tests that require permutations
MMAX = 10
alpha = 0.05
RBENCHMARK_REPETITION = 100

#MODES:
MODE_SUBSECTION_C_PLOT_SETTINGS = TRUE
MODE_SUBSECTION_D_MEASURE_TIMES = FALSE
MODE_SUBSECTION_E_MIC_NULL_TABLE = FALSE
MODE_SUBSECTION_F_GENERATE_NULL_TABLES = FALSE
MODE_SUBSECTION_G_RUN_SCENARIOS = FALSE
MODE_SUBSECTION_H_ANALYZE_RESULTS = TRUE
MODE_SUBSECTION_I_PLOT_POWER_BY_N = TRUE

#functions for filenames:

#function for time measurements by N:
get_filename_timemeasurements = function(N){
  return(paste0('SIMULATION_RUN_TIMES_N_',N,'.RData'))
}

#function for MIC table by N:
get_filename_MIC_Threshold = function(N){
  return(paste0('MIC_THRESHOLD_VALUE_N_',N,'.RData'))
}

#function for ADP Null Tables by N:
get_filename_ADP_Null_Tables = function(N){
  return(paste0('NULL_TABLES_N_',N,'.RData'))
}

#function for results_filename
get_filename_Power_Results = function(N){
  return(paste0('Ind_Power_Results_N_',N,'.RData'))
}

get_filename_Power_Results_Graph = function(N){
  return(paste0('IndPowerResults_N_',N,'.PDF'))
}

#Subsection C: Plot Settings
if(MODE_SUBSECTION_C_PLOT_SETTINGS){
  set.seed(1)
  n.plot = 5000
  plts = list()
  for (i in 1:NR_SCENARIOS) {
    current_N = N_vec[length(N_vec)]
    dat = data.frame(t(Scenario_list[[i]](n.plot, 0)))
    names(dat) = c('x', 'y')
    dat.noisy = data.frame(t(Scenario_list[[i]](current_N)))
    names(dat.noisy) = c('x', 'y')
    dat = rbind(cbind(dat, grp = 'clean'),cbind(dat.noisy, grp = 'noisy'))
    dat = dat[c((1:current_N)+n.plot,1:n.plot),]
    plts[[i]] = ggplot(data = dat) + 
      geom_point(aes(x = x, y = y, colour = as.factor(grp), size = grp)) + #, alpha = 0.8
      theme_minimal(base_size = 10) + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
      ggtitle(Scenario_names[i]) + scale_color_discrete(guide = 'none') + scale_size_manual(values = c(0.1,0.2),guide = 'none')
    scale_size_discrete(range = c(1, 2), guide = 'none')
  }
  
  pdf('IndScenarios.PDF', width = 8, height = 2.2)
  
  grid.newpage()
  nr.cols = 4
  nr.rows = 1
  pushViewport(viewport(layout = grid.layout(nr.rows, nr.cols, widths = c(300, 300))))
  
  ix = 1
  iy = 1
  for (i in 1:NR_SCENARIOS){
    print(plts[[i]], vp = viewport(layout.pos.row = iy, layout.pos.col = ix))
    if (NR_SCENARIOS <= 6) {
      iy = iy + 1
      if (iy > nr.rows) {
        iy = 1
        ix = ix + 1
      }
    } else {
      ix = ix + 1
      if (ix > nr.cols) {
        ix = 1
        iy = iy + 1
      }
    }
  }
  dev.off()
}


#Subsection D: Time Measurements & Constants (e.g. MIC Threshold value)
if(MODE_SUBSECTION_D_MEASURE_TIMES){
  for(current_N_index in 1:length(N_vec)){
    current_N = N_vec[current_N_index]
    print(paste0('Measuring times for N: ',current_N))
    
    set.seed(1)
    x = rnorm(current_N)
    y = x + rnorm(current_N)
    
    microbenchmark_wrapper_Fast_NA_15_ML = function(){
      res = Fast.independence.test(x,y,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 15)
    }
    
    microbenchmark_wrapper_Fast_NA_15_MM = function(){
      res = Fast.independence.test(x,y,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 15, variant = 'ADP-EQP')
    }
    
    microbenchmark_wrapper_Fast_NA_30_ML = function(){
      res = Fast.independence.test(x,y,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 30)
    }
    
    microbenchmark_wrapper_Fast_NA_30_MM = function(){
      res = Fast.independence.test(x,y,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 30, variant = 'ADP-EQP')
    }
    
    microbenchmark_wrapper_Fast_NA_45_ML = function(){
      res = Fast.independence.test(x,y,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 45)
    }
    
    microbenchmark_wrapper_Fast_NA_45_MM = function(){
      res = Fast.independence.test(x,y,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 45, variant = 'ADP-EQP')
    }
    
    microbenchmark_wrapper_Fast_NA_60_ML = function(){
      res = Fast.independence.test(x,y,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 60)
    }
    
    microbenchmark_wrapper_Fast_NA_60_MM = function(){
      res = Fast.independence.test(x,y,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 60, variant = 'ADP-EQP')
    }
    
    microbenchmark_wrapper_dcov = function(){
      res = energy::dcov.test(x,y,R = PERMUTATIONS_FOR_TEST)
    }
    
    microbenchmark_wrapper_dHSIC = function(){
      res = dHSIC::dhsic.test(x, y, method="permutation", B = PERMUTATIONS_FOR_TEST)
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
    
    
    print('Measureing Time for 15 atoms')
    mcb_F15_ML = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_15_ML(),replications = RBENCHMARK_REPETITION)
    mcb_F15_MM = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_15_MM(),replications = RBENCHMARK_REPETITION)
    print('Measureing Time for 30 atoms')
    mcb_F30_ML = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_30_ML(),replications = RBENCHMARK_REPETITION)
    mcb_F30_MM = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_30_MM(),replications = RBENCHMARK_REPETITION)
    print('Measureing Time for 45 atoms')
    mcb_F45_ML = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_45_ML(),replications = RBENCHMARK_REPETITION)
    mcb_F45_MM = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_45_MM(),replications = RBENCHMARK_REPETITION)
    print('Measureing Time for 60 atoms')
    mcb_F60_ML = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_60_ML(),replications = RBENCHMARK_REPETITION)
    mcb_F60_MM = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_60_MM(),replications = RBENCHMARK_REPETITION)
    print('Measureing Time for dCOV')
    mcb_dcov = rbenchmark::benchmark(microbenchmark_wrapper_dcov(),replications = RBENCHMARK_REPETITION)
    print('Measureing Time for dHSIC')
    mcb_dHSIC = rbenchmark::benchmark(microbenchmark_wrapper_dHSIC(),replications = RBENCHMARK_REPETITION)
    print('Measureing Time for MIC')
    mcb_MIC = rbenchmark::benchmark(microbenchmark_wrapper_MIC(),replications = RBENCHMARK_REPETITION)
    
    
    run_times = c(
      ML_15Atoms = mcb_F15_ML$elapsed/RBENCHMARK_REPETITION,
      MM_15Atoms = mcb_F15_MM$elapsed/RBENCHMARK_REPETITION,
      ML_30Atoms = mcb_F30_ML$elapsed/RBENCHMARK_REPETITION,
      MM_30Atoms = mcb_F30_MM$elapsed/RBENCHMARK_REPETITION,
      ML_45Atoms = mcb_F45_ML$elapsed/RBENCHMARK_REPETITION,
      MM_45Atoms = mcb_F45_MM$elapsed/RBENCHMARK_REPETITION,
      ML_60Atoms = mcb_F60_ML$elapsed/RBENCHMARK_REPETITION,
      MM_60Atoms = mcb_F60_MM$elapsed/RBENCHMARK_REPETITION,
      dCOV = mcb_dcov$elapsed/RBENCHMARK_REPETITION,
      dHSIC = mcb_dHSIC$elapsed/RBENCHMARK_REPETITION,
      MIC = mcb_MIC$elapsed/RBENCHMARK_REPETITION * (NULL_TABLE_SIZE + 1)
    )
    filename = get_filename_timemeasurements(current_N)
    save(run_times,file = filename)
    
  }
}


if(MODE_SUBSECTION_E_MIC_NULL_TABLE){
  for(current_N_i in 1:length(N_vec)){
    current_N = N_vec[current_N_i]
    print(paste0("MIC null table for N: ",current_N))
    NR.WORKERS = NR_CORES
    cl <- makeCluster(NR.WORKERS)
    registerDoParallel(cl)
    
    
    MIC_worker_function = function(){
      
      NR_REPS_PER_WORKER = ceiling(NULL_TABLE_SIZE/NR.WORKERS)
      NULL_DIST = rep(NA,NR_REPS_PER_WORKER)
      for(i in 1:NR_REPS_PER_WORKER){
        NULL_DIST[i] = minerva::mine(1:current_N,sample(1:current_N))$MIC
      }
      return(NULL_DIST)
    }
    
    MIC_null_dist_res <- foreach(i=1:NR.WORKERS, .options.RNG=1234) %dorng% { MIC_worker_function() }
    MIC_null_dist_as_array = unlist(MIC_null_dist_res)
    MIC_THRESHOLD_VALUE = as.numeric(quantile(MIC_null_dist_as_array,probs = 1-alpha))
    filename_MIC = get_filename_MIC_Threshold(current_N)
    save(MIC_THRESHOLD_VALUE,file = filename_MIC)
    stopCluster(cl)
  }
  
}

if(MODE_SUBSECTION_F_GENERATE_NULL_TABLES){
    for(current_N_i in 1:length(N_vec)){
      current_N = N_vec[current_N_i]
      print(paste0('Building Null Tables for N: ',current_N))
      set.seed(1)
      NULL_TABLES  = list()
      NULL_TABLES[[1]] = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 15)
      NULL_TABLES[[2]] = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 15, variant = 'ADP-EQP')
      NULL_TABLES[[3]] = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 30)
      NULL_TABLES[[4]] = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 30, variant = 'ADP-EQP')
      NULL_TABLES[[5]] = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 45)
      NULL_TABLES[[6]] = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 45, variant = 'ADP-EQP')
      NULL_TABLES[[7]] = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 60)
      NULL_TABLES[[8]] = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 60, variant = 'ADP-EQP')
      
      PROPOSED_TESTS = c('ML_15ATOMS','MM_15ATOMS','ML_30ATOMS','MM_30ATOMS','ML_45ATOMS','MM_45ATOMS','ML_60ATOMS','MM_60ATOMS')
      names(NULL_TABLES) = PROPOSED_TESTS
      filename = get_filename_ADP_Null_Tables(current_N)
      save(NULL_TABLES,file = filename)
    }
}


#Subsection G: Run Scenarios (multiple core)
if(MODE_SUBSECTION_G_RUN_SCENARIOS){
  
  Simulation_results = NULL
  
  POWER_SIM_WORKER_FUNCTION = function(SCENARIO_ID,sample_size_index){
    
    library(HHG)
    library(dHSIC)
    library(minerva)
    library(energy)
    current_N = N_vec[sample_size_index]
    load(file = get_filename_MIC_Threshold(current_N)) #=> MIC_THRESHOLD_VALUE
    load(file = get_filename_ADP_Null_Tables(current_N)) #=> NULL_TABLES
    
    NR_REPS_PER_WORKER = ceiling(POWER_REPETITIONS/NR.WORKERS)
    result_names = c('REPS',names(NULL_TABLES),'dCOV','dHSIC','MIC')
    results = matrix(0,nrow = 1,ncol = length(result_names))     
    colnames(results) = result_names
    for(i in 1:NR_REPS_PER_WORKER){
      data = Scenario_list[[SCENARIO_ID]](current_N,N_vec[length(N_vec)])
      x = data[1,]
      y = data[2,]
      
      results[1,1] = results[1,1] + 1
      for(n in 1:length(NULL_TABLES)){
        MinP_res = Fast.independence.test(x,y,NullTable = NULL_TABLES[[n]])
        results[1,n + 1] = results[1,n + 1] + 1 *(MinP_res$MinP.pvalue  <= alpha) 
      }
      
      dCov_res = energy::dcov.test(x,y,R = PERMUTATIONS_FOR_TEST)
      results[1,length(NULL_TABLES) + 2] = results[1,length(NULL_TABLES) + 2] + 1 * (dCov_res$p.value <= alpha)
      
      dHSIC_res = dHSIC::dhsic.test(x, y, method="permutation", B = PERMUTATIONS_FOR_TEST)
      results[1,length(NULL_TABLES) + 3] = results[1,length(NULL_TABLES) + 3] + 1 * (dHSIC_res$p.value <= alpha)
      
      MIC_res = mic_res = minerva::mine(x,y)
      results[1,length(NULL_TABLES) + 4] = results[1,length(NULL_TABLES) + 4] + 1 * (MIC_res$MIC >= MIC_THRESHOLD_VALUE)
    }
    return(results)
  }
  
  Scenarios_to_run = 1:NR_SCENARIOS
  
  for(current_N_ind in 1:length(N_vec)){
    Simulation_results = NULL
    current_N = N_vec[current_N_ind]
    for(SCENARIO_ID_FOR_SIM  in Scenarios_to_run){
      print(paste0('Running Scenario ID ',SCENARIO_ID_FOR_SIM,' N index :',current_N_ind,' Sample size: ',current_N))
      NR.WORKERS = NR_CORES
      cl <- makeCluster(NR.WORKERS)
      registerDoParallel(cl)
      SCENARIO_SIM_RESULTS <- foreach(i=1:NR.WORKERS, .options.RNG=1234,.combine = 'rbind') %dorng% { POWER_SIM_WORKER_FUNCTION(SCENARIO_ID_FOR_SIM,current_N_ind) }
      stopCluster(cl)
      sim_results_as_row = apply(SCENARIO_SIM_RESULTS,2,sum)
      save(sim_results_as_row,file = paste0('ind_sim_results_as_row_',SCENARIO_ID_FOR_SIM,'_N_',current_N,'.RData'))
      Simulation_results = rbind(Simulation_results,sim_results_as_row) 
    }
    
    Power_results = Simulation_results
    for(i in 1:nrow(Power_results))
      Power_results[i,] = Power_results[i,] / Power_results[i,1]
    rownames(Power_results) = Scenario_names[Scenarios_to_run]
    save(Power_results,file = get_filename_Power_Results(current_N))
  }
}

#Subsection H: Analyze Results
if(MODE_SUBSECTION_H_ANALYZE_RESULTS){
  
  for(current_N_index in 1:length(N_vec)){
    current_N = N_vec[current_N_index]
    
    load(file = get_filename_Power_Results(current_N)) # => Power_results
    load(file = get_filename_timemeasurements(current_N)) #=> run_times
    
    colnames(Power_results) = c("REPS","mXl, 15 Atoms", "mXm, 15 Atoms", "mXl, 30 Atoms", "mXm, 30 Atoms", "mXl, 45 Atoms", "mXm, 45 Atoms", "mXl, 60 Atoms", "mXm, 60 Atoms", "dCOV", "dHSIC" , "MIC" )
    #run_times
    power_res_plot = data.frame(Test = NA, Scenario = NA, RunTime = NA, Power = NA)
    current_row = 1
    for(i in 1:nrow(Power_results)){
      for(j in 2:ncol(Power_results)){
        current_test = colnames(Power_results)[j]
        current_scenario = rownames(Power_results)[i]
        current_power = Power_results[i,j]
        current_Runtime = run_times[j-1]
        power_res_plot[current_row,] = c(current_test,current_scenario,current_Runtime,current_power)
        current_row = current_row + 1
      }
    }
    power_res_plot$RunTime = as.numeric(power_res_plot$RunTime)
    power_res_plot$Power = as.numeric(power_res_plot$Power)
    power_res_plot$Scenario = factor(power_res_plot$Scenario ,levels = c('Line','Exp2x','Circles','Sine'))
    Test_names_vec = colnames(Power_results)[c(10,11,12,2,4,6,8,3,5,7,9)]
    
    pdf(get_filename_Power_Results_Graph(current_N), width = 8, height = 3)
    
    print(ggplot(power_res_plot) + geom_point(aes(x = log(RunTime) , y = Power ,shape = Test,color = Test))+facet_wrap(~Scenario,nrow = 1,ncol = 4) +
      scale_color_manual(breaks = Test_names_vec,values = c(1,1,1,2,3,4,1,2,3,4,1)) +
      scale_shape_manual(breaks = Test_names_vec,values = c(0,1,2,16,16,16,16,17,17,17,17)) +
      xlab("log(RunTime[Seconds])") + ylim(c(0,1)) + theme_bw())
    
    dev.off()
  }
 
}

if(MODE_SUBSECTION_I_PLOT_POWER_BY_N){
  cols_ind_to_graph = c(6,10,11,12)
  power_list = list()
  for(current_N_ind in 1:length(N_vec)){
    current_N = N_vec[current_N_ind]  
    load(file = get_filename_Power_Results(current_N)) # => Power_results
    power_list[[current_N_ind]] = Power_results
  }
  power_matrix = data.frame(N = NA,Test = NA, Power = NA,Scenario = NA)
  test_names = c("mxl, 45 Atoms", "dCOV","dHSIC","MIC" )
  power_matrix_pointer = 1
  for(i in 1:length(power_list)){
    for(j in 1:length(cols_ind_to_graph)){
      for(Scenario_ID in 1:length(Scenario_names)){
        power_matrix[power_matrix_pointer , ] = c(N_vec[i],test_names[j],power_list[[i]][Scenario_ID, cols_ind_to_graph[j]],Scenario_names[Scenario_ID])
        power_matrix_pointer = power_matrix_pointer + 1
      }
    }
  }
  power_matrix$N = as.numeric(power_matrix$N)
  power_matrix$Power = as.numeric(power_matrix$Power)
  power_matrix$Scenario = factor(power_matrix$Scenario,levels = Scenario_names)
  
  
  pdf("IndPower_by_N.pdf", width = 8, height = 2.5)
  
  ggplot(data = power_matrix,mapping = aes(x= N,y = Power,color = Test)) +geom_point()+geom_line() +
    facet_wrap(~Scenario,nrow =1, ncol = 4) + ylim(c(0,1)) +  theme_bw() +
    theme(axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 6)) 
  
  dev.off()
}