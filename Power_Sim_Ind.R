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
N = 2500 #Sample Size for simulation
NULL_TABLE_SIZE = 1000
POWER_REPETITIONS = 2000 #Number of realizations for power evaluation
PERMUTATIONS_FOR_TEST = NULL_TABLE_SIZE #Number of permutations for tests that require permutations
MMAX = 8
alpha = 0.05

#MODES:
MODE_SUBSECTION_C_PLOT_SETTINGS = FALSE
MODE_SUBSECTION_D_MEASURE_TIMES = FALSE
MODE_SUBSECTION_E_MIC_NULL_TABLE = TRUE
MODE_SUBSECTION_F_GENERATE_NULL_TABLES = TRUE
MODE_SUBSECTION_G_RUN_SCENARIOS = TRUE
MODE_SUBSECTION_H_ANALYZE_RESULTS = TRUE



#Subsection C: Plot Settings
if(MODE_SUBSECTION_C_PLOT_SETTINGS){
  set.seed(1)
  n.plot = 5000
  plts = list()
  for (i in 1:NR_SCENARIOS) {
    dat = data.frame(t(Scenario_list[[i]](n.plot, 0)))
    names(dat) = c('x', 'y')
    dat.noisy = data.frame(t(Scenario_list[[i]](N)))
    names(dat.noisy) = c('x', 'y')
    dat = rbind(cbind(dat, grp = 'clean'),cbind(dat.noisy, grp = 'noisy'))
    dat = dat[c((1:N)+n.plot,1:n.plot),]
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
    set.seed(1)
    x = rnorm(N)
    y = x + rnorm(N)
    
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
    
    MICROBENCHMARK_REPETITION = 100
    
    mcb_F15_ML = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_15_ML(),replications = MICROBENCHMARK_REPETITION)
    mcb_F15_MM = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_15_MM(),replications = MICROBENCHMARK_REPETITION)
    mcb_F30_ML = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_30_ML(),replications = MICROBENCHMARK_REPETITION)
    mcb_F30_MM = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_30_MM(),replications = MICROBENCHMARK_REPETITION)
    mcb_F45_ML = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_45_ML(),replications = MICROBENCHMARK_REPETITION)
    mcb_F45_MM = rbenchmark::benchmark(microbenchmark_wrapper_Fast_NA_45_MM(),replications = MICROBENCHMARK_REPETITION)
    mcb_dcov = rbenchmark::benchmark(microbenchmark_wrapper_dcov(),replications = MICROBENCHMARK_REPETITION)
    mcb_dHSIC = rbenchmark::benchmark(microbenchmark_wrapper_dHSIC(),replications = MICROBENCHMARK_REPETITION)
    mcb_MIC = rbenchmark::benchmark(microbenchmark_wrapper_MIC(),replications = MICROBENCHMARK_REPETITION)
    
    
    run_times = c(
      ML_15Atoms = mcb_F15_ML$elapsed/MICROBENCHMARK_REPETITION,
      MM_15Atoms = mcb_F15_MM$elapsed/MICROBENCHMARK_REPETITION,
      ML_30Atoms = mcb_F30_ML$elapsed/MICROBENCHMARK_REPETITION,
      MM_30Atoms = mcb_F30_MM$elapsed/MICROBENCHMARK_REPETITION,
      ML_45Atoms = mcb_F45_ML$elapsed/MICROBENCHMARK_REPETITION,
      MM_45Atoms = mcb_F45_MM$elapsed/MICROBENCHMARK_REPETITION,
      dCOV = mcb_dcov$elapsed/MICROBENCHMARK_REPETITION,
      dHSIC = mcb_dHSIC$elapsed/MICROBENCHMARK_REPETITION,
      MIC = mcb_MIC$elapsed/MICROBENCHMARK_REPETITION * (NULL_TABLE_SIZE + 1)
    )
    
    save(run_times,file = 'SIMULATION_RUN_TIMES.RData')
    
}


if(MODE_SUBSECTION_E_MIC_NULL_TABLE){
  
  NR.WORKERS = NR_CORES
  cl <- makeCluster(NR.WORKERS)
  registerDoParallel(cl)
  
  
  MIC_worker_function = function(){
    
    NR_REPS_PER_WORKER = ceiling(NULL_TABLE_SIZE/NR.WORKERS)
    NULL_DIST = rep(NA,NR_REPS_PER_WORKER)
    for(i in 1:NR_REPS_PER_WORKER){
      NULL_DIST[i] = minerva::mine(1:N,sample(1:N))$MIC
    }
    return(NULL_DIST)
  }
  
  MIC_null_dist_res <- foreach(i=1:NR.WORKERS, .options.RNG=1234) %dorng% { MIC_worker_function() }
  MIC_null_dist_as_array = unlist(MIC_null_dist_res)
  MIC_THRESHOLD_VALUE = as.numeric(quantile(MIC_null_dist_as_array,probs = 1-alpha))
  save(MIC_THRESHOLD_VALUE,file = 'MIC_THRESHOLD_VALUE.RData')
  stopCluster(cl)
}

if(MODE_SUBSECTION_F_GENERATE_NULL_TABLES){
  
    set.seed(1)
    NULL_TABLES  = list()
    NULL_TABLES[[1]] = Fast.independence.test.nulltable(N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 15)
    NULL_TABLES[[2]] = Fast.independence.test.nulltable(N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 15, variant = 'ADP-EQP')
    NULL_TABLES[[3]] = Fast.independence.test.nulltable(N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 30)
    NULL_TABLES[[4]] = Fast.independence.test.nulltable(N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 30, variant = 'ADP-EQP')
    NULL_TABLES[[5]] = Fast.independence.test.nulltable(N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 45)
    NULL_TABLES[[6]] = Fast.independence.test.nulltable(N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 45, variant = 'ADP-EQP')
    
    PROPOSED_TESTS = c('ML_15ATOMS','MM_15ATOMS','ML_30ATOMS','MM_30ATOMS','ML_45ATOMS','MM_45ATOMS')
    names(NULL_TABLES) = PROPOSED_TESTS
    save(NULL_TABLES,file = 'NULL_TABLES.RData')     
}


#Subsection G: Run Scenarios (multiple core)
if(MODE_SUBSECTION_G_RUN_SCENARIOS){
  
  Simulation_results = NULL
  
  POWER_SIM_WORKER_FUNCTION = function(SCENARIO_ID){
    
    library(HHG)
    library(dHSIC)
    library(minerva)
    library(energy)
    
    load(file = 'MIC_THRESHOLD_VALUE.RData') #=> MIC_THRESHOLD_VALUE
    load(file = 'NULL_TABLES.RData') #=> NULL_TABLES
    
    NR_REPS_PER_WORKER = ceiling(POWER_REPETITIONS/NR.WORKERS)
    result_names = c('REPS',names(NULL_TABLES),'dCOV','dHSIC','MIC')
    results = matrix(0,nrow = 1,ncol = length(result_names))     
    colnames(results) = result_names
    for(i in 1:NR_REPS_PER_WORKER){
      data = Scenario_list[[SCENARIO_ID]](N)
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
  for(SCENARIO_ID_FOR_SIM  in Scenarios_to_run){
    print(paste0('Running Scenario ID ',SCENARIO_ID_FOR_SIM))
    NR.WORKERS = NR_CORES
    cl <- makeCluster(NR.WORKERS)
    registerDoParallel(cl)
    SCENARIO_SIM_RESULTS <- foreach(i=1:NR.WORKERS, .options.RNG=1234,.combine = 'rbind') %dorng% { POWER_SIM_WORKER_FUNCTION(SCENARIO_ID_FOR_SIM) }
    stopCluster(cl)
    sim_results_as_row = apply(SCENARIO_SIM_RESULTS,2,sum)
    save(sim_results_as_row,file = paste0('ind_sim_results_as_row_',SCENARIO_ID_FOR_SIM,'.RData'))
    Simulation_results = rbind(Simulation_results,sim_results_as_row) 
  }
  
  Power_results = Simulation_results
  for(i in 1:nrow(Power_results))
    Power_results[i,] = Power_results[i,] / Power_results[i,1]
  rownames(Power_results) = Scenario_names[Scenarios_to_run]
  save(Power_results,file = paste0('Ind_Power_Results.RData'))
}

#Subsection H: Analyze Results
if(MODE_SUBSECTION_H_ANALYZE_RESULTS){
  load(file = 'Ind_Power_Results.RData') # => Power_results
  load(file = 'SIMULATION_RUN_TIMES.RData') #=> run_times
  
  colnames(Power_results) = c("REPS","MXL, 15 Atoms", "MXM, 15 Atoms", "MXL, 30 Atoms", "MXM, 30 Atoms", "MXL, 45 Atoms", "MXM, 45 Atoms", "dCOV", "dHSIC" , "MIC" )
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
  Test_names_vec = colnames(Power_results)[c(8,9,10,2,4,6,3,5,7)]
  
  pdf('IndPowerResults.PDF', width = 8, height = 3)
  
  ggplot(power_res_plot) + geom_point(aes(x = log(RunTime) , y = Power ,shape = Test,color = Test))+facet_wrap(~Scenario,nrow = 1,ncol = 4) +
    scale_color_manual(breaks = Test_names_vec,values = c(1,1,1,2,3,4,2,3,4)) +
    scale_shape_manual(breaks = Test_names_vec,values = c(0,1,2,16,16,16,17,17,17)) +
    xlab("log(RunTime[Seconds])") + ylim(c(0,1)) + theme_bw()
  
  dev.off()
}