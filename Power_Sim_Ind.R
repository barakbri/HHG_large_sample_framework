# Power Simulation for independence settings 
# See brill-heller-heller.R for an index used for running all scripts.

#setwd('~') #set in main script, disabled here
NR_CORES = detectCores() - 1 # The simulation uses the maximum number of cores available.
# For reproducability, please run on an Amazon C5.18XL machine with 72 cores.
# On any machine with 72 cores (or by setting NR_CORES = 71), results would be reproducible with paper.
# On an amazon C5 machine, run times will be reproducible as well.


#load libraries required
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
#***********************
# Here we define the scenarios for the simulation.
#Each scenario has a parameter n, setting the sample size, and a parameter m, setting the noise level.
#Several sample sizes have a default noise level.


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

#This is a list of the scenarios, for iterating over them, along with their names
Scenario_list = list(datagenLine,datagenExp2x,datagenCircles,datagenSine)
Scenario_names = c('Line','Exp2x','Circles','Sine')

#Subsection B: Declarations
#***********************

NR_SCENARIOS = length(Scenario_list) # the number of scenarios in the simulation
N_vec = c(500,1000,1500,2000,2500) #sample sizes for simulation. The noise level is set by the last one.
NULL_TABLE_SIZE = 1000 # Size of null table for ADP and MIC tests
POWER_REPETITIONS = 2000 #Number of realizations for power evaluation, 
PERMUTATIONS_FOR_TEST = NULL_TABLE_SIZE #Number of permutations for tests that require permutations
MMAX = 10 # parameter for m.max is set to be the package default
alpha = 0.05 #requried type I error rate
RBENCHMARK_REPETITION = 100 #number of repetitions, for evaluating time
MULTIPLIER_FOR_MINP = 5 #Multiplier for the number of repetitions for power evaluation, for the MINP procedure

#MODES:
MODE_SUBSECTION_C_PLOT_SETTINGS = FALSE # subsection used for plotting the settings to file
MODE_SUBSECTION_D_MEASURE_TIMES = FALSE # subsection used for measuring times, for each of the tests
MODE_SUBSECTION_E_MIC_NULL_TABLE = FALSE # subsection for computing the MIC critical value for rejection. (multiple cores)
MODE_SUBSECTION_F_GENERATE_NULL_TABLES = FALSE #subsection for generating null tables for MinP ADP.
MODE_SUBSECTION_G_RUN_SCENARIOS = FALSE #subsection for running the actual simulations. (multiple cores)
MODE_SUBSECTION_H_ANALYZE_RESULTS = TRUE #subsection for plotting results, by N
MODE_SUBSECTION_I_PLOT_POWER_BY_N = TRUE # subsection for plotting results, comparinig power over several N

MODE_ADD_ATOMS = TRUE # This mode is used for computing results also for N_Atoms = 5,10, as an extension of the original simulation. With this flag on, run times for all other statistics will be taken from previous results

#functions for filenames:
#these are used to access files, according to their sample sizes.

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
#***********************
if(MODE_SUBSECTION_C_PLOT_SETTINGS){
  set.seed(1)
  n.plot = 5000 #number of points used to plot the noiseless curve
  plts = list()
  for (i in 1:NR_SCENARIOS) {
    current_N = N_vec[length(N_vec)] #generate data with no noise
    dat = data.frame(t(Scenario_list[[i]](n.plot, 0)))
    names(dat) = c('x', 'y')
    dat.noisy = data.frame(t(Scenario_list[[i]](current_N))) #generate data with nois
    names(dat.noisy) = c('x', 'y')
    #combine the datasets
    dat = rbind(cbind(dat, grp = 'clean'),cbind(dat.noisy, grp = 'noisy'))
    dat = dat[c((1:current_N)+n.plot,1:n.plot),]
    #plot both datasets, in different colors
    plts[[i]] = ggplot(data = dat) + 
      geom_point(aes(x = x, y = y, colour = as.factor(grp), size = grp)) + #, alpha = 0.8
      theme_minimal(base_size = 10) + 
      theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.title = element_text(hjust = 0.5)) +
      ggtitle(Scenario_names[i]) + scale_color_discrete(guide = 'none') + scale_size_manual(values = c(0.1,0.2),guide = 'none')
    scale_size_discrete(range = c(1, 2), guide = 'none')
  }
  
  # plot to file, over a grid
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
#***********************
if(MODE_SUBSECTION_D_MEASURE_TIMES){
  #we measure times for each of the sample sizes
  for(current_N_index in 1:length(N_vec)){
    current_N = N_vec[current_N_index]
    print(paste0('Measuring times for N: ',current_N))
    
    Prev_run_times = NULL
    if(MODE_ADD_ATOMS){
      filename = get_filename_timemeasurements(current_N)
      load(filename) # => run_times
      Prev_run_times = run_times
    }
    
    
    ind_in_prev = function(s){return(which(names(Prev_run_times)==s))} # we assume there is only one, used for MODE_ADD_ATOMS
    
    #generate data
    set.seed(1)
    x = rnorm(current_N)
    y = x + rnorm(current_N)
    
    #wrappers for ADP/MinP, with all steps included - null table generation, statistic computation, from start to end.
    #wrappers include 15,30,45,60 atoms, for mXl and mXm variants.
    
    microbenchmark_wrapper_Fast_NA_ML = function(N_A){
      res = Fast.independence.test(x,y,mmax = min(MMAX,N_A), nr.perm = NULL_TABLE_SIZE,nr.atoms = N_A)
    }
    
    microbenchmark_wrapper_Fast_NA_MM = function(N_A){
      res = Fast.independence.test(x,y,mmax = min(MMAX,N_A), nr.perm = NULL_TABLE_SIZE,nr.atoms = N_A, variant = 'ADP-EQP')
    }
    
    #wrapper for measuring time, for dCov
    microbenchmark_wrapper_dcov = function(){
      res = energy::dcov.test(x,y,R = PERMUTATIONS_FOR_TEST)
    }
    
    #wrapper for measuring time, for dHSIC
    microbenchmark_wrapper_dHSIC = function(){
      res = dHSIC::dhsic.test(x, y, method="permutation", B = PERMUTATIONS_FOR_TEST)
    }
    
    #wrapper for measurign time, for MIC
    microbenchmark_wrapper_MIC = function(){
      mic_res = minerva::mine(x,y)
      #MIC is very time expensive for full 1000 permutation test and time repetitions.
      #We will compute a single statistic and multiply time by 1001, to imitate a full null table.
      #mic_res = minerva::mine(x,y)
      #mine_permutations = rep(NA,1000)
      #for(b in 1:1000){print(b);mine_permutations[b] = minerva::mine(x,sample(y))$MIC}
      #mic_pvalue = (1+sum(mine_permutations>=mic_res))/1001
      #mic_pvalue 
    }
    
    # call the time measurement fucntions
    
    mcb_F15_ML = NULL
    mcb_F15_MM = NULL
    mcb_F30_ML = NULL
    mcb_F30_MM = NULL
    mcb_F45_ML = NULL
    mcb_F45_MM = NULL
    mcb_F60_ML = NULL
    mcb_F60_MM = NULL
    if(!MODE_ADD_ATOMS){
      print('Measureing Time for 15 atoms')
      mcb_F15_ML = rbenchmark::benchmark({function(){microbenchmark_wrapper_Fast_NA_ML(15)}}(),replications = RBENCHMARK_REPETITION)
      mcb_F15_MM = rbenchmark::benchmark({function(){microbenchmark_wrapper_Fast_NA_MM(15)}}(),replications = RBENCHMARK_REPETITION)
      print('Measureing Time for 30 atoms')
      mcb_F30_ML = rbenchmark::benchmark({function(){microbenchmark_wrapper_Fast_NA_ML(30)}}(),replications = RBENCHMARK_REPETITION)
      mcb_F30_MM = rbenchmark::benchmark({function(){microbenchmark_wrapper_Fast_NA_MM(30)}}(),replications = RBENCHMARK_REPETITION)
      print('Measureing Time for 45 atoms')
      mcb_F45_ML = rbenchmark::benchmark({function(){microbenchmark_wrapper_Fast_NA_ML(45)}}(),replications = RBENCHMARK_REPETITION)
      mcb_F45_MM = rbenchmark::benchmark({function(){microbenchmark_wrapper_Fast_NA_MM(45)}}(),replications = RBENCHMARK_REPETITION)
      print('Measureing Time for 60 atoms')
      mcb_F60_ML = rbenchmark::benchmark({function(){microbenchmark_wrapper_Fast_NA_ML(60)}}(),replications = RBENCHMARK_REPETITION)
      mcb_F60_MM = rbenchmark::benchmark({function(){microbenchmark_wrapper_Fast_NA_MM(60)}}(),replications = RBENCHMARK_REPETITION)  
    }else{ #we are running in add atom modes - long run times (from previous simulation) are taken from result files
      print('Getting Times for Atoms 15,30,45,60')
      
      mcb_F15_ML = list(elapsed = Prev_run_times[ind_in_prev('ML_15Atoms')] * RBENCHMARK_REPETITION)
      mcb_F15_MM = list(elapsed = Prev_run_times[ind_in_prev('MM_15Atoms')] * RBENCHMARK_REPETITION)
      mcb_F30_ML = list(elapsed = Prev_run_times[ind_in_prev('ML_30Atoms')] * RBENCHMARK_REPETITION)
      mcb_F30_MM = list(elapsed = Prev_run_times[ind_in_prev('MM_30Atoms')] * RBENCHMARK_REPETITION)
      mcb_F45_ML = list(elapsed = Prev_run_times[ind_in_prev('ML_45Atoms')] * RBENCHMARK_REPETITION)
      mcb_F45_MM = list(elapsed = Prev_run_times[ind_in_prev('MM_45Atoms')] * RBENCHMARK_REPETITION)
      mcb_F60_ML = list(elapsed = Prev_run_times[ind_in_prev('ML_60Atoms')] * RBENCHMARK_REPETITION)
      mcb_F60_MM = list(elapsed = Prev_run_times[ind_in_prev('MM_60Atoms')] * RBENCHMARK_REPETITION)
      
    }
    #times for 5,10 atoms are always measured
    mcb_F5_ML = NULL
    mcb_F5_MM = NULL
    mcb_F10_ML = NULL
    mcb_F10_MM = NULL
    
    
    print('Measuring Time for mXl, 5 Atoms'); mcb_F5_ML  = rbenchmark::benchmark({function(){microbenchmark_wrapper_Fast_NA_ML(5)}}() ,replications = RBENCHMARK_REPETITION)
    print('Measuring Time for mXm, 5 Atoms'); mcb_F5_MM  = rbenchmark::benchmark({function(){microbenchmark_wrapper_Fast_NA_MM(5)}}() ,replications = RBENCHMARK_REPETITION)
    print('Measuring Time for mXl, 10 Atoms'); mcb_F10_ML = rbenchmark::benchmark({function(){microbenchmark_wrapper_Fast_NA_ML(10)}}(),replications = RBENCHMARK_REPETITION)
    print('Measuring Time for mXm, 10 Atoms'); mcb_F10_MM = rbenchmark::benchmark({function(){microbenchmark_wrapper_Fast_NA_MM(10)}}(),replications = RBENCHMARK_REPETITION)
    
    
    print('Measuring Time for dCOV')
    mcb_dcov = NULL
    if(MODE_ADD_ATOMS){
      mcb_dcov = list(elapsed = Prev_run_times[ind_in_prev("dCOV")] * RBENCHMARK_REPETITION)
    }else{
      mcb_dcov = rbenchmark::benchmark(microbenchmark_wrapper_dcov(),replications = RBENCHMARK_REPETITION)  
    }
    
    
    print('Measuring Time for dHSIC')
    mcb_dHSIC = NULL
    if(MODE_ADD_ATOMS){
      mcb_dHSIC = list(elapsed = Prev_run_times[ind_in_prev("dHSIC")] * RBENCHMARK_REPETITION)
    }else{
      mcb_dHSIC = rbenchmark::benchmark(microbenchmark_wrapper_dHSIC(),replications = RBENCHMARK_REPETITION)  
    }
    
    
    print('Measuring Time for MIC')
    mcb_MIC = NULL
    if(MODE_ADD_ATOMS){
      mcb_MIC = list(elapsed = Prev_run_times[ind_in_prev("MIC")] * RBENCHMARK_REPETITION)
    }else{
      mcb_MIC = rbenchmark::benchmark(microbenchmark_wrapper_MIC(),replications = RBENCHMARK_REPETITION)  
    }
    
    
    #organize results and save to file
    run_times = c(
      ML_5Atoms  = mcb_F5_ML$elapsed  /(RBENCHMARK_REPETITION  * SMALL_SAMPLE_MULTIPLIER),
      MM_5Atoms  = mcb_F5_MM$elapsed  /(RBENCHMARK_REPETITION  * SMALL_SAMPLE_MULTIPLIER),
      ML_10Atoms = mcb_F10_ML$elapsed /(RBENCHMARK_REPETITION  * SMALL_SAMPLE_MULTIPLIER),
      MM_10Atoms = mcb_F10_MM$elapsed /(RBENCHMARK_REPETITION  * SMALL_SAMPLE_MULTIPLIER),
      
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

#Generate thresholds for MIC
#***********************
if(MODE_SUBSECTION_E_MIC_NULL_TABLE){
  #for each sample size
  for(current_N_i in 1:length(N_vec)){
    current_N = N_vec[current_N_i]
    print(paste0("MIC null table for N: ",current_N))
    NR.WORKERS = NR_CORES
    
    #we define a cluster of workers
    cl <- makeCluster(NR.WORKERS)
    registerDoParallel(cl)
    
    # this is the base worker funcion in each cluster
    MIC_worker_function = function(){
      
      NR_REPS_PER_WORKER = ceiling(NULL_TABLE_SIZE/NR.WORKERS)
      NULL_DIST = rep(NA,NR_REPS_PER_WORKER)
      for(i in 1:NR_REPS_PER_WORKER){
        NULL_DIST[i] = minerva::mine(1:current_N,sample(1:current_N))$MIC
      }
      return(NULL_DIST)
    }
    
    # run to get results
    MIC_null_dist_res <- foreach(i=1:NR.WORKERS, .options.RNG=1234) %dorng% { MIC_worker_function() }
    # get this as array
    MIC_null_dist_as_array = unlist(MIC_null_dist_res)
    #compute threshold
    MIC_THRESHOLD_VALUE = as.numeric(quantile(MIC_null_dist_as_array,probs = 1-alpha))
    #save threshold to file
    filename_MIC = get_filename_MIC_Threshold(current_N)
    save(MIC_THRESHOLD_VALUE,file = filename_MIC)
    stopCluster(cl)
  }
  
}

#Generate Null tables for MinP/ADP
#***********************
if(MODE_SUBSECTION_F_GENERATE_NULL_TABLES){
   #for each sample size
    for(current_N_i in 1:length(N_vec)){
      current_N = N_vec[current_N_i]
      print(paste0('Building Null Tables for N: ',current_N))
      #generate null tables for 5,10,15,30,45,60 atoms. Variants are mXl and mXm
      
      set.seed(1)
      NULL_TABLES  = list()
      NULL_TABLES[[1]]  = Fast.independence.test.nulltable(current_N,mmax = min(MMAX,5) ,nr.perm = NULL_TABLE_SIZE,nr.atoms = 5 )
      NULL_TABLES[[2]]  = Fast.independence.test.nulltable(current_N,mmax = min(MMAX,5) ,nr.perm = NULL_TABLE_SIZE,nr.atoms = 5, variant = 'ADP-EQP' )
      NULL_TABLES[[3]]  = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 10)
      NULL_TABLES[[4]]  = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 10, variant = 'ADP-EQP')
      NULL_TABLES[[5]]  = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 15)
      NULL_TABLES[[6]]  = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 15, variant = 'ADP-EQP')
      NULL_TABLES[[7]]  = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 30)
      NULL_TABLES[[8]]  = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 30, variant = 'ADP-EQP')
      NULL_TABLES[[9]]  = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 45)
      NULL_TABLES[[10]] = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 45, variant = 'ADP-EQP')
      NULL_TABLES[[11]] = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 60)
      NULL_TABLES[[12]] = Fast.independence.test.nulltable(current_N,mmax = MMAX,nr.perm = NULL_TABLE_SIZE,nr.atoms = 60, variant = 'ADP-EQP')
      
      #save results to null table file, for specific N
      PROPOSED_TESTS = c('ML_5ATOMS','MM_5ATOMS','ML_10ATOMS','MM_10ATOMS','ML_15ATOMS','MM_15ATOMS','ML_30ATOMS','MM_30ATOMS','ML_45ATOMS','MM_45ATOMS','ML_60ATOMS','MM_60ATOMS')
      names(NULL_TABLES) = PROPOSED_TESTS
      filename = get_filename_ADP_Null_Tables(current_N)
      save(NULL_TABLES,file = filename)
    }
}


#Subsection G: Run Scenarios (multiple core)
#***********************
if(MODE_SUBSECTION_G_RUN_SCENARIOS){
  #This will store the simulation results
  Simulation_results = NULL
  
  #This is the worker function, for each worker node in the cluster.
  #SCENARIO_ID - 1,2,3,4 - this is the scenario being run
  #sample_size_index - the index of the sample size used, out of N_vec
  POWER_SIM_WORKER_FUNCTION = function(SCENARIO_ID,sample_size_index){
    
    #this is a different R session, we load all the parameters we need, along with results for MIC threshold and null tables.
    # we also load the packages
    library(HHG)
    library(dHSIC)
    library(minerva)
    library(energy)
    
    current_N = N_vec[sample_size_index]
    load(file = get_filename_MIC_Threshold(current_N)) #=> MIC_THRESHOLD_VALUE
    load(file = get_filename_ADP_Null_Tables(current_N)) #=> NULL_TABLES
    
    NR_REPS_PER_WORKER = ceiling(POWER_REPETITIONS/NR.WORKERS)
    result_names = c('REPS',names(NULL_TABLES),'dCOV','dHSIC','MIC')
    results = matrix(0,nrow = 1,ncol = length(result_names))      #results for the competition
    results_MinP = matrix(0,nrow = 1,ncol = length(result_names)) #this stores results for the MinP procedure. We use two different arrays as we have a multiplier for the repetitions of data generations.     
    
    colnames(results) = result_names
    colnames(results_MinP) = result_names
    
    
    #repeat by the required number of repetitions PER WORKER
    
    for(i in 1:(NR_REPS_PER_WORKER * MULTIPLIER_FOR_MINP )){
      # generate data
      data = Scenario_list[[SCENARIO_ID]](current_N,N_vec[length(N_vec)])
      x = data[1,]
      y = data[2,]
      
      
      results_MinP[1,1] = results_MinP[1,1] + 1 # this counts the number of repetitions
      # run our tests. Remember that parameters are defined by the null tables, so we are just iterating over null tables
      for(n in 1:length(NULL_TABLES)){
        MinP_res = Fast.independence.test(x,y,NullTable = NULL_TABLES[[n]])
        results_MinP[1,n + 1] = results_MinP[1,n + 1] + 1 *(MinP_res$MinP.pvalue  <= alpha) 
      }
      
      #Run tests for competition
      if(i %% MULTIPLIER_FOR_MINP == 0){
        results[1,1] = results[1,1] + 1 # counts the number of repetitions for the competition
        
        #run dCOV
        dCov_res = energy::dcov.test(x,y,R = PERMUTATIONS_FOR_TEST)
        results[1,length(NULL_TABLES) + 2] = results[1,length(NULL_TABLES) + 2] + 1 * (dCov_res$p.value <= alpha)
        
        #run dHSIC
        dHSIC_res = dHSIC::dhsic.test(x, y, method="permutation", B = PERMUTATIONS_FOR_TEST)
        results[1,length(NULL_TABLES) + 3] = results[1,length(NULL_TABLES) + 3] + 1 * (dHSIC_res$p.value <= alpha)
        
        #run MIC
        MIC_res = mic_res = minerva::mine(x,y)
        results[1,length(NULL_TABLES) + 4] = results[1,length(NULL_TABLES) + 4] + 1 * (MIC_res$MIC >= MIC_THRESHOLD_VALUE)
      }
      
    }
    
    #normalize the tables
    results = results / results[1,1]
    results_MinP = results_MinP / results_MinP[1,1]
    ret = results + results_MinP # return the results from both rows
    ret[1,1] = 1 #normalized back to 1
    return(ret) 
  }
  
  #these are the scenarios to be run in the simulation. By default, it is all of them.
  Scenarios_to_run = 1:NR_SCENARIOS
  
  # for each of the sample sizes
  for(current_N_ind in 1:length(N_vec)){
    Simulation_results = NULL
    current_N = N_vec[current_N_ind]
    # we iterate over the scenarios
    for(SCENARIO_ID_FOR_SIM  in Scenarios_to_run){
      
      print(paste0('Running Scenario ID ',SCENARIO_ID_FOR_SIM,' N index :',current_N_ind,' Sample size: ',current_N))
      
      # we start a cluster of the required number of workers
      NR.WORKERS = NR_CORES
      cl <- makeCluster(NR.WORKERS)
      registerDoParallel(cl)
      SCENARIO_SIM_RESULTS <- foreach(i=1:NR.WORKERS, .options.RNG=1234,.combine = 'rbind') %dorng% { POWER_SIM_WORKER_FUNCTION(SCENARIO_ID_FOR_SIM,current_N_ind) }
      stopCluster(cl)
      
      #gather the results from the different workers, and save them to file
      sim_results_as_row = apply(SCENARIO_SIM_RESULTS,2,sum)
      save(sim_results_as_row,file = paste0('ind_sim_results_as_row_',SCENARIO_ID_FOR_SIM,'_N_',current_N,'.RData'))
      Simulation_results = rbind(Simulation_results,sim_results_as_row) 
    }
    # power is obtained by dividing the number of rejecetions, by the number of repetitions.
    # (note that the number of repetitions may be a bit higher than required, as it is a complete multiple of the number of cores)
    Power_results = Simulation_results
    for(i in 1:nrow(Power_results))
      Power_results[i,] = Power_results[i,] / Power_results[i,1]
    rownames(Power_results) = Scenario_names[Scenarios_to_run]
    save(Power_results,file = get_filename_Power_Results(current_N))
  }
}

#Subsection H: Analyze Results - plot power and runtime by different N
#***********************
if(MODE_SUBSECTION_H_ANALYZE_RESULTS){
  # for each of the sample sizes, we load the data
  for(current_N_index in 1:length(N_vec)){
    current_N = N_vec[current_N_index]
    
    load(file = get_filename_Power_Results(current_N)) # => Power_results
    load(file = get_filename_timemeasurements(current_N)) #=> run_times
    
    # organize the reults in a long format, readable by ggplot 2
    colnames(Power_results) = c("REPS","mXl, 5 Atoms", "mXm, 5 Atoms","mXl, 10 Atoms", "mXm, 10 Atoms","mXl, 15 Atoms", "mXm, 15 Atoms", "mXl, 30 Atoms", "mXm, 30 Atoms", "mXl, 45 Atoms", "mXm, 45 Atoms", "mXl, 60 Atoms", "mXm, 60 Atoms", "dCOV", "dHSIC" , "MIC" )
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
    power_res_plot$SE = qnorm(1-0.05/2)* sqrt(0.5*0.5/ (POWER_REPETITIONS))
    ind_MINP = which(!(power_res_plot$Test %in% c('dCOV','MIC','dHSIC')))
    power_res_plot$SE[ind_MINP] = qnorm(1-0.05/2)* sqrt(0.5*0.5/ (POWER_REPETITIONS * MULTIPLIER_FOR_MINP))
    Test_names_vec = colnames(Power_results)[c(14,15,16,2,4,6,8,10,12,3,5,7,9,11,13)]
    
    #plot to file
    pdf(get_filename_Power_Results_Graph(current_N), width = 8, height = 5)
    
    print(
      ggplot(power_res_plot) + geom_point(aes(x = log(RunTime) , y = Power ,shape = Test,color = Test), size = 0.85)+facet_wrap(~Scenario,nrow = 2,ncol = 2) +
      scale_color_manual(breaks = Test_names_vec,values = c(1,1,1,2,3,4,1,8,6,2,3,4,1,8,6)) +
      scale_shape_manual(breaks = Test_names_vec,values = c(0,1,2,16,16,16,16,16,16,17,17,17,17,17,17)) + geom_errorbar(aes(x = log(RunTime), ymin=Power-SE, ymax=Power+SE), color = 'gray56', width = 0.3,lwd = 0.25) + #'gray66' is also an option for error bar colors
      xlab("log(RunTime[Seconds])") + ylim(c(0,1)) + theme_bw() + geom_point(aes(x = log(RunTime) , y = Power ,shape = Test,color = Test), size = 0.7) + theme(legend.text=element_text(size=10) , legend.key.height=unit(1.0,"line"))
      ) 
    
    dev.off()
  }
 
}

#Subsection I: Analyze Results - plot power for all N.
#***********************
if(MODE_SUBSECTION_I_PLOT_POWER_BY_N){
  #we load data for all N's
  cols_ind_to_graph = c(10,14,15,16)
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
  power_matrix$SE = qnorm(1-0.05/2)* sqrt(0.5*0.5/ (POWER_REPETITIONS))
  ind_MINP = which(!(power_matrix$Test %in% c('dCOV','MIC','dHSIC')))
  power_matrix$SE[ind_MINP] = qnorm(1-0.05/2)* sqrt(0.5*0.5/ (POWER_REPETITIONS * MULTIPLIER_FOR_MINP))
  #plot to file
  pdf("IndPower_by_N.pdf", width = 8, height = 2.5)
  
  ggplot(data = power_matrix,mapping = aes(x= N,y = Power,color = Test)) +geom_point(size=0.25)+geom_line(lwd = 0.2) +
    facet_wrap(~Scenario,nrow =1, ncol = 4) + ylim(c(0,1)) +  theme_bw() + geom_errorbar(aes(ymin=Power-SE, ymax=Power+SE), width = 60,lwd = 0.25) +
    theme(axis.text.x = element_text(size = 6),axis.text.y = element_text(size = 6)) 
  
  dev.off()
}