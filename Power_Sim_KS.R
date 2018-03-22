# Script for K-sample simulation study for power
# See brill-heller-heller.R for an index used for running all scripts.

#setwd('~') #set in main script, disabled here

library(kernlab)
library(energy)
library(HHG)
library(ggplot2)
library(parallel)
library(doParallel)
library(foreach)
library(doRNG)
NR_CORES = detectCores() - 1 # The simulation uses the maximum number of cores available.
# For reproducability, please run on an Amazon C5.18XL machine with 72 cores.
# On any machine with 72 cores (or by setting NR_CORES = 71), results would be reproducible with paper.
# On an amazon C5 machine, run times will be reproducible as well.

set.seed(1)

#Subsection A: Functions for data generation
#***********************
# Here we define the scenarios for the simulation.
#Each scenario has a parameter n, setting the sample size, and a parameter m, setting the noise level.
#Several sample sizes have a default noise level.



#mixture - 3 components
Scenario_KS_Mixture = function(n1=500,n2=500){
  delta_1 = sample(x = c(-4,0,4),replace = T,size = n1)
  delta_2 = sample(x = c(-4,0,4),replace = T,size = n2)
  sigma1 = 1
  sigma2 = 0.8
  x1 = rnorm(n1,0,sigma1) + delta_1
  x2 = rnorm(n2,0,sigma2) + delta_2
  x = c(x1,x2)
  y = c(rep(0,n1),rep(1,n2))
  return(list(x=x,y=y))
}

#normals with shift, scenario
Scenario_KS_Normal = function(n1=500,n2=500){
  y = c(rep(0,n1),rep(1,n2))
  x = 0.12 * y + rnorm(n1+n2)
  return(list(x=x,y=y))
}

#mixture, 2 components
Scenario_KS_Normal_Mixture_2_Components = function(n1=500,n2=500,contamination = 0.6){
  y = c(rep(0,n1),rep(1,n2))
  u = sample(c(0,1),size = n1+n2,prob = c(1-contamination,contamination),replace = T)
  x = (1-u) *(0.35 * y + rnorm(n1+n2)) + u * (rnorm(n1+n2,sd = 8))
  return(list(x=x,y=y))
}

Scenario_list = list(Scenario_KS_Normal,Scenario_KS_Normal_Mixture_2_Components,Scenario_KS_Mixture)
Scenario_names = c('Normal, Shift', 'Mixture, 2 Components', 'Mixture, 3 Components')

#Subsection B: Declarations
#***********************
NR_SCENARIOS = length(Scenario_list) # number of scenarios in the simulation
NULL_TABLE_SIZE = 1000 # the null table size
POWER_REPETITIONS = 2000 #Number of realizations for power evaluation
PERMUTATIONS_FOR_TEST = NULL_TABLE_SIZE #Number of permutations for tests that require permutations
MULTIPLIER_FOR_MINP = 10 #Multiplier for the number of repetitions for power evaluation, for the MINP procedure
MMAX_FOR_KS = 10 #parameter for mmax
alpha = 0.05 #Type I error level for rejection
N1=N2 = 1000 # sample size is two groups, each of size 1000
N = N1+N2 #total sample size
BENCHMARK_REPETITION = 100 #number of repetitions for measurement of running times

#MODES:
MODE_SUBSECTION_KS_C_PLOT_SETTINGS = FALSE # Subsection is used for plotting settings
MODE_SUBSECTION_KS_D_GENERATE_NULL_TABLES = FALSE # Subsection used for generating null tables
MODE_SUBSECTION_KS_E_MEASURE_TIMES = FALSE #Subsection used for meauring running times
MODE_SUBSECTION_KS_F_RUN_SCENARIOS = TRUE #Subsection used for running the power simulation, for the different scenraio (multiple cores).
MODE_SUBSECTION_KS_G_ANALYZE_RESULTS = TRUE #Subsection used for analyzing rersults and plotting them

#Subsection C: Plot Settings 
#***********************
if(MODE_SUBSECTION_KS_C_PLOT_SETTINGS){
  # we plot the scenarios by sampling two milltion obs from each setting, and then doing a kernel estimation of building the density.
  
  #sample first setting
  plot_data = Scenario_KS_Normal(10^6,10^6) #sample 2*10^6 points
  dt_gg_normal = data.frame(X = plot_data$x,Group = factor(plot_data$y,labels = c('Group 1','Group 2')),Setting = Scenario_names[1])
  #sample second setting
  plot_data = Scenario_KS_Normal_Mixture_2_Components(10^6,10^6) #sample 2*10^6 points
  dt_gg_normalcont = data.frame(X = plot_data$x,Group = factor(plot_data$y,labels = c('Group 1','Group 2')),Setting = Scenario_names[2])
  #sample third setting
  plot_data = Scenario_KS_Mixture(10^6,10^6) #sample 2*10^6 points
  dt_gg_mixture = data.frame(X = plot_data$x,Group = factor(plot_data$y,labels = c('Group 1','Group 2')),Setting = Scenario_names[3])
  #combine results
  dt_gg = rbind(dt_gg_normal,dt_gg_normalcont,dt_gg_mixture)
  
  #x_plot_min = -10
  #x_plot_max = 7.5
  #dt_gg = dt_gg[dt_gg$X >= x_plot_min & dt_gg$X <= x_plot_max,]
  
  #plot the graph, different facet for each scenario
  PLOT_TO_FILE = T
  if(PLOT_TO_FILE){
    pdf('TwoSample-Settings.pdf',height = 2.7,width = 8)
  }
 
  ggplot(dt_gg) +  geom_density(aes(X, fill = Group, colour = Group),alpha = 0.1,trim = TRUE,lwd = 0.2) + theme_classic() + facet_wrap(~Setting,nrow = 1,ncol = 3) +xlim(c(-10,10))
  
  if(PLOT_TO_FILE){
    dev.off() 
  }
}

#Subsection D: Generate NUll Tables
#***********************
if(MODE_SUBSECTION_KS_D_GENERATE_NULL_TABLES){
  set.seed(1)
  #we generate null tables for the MinP procedure, for the Sm statistics. null tables differ by number of atoms, 25,50,100,150
  null_tables = list()
  null_tables[[1]] = HHG::hhg.univariate.ks.nulltable(c(N1,N2),variant = 'KSample-Equipartition', nr.atoms = 25,mmax = MMAX_FOR_KS, nr.replicates = NULL_TABLE_SIZE)
  null_tables[[2]] = HHG::hhg.univariate.ks.nulltable(c(N1,N2),variant = 'KSample-Equipartition', nr.atoms = 50,mmax = MMAX_FOR_KS, nr.replicates = NULL_TABLE_SIZE)
  null_tables[[3]] = HHG::hhg.univariate.ks.nulltable(c(N1,N2),variant = 'KSample-Equipartition', nr.atoms = 100,mmax = MMAX_FOR_KS, nr.replicates = NULL_TABLE_SIZE)
  null_tables[[4]] = HHG::hhg.univariate.ks.nulltable(c(N1,N2),variant = 'KSample-Equipartition', nr.atoms = 150,mmax = MMAX_FOR_KS, nr.replicates = NULL_TABLE_SIZE)
  #save results to file
  names(null_tables) = c('Sm 25 Atoms','Sm 50 Atoms', 'Sm 100 Atoms','Sm 150 Atoms')
  save(null_tables,file = 'KS_NULL_TABLES.RData')
}


#Subsection E: Measure Times
#***********************
if(MODE_SUBSECTION_KS_E_MEASURE_TIMES){
  # we generate data, used for testing
  set.seed(1)
  data = Scenario_KS_Mixture(N1,N2)
  
  
  #This is a wrapper for measuring times for the MinP procedure. It receives as a parameter the number of atoms required.
  #time measurement is for the procedure as a whole: Null table construction, including permuatations and tabulation, computation of test statistic and P.value
  microbenchmark_wrapper_Sm = function(n_atoms){
    res = HHG::hhg.univariate.ks.combined.test(data$x,data$y,variant = 'KSample-Equipartition',nr.atoms = n_atoms,mmax = MMAX_FOR_KS,nr.perm = NULL_TABLE_SIZE)
  }

  #wrapper for the KMMD test  
  microbenchmark_wrapper_KMMD = function(){
    kmmd_res = kernlab::kmmd(matrix(data$x[1:N1],ncol = 1)
                             ,matrix(data$x[(N1+1):N],ncol = 1),asymptotic = TRUE,ntimes = PERMUTATIONS_FOR_TEST)
  }
  
  #wrapper for the energy test
  microbenchmark_wrapper_ENERGY = function(){
    energy_res = energy::eqdist.etest(data$x,sizes = c(N1,N2),distance = FALSE,R=PERMUTATIONS_FOR_TEST)
  }
  
  
  #run the different wrappers, and measure times.
  print('Measuring times for Sm - 25 Atoms')
  mcb_Sm_25 = rbenchmark::benchmark(microbenchmark_wrapper_Sm(25),replications = BENCHMARK_REPETITION)
  print('Measuring times for Sm - 50 Atoms')
  mcb_Sm_50 = rbenchmark::benchmark(microbenchmark_wrapper_Sm(50),replications = BENCHMARK_REPETITION)
  print('Measuring times for Sm - 100 Atoms')
  mcb_Sm_100 = rbenchmark::benchmark(microbenchmark_wrapper_Sm(100),replications = BENCHMARK_REPETITION)
  print('Measuring times for Sm - 150 Atoms')
  mcb_Sm_150 = rbenchmark::benchmark(microbenchmark_wrapper_Sm(150),replications = BENCHMARK_REPETITION)
  print('Measuring times for Energy')
  mcb_ENERGY = rbenchmark::benchmark(microbenchmark_wrapper_ENERGY(),replications = BENCHMARK_REPETITION)
  print('Measuring times for kmmd')
  mcb_KMMD = rbenchmark::benchmark(microbenchmark_wrapper_KMMD(),replications = BENCHMARK_REPETITION)
  
  #save run times to disk
  run_times = c(
    Sm_25 = mcb_Sm_25$elapsed/BENCHMARK_REPETITION,
    Sm_50 = mcb_Sm_50$elapsed/BENCHMARK_REPETITION,
    Sm_100 = mcb_Sm_100$elapsed/BENCHMARK_REPETITION,
    Sm_150 = mcb_Sm_150$elapsed/BENCHMARK_REPETITION,
    ENERGY = mcb_ENERGY$elapsed/BENCHMARK_REPETITION,
    KMMD = mcb_KMMD$elapsed/BENCHMARK_REPETITION
  )
  
  run_times
  save(run_times,file = 'SIMULATION_RUN_TIMES_KS.RData')
  
}


#Subsection F: Run Scenarios (multiple cores)
#***********************
if(MODE_SUBSECTION_KS_F_RUN_SCENARIOS){
  #Place holder for simulation results
  Simulation_results = NULL
  
  #This is the worker function run on each of the worker nodes
  # is receives at a parameter the scenario ID to be run
  POWER_SIM_WORKER_FUNCTION = function(SCENARIO_ID){
    #this is a new R session, we load all packages and null tables that we need
    library(kernlab)
    library(energy)
    library(HHG)
    
    #load null tables
    load(file = 'KS_NULL_TABLES.RData') #=> null_tables

    NR_REPS_PER_WORKER = ceiling(POWER_REPETITIONS/NR.WORKERS)
    result_names = c('REPS',names(null_tables),'ENERGY','KMMD')
    results = matrix(0,nrow = 1,ncol = length(result_names))      #this stores results for the competition
    results_MinP = matrix(0,nrow = 1,ncol = length(result_names)) #this stores results for the MinP procedure. We use two different arrays as we have a multiplier for the repetitions of data generations.     
    colnames(results) = result_names
    colnames(results_MinP) = result_names
    
    #for each of the required repetitions, PER WORKER
    for(i in 1:( NR_REPS_PER_WORKER * MULTIPLIER_FOR_MINP )){
      #we generate data for the required scenario
      data = Scenario_list[[SCENARIO_ID]](N1,N2)
      x = data$x
      y = data$y
      
      #we increment the number of tries by 1
      results_MinP[1,1] = results_MinP[1,1] + 1
      
      #We iterate over the null tables we generated, and test for rejection. note that all test parameters are stored in the null table.
      for(j in 1:length(null_tables)){
        MinP_res = hhg.univariate.ks.combined.test(x,y,null_tables[[j]])
        results_MinP[1,j+1] = results_MinP[1,j+1] + 1 *( MinP_res$MinP.pvalue<=alpha)
      }
      if(i %% MULTIPLIER_FOR_MINP == 0){
        results[1,1] = results[1,1] + 1
        # We test using the energy statistic
        energy_res = energy::eqdist.etest(x,sizes = c(N1,N2),distance = FALSE,R=PERMUTATIONS_FOR_TEST)
        results[1,length(result_names)-1] = results[1,length(result_names)-1] + 1 *( energy_res$p.value<=alpha)
        
        
        # we test using the kmmd test
        kmmd_res = kernlab::kmmd(matrix(x[y==0],ncol = 1)
                                 ,matrix(x[y==1],ncol = 1),asymptotic = TRUE,ntimes = PERMUTATIONS_FOR_TEST)
        results[1,length(result_names)] = results[1,length(result_names)] + 1 *( kmmd_res@AsympH0)  
      }
      
    }
    #normalize the tables
    results = results / results[1,1]
    results_MinP = results_MinP / results_MinP[1,1]
    ret = results + results_MinP # return the results from both rows
    ret[1,1] = 1 #normalized back to 1
    return(ret) 
  }
  
  #these are the scenarios to be run. By default, we run all scenarios
  Scenarios_to_run = 1:NR_SCENARIOS
  
  #for each of the required scenarios:
  for(SCENARIO_ID_FOR_SIM  in Scenarios_to_run){
    print(paste0('Running Scenario ID ',SCENARIO_ID_FOR_SIM))
    
    #we start a cluster of R workers on the machine
    NR.WORKERS = NR_CORES
    cl <- makeCluster(NR.WORKERS)
    registerDoParallel(cl)
    SCENARIO_SIM_RESULTS <- foreach(i=1:NR.WORKERS, .options.RNG=1234,.combine = 'rbind') %dorng% { POWER_SIM_WORKER_FUNCTION(SCENARIO_ID_FOR_SIM)}
    stopCluster(cl)
    sim_results_as_row = apply(SCENARIO_SIM_RESULTS,2,sum)
    #we save results
    save(sim_results_as_row,file = paste0('sim_results_as_row_',SCENARIO_ID_FOR_SIM,'.RData'))
    Simulation_results = rbind(Simulation_results,sim_results_as_row) 
  }
  #we convert the number of rejections to power.
  # the number of actial repetitions performed might be larger than required, since number of actual repetitions is a complete multiple of number of workers
  Power_results = Simulation_results
  for(i in 1:nrow(Power_results))
    Power_results[i,] = Power_results[i,] / Power_results[i,1]
  rownames(Power_results) = Scenario_names[Scenarios_to_run]
  
  #we save the power results to disk
  save(Power_results,file = paste0('KS_Power_Results.RData'))
  
}

#Subsection G: Analyze Results
#***********************
if(MODE_SUBSECTION_KS_G_ANALYZE_RESULTS){
  #we load power results and running times
  load(file = 'KS_Power_Results.RData') # => Power_results
  load(file = 'SIMULATION_RUN_TIMES_KS.RData') #=> run_times
  
  #we organize data in a long format, as preffered by GGPLOT2
  colnames(Power_results) = c("REPS","25 Atoms","50 Atoms", "100 Atoms","150 Atoms","ENERGY","KMMD")
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
  Test_names_vec = colnames(Power_results)[c(2,3,4)]
  
  power_res_plot$Scenario
  
  power_res_plot$Scenario <- factor(as.character(power_res_plot$Scenario), levels = Scenario_names)
  power_res_plot$Test <- factor(as.character(power_res_plot$Test), levels = colnames(Power_results)[-1])
  #we plot the result, log run time by power. Each scenario is given a different facet.
  pdf('KS_PowerResults.PDF', width = 8.5, height = 2.40)
  
  ggplot(power_res_plot) + geom_point(aes(x = log(RunTime) , y = Power ,shape = Test,color = Test)) + facet_wrap(~Scenario) +
    xlab("log(RunTime[Seconds])") + ylim(c(0,1)) + theme_bw()
  
  dev.off()
  
}