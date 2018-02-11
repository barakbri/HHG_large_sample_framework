setwd('~')
NR_CORES = detectCores() - 1 # currently set to 7, for reproducability


library(kernlab)
library(energy)
library(HHG)
library(ggplot2)
library(parallel)
library(doParallel)
library(foreach)
library(doRNG)

set.seed(1)

#function for generating data
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

Scenario_KS_Normal = function(n1=500,n2=500){
  y = c(rep(0,n1),rep(1,n2))
  x = 0.075 * y + rnorm(n1+n2)
  return(list(x=x,y=y))
}


Scenario_KS_Normal_Contaminated = function(n1=500,n2=500,contamination = 0.3){
  y = c(rep(0,n1),rep(1,n2))
  u = sample(c(0,1),size = n1+n2,prob = c(1-contamination,contamination),replace = T)
  x = (1-u) *(0.15 * y + rnorm(n1+n2)) + u * (rnorm(n1+n2,sd = 8))
  return(list(x=x,y=y))
}

Scenario_list = list(Scenario_KS_Normal,Scenario_KS_Normal_Contaminated,Scenario_KS_Mixture)
Scenario_names = c('Normal, Shift', 'Normal, Shift with Contamination', 'Mixture, 3 Components')

#Subsection B: Declarations

NR_SCENARIOS = length(Scenario_list)
NULL_TABLE_SIZE = 1000
POWER_REPETITIONS = 2000 #Number of realizations for power evaluation
PERMUTATIONS_FOR_TEST = NULL_TABLE_SIZE #Number of permutations for tests that require permutations
MMAX_FOR_KS = 10
alpha = 0.05
N1=N2 = 2500
N = N1+N2

#MODES:
MODE_SUBSECTION_KS_C_PLOT_SETTINGS = FALSE
MODE_SUBSECTION_KS_D_GENERATE_NULL_TABLES = FALSE
MODE_SUBSECTION_KS_E_MEASURE_TIMES = FALSE
MODE_SUBSECTION_KS_F_RUN_SCENARIOS = TRUE
MODE_SUBSECTION_KS_G_ANALYZE_RESULTS = TRUE

#Subsection C: Plot Settings 

if(MODE_SUBSECTION_KS_C_PLOT_SETTINGS){
  PLOT_TO_FILE = T
  if(PLOT_TO_FILE){
    pdf('TwoSample-Settings.pdf',height = 6)
  }
  
  plot_data = Scenario_KS_Normal(10^6,10^6) #sample 2*10^6 points
  dt_gg_normal = data.frame(X = plot_data$x,Group = factor(plot_data$y,labels = c('Group 1','Group 2')),Setting = Scenario_names[1])
  plot_data = Scenario_KS_Normal_Contaminated(10^6,10^6) #sample 2*10^6 points
  dt_gg_normalcont = data.frame(X = plot_data$x,Group = factor(plot_data$y,labels = c('Group 1','Group 2')),Setting = Scenario_names[2])
  plot_data = Scenario_KS_Mixture(10^6,10^6) #sample 2*10^6 points
  dt_gg_mixture = data.frame(X = plot_data$x,Group = factor(plot_data$y,labels = c('Group 1','Group 2')),Setting = Scenario_names[3])
  dt_gg = rbind(dt_gg_normal,dt_gg_normalcont,dt_gg_mixture)
  x_plot_min = -7.5
  x_plot_max = 7.5
  dt_gg = dt_gg[dt_gg$X >= x_plot_min & dt_gg$X <= x_plot_max,]
  ggplot(dt_gg) +  geom_density(aes(X, fill = Group, colour = Group),alpha = 0.1,trim = TRUE) + theme_classic() + facet_wrap(~Setting,nrow = 3,ncol = 1) +xlim(c(-7.5,7.5))
  
  
  
  if(PLOT_TO_FILE){
    dev.off() 
  }
}

#Subsection D: Generate NUll Tables
if(MODE_SUBSECTION_KS_D_GENERATE_NULL_TABLES){
  set.seed(1)
  
  null_tables = list()
  
  null_tables[[1]] = HHG::hhg.univariate.ks.nulltable(c(N1,N2),variant = 'KSample-Equipartition', nr.atoms = 50,mmax = MMAX_FOR_KS, nr.replicates = NULL_TABLE_SIZE)
  null_tables[[2]] = HHG::hhg.univariate.ks.nulltable(c(N1,N2),variant = 'KSample-Equipartition', nr.atoms = 100,mmax = MMAX_FOR_KS, nr.replicates = NULL_TABLE_SIZE)
  
  names(null_tables) = c('Sm, 50 Atoms', 'Sm, 100 Atoms' )
  save(null_tables,file = 'KS_NULL_TABLES.RData')
}


#Subsection E: Measure Times
if(MODE_SUBSECTION_KS_E_MEASURE_TIMES){
  
  set.seed(1)
  data = Scenario_KS_Mixture(N1,N2)
  
  microbenchmark_wrapper_Sm = function(n_atoms){
    res = HHG::hhg.univariate.ks.combined.test(data$x,data$y,variant = 'KSample-Equipartition',nr.atoms = n_atoms,mmax = MMAX_FOR_KS,nr.perm = 1000)
  }
  
  microbenchmark_wrapper_KMMD = function(){
    kmmd_res = kernlab::kmmd(matrix(data$x[1:N1],ncol = 1)
                             ,matrix(data$x[(N1+1):N],ncol = 1),asymptotic = FALSE,ntimes = PERMUTATIONS_FOR_TEST)
  }
  
  microbenchmark_wrapper_ENERGY = function(){
    energy_res = energy::eqdist.etest(data$x,sizes = c(N1,N2),distance = FALSE,R=PERMUTATIONS_FOR_TEST)
  }
  
  MICROBENCHMARK_REPETITION = 100
  
  mcb_Sm_50 = rbenchmark::benchmark(microbenchmark_wrapper_Sm(50),replications = MICROBENCHMARK_REPETITION)
  mcb_Sm_100 = rbenchmark::benchmark(microbenchmark_wrapper_Sm(100),replications = MICROBENCHMARK_REPETITION)
  mcb_ENERGY = rbenchmark::benchmark(microbenchmark_wrapper_ENERGY(),replications = MICROBENCHMARK_REPETITION)
  mcb_KMMD = rbenchmark::benchmark(microbenchmark_wrapper_KMMD(),replications = MICROBENCHMARK_REPETITION)
  
  run_times = c(
    Sm_50 = mcb_Sm_50$elapsed/MICROBENCHMARK_REPETITION,
    Sm_100 = mcb_Sm_100$elapsed/MICROBENCHMARK_REPETITION,
    ENERGY = mcb_ENERGY$elapsed/MICROBENCHMARK_REPETITION,
    KMMD = mcb_KMMD$elapsed/MICROBENCHMARK_REPETITION
  )
  
  run_times
  save(run_times,file = 'SIMULATION_RUN_TIMES_KS.RData')
  
}
#Subsection F: Run Scenarios
if(MODE_SUBSECTION_KS_F_RUN_SCENARIOS){
  
  Simulation_results = NULL
  
  POWER_SIM_WORKER_FUNCTION = function(SCENARIO_ID){
    library(kernlab)
    library(energy)
    library(HHG)
    
    load(file = 'KS_NULL_TABLES.RData') #=> null_tables

    NR_REPS_PER_WORKER = ceiling(POWER_REPETITIONS/NR.WORKERS)
    result_names = c('REPS',names(null_tables),'ENERGY','KMMD')
    results = matrix(0,nrow = 1,ncol = length(result_names))     
    colnames(results) = result_names
    for(i in 1:NR_REPS_PER_WORKER){
      data = Scenario_list[[SCENARIO_ID]](N1,N2)
      x = data$x
      y = data$y
      
      results[1,1] = results[1,1] + 1
      
      for(j in 1:length(null_tables)){
        MinP_res = hhg.univariate.ks.combined.test(x,y,null_tables[[j]])
        results[1,j+1] = results[1,j+1] + 1 *( MinP_res$MinP.pvalue<=alpha)
      }

      energy_res = energy::eqdist.etest(x,sizes = c(N1,N2),distance = FALSE,R=PERMUTATIONS_FOR_TEST)
      results[1,4] = results[1,4] + 1 *( energy_res$p.value<=alpha)
      
      kmmd_res = kernlab::kmmd(matrix(x[y==0],ncol = 1)
                               ,matrix(x[y==1],ncol = 1),asymptotic = FALSE,ntimes = PERMUTATIONS_FOR_TEST)
      results[1,5] = results[1,5] + 1 *( kmmd_res@H0)
    }
    return(results)
  }
  
  Scenarios_to_run = 1:NR_SCENARIOS
  for(SCENARIO_ID_FOR_SIM  in Scenarios_to_run){
    print(paste0('Running Scenario ID ',SCENARIO_ID_FOR_SIM))
    NR.WORKERS = NR_CORES
    cl <- makeCluster(NR.WORKERS)
    registerDoParallel(cl)
    SCENARIO_SIM_RESULTS <- foreach(i=1:NR.WORKERS, .options.RNG=1234,.combine = 'rbind') %dorng% { POWER_SIM_WORKER_FUNCTION(SCENARIO_ID_FOR_SIM)}
    stopCluster(cl)
    sim_results_as_row = apply(SCENARIO_SIM_RESULTS,2,sum)
    save(sim_results_as_row,file = paste0('sim_results_as_row_',SCENARIO_ID_FOR_SIM,'.RData'))
    Simulation_results = rbind(Simulation_results,sim_results_as_row) 
  }
  
  Power_results = Simulation_results
  for(i in 1:nrow(Power_results))
    Power_results[i,] = Power_results[i,] / Power_results[i,1]
  rownames(Power_results) = Scenario_names[Scenarios_to_run]
  save(Power_results,file = paste0('KS_Power_Results.RData'))
  
}

#Subsection G: Analyze Results
if(MODE_SUBSECTION_KS_G_ANALYZE_RESULTS){
  load(file = 'KS_Power_Results.RData') # => Power_results
  load(file = 'SIMULATION_RUN_TIMES_KS.RData') #=> run_times
  
  colnames(Power_results) = c("REPS","Sm, 50 Atoms", "Sm, 100 Atoms","ENERGY","KMMD" )
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
  
  pdf('KS_PowerResults.PDF', width = 8.5, height = 2.40)
  
  ggplot(power_res_plot) + geom_point(aes(x = log(RunTime) , y = Power ,shape = Test,color = Test))+facet_wrap(~Scenario) +
    xlab("log(RunTime[Seconds])") + ylim(c(0,1)) + theme_bw()
  
  dev.off()
  
}