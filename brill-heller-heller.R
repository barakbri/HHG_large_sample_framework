# This file contains the code for the paper by Brill, Heller, Heller
# This file is the main index file, it is used to run all other files. 
# See the flags section below, on what each part of the code base does
# date 2018-02-22
# flags for running the different sections

#################
# flags:
#################
PART_1_USAGE_EXAMPLES = TRUE #Independence tests usage cases and plots. This is the example used to create the example output.
PART_2_KSAMPLE_USAGE_EXAMPLES = TRUE #K-Sample tests usage cases and plots  . This was used to create the example for the K-Sample test.
PART_3_PLOT_ADP_BOARD = TRUE #Create the demo plot for the methods section.

PART_4_INDEPENDENCE_TESTING_SIMULATION = TRUE # Simualtion for comparing running times and power to competitors, over scenarios.
PART_5_KSAMPLE_TESTING_SIMULATION = TRUE # Simulation for comparing running times and power to competitors, over K-sample scenarios
PART_6_RUNNING_TIME_SIMULATION = TRUE #simulation for running times for different N and N.Atoms
PART_7_SMALL_N_INDEPENDENCE_SIMULATION = TRUE #Testing power for small N

#set the working directory
#setwd('~')
# it is assumed that all source files are in the same directory. Plots and results will be saved to the same directory.
# The GitHub repository saves this directory, with all results and graphs.

#Load packages required
library(HHG)
library(energy)
library(dHSIC)
library(minerva)
library(ggplot2)
library(kernlab)

# This part shows the usage examples for independence testing
if(PART_1_USAGE_EXAMPLES){
  
  #first we generate data, by noising data, from a known curve
  
  set.seed(1)
  require(ggplot2)
  #Code for the Batman curve taken from:
  #https://www.r-bloggers.com/the-batman-equation/
  #http://ygc.name/2011/08/14/bat-man/
  #Full citation also given in paper
  
  batman_curve = function(){
    
    f1 <- function(x) {
      y1 <- 3*sqrt(1-(x/7)^2)
      y2 <- -3*sqrt(1-(x/7)^2)
      y <- c(y1,y2)
      d <- data.frame(x=x,y=y)
      d <- d[d$y > -3*sqrt(33)/7,]
      return(d)
    }
    
    x1 <- c(seq(3, 7, 0.001), seq(-7, -3, 0.001))
    d1 <- f1(x1)
    x1_plot = d1$x
    y1_plot = d1$y
    p1 <- ggplot(d1,aes(x,y)) + geom_point(color="red")
    
    x2 <- seq(-4,4, 0.001)
    y2 <- abs(x2/2)-(3*sqrt(33)-7)*x2^2/112-3 + sqrt(1-(abs(abs(x2)-2)-1)^2)
    
    x3 <- c(seq(0.75,1,0.001), seq(-1,-0.75,0.001))
    y3 <- 9-8*abs(x3)
    
    x4 <- c(seq(0.5,0.75,0.001), seq(-0.75,-0.5,0.001))
    y4 <- 3*abs(x4)+0.75
    
    
    x5 <- seq(-0.5,0.5,0.001)
    y5 <- rep(2.25,length(x5))
    
    
    x6 <- c(seq(-3,-1,0.001), seq(1,3,0.001))
    y6 <- 6 * sqrt(10)/7 +
      (1.5 - 0.5 * abs(x6)) * sqrt(abs(abs(x6)-1)/(abs(x6)-1)) -
      6 * sqrt(10) * sqrt(4-(abs(x6)-1)^2)/14
    
    
    x_plot = c(x1_plot,x2,x3,x4,x5,x6)
    y_plot = c(y1_plot,y2,y3,y4,y5,y6)
    ret = list()
    ret$x = x_plot
    ret$y = y_plot
    return(ret)
  }
  
  #Create data sample
  batman = batman_curve()
  n = 1200
  sigma = 1
  ind_sample = sample(1:length(batman$x), size = n, replace = T)
  x = batman$x[ind_sample] + rnorm(n,0,sigma)
  y = batman$y[ind_sample] + rnorm(n,0,sigma)
  
  #plot to file
  PLOT_TO_FILE = T
  if(PLOT_TO_FILE){
    
    png('BatMan_plot.png',
        width     = 3.25,
        height    = 3.25,
        units     = "in",
        res       = 1200,
        pointsize = 4)  
  }
  par(bg = "black",col.axis	= "white", col.lab = "white", mar = c(2,2,1.5,1.5))
  plot(x,y,col = 'gray',pch=16,cex=0.5)
  box(col = 'white')
  points(batman$x,batman$y,col = 'yellow',cex=0.5)
  par(bg = "white",col.axis	= "black", col.lab = "black")
  
  if(PLOT_TO_FILE){
   dev.off() 
  }
  
  #generate null table, using base paramters, m.max = 10, nr.atoms = 40.
  nt                      = Fast.independence.test.nulltable(n)
  
  #run test, using null table
  res = Fast.independence.test(x,y,nt)
  res$MinP.pvalue #0.009950249
  
  #run test, same sample size, same null table
  x.2 = rnorm(1200)
  y.2 = x.2 + rnorm(1200)
  res2 = Fast.independence.test(x.2,y.2,nt)
  res2$MinP.pvalue #0.004975124
  
  #run Batman with dcov
  dcov.res = energy::dcov.test(x,y,R=1000)
  dcov.res$p.value #0.5294705
  
  #dHSIC
  dHSIC.res = dHSIC::dhsic.test(x, y, method="permutation", B = 1000)
  dHSIC.res$p.value #0.3786214
  
  #mic
  mic_res = minerva::mine(x,y)$MIC
  mine_permutations = rep(NA,1000)
  for(b in 1:1000){print(b);mine_permutations[b] = minerva::mine(x,sample(y))$MIC}
  mic_pvalue = (1+sum(mine_permutations>mic_res))/1001
  mic_pvalue #0.5324675
  
}

#Usagee examples for the KSample case
if(PART_2_KSAMPLE_USAGE_EXAMPLES){
  
  set.seed(1)
  
  #function for generating data
  sample_data = function(n1=500,n2=500){
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
  
  #plot two sample density to file
  
  PLOT_TO_FILE = T
  if(PLOT_TO_FILE){
    pdf('TwoSample-Example.pdf',height = 2)
  }
  
  plot_data = sample_data(10^6,10^6) #sample 2*10^6 points, actually show data
  dt_gg = data.frame(X = plot_data$x,Group = factor(plot_data$y,labels = c('Group 1','Group 2')))
  ggplot(dt_gg) +  geom_density(aes(X, fill = Group, colour = Group),alpha = 0.1) + theme_classic()
  
  if(PLOT_TO_FILE){
    dev.off() 
  }
  
  #sample data
  test_data = sample_data()
  X = test_data$x
  Y = test_data$y
  #generate null table
  nt = hhg.univariate.ks.nulltable(c(500,500), # This is the sample size
                                   variant = "KSample-Equipartition", #there are two variants: "KSample-Variant" for small data samples and 'KSample-Equipartition' for large data samples
                                   mmax = 10, # parameter for m.max
                                   nr.atoms = 30, #number of atoms
                                   nr.replicates = 10^4) #number of replicates in tables. This gives the minimum P-value attainable by the permutation test
  
  #run test - parameters are taken from the null table
  res = hhg.univariate.ks.combined.test(X,Y,nt)  
  res
  res$MinP.pvalue #0.0065
  
  #run kmmd, does not reject
  kmmd_res = kernlab::kmmd(matrix(X[Y==0],ncol = 1)
                                     ,matrix(X[Y==1],ncol = 1),asymptotic = TRUE,ntimes = 1000)
  kmmd_res@H0 #False
  kmmd_res@AsympH0 #False
  kmmd_res@Asymbound #0.004048766
  kmmd_res@Radbound #0.2442518
  kmmd_res@mmdstats #0.033663244 -0.001152192
  
  #run energy, does not reject
  energy_res = energy::eqdist.etest(X,sizes = c(500,500),distance = FALSE,R=1000)
  energy_res #0.4535
  
  #distribution freedom - null table can be used for any test of required size:
  set.seed(1)
  X2 = rnorm(length(Y),0,1)+1*Y
  res2 = hhg.univariate.ks.combined.test(X2,Y,nt)  
  res2
  res2$MinP.pvalue # 10^(-4)
  
}

# This Section draws the ADP board for Section 2 of the paper. The code deals mainyl with graphics, colors and lines.
# See the paper for full explanation regarding the plot.
if(PART_3_PLOT_ADP_BOARD){
  #draw the ADP board, for the methods section
  #sample data
  set.seed(1)
  N = 50
  
  x = (rnorm(N))
  y = rank(3*cos(2*pi*x) -1.5*x + rnorm(N))
  x = rank(x)
  
  PLOT_TO_FILE = T
  if(PLOT_TO_FILE){
    pdf('Atoms_Scheme.pdf',width     = 5,
            height    = 5,
            pointsize = 2)
  }
  
  #plot dots
  plot_cex = 1.3
  text_cex = 1
  plot(x,y,xaxt = 'n',yaxt='n',pch=16,cex=plot_cex,xlab='rank(x)',ylab = 'rank(y)')
  
  #plots grid lines, for atoms, partition and selected cell
  axis(1, at = c(1,(1:20)*5), las=1)
  axis(2, at = c(1,(1:20)*5), las=2)
  Atom_locations = c(1:9)*5+0.5
  Atom_gridLine_color = '#89898c'
  Atom_gridLine_selected_color = '#4466ce'#'#db3f3f'
  Atom_gridLine_selected_lwd = 1.3
  Atom_gridLine_partition_lwd = 1.3
  Atom_gridLine_partition_color = '#4466ce'#'#4466ce'#'#db8f3f'
  Atom_ind_x_left = 3
  Atom_ind_y_bottom = 6
  Atom_ind_x_right = 5
  Atom_ind_y_top = 9
  Atom_ind_additional_h_marker = c(2)
  Atom_ind_additional_v_marker = c(1,7)
  
  rect(xleft = Atom_locations[Atom_ind_x_left],
       ybottom = Atom_locations[Atom_ind_y_bottom],
       xright = Atom_locations[Atom_ind_x_right],
       ytop = Atom_locations[Atom_ind_y_top],
       col='lightblue',border = NA)#'#efddda'
  
  points(x,y,pch=16,cex = plot_cex)
  for(i in 1:length(Atom_locations)){
    if(!( i %in% c(Atom_ind_y_bottom,Atom_ind_y_top,Atom_ind_additional_h_marker)))
      abline(h=Atom_locations[i],lty=2,col=Atom_gridLine_color)
    if(!( i %in% c(Atom_ind_x_left,Atom_ind_x_right,Atom_ind_additional_v_marker)))
      abline(v=Atom_locations[i],lty=2,col=Atom_gridLine_color)
  }
  
  
  abline(h=Atom_locations[Atom_ind_y_bottom],lty=1,col=Atom_gridLine_selected_color,lwd=Atom_gridLine_selected_lwd)
  abline(h=Atom_locations[Atom_ind_y_top],lty=1,col=Atom_gridLine_selected_color,lwd=Atom_gridLine_selected_lwd)
  abline(v=Atom_locations[Atom_ind_x_right],lty=1,col=Atom_gridLine_selected_color,lwd=Atom_gridLine_selected_lwd)
  abline(v=Atom_locations[Atom_ind_x_left],lty=1,col=Atom_gridLine_selected_color,lwd=Atom_gridLine_selected_lwd)
  for(i in 1:length(Atom_ind_additional_h_marker))
    abline(h = Atom_locations[Atom_ind_additional_h_marker[i]], lty=1, col=Atom_gridLine_partition_color,lwd=Atom_gridLine_partition_lwd)
  for(i in 1:length(Atom_ind_additional_v_marker))
    abline(v = Atom_locations[Atom_ind_additional_v_marker[i]], lty=1, col=Atom_gridLine_partition_color,lwd=Atom_gridLine_partition_lwd)
  
  
  
  
  if(PLOT_TO_FILE){
    dev.off() 
  }
}

#This part deals with the power simulation for the tests of independence.
# See comments in script.
if(PART_4_INDEPENDENCE_TESTING_SIMULATION){
  source('Power_Sim_Ind.R')
}

#This part deals with the power simulation for the tests of equality of distributions (K-sample).
# See comments in script.
if(PART_5_KSAMPLE_TESTING_SIMULATION){
  source('Power_Sim_KS.R')
}

#Simulation script for measuring times by N, and by Nr.atoms
if(PART_6_RUNNING_TIME_SIMULATION){
  source('TimeMeasurements_Comparison_To_Other_Methods.R')
}

#Simulation script for testing power for small N.
if(PART_7_SMALL_N_INDEPENDENCE_SIMULATION){
  source('PowerSimulation_SmallN.R')
}