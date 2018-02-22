# HHG_large_sample_framework
This Git repository houses the Code and framework for the simulations and results of the paper regarding the nonparametric test of independence  given in the function Fast.independence.test(...). This test is computationally efficient for large sample sizes.

## Directory structure

The Script files are as follows:

* brill-heller-heller.R - This is the index file running all other files. This script also shows the usage examples detailed in the paper and used to produce the example outputs for the test. The script also graphs the ADP plot used for the description of the algorithm used in the paper. 

* Power_Sim_Ind.R - This simulation script runs the power simulations for large sample sizes, comparing the $MinP$/$S^{m\times l}$ test to other tests of independence available on CRAN. The script produces graphs for the different setting, along with run time and power measurements for different sample sizes. Plots produced for results include power by logarithm of runtime, for the different scenarios (and sample sizes)  and power by different sample sizes.

* Power_Sim_KS.R - This simulation script runs the power simulations for large sample sizes, for the K-sample tets. The script plots the different scenarios, along with a graph for power by logarithm of runtime.

* TimeMeasurements_Comparison_To_Other_Methods.R - The script runs simulations measuring times for the different tests for different sample sizes. An additional simulation measures the run time of the test of a different number of Atoms (method resolution parameter, setting a tradeoff between runtime and computational complexity)

* PowerSimulation_SmallN.R - This script performs a power simulation for small sample sizes ($N = 10, 20, 30, 40$), comparing the power of the method to the spearman test.

Other files in the directory are the simulation results and graphs

## Reproducability
All simulation results and graphs used in producing the paper, are available in this repository.
Simulation flags (in each file) are configured such that if scripts are run, long simulations (taking hours/ days) are loaded from memory (and not rerun), All graphs and figures are fully reproduced.
To rerun simualtions, change the 'MODE_XXX' flags at the top of each file to TRUE.

Power simulations were run on a C5.18XL amazon machine with 72 cores. The number of cores in each simulation is set according to the machine parameters it is run on. To fully reproduce paper result, please set number of cores to 72. Running times were measured on a single core.
On a different machine, running times may differ.

## Further Reading

* A detailed vignette is available in our CRAN packge. Install the latest packge using 'install.packages('HHG')' in R.
The vignette is available using 'vignette(HHG")'

* Further simulations and theory regarding the method for nonparametric independence testing, is available in papers published:

1. Heller, R., Heller, Y., and Gorfine, M. (2013). A consistent multivariate test of association based on ranks of distances. Biometrika, 100(2), 503-510.

2. Heller, R., Heller, Y., Kaufman S., Brill B, & Gorfine, M. (2016). Consistent Distribution-Free K-Sample and Independence Tests for Univariate Random Variables, JMLR 17(29):1-54 http://www.jmlr.org/papers/volume17/14-441/14-441.pdf

3. Brill B. (2016) Scalable Non-Parametric Tests of Independence (master's thesis). http://primage.tau.ac.il/libraries/theses/exeng/free/2899741.pdf
