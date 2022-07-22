# L0TF
This folder contains all the necesaary codes for duplicating the results in the  paper: 
**L_0 Trend Filtering
**Canhong Wen, Xueqin Wang, Aijun Zhang

Last update: 2022-07-11

#### __code__ folder contains the following files in R language:
* SimuL0TF.Rmd: generate Figures 1-6 and include some more illustrative simulated examples.
* AlgoAnalysis.Rmd: generate Figures 7-14 in Section 4.1.
* utils.R and amiasutils.R: source codes used in AlgoAnalysis.Rmd.
* simulations.R: replicate the results in Section 4.2 except for the l0-MIP.
* nsimul0tfc.R and tsimul0tfc.R: replicate the results of l0-MIP in Section 4.2.
* RealData.R: replicate the results and generate all the graphs in Section 4.3.

#### __figs__ folder contains all the figures we generate by running the R code in the __code__ folder.

#### __Rda__ folder contains the output RData by runing the R code in the __code__ folder.

#### __data__ folder contains the data in the real data application in Section 4.3.

#### AMIAS_1.0.3.tar.gz contains the source file for implementing the algorithm. After downloading it, you need to run the following code in R to install it.

install.packages("Your_download_path/AMIAS_1.0.3.tar.gz", repos = NULL)
