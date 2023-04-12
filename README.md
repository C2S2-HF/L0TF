[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# $\ell_0$ Trend Filtering

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and R scripts in this repository are a snapshot of the software and code that were used in the research reported on in the paper 
[$\ell_0$ Trend Filtering](https://doi.org/10.1287/ijoc.2019.0000) by C. Wen and X. Wang and A. Zhang. 

**Important: This code is being developed on an on-going basis at 
https://github.com/C2S2-HF/L0TF. Please go there if you would like to
get a more recent version or would like support**

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2019.0000

https://doi.org/10.1287/ijoc.2019.0000.cd

Below is the BibTex for citing this snapshot of the respoitory.

```
@article{wenl0trend,
  author =        {Canhong Wen, Xueqin Wang, Aijun Zhang},
  publisher =     {INFORMS Journal on Computing},
  title =         {$\ell_0$ Trend Filtering},
  year =          {2023},
  doi =           {10.1287/ijoc.2019.0000.cd},
  url =           {https://github.com/INFORMSJoC/2021.0313},
}  
```

## Description

The goal of this repository is to share softare and R scripts of our paper *$\ell_0$ Trend Filtering*. Our motivation is to present our code and results in a reproducible way and facilitate the coding effort of thos who want to run further experiments or improve our model.

## Repository Structure

#### [code](code) folder contains the following files in R language:
* SimuL0TF.Rmd: generate Figures 2-7 and include some more illustrative simulated examples.
* AlgoAnalysis.Rmd: replicate the reuslts and generate Figures 8-11 in Section 4.1.
* utils.R and amiasutils.R: source codes used in AlgoAnalysis.Rmd.
* RealData.R: replicate the results and generate all the graphs in Section 4.3.
* AlgoAnalysis_APP.Rmd: replicate the reuslts and generate Figures B.1-B.7 in Appendix B.1.
* [simu](simu) folder contains R scripts used in Appendix B.2: 
    * nsimu.R and tsimu.R: replicate the results for all methods except for the l0-MIP with large sample size.
    * nsimul0tfc.R and tsimul0tfc.R: replicate the results of the l0-MIP method when sample size is large.
* [simu_plots](simu_plots) folder contains the R Scripts used to generate Figures B.9-B.20 in Appendix B.2.
    * post_plot.R: generate Figures B.9, B.11 and B.13. 
    * pre_plot.R: generate Figures B.10, B.12 and B.14.
    * post_tplot.R: generate Figures B.15, B.17 and B.19. 
    * pre_tplot.R: generate Figures B.16, B.18 and B.20. 
    * combine_RData.R: combine the RData and needed to be run before generating all the figures.
    
#### [data](data) folder contains the data in the real data application in Section 4.2. Please see [spreadsheet file](data/air_hourly.csv) to view the data.

#### [AMIAS_1.0.3.tar.gz](AMIAS_1.0.3.tar.gz) contains the source file of the R package for implementing the AMIAS algorithm proposed in our paper. After downloading it, you need to run the following code in R to install it.

    install.packages("Your_download_path/AMIAS_1.0.3.tar.gz", repos = NULL)


## Support

For support in using the scripts, you can reach the authors by email wench@ustc.edu.cn.




