# AMIAS: an R package for solving the generalized $\ell_0$ problem

`AMIAS` (Alternating Minimization Induced Active Set) is an R package that aims to solve the generalized $\ell_0$ problem. That is, given an observed vector $y$, the generalized $\ell_0$ problem minimizer the sum of squared residuals with $\ell_0$-regularization on some linear combination of the estimator. The AMIAS method is based on the necessary optimality conditions derived from an augmented Lagrangian framework. The proposed method takes full advantage of the primal and dual variables with complementary supports, and decouples the high-dimensional problem into two sub-systems on the active and inactive sets, respectively. A sequential AMIAS algorithm with warm start initialization is developed for efficient determination of the cardinality parameter, along with the output of solution paths.

The most typical example is the $\ell_0$ trend filtering, which is an effective tool for nonparametric regression with the power of automatic knot detection
in function values or derivatives. For more details, please see the paper **$\ell_0$ Trend Filtering** by C. Wen and X. Wang and A. Zhang. 


## Installation

To install the `AMIAS` R package from Github, just run:

```r
if(!require(devtools)) install.packages('devtools')
devtools::install_github("wencanhong/AMIAS")
```

## References

- Canhong Wen, Xueqin Wang, and Aijun Zhang (2022). $\ell_0$ Trend Filtering. 
