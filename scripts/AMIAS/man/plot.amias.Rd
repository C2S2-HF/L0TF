\name{plot.amias}
\alias{plot.amias}
\alias{plot.samias}
\title{Plotting method for amias and samias objects}
\description{
   The function \code{plot.amias} plots the coefficient from an \code{"amias"} object with a user-specified number of knots as a function of the input positions (which are assumed to be evenly spaced if not specified). The function \code{plot.samias} produces either a profile plot of the solution path for the primal variable or the optimal fitted coefficients (which, recall, are determined by minimizing the Bayesian Information Criterion) from a fitted \code{"samias"} object.
}
\usage{
\method{plot}{amias}(object, add.knots = TRUE,...)
\method{plot}{samias}(object, type = c("coef", "vpath"), k, add.label = TRUE, add.knots = TRUE,...)

}
\arguments{
  \item{object}{A list with class attribute 'amias' or 'samias'.}
  \item{add.knots}{Whether to add the locations of the detected knots at the top of the figure. It only works when \code{type = "opt"} for the function \code{plot.samias}.}
  \item{type}{Types of the plotting figure, either "v" or "opt".  For \code{type = "v"}, we plot the solution path of the primal variable \eqn{v} versus the cardinality \eqn{k} of knots. For \code{type = "opt"}, we plot the optimal \eqn{l_0} fitted coefficients as well as the raw data.}
  \item{k}{The number of knots. The user can choose k from \eqn{1,2,\dots, kmax}; if none are specified, then it plots the optimal coefficient determined by miniming the BIC value.}
  \item{add.label}{Whether to add labels at the right side to indicate the locations of the detected knots. It only works when \code{type = "vpath"} for the function \code{plot.samias}.}
  \item{\dots}{Other arguments,see \code{\link[graphics]{plot.default}} for details.}
}

\seealso{
  \code{\link{amias}}, \code{\link{samias}}.
}

\author{
  Canhong Wen, Xueqin Wang, Shijie Quan, Zelin Hong and Aijun Zhang.

  Maintainer: Canhong Wen <wench@ustc.edu.cn>
}
\references{
  Wen, C., Wang, X., and Zhang, A. (2022) \emph{\eqn{l_0} trend filtering}, technique report.
}


\examples{
  # A simple toy example
  set.seed(0)
  n <- 100
  alpha0 <- rep(sample(1:10,5), each = n/5)
  y = alpha0 + rnorm(n, sd = 0.8)

  fit <- amias(y, k = 4)
  plot(fit)

  fit <- samias(y, kmax = 5)
  plot(fit, type = "vpath")

  # Piecewise constant trend filtering example in Wen et al.(2018).
  set.seed(1)
  data <- sim.blocks(2048)

  fit <- samias(data$y, kmax = 12)

  plot(fit, main = "Blocks")
  plot(fit, type= "vpath", main = "Blocks")
  plot(fit, type= "vpath", add.knots = FALSE, main = "Blocks")

}
