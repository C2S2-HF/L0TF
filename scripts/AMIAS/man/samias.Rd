\name{samias}
\alias{samias}
\title{Solve Generalized \eqn{l_0} Problem with SAMIAS Method}
\description{
  This function solves the generalized \eqn{l_0} problem with a sequential list of numbers of knots. The generalized \eqn{l_0} coefficient is computed via the sequential alternating minimization induced active set (SAMIAS) algorithm. Can deal with any polynomial order or even a general penalty matrix for structural filtering.
}
\usage{
samias(y, D = NULL, D_type = c("tf0", "tfq", "user"), q = 0, kmax = min(20,m-1),
    rho = n^(q+1), tmax = 10, eps = 0, smooth = FALSE, h = 5, adjust = FALSE,
    delta = 10, adjust.max = 10, ...)
}
\arguments{
  \item{y}{Observed data, of length n.}
  \item{D}{Penalty matrix of size m x n (see Details). Default is \code{NULL}. }
  \item{D_type}{Types of \eqn{D}. Either "tf0", "tfq" or "user", depending on what types of contraint that user wants to impose on the data.  For \code{D_type = "tfq"} and \code{D = NULL}, we solve the \eqn{l_0} trend filtering of order \eqn{q+1}, where \eqn{q} is determined by the argument \code{q}. For \code{D_type = "tf0"} and \code{D = NULL}, we solve the \eqn{l_0} trend filtering of order 0, i.e., with piecewise constant constraint on the data. For \code{D_type = "user"}, the penalty matrix \code{D} must be specified.}
  \item{q}{Nonnegative integer used to specify the order of penalty matrix. See \code{\link{genDtf1d}} for Details.}
  \item{kmax}{The maximum number of knots. The number of knots is thus \eqn{{1,2,\dots,kmax}}. Default is \eqn{min(20,m-1)}.}
   \item{rho}{The hyperparameter \eqn{\rho} in the augmented Lagrangian of the \eqn{l_0} trend filtering problem. Default is \eqn{n^{q+1}}.}
  \item{tmax}{The maximum number of iterations in the AMIAS algorithm.}
  \item{eps}{The tolerance \eqn{\epsilon} for an early stoping rule in the SAMIAS algorithm. The early stopping rule is defined as \eqn{\|y- \alpha\|_2^2/n \le \epsilon}, where \eqn{\alpha} is the fitting coefficient.}
  \item{smooth}{Whether to smooth the data, if \code{TRUE}, it smoothes the input data. Default is FALSE.}
  \item{h}{Bandwidth in smoothing data. See \code{\link{my.rollmean}} for details.}
  \item{adjust}{Whether to adjust the indexes of the active set in the AMIAS algorithm. If \code{TRUE}, it implements the adjustment when the indexes in the active set are not well separated. Default is FALSE.}
  \item{delta}{The minimum gap etween the adjacent knots. Only used when \code{adjust = TRUE}.}
  \item{adjust.max}{The number of iterations in the adjustment step. Only used when \code{adjust = TRUE}.}
  \item{\dots}{Other arguments.}
}
\details{
  The generalized \eqn{l_0} problem with a maximal number of possible knots \eqn{kmax} is
  \deqn{\min_\alpha 1/2 \|y-\alpha\|_2^2 \;\;{\rm s.t.}\;\; \|D\alpha\|_0 \le kmax.}
  The  penalty matrix \eqn{D} is either the discrete difference operator of order \eqn{q + 1} or a general matrix specified by users to enforce some special structure.  We solve this problem by sequentially considering the following sub-problem.
    \deqn{\min_\alpha 1/2 \|y-\alpha\|_2^2 \;\;{\rm s.t.}\;\; \|D\alpha\|_0 = k}
for \eqn{k} ranging from 1 to \eqn{kmax}. This sub-problem can be solved via the AMIAS method, see \code{\link{amias}}. Thus a sequential AMIAS algorithm named SAMIAS is proposed, which is a simple extension of AMIAS by adopting the warm start strategy.

Since it outputs the whole sequence of results for k = 1, . . . , kmax, we can naturally draw the solution path for increasing k and apply the Bayesian information criterion (BIC) for hyperparameter determination. The BIC is defined by \deqn{n log(mse) + 2*log(n)*df,}
where \eqn{mse} is the mean squared error, i.e., \eqn{\|y-\alpha\|_2^2/n}, and \eqn{df = k+q+1}.



}
\value{
A list with class attribute 'samias' and named components:
  \item{call}{The call that produces this object.}
  \item{y}{Observe sequence, if smooth, the smooth y will be returned.}
  \item{D_type}{Types of \eqn{D}.}
  \item{q}{The order of penalty matrix.}
  \item{k}{The sequential list of number of knots, i.e., \eqn{{1,2,\dots,kmax}}}
  \item{alpha}{The optimal coefficients \eqn{\alpha} determined by minimizing BIC.}
  \item{v}{The optimal primal variable or split variable of the argumented lagrange form in \eqn{D\alpha}.}
  \item{u}{The optimal dual variable or lagrange operator of the argumented lagrange form in \eqn{D\alpha} for linear item.}
  \item{A}{The optimal estimate of active sets, i.e., set of detected knots.}
  \item{df}{Degree of freedom for all the candidate models.}
  \item{alpha.all}{Coefficients Matrix with \eqn{k} ranging from 1 to \code{kmax}, of dimension \eqn{n x kmax}.}
  \item{v.all}{Matrix of primal variables with \eqn{k} ranging from 1 to \code{kmax}, of dimension \eqn{(n-q-1) x kmax}.}
  \item{u.all}{Matrix of dual variables with \eqn{k} ranging from 1 to \code{kmax}, of dimension \eqn{(n-q-1) x kmax}.}
  \item{A.all}{List of active sets, of length \code{kmax}.}
  \item{bic}{The BIC value for each value of \eqn{k}.}
  \item{iter}{A vector of length \code{kmax}, record the iterations in the AMIAS algorithm.}
  \item{smooth}{Whether to smooth the data.}
}

\author{
  Canhong Wen, Xueqin Wang, Shijie Quan, Zelin Hong and Aijun Zhang.

  Maintainer: Canhong Wen <wench@ustc.edu.cn>
}
\references{
  Wen, C., Wang, X., and Zhang, A. (2022) \emph{\eqn{l_0} trend filtering}, technique report.
}
\seealso{
   \code{\link{print.samias}}, \code{\link{coef.samias}} and \code{\link{plot.samias}} method, and the \code{\link{amias}} function.
}

\examples{
  ##----- A simple toy example ------
  set.seed(0)
  n <- 100
  alpha0 <- rep(sample(1:10,5), each = n/5)
  y = alpha0 + rnorm(n, sd = 0.8)

  # With default input argument
  fit <- samias(y, kmax = 6)
  print(fit)
  plot(fit)
  plot(fit, k = 6)
  plot(fit, type="vpath")

  # Piecewise constant example
  set.seed(1)
  data <- sim.blocks(2048)

  # L_0 trend filtering with a series of number of knots.
  fit <- samias(data$y, kmax = 12)
  plot(fit, type="coef", main = "Blocks")
  plot(fit, type= "vpath", main = "Blocks")

  # Linear trend filtering example
  set.seed(1)
  data <- sim.wave(2048)

  # L_0 trend filtering with a series of number of knots.
  fit <- samias(data$y, kmax = 20, q = 1, D_type = "tfq", adjust = TRUE, delta = 30)
  plot(fit, main = "Wave")
  plot(fit, type= "vpath", main = "Wave")


  # Doppler example
  set.seed(1)
  data <- sim.doppler(2048)

  # L_0 trend filtering with a series of number of knots.
  fit <- samias(data$y, kmax = 30, q = 2, D_type = "tfq")
  plot(fit, type= "coef", main = "Doppler")
}
