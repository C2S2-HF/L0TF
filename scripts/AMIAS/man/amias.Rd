\name{amias}
\alias{amias}
\title{Solve Generalized \eqn{l_0} Problem with AMIAS Method}
\description{
  The function solves the Generalized \eqn{l_0} problem with a user-specified number of knots. The generalized \eqn{l_0} coefficient is computed via the alternating minimization induced active set (AMIAS) algorithm. Can deal with any polynomial order or even a general penalty matrix for structural filtering.
}
\usage{
amias(y, D = NULL, D_type = c("tf0", "tfq", "user"), q = 0, k = 3,
    rho = n^(q+1), tmax = 10, A = NULL, smooth = FALSE, h = 5, adjust = FALSE,
    delta = 10, adjust.max = 10,...)
}
\arguments{
  \item{y}{Observed data, of length n.}
  \item{D}{Penalty matrix of size m x n (see Details). Default is \code{NULL}. }
  \item{D_type}{Types of \eqn{D}. Either "tf0", "tfq" or "user", depending on what types of contraint that user wants to impose on the data.  For \code{D_type = "tfq"} and \code{D = NULL}, we solve the \eqn{l_0} trend filtering of order \eqn{q+1}, where \eqn{q} is determined by the argument \code{q}. For \code{D_type = "tf0"} and \code{D = NULL}, we solve the \eqn{l_0} trend filtering of order 0, i.e., with piecewise constant constraint on the data. For \code{D_type = "user"}, the penalty matrix \code{D} must be specified.}
  \item{q}{Nonnegative integer used to specify the order of penalty matrix. See \code{\link{genDtf1d}} for Details.}
  \item{k}{The number of knots. It is also the size of active set \code{A} or the number of equality satisfied in the contraint \eqn{D\beta = 0}. Default is 3.}
  \item{rho}{The hyperparameter \eqn{\rho} in the augmented Lagrangian of the \eqn{l_0} trend filtering problem. Default is \eqn{n^{q+1}}.}
  \item{tmax}{The maximum number of iterations in the AMIAS algorithm.}
  \item{smooth}{Whether to smooth the data, if \code{TRUE}, it smoothes the input data. Default is FALSE.}
  \item{h}{Bandwidth in smoothing data. See \code{\link{my.rollmean}} for details.}
  \item{A}{Initialization for the active set. Default is NULL, which corresponds to an empty set.}
  \item{adjust}{Whether to adjust the indexes of the active set in the AMIAS algorithm. If \code{TRUE}, it implements the adjustment when the indexes in the active set are not well separated. Default is FALSE.}
  \item{delta}{The minimum gap etween the adjacent knots. Only used when \code{adjust = TRUE}.}
  \item{adjust.max}{The number of iterations in the adjustment step. Only used when \code{adjust = TRUE}.}
  \item{\dots}{Other arguments.}
}
\details{
  The generalized \eqn{l_0} problem with a user-specified number of knots \eqn{k} is
  \deqn{\min_\alpha 1/2 \|y-\alpha\|_2^2 \;\;{\rm s.t.}\;\; \|D\alpha\|_0 = k.}
  The  penalty matrix \eqn{D} is either the discrete difference operator of order \eqn{q + 1} or a general matrix specified by users to enforce some special structure.

  The generalized \eqn{l_0} problem is solved via the alternating minimization induced active set (AMIAS) method, which is proposed by Wen et al. (2018). In AMIAS, by the augmented Lagrangian formulation with a variable splitting technique, we consider the associated augmented Lagrangian:
  \deqn{\min_\alpha 1/2 \|y-\alpha\|_2^2 + u^T(D\alpha-v) + \rho/2\|D\alpha-v\|_2^2;\;{\rm s.t.}\;\; \|v\|_0 = k.}
Based on this formulation, we derive the necessary optimality conditions through alternating minimization and coordinatewise hard thresholding. The conditions are based on primal variable \eqn{v} and dual variable \eqn{u} with complementary supports, namely active and inactive sets (denoted by \eqn{A} and \eqn{I}). Given \eqn{{A, I}}, the KKT conditions on the high dimensional variables are decoupled into two sub-systems with elegant linear algebraic solutions. For more details, please see Wen et al. (2019).


 When \eqn{D} is the discrete difference operator of order \eqn{q + 1}, it reduces to the \eqn{l_0} trend filtering problem, which will produce a piecewise \eqn{q}-th polynomial curve with automatically identified knots. In this case, the penalty matrix is a banded matrix and thus banded Cholesky decomposition is applied to accelerate the computation. In particular, when \eqn{q = 0}, it produces a piecewise constant estimate with automatically detected change points. We use a fast and memory-saving strategy to further improve the computational efficiency.

 To prevent the adjacent detected knots not to be too close, an adjustment step can be added in each step of the AMIAS method. When an adjustment is used, users need to specify the minimum possible gap between the adjacent knots and the maximum number of iterations in one adjustment step.

}
\value{
A list with class attribute 'amias' and named components:
  \item{call}{The call that produces this object.}
  \item{y}{Observed sequence, if smooth, the smooth y will be returned.}
  \item{D_type}{Types of \eqn{D}.}
  \item{q}{The order of penalty matrix.}
  \item{k}{The number of knots.}
  \item{alpha}{The fitting coefficients \eqn{\alpha}.}
  \item{v}{The primal variable or splitting variable of the argumented lagrangian form in \eqn{D\alpha}.}
  \item{u}{The dual variable or lagrangian operator of the argumented lagrangian form in \eqn{D\alpha} for linear item.}
  \item{A}{The final estimate of active set, i.e., set of the detected knots.}
  \item{df}{Degree of freedom of the fitting model, which is defined as \eqn{df = k+q+1}.}
  \item{iter}{The iterations used.}
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
   \code{\link{print.amias}}, \code{\link{coef.amias}} and \code{\link{plot.amias}} method, and the \code{\link{samias}} function.
}

\examples{

  ##----- A simple toy example ------
  set.seed(0)
  n <- 100
  alpha0 <- rep(sample(1:10,5), each = n/5)
  y = alpha0 + rnorm(n, sd = 0.8)

  # With default input argument
  fit <- amias(y, k = 4)
  plot(fit)
  coef(fit)                # extract the fitting coefficients
  coef(fit, type="active") # extract the detected knots

  # with user-specified initialized active set
  fit <- amias(y, k = 4, A = seq(20, 80, length.out = 4))
  plot(fit)

  # With default input argument
  fit <- amias(y, k = 4, adjust = TRUE)
  plot(fit)


  ##------ Piecewise constant trend filtering example in Wen et al.(2018).-----
  set.seed(1)
  data <- sim.blocks(2048, sigma=0.5)
  fit <- amias(data$y, k = 11, adjust = TRUE, delta = 20)
  plot(fit)
  lines(seq(1/2048,1,length.out = 2048), data$y0, type="s") # Add the true signal for reference


  ##------ Piecewise linear trend filtering example in Wen et al.(2018).-----
  set.seed(1)
  n <- 512
  data <- sim.wave(n)
  fit <- amias(data$y, k = 8, q = 1, D_type = "tfq", A = sample(n,8), adjust = TRUE)
  plot(fit, main = "Wave")


  ##------ Doppler example in Wen et al.(2018).-----
  set.seed(1)
  n <- 2048
  data <- sim.doppler(n)
  fit <- amias(data$y, k = 30, q = 2, D_type = "tfq", A = sample(n,30), adjust = TRUE)
  plot(fit, main = "Doppler")

}
