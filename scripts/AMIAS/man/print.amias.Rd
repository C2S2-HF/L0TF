\name{print.amias}
\alias{print.amias}
\alias{print.samias}
\title{Print an amias or a samias object}
\description{
  Print a summary of the fitted \code{"amias"} or \code{"samias"} object.
}
\usage{
\method{print}{amias}(object, ...)
\method{print}{samias}(object, digits = max(3, getOption("digits") - 3), ...)

}
\arguments{
  \item{object}{A list with class attribute 'amias' or 'samias'.}
  \item{digits}{significant digits in printout}
  \item{\dots}{Additional print arguments.}
}

\seealso{
  \code{\link{amias}}, \code{\link{samias}}.
}
\details{
 For an \code{"amias"} object, the call that produced the object \code{x} is printed, followed by a summary of the output.

  For a \code{"samias"} object, the call that produced the object \code{x} is printed, followed by a three-column matrix with columns \code{k}, \code{Df} and \code{\%BIC}. The \code{k} column is the number of knots and the \code{Df} column is the degree of freedom defined as \eqn{Df = k+q+1}. \code{BIC} is the BIC value.

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
  print(fit)

  fit <- samias(y, kmax = 6)
  print(fit)


}
