\name{coef.amias}
\alias{coef.amias}
\alias{coef.samias}
\title{
  Extract cefficients from a amias or samias object
}
\description{
  This function extracts coefficients from a amias or samias object.
}
\usage{
\method{coef}{amias}(object, type=c("coef","primal", "dual", "active"), ...)
\method{coef}{samias}(object, type = c("coef","primal", "dual", "active"), k, ...)
}
\arguments{
  \item{object}{A list with class attribute 'amias' or 'samias'.}
  \item{type}{A character string, one of \code{"coef"}, \code{"primal"}, \code{"dual"}, or \code{"active"}, indicating whether generalized \eqn{l_0} coefficients, primal coefficients, dual coefficients, or active sets, should be returned. Default is \code{"coef"}, which corresponds to the solution of the original problem.
  }
  \item{k}{An integer vector of numbers of knots at which coefficients should be calculated. The user can choose to a subset of index from \eqn{1,2,\dots, kmax}; none are specified, then it returns the optimal coefficient determined by miniming the BIC value.}
  \item{...}{
    additional arguments passed to coef.
  }
}
\value{
  Returns a list with the following components:
    \item{alpha}{If the type is \code{"coef"}, a matrix containing the generalized \eqn{l_0} coefficients, each column corresponding to a value of \code{k}.}
    \item{v}{If the type is \code{"primal"}, a matrix containing the primal variables, each column corresponding to a value of \code{k}.}
   \item{u}{If the type is \code{"dual"}, a matrix containing the dual variables, each column corresponding to a value of \code{k}.}
   \item{A}{If the type is \code{"active"}, a list containing the active sets, each sublist corresponding to a value of \code{k}.}
  \item{k}{
    An integer that specifies the optimal \code{k} or an integer vector containing the sequence of \code{k}
    values corresponding to the columns of \code{beta}, \code{v}, and \code{u}. 
  }
}
\seealso{
  \code{\link{amias}}, \code{\link{samias}}.
}
\examples{
  # A simple toy example
  set.seed(0)
  n <- 100
  alpha0 <- rep(sample(1:10,5), each = n/5)
  y <- alpha0 + rnorm(n, sd = 0.8)
  
  fit <- amias(y, k = 4) 
  coef(fit)
  
  fit <- samias(y, kmax = 6)
  coef(fit)                         # get the coefficients with the minimum BIC value
  coef(fit, k = 3:5)                # get the coefficients that number of knots being 3, 4, and 5 
  
  coef(fit, type="active")          # get the active set with minimum BIC value
  coef(fit, type="active", k = 3:5) # get the active sets that number of knots being 3, 4, and 5 
}
\keyword{methods}
