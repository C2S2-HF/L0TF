\name{my.rollmean}
\alias{my.rollmean}
\title{Specil Rolling Means}
\description{
  Generic functions for computing rolling means of ordered observations.
}
\usage{
   my.rollmean(y,h=5,...)
}
\arguments{
  \item{y}{Observe sequence, with length n.}
  \item{h}{integer width of the rolling window, \code{\link[zoo]{rollmean}} would be used with \code{k = h*2+1}.}
  \item{\dots}{Other arguments, see \code{\link[zoo]{rollmean}} for details.}
}
\value{
  \item{vector}{Vector: smooth observe sequence, with length n.}
}

\author{
  Canhong Wen, Xueqin Wang, Shijie Quan, Zelin Hong and Aijun Zhang.
}
\references{
  Wen, C., Wang, X., and Zhang, A. (2022) \emph{\eqn{l_0} trend filtering}, technique report.
}

\examples{
  library(AMIAS)

  set.seed(12580)
  n <- 100
  sigma <- 0.3
  y0 <- rep(0,n)
  y0[10:15] <- 2
  y0[40:60] <- -1
  y0[80:82] <- 4
  y <- y0 + sigma*rnorm(n)
  y[80:82] <- y0[80:82] + sigma*rnorm(3)
  my.rollmean(y,1)
}
