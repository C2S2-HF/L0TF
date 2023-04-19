\name{sim.blocks}
\alias{sim.blocks}
\alias{sim.wave}
\alias{sim.doppler}
\title{Generate simulated data}
\description{
  Generate data for simulations.
}
\usage{
  sim.blocks(n, sigma = 1)
  sim.wave(n, sigma = 1)
  sim.doppler(n, sigma = 1)

}
\arguments{
  \item{n}{The number of observations. }
  \item{sigma}{A parameter used to control the signal-to-noise ratio.}
}

\details{
da
}
\value{
A list with the following components:
  \item{y}{The simulated data of length \eqn{n}.}
  \item{y0}{The true signal of length \eqn{n}.}
  \item{knots}{The true knots in the simulated data.}
}


\author{
  Canhong Wen, Xueqin Wang, Shijie Quan, Zelin Hong and Aijun Zhang.
}
\references{
  Wen, C., Wang, X., and Zhang, A. (2022) \emph{\eqn{l_0} trend filtering}, technique report.
}

\examples{
  # Generate simulated data
  set.seed(1)
  data <- sim.blocks(2048)

  # L_0 trend filtering with given number of knots
  fit <- amias(data$y, k = 11)
  plot(fit)

}
