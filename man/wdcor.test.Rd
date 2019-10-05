\name{wdcor.test}
\alias{dcor.test}
\alias{wdcor.test}
\title{ Weighted Distance Correlation test}
\description{
  Distance correlation test of multivariate independence.
  Distance correlation are multivariate measures of dependence.}
\usage{
  wdcor.test(x, y, G.list=NULL, R=NULL)
}
\arguments{
  \item{x}{ data of first sample}
  \item{y}{ data of second sample which needed to be weighted}
  \item{G.list}{positive integer set of the exponent on the weight }
  \item{R}{number of replicates}
}
\details{
  \code{wdcor.test} are nonparametric tests of multivariate independence.
  The test decision is obtained via permutation bootstrap, with \code{R} replicates.
  %%. add more
  The sample sizes (number of rows) of the two samples must
  agree, and samples must not contain missing values. Arguments
  \code{x}, \code{y} can optionally be \code{\link{dist}} objects;
  otherwise these arguments are treated as data.
  %%---new

  Weighted distance correlation is a new measure of association analysis
  between two random vectors introduced by Wen, Yang, Xiao, Huang and Pan (2020).
  The weighted distance correlation was motivated by distance correlation, and
  generalized to solve the problem of power declining rapidly of distance correlation in
  ultra-high dimensional data.

  See more details from the manual page for \code{wdcor}.
}
\value{
  \code{wdcor.test} returns a list with class \code{wdcor test} containing
  \item{     method}{ description of test}
  \item{  statistic}{ observed value of the test statistic, wdcor(x,y)}
  \item{    p.value}{ approximate p-value of the test}
  \item{ replicates}{ replicates of the test statistic}
  \item{          n}{ sample size}
  \item{  data.name}{ description of data}
}

\seealso{
  \code{\link{wdcor} }
}

\examples{
  x <- matrix(c(1:9),3,3)
  y <- matrix(rnorm(9),3,3)
  set.seed(1)
  wdcor.test(x, y, G.list=seq(1,5,2), R=199)

  ## Compare the power of wdCor and dCor in ultra-high dimension situations.
  library(MASS)
  library(energy)
  alpha = 0.05  ## Significant level
  set.seed(1)
  M = 100  ## Monte Carlo realizations times
  n = 20 ## Sample size
  q = 500  ## the dimension of Y
  q1 = 5 ## the true dependent dimension of Y
  pow1 <- pow2 <- numeric()
  for (m in 1:M){
    X <- matrix(runif(n),n,1)
    Beta <- c(rep(4,q1),rep(0,q-q1))
    E <- floor(mvrnorm(n,rep(0,q),diag(q)))
    Y <- X%*%Beta + rep(10,q) + E
    pow1[m] <- dcov.test(X,Y,R=199)$p.value<alpha
    pow2[m] <- wdcor.test(X,Y,R=199,G.list=seq(1,15,2))$p.value<alpha
  }
  print(c(power1 = mean(pow1), power2=mean(pow2)))  ## the power of dCor and wdCor
}
\keyword{ weighted distance correlation test}
\keyword{ multivariate }
\keyword{ nonparametric }
\concept{ independence }
\concept{ multivariate }
\concept{ distance correlation }



