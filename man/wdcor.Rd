\name{weighted distance correlation}
\alias{wdcor}
\alias{dcor}
\title{ Weighted Distance Correlation Statistic}
\description{
 Computes weighted distance correlation statistic,
 which is a multivariate measure of dependence.
 }
\usage{
wdcor(x, y, G.list= NULL)
}
\arguments{
  \item{x}{ data of first sample}
  \item{y}{ data of second sample which needed to be weighted}
  \item{G.list}{positive integer set of the exponent on the weight}
}
\details{
\code{wdcor} compute weighted distance correlation statistics.

  The sample sizes (number of rows) of the two samples must
  agree, and samples must not contain missing values. Arguments
  \code{x}, \code{y} can optionally be \code{\link{dist}} objects;
 otherwise these arguments are treated as data.

%%---new
Weighted distance correlation is a new measure of association analysis
between two random variables introduced by Wen, Yang, Xiao, Huang and Pan (2019).
The weighted distance correlation was motivated by distance correlation, and
generalized to solve the problem of power declining rapidly of distance correlation in
ultra-high dimensional cases.

Weighted distance correlation (\code{wdCor}) is conducted by assigning positive weights to the true dependent parts and negligible weights to the remaining independent parts, and evaluating distance correlation based on the weighted variables. In this way, it can alleviate the power loss of \code{dCor} caused by the curse of dimensionality.

The empirical weighted distance correlation is defined as

\deqn{
\widehat{wdCor}(\mathbf{X},\mathbf{Y};\widehat{\omega}_{opt})=
\widehat{dCor}(\mathbf{X}, (\hat{\omega}_{opt,1}\mathbf{Y_1},\ldots,\hat{\omega}_{opt,q}\mathbf{Y_q}))
}{
\widehat wdCor (X,Y;\hat \omega_opt)=
\widehat dCor (X, (\hat \omega_{opt,1}Y_1,\ldots,\hat \omega_{opt,q}Y_q))
}

Where the empirical optimal weight can be obtained as follows.
Consider a function of \eqn{\beta_j = dCor(\mathbf{X},\mathbf{Y_j})}{\beta_j=dCor(X,Y_j)} as the weight function, and define
\eqn{\omega_j = \frac{{\beta_j}^{\gamma}}{\parallel \beta^{\gamma} \parallel}}{\omega_j = (\beta_j ^\gamma) / (||\beta^\gamma||)}. When \eqn{\gamma} traverses the positive integer \eqn{\Gamma}, the eimprical optimal empirical solution can be obtained:
\deqn{
\hat{\omega}_{opt} = \arg\mathop{\max}_{\gamma \in \Gamma} \widehat{ wdCor}(\mathbf{X},\mathbf{Y};\omega)=\arg\mathop{\max}_{\gamma \in \Gamma} \widehat{ dCor}(\mathbf{X},(\omega_1 \mathbf{Y_1},\ldots,\omega_1 \mathbf{Y_q}))
}{
\hat \omega_opt=argmax_{\gamma \in \Gamma} \widehat wdCor(X,Y;\omega)=argmax_{\gamma \in \Gamma}\widehat dCor(X,(\omega_1Y_1,\ldots,\omega_1Y_1))
}

See \code{\link{wdcor.test}} for a test of multivariate independence
based on the weighted distance covariance statistic.}

\value{
\code{wdcor} returns the sample weighted distance correlation
}

\seealso{
\code{\link{wdcor.test}}
}
\author{ Wen C.H. \email{wench@ustc.edu.cn} and
Yang Y.H., Xiao Q., Huang M.Y. and Pan W.L.
}
\examples{
 x <- matrix(c(1:9),3,3)
 y <- matrix(rnorm(9),3,3)
 set.seed(1)
 wdcor(x, y, G.list=seq(1,15,2))
 }
\keyword{ weighted distance correlation}
\keyword{ multivariate }
\concept{ independence }
\concept{ distance correlation }



