\name{Spgr}
\alias{Spgr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate subgroups without repeated measures with both z and x}
\description{
Estimate subgroups based on the alternating direction method of multipliers (ADMM) algorithm and smoothly clipped absolute deviation (SCAD) penalty in pairwise penalty.
}
\usage{
Spgr(y, z, x, weights, betam0, ...)
}

\arguments{
  \item{y}{ response vector with n observations .}
  \item{z}{ explanatory variables matrix which have the same coefficients among individuals.}
  \item{x}{ explanatory variables matrix with p columns which have individual coefficients.}
  \item{weights}{weights vector associated with pairwise penalties with length n(n-1)/2 .}
  \item{betam0}{initial values matrix with n rows and p columns for beta.}
  \item{nu}{penalty parameter in, default value is 1.}
  \item{gam}{parameter in SCAD penalty, default value is 3.}
  \item{lam}{tuning parameter, default value is 0.5.}
  \item{maxiter}{maximum number of iterations, default value is 1000.}
  \item{tolabs}{absolute tolerance, default value is 1e-4. }
  \item{tolrel}{relative tolerance, default value is 1e-2. }
}
\details{
The model is \deqn{y_i = z_i \eta + x_i \beta_i + \epsilon_i.} The details of absolute tolerance and relative tolerance can be found in Boyd, S., Parikh, N., Chu, E., Peleato, B., and Eckstein, J. (2011). Distributed optimization and statistical learning via the alternating direction method of multipliers.
}
\value{
\item{beta}{estimated individual coefficients.}
\item{betaest}{estimated group coefficients.}
\item{eta}{estimated common coefficients for z.}
\item{sig2}{estimated varaince.}
\item{group}{estimated group.}
\item{deltam}{coefficient differences.}
\item{flag}{code for convergence, 0 means converged.}
\item{rm}{primal residuals.}
\item{sm}{dual residuals.}
\item{tolpri}{primal tolerance.}
\item{toldual}{dual tolerance.}
\item{niteration}{number of iterations.}
}
\references{
Ma, S. and Huang, J. (2016). Estimating subgroup-specific treatment effects via concave fusion. arXiv preprint arXiv:1607.03717.

Ma, S. and Huang, J. (2017). A concave pairwise fusion approach to subgroup analysis. Journal of the American Statistical Association, 112(517):410–423.

Boyd, S., Parikh, N., Chu, E., Peleato, B., and Eckstein, J. (2011). Distributed opti- mization and statistical learning via the alternating direction method of multipliers. Foundations and Trends in Machine Learning, 3(1):1–122.

}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cal_initial}}, \code{\link{getgroup}}, \code{\link{BIClog}},  \code{\link{BICc}}
}
\examples{
data(simdat1)
y <- simdat1$y
z <- as.matrix(simdat1[,paste("z",1:5,sep="")])
x <- as.matrix(simdat1[,paste("x",1:2,sep="")])
betam0 <- cal_initial(y = y, z = z, x = x)
nr <- length(y)

### equal weight
wt <- rep(1, nr*(nr-1)/2)
res <- Spgr(y, z, x, wt, betam0, lam = 0.7,maxiter = 1000)
groupest <- getgroup(res$deltam, nr) ## the order of estimated group may not be the same as true group
BICc(res,y,z,x,c0 = 0.2)

## the same as res$group

### spatial weight
data(Cmat)
ordervalue <- getorder(Cmat)
wts <- exp(0.7*(1-ordervalue))
res2 <- Spgr(y, z, x, wts, betam0, lam = 1, maxiter = 3000)
groupest2 <- getgroup(res2$deltam, nr) ### the same as res2$group
BICc(res2,y,z,x,c0 = 0.2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
