\name{Spgrr}
\alias{Spgrr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use Spgr to estimate subgroups with repeated measures with both z
 and x}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Spgrr(indexy,y, z, x, weights, betam0, ...)
}

\arguments{
  \item{indexy}{numeric numbers, index of repeated meausre sort assending}
  \item{y}{ response variable.}
  \item{z}{ explanatory variables which have the same coefficients among individuals.}
  \item{x}{ explanatory variables which have individual coefficients.}
  \item{weights}{weights associated with pairwise penalties.}
  \item{betam0}{initial values for betas.}
  \item{nu}{penalty parameter in, default value is 1.}
  \item{gam}{parameter in SCAD penalty, default value is 1.}
  \item{lam}{tuning parameter, default value is 0.5.}
  \item{maxiter}{maximum number of iterations.}
  \item{tolabs}{absolute tolerance, default value is 1e-4. }
  \item{tolrel}{relative tolerance, default value is 1e-2. }
}
\details{
The details of absolute tolerance and relative tolerance can be found in Boyd, S., Parikh, N., Chu, E., Peleato, B., and Eckstein, J. (2011). Distributed opti- mization and statistical learning via the alternating direction method of multipliers.
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
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  data(simdat2)
  y <- simdat2$y
  z <- as.matrix(simdat2[,paste("z",1:5,sep="")])
  x <- as.matrix(simdat2[,paste("x",1:2,sep="")])
  
  nr <- 49
  betam0 <- cal_initialr(indexy = simdat2$indexy,y = simdat2$y,z = z,x = x )
  wt <- rep(1, nr*(nr-1)/2)
  res1 <- Spgrr(indexy = simdat2$indexy, y, z, x, wt, betam0, lam = 0.1,maxiter = 1000)
  groupest1 <- getgroup(res1$deltam, nr) 
  BICcr(obj = res1,indexy = simdat2$indexy,y,z,x,c0 = 0.2)
  
  data(Cmat2)
  ordervalue <- getorder(Cmat2)
  wts <- exp(0.7*(1-ordervalue))
  res2 <- Spgrr(indexy = simdat2$indexy,y, z, x, wts, betam0, lam = 0.2, maxiter = 1000)
  groupest2 <- getgroup(res2$deltam, nr) ### the same as res2$group
  BICcr(obj = res2,indexy = simdat2$indexy,y,z,x,c0 = 0.2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
