\name{Spgrrx}
\alias{Spgrrx}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Use Spgr to estimate subgroups with repeated measures with individual covariates x}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Spgrrx(y, indexy, x, weights, betam0, ...)
}

\arguments{
  \item{indexy}{numeric numbers, index of repeated meausre sort assending}
  \item{y}{ response variable.}
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
