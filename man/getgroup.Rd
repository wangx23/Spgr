\name{getgroup}
\alias{getgroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Obtain estimated group
}
\description{
A function to obtain the estimated group based on the results from \code{\link{concavefusion}}. 
}
\usage{
getgroup(deltam, n, tol = 1e-4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{deltam}{A matrix with p rows (dimension of x) and  n(n-1)/2 columns, the matrix of the difference among \eqn{\beta}s.}
  \item{n}{number of observations.}
  \item{tol}{numeric value to define 0s for the difference between \eqn{\beta_i - \beta_j}, default value is 1e-4.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
estimated group
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
## see the example in concavefusion}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
