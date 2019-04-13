\name{BICog}
\alias{BICog}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BIC function 
}
\description{
A function gives the BIC value based on estimated results from \code{\link{Spgr}}. Cn = log(np+q)
}
\usage{
BICog(obj, y, z, x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{output of Spgr}
  \item{y}{ response vector with n observations .}
  \item{z}{ explanatory variables matrix which have the same coefficients among individuals.}
  \item{x}{ explanatory variables matrix with p columns which have individual coefficients.}
}
\details{
BIClog is a specical case of BIClogr without repeated measures. Cn in BIC is log(np+q).
}
\value{
a numeric BIC value
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
\code{\link{concavefusion}}, \code{\link{BIClogr}}
}
\examples{
## see the example in concavefusion

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
