\name{cal_initial}
\alias{cal_initial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to calculate initial values in \code{\link{Spgrr}}. 
}
\description{
Calculate initial values for individual coefficients with repeated measures.
}
\usage{
cal_initialr(indexy, y, z, x,lam0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{indexy}{numeric numbers, index of repeated meausre sort assending}
  \item{y}{response vector}
  \item{z}{ explanatory variables matrix which have the same coefficients among individuals.}
  \item{x}{ explanatory variables matrix with p columns which have individual coefficients.}
  \item{lam0}{a constant, default value is 0.0001}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a numeric matrix with the same dimension as x.
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
\code{\link{Spgrr}}
}
\examples{
## calculate initial values 
data(simdat2)
  y <- simdat2$y
  z <- as.matrix(simdat2[,paste("z",1:5,sep="")])
  x <- as.matrix(simdat2[,paste("x",1:2,sep="")])
  betam0 <- cal_initialr(indexy = simdat2$indexy,y = simdat2$y,z = z,x = x )
  plot(betam0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
