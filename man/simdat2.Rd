\name{simdat1}
\alias{simdat1}
\docType{data}
\title{
A simulated data set for \code{\link{concavefusion}}.
}
\description{
This is a simulated data set to fit model \eqn{y_{ih} = z_{ih}\eta + x_{ih}\beta_i + e_{ih}}.
}
\usage{data("simdat2")}
\format{
  A data frame with 49 locations (7 squared grid) and 30 repeated measures for each location on the following 10 variables.
  \describe{
    \item{index}{a numeric vector to index the observation.}
    \item{y}{a numeric response varaible.}
    \item{z1,z2,z3,z4,z5}{z1 is the intercept column, z2,z3,z4,z5 are numeric vectors.}
    \item{x1,x2}{numeric vectors }
    \item{group}{true group information.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(simdat2)
summary(simdat2$y)
}
\keyword{datasets}
