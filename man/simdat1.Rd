\name{simdat1}
\alias{simdat1}
\docType{data}
\title{
A simulated data set for \code{\link{concavefusion}}.
}
\description{
This is a simulated data set to fit model \eqn{y_i = z_i\eta + x_i\beta_i}.
}
\usage{data("simdat1")}
\format{
  A data frame with 400 observations on the following 10 variables.
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
data(simdat1)
summary(simdat1$y)
}
\keyword{datasets}
