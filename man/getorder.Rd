\name{getorder}
\alias{getorder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Neighbor order vecotr
}
\description{
A function gives the neighbor order vector which can be used to define weights in pairwise penalties.
}
\usage{
getorder(Cmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Cmat}{Sparse adjacency matrix with diagnoal zero values}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

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
\code{\link{concavefusion}}, \code{\link{Cmat}}
}
\examples{
data(Cmat)
ordervalue <- getorder(Cmat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
