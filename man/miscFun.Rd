\name{Fitting Functions}
\alias{majVote}
\alias{meanFit}
% - Also NEED an '\alias' for EACH other topic documented here.
\title{
    Prediction Routines
}
\description{
    Predicts from the Clustering labels using majority vote or cluster mean
}
\usage{
meanFit(Y, lab)
majVote(Y, lab)
}
% - maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
       Outcome variable
}
  \item{lab}{
       cluster labels.
}
}
\details{
    ~~ If necessary, more details than the description above ~~
}
\value{
    The average outcome for each cluster
}
\references{
   ~put references to the literature/web site here ~
}
\author{
    Luca Pozzi, p.luc@stat.berkeley.edu
}
\note{
    ~~further notes~~
}

%   ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
   ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
XY <- simHop(n=50,p=2,mu=0, muk=c(0,5,10),sX=0.5,sY=1,scale=10)
meanFit(XY$Y, XY$lab)
}
%  Add one or more standard keywords, see file 'KEYWORDS' in the
%  R documentation directory.
\keyword{ Majority Vote }
\keyword{ Clustering }  
\keyword{ Prediction }  
