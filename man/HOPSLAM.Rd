\name{HOPSLAM}
\alias{HOPSLAM}
\alias{fitHop}
\alias{crossHoptimum}

% - Also NEED an '\alias' for EACH other topic documented here.
\title{
   Supervised Clustering 
}
\description{
    Fit a Supervised Clustering Regression
}
\usage{
HOPSLAM <- function(X,Y, d ="euclid",dd=NULL,fold=5L,method=c("regression","classification"),L=L2,Fit=meanFit,pi=1,lambda=1,...)

fitHoptimum(X, Y, d = "euclid", dd = NULL, fold = 5L, method = c("regression", "classification"), L = L2, Fit = meanFit, pi = 1, lambda = 1, ...)

crossHop(fold,Y,X,L=L2,Fit=meanFit,dmat, d, K,...)
}
 - maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
       Feature matrix
}
  \item{Y}{
       Outcomes
}
  \item{d}{
       distance
}
  \item{dd}{
       Distance Matrix
}
  \item{fold}{
       Cross Validation Fold
}
  \item{method}{
       Regression or Classification
}
  \item{L}{
       Loss Function
}
  \item{Fit}{
       Fit Function
}
  \item{pi}{
       Contribution of the outcome in the silhouette score (\code{pi}=1 iff unsupervised)
}
  \item{lambda}{
       Contribution of the outcome in the distance matrix score (\code{lambda}=1 iff unsupervised)
}
  \item{\dots}{
       Arguments to be passed to \code{hopach}
}
}
\details{
    \code{HOPSLAM} call \code{fitHop} and creates a \code{HIP} object.
}
\value{
    an Object of class \code{HIP}
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
res <- HOPSLAM(XY$X,XY$Y,fold=10)
}
%  Add one or more standard keywords, see file 'KEYWORDS' in the
%  R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }  __ONLY ONE__ keyword per line
