\name{simHop}
\alias{simHop}

\title{
Simulation Routine    
}
\description{
    Function for simulating features $X$ in $R^p$ and outcome $Y$ according to a partition of the feature space.
}
\usage{
simHop(n,p,q=p,mu=0, muk,sX=1,sY=1,scale=1,seed=sample(.Random.seed,1),disj=FALSE,easy=TRUE)
}

\arguments{
  \item{n}{
       Number of observations
}
  \item{p}{
       Size of the feature space
}
  \item{q}{
       "real" size of the feature space, i.e. size of the informative subset
}
  \item{mu}{
       non-random effect for the $Y$
}
  \item{muk}{
       Vector of group mean for the $Y$ in each cluster
}
  \item{sX}{
       feature variance
}
  \item{sY}{
       outcome variance
}
  \item{scale}{
       radius of the hypersphere the medoids are on
}
  \item{seed}{
       random seed
}
  \item{easy}{
       how fuzzy the real cluster are
}
}
\value{
    A list with components
    \item{X }{The features}
    \item{Y }{The outcome}
\item{M }{The Medoids}
\item{labels}{The cluster labels}
\item{seed}{The random seed}
}
\author{
    Luca Pozzi, p.luc@stat.berkeley.edu
}
\examples{
XY <- simHop(n=50,p=2,mu=0, muk=c(0,5,10),sX=0.5,sY=1,scale=10)
plot(XY$X,col=XY$lab,xlab=expression(X[1]),ylab=expression(X[2]),main="Feature Space")


}

\keyword{ simulation }
\keyword{ linear mixed models }  
