simHop <-
function(n,p,q=0,mu=0, muk,sX=1,sY=1,scale=1,seed=NULL){#,easy=TRUE
	
	if(is.null(seed)){
		seed <- runif(1)
		seed <- sample(.Random.seed,1)
		}else{
			set.seed(seed)
			}
	
	K <- length(muk)
	
	mm <- sample(1:(p-q),p-q,replace=FALSE)%%K + 1
	MM <- do.call(rbind,lapply(1:K,function(k) as.numeric(mm==k)))
	MM <- MM*scale
	MM <- cbind(MM,diag(0,nrow=K,ncol=q))
	MM <- MM[,sample(1:p,p,replace=FALSE)]
	
	lab<- sample(1:K,n,replace=TRUE)

	XX <- do.call(rbind,lapply(lab,function(l)rnorm(p,mean=MM[l,],sd=sX)))
	
	ZZ <- t(apply(XX,1,function(x){
		z <- apply(MM,1,function(m){
			distancematrix(rbind(x,m), d="euclid")[2,1]
			})
			which.min(z)
		}))
		
	#lab <- ZZ
			
	YY <- mu + muk[lab]
	YY <- rnorm(n,mean=YY,sd=sY)
	
	lab <- as.factor(lab)
	#levels(ZZ) <- as.character(levels(lab))
	
	#return(list(X=XX,Y=YY,M=MM,labels=lab,seed=seed))
	return(list(X=XX,Y=as.matrix(YY),Z=ZZ,M=MM,labels=lab,seed=seed))
	
	}
