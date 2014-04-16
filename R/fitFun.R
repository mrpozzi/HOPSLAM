.build <- function(X,Y, d ="euclid", dMat =NULL,fold=5L,L=L2,Fit=meanFit,balanced=FALSE,verbose=FALSE,...){
	
	Hop <- hopach(X,dmat=dMat, d=d,...)
	nLev <- length(strsplit(format(Hop$final$labels[1],scientific=FALSE),split="")[[1]])

	if(verbose)cat("nLev=",nLev,"\n")
	

	if(balanced){
		if(is.numeric(Y)|length(unique(Y))==2){
			N <- lapply(unique(Y),function(y)which(Y==y))
			V <- unlist(lapply(N,function(n)sample(n,length(n),replace=FALSE)%%fold + 1))
			}else{
				n <- length(Y)
				V <- sample(1:n,n,replace=FALSE)
				V <- V%%fold + 1
				}
			}else{
				n <- length(Y)
				V <- sample(1:n,n,replace=FALSE)
				V <- V%%fold + 1
				}
	
	
	hopFold <- lapply(1:fold,function(v){
		if(verbose) cat("*")
		tryCatch(hopach(X[V!=v,],dmat= dMat[V!=v,V!=v], d=d,...),error=function(e)NA)
		})
	cat("\n")
	
	cvR <- lapply(1:fold,function(v,...){
		# options(warn=-1)
		# if(is.na(hopFold[[v]])){
			# hop <- hopach(X[V!=v,],dmat=as.matrix(dMat[V!=v,V!=v]), d=d,...)
			# } else 
		hop <- hopFold[[v]]
		# options(warn=0)
		
		labels <- format(hop$final$labels,scientific=FALSE)
		alphas <- .findAlphas(labels,Y[V!=v],L=L,Fit=Fit)
		
		err <- do.call(rbind,lapply(alphas,function(alpha){
			lab <- .pruneTree(labels,Y[V!=v],L=L,Fit=Fit,alpha)
			med <- hop$final$medoids[as.character(hop$final$medoids[,"label"])%in%lab,"medoid"]
			M <- X[V!=v,][med,]
			
			fit <- Fit(Y[V!=v],lab)
			Z <- t(predictHop(M=M,X=X[V==v,],d=d,...))
			YHat <- fit[Z]
			
			return(c(alpha,mean(L(Y[V==v],YHat),na.rm=TRUE),v))
			}))		
			stepfun(err[,1],c(err[1,2],err[,2]))
			})
	
	Risk <- Vectorize(function(alpha){mean(unlist(lapply(cvR,function(cv)cv(alpha))))})	
	RiskSq <- Vectorize(function(alpha){mean(unlist(lapply(cvR,function(cv)cv(alpha)^2)))})	
	#plot(Risk,xlim=c(0,5))
	
	best  <- optimize(Risk,c(0,10))
	alphaOpt <- best$minimum
	k <- 1:ncol(X)
		
	labels <- format(Hop$final$labels,scientific=FALSE)
	lab <- .pruneTree(labels,Y,L=L,Fit=Fit,alpha=alphaOpt)
	#table(lab)
	
	medoids <- Hop$final$medoids
	rownames(medoids) <- gsub(" ","",format(Hop$final$medoids[,"label"],scientific=FALSE))
	medoids <- medoids[unique(lab),2]
	
	err <- best$objective		
	err <- c(best$objective,sqrt(RiskSq(alphaOpt)-best$objective^2))
			
	list(alpha=alphaOpt,k=k,labels=lab,medoids=medoids,dMat=dMat,hop=Hop,err=err,fold=V)
			
	}



.buildNselect <- function(X,Y, d ="euclid", dMat =NULL,fold=5L,L=L2,Fit=meanFit,balanced=FALSE,verbose=FALSE,bNb=FALSE,forward=TRUE,maxVars=ncol(X),...){
	
	Hop <- hopach(X,dmat=dMat, d=d,...)
	# nLev <- length(strsplit(format(Hop$final$labels[1],scientific=FALSE),split="")[[1]])
	# if(verbose)cat("nLev=",nLev,"\n")
	
	if(balanced&(is.numeric(Y)|length(unique(Y))==2)){
		N <- lapply(unique(Y),function(y)which(Y==y))
		V <- unlist(lapply(N,function(n)sample(n,length(n),replace=FALSE)%%fold + 1))
		}else{
			n <- length(Y)
			V <- sample(1:n,n,replace=FALSE)
			V <- V%%fold + 1
			}
	
			
	cvR <- lapply(1:fold,function(v,...){
		
		cat("fold: ",v,"\n")
		varTree<-list()
		vars <- 1:ncol(X)
		if(forward) vars <- -vars
		for(k in 1:maxVars){
			
			print(vars)
			levFun <- lapply(1:length(vars),function(i){
				
				XT<-matrix(X[V!=v,vars[-i]],nrow=sum(V!=v))
				XV <- matrix(X[V==v,vars[-i]],nrow=sum(V==v))
				dMat <- as.matrix(distancematrix(XT, d=d))
				# hop <- tryCatch(hopach(XT,dmat=dMat, d=d,...),error=function(e)NA)
				# options(warn=-1)
				# if(is.na(hop)){
				hop <- hopach(XT,dmat=dMat, d=d,...)
					# }
				# options(warn=0)
				
				labels <- format(hop$final$labels,scientific=FALSE)
				alphas <- .findAlphas(labels,Y[V!=v],L=L,Fit=Fit)
				
				err <- do.call(rbind,lapply(alphas,function(alpha){
					lab <- .pruneTree(labels,Y[V!=v],L=L,Fit=Fit,alpha)
					med <- hop$final$medoids[as.character(hop$final$medoids[,"label"])%in%lab,"medoid"]
					M <- XT[med,]
					fit <- Fit(Y[V!=v],lab)
					names(fit)<-unique(lab)
					Z <- t(predictHop(M=M,X=XV,d=d,...))
					YHat <- fit[Z]
					return(c(alpha,mean(L(Y[V==v],c(YHat)),na.rm=TRUE),
					               mean(L(Y[V!=v],c(fit[lab])),na.rm=TRUE)))
					}))
				### return two functions and base the choice on the error
				list("Risk"=stepfun(err[,1],c(err[1,2],err[,2])),
				"Err"=stepfun(err[,1],c(err[1,3],err[,3])))
				
				})
			
			#Integrated Error
			IErr <- unlist(lapply(levFun,function(fun){
				sum(get("y",envir=environment(fun[[2]]))[-1]*diff(get("x",envir=environment(fun[[2]]))))
				}))
			varTree[[k]] <- levFun[which.min(IErr)][[1]]$Risk
			vars <- vars[-which.min(IErr)]
			
			}		
		
		return(varTree)
		
		})
	#browser()
	
	Risk <- lapply(1:length(cvR[[1]]),function(k){
		k<-k
		function(alpha){
			return(mean(unlist(lapply(cvR,function(cv)cv[[k]](alpha)))))
			}
		})
	
	RiskSq <- lapply(1:length(cvR[[1]]),function(k){
		k<-k
		function(alpha){
			return(mean(unlist(lapply(cvR,function(cv)cv[[k]](alpha)^2))))
			}
		})
	
	
	best  <- lapply(Risk,function(R) optimize(R,c(0,10)))
	
	# browser()
	# R <- unlist(lapply(best,function(b)b$objective))
	# R2 <- unlist(lapply(1:length(RiskSq),function(k)RiskSq[[k]](best[[k]]$minimum)))
	# sdR <- sqrt(R2-R^2)
	# plot(R,pch=16,ylim=c(min(R-1.96*sdR/sqrt(nrow(X))),max(R+1.96*sdR/sqrt(nrow(X)))))
	# segments(1:length(R),y0=R-1.96*sdR/sqrt(nrow(X)),y1=R+1.96*sdR/sqrt(nrow(X)))
	
	# maxVars <- which.min(unlist(lapply(best,function(b)b$objective)))
	maxVars <- which.min(unlist(lapply(best,function(b)b$objective)))
	best<-best[[maxVars]]
	#best  <- nlminb(start=rep(1,2),objective=Risk)
	alphaOpt <- best$minimum
	
	# browser()
	vars <- 1:ncol(X)
	if(forward) vars <- -vars
	vars <- 1:ncol(X)
	if(forward) vars <- -vars
	for(k in 1:maxVars){
		
		hopSel <- lapply(1:length(vars),function(i){
			
			XT<-matrix(X[,vars[-i]],nrow=length(Y))

			dMat <- distancematrix(XT, d=d)
			hop <- tryCatch(hopach(XT,dmat= dMat, d=d,...),error=function(e)NA)
			if(is.na(hop)){
				hop <- hopach(XT,dmat= dMat, d=d,...)
				}
				
			labels <- format(hop$final$labels,scientific=FALSE)
			
			lab <- .pruneTree(labels,Y,L=L,Fit=Fit,alphaOpt)
			fit <- Fit(Y,lab)
					
			return(list(mean(L(Y,c(fit[lab])),na.rm=TRUE),hop,lab))
			
			})
		
		Err <- unlist(lapply(hopSel,function(hop)hop[[1]]))
		vars <- vars[-which.min(Err)]
		Hop <- hopSel[[which.min(Err)]][[2]]
		lab <- hopSel[[which.min(Err)]][[3]]
		
		}	
	
	medoids <- Hop$final$medoids
	rownames(medoids) <- gsub(" ","",format(Hop$final$medoids[,"label"],scientific=FALSE))
	medoids <- medoids[unique(lab),2]
	
	err <- best$objective		
	err <- c(best$objective,sqrt(RiskSq[[maxVars]](alphaOpt)-best$objective^2))
	#err <- c(CVerr[best,1],sqrt(CVerr[best,2]-CVerr[best,1]^2))
			
	list(alpha=alphaOpt,k= vars,labels=lab,medoids=medoids,dMat=dMat,hop=Hop,err=err,fold=V)
			
	}



predictHop <- function(M,X,d,...){
	unlist(apply(X,1,function(x)which.min(distancevector(X=M,y=as.vector(x), d=d,na.rm=TRUE))))
	}


## X is the new data
## M is an HIP object
setMethod("predict", signature = "HIP", definition = function(object,X,...) {#,geo=FALSE
	 if(is.null(object@args$Fit)){
	 	if(is.null(object@args$method)){
	 		Fit <- meanFit
	 		}else{
	 			if(object@args$method=="regression") Fit <- meanFit
	 			if(object@args$method=="classification") Fit <- majVote
	 			}
	 	}else Fit <- object@args$Fit
	 	
	 if(is.null(object@args$d)){
	 	d <- "euclid"
	 	}else d <- object@args$d
	 	
	 if(missing(X))X<-object@args$X
	 
	 predictHop(M=object@args$X[object@M,object@k],X=X[,object@k],d=d,...)
	 
	 })

