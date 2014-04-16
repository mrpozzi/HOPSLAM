HOPSLAM <- function(X,Y, d ="euclid", dMat =NULL,fold=5L,method=c("regression","classification"),L=L2,Fit=meanFit,balanced=FALSE,verbose=FALSE,bNb=FALSE,forward=TRUE,maxVars=ncol(X)-1,...){
	
	
	method <- match.arg(method)
	call <- match.call()
	
	if(missing(X)){
		if(is.null(dMat)){
			error("You Must provide either a distance matrix or data!")
			}else{
				X <- matrix(NA,nrow(dMat),nrow(dMat))
				}
		}
	
	if((method=="classification")|is.factor(Y)){
		if(missing(L)) L <- L0
		if(missing(Fit)) Fit <- majVote
		Y <- as.character(Y)
		}

	
	argz <- lapply(2:length(call),function(i)eval(call[[i]]))
	names(argz) <- names(call)[2:length(names(call))]
		
	if(is.null(dMat)) dMat <- distancematrix(X, d=d)
	dMat <- as.matrix(dMat)
	
	if(!bNb){
		fit <- .build(X,Y,d,dMat=dMat,fold,method,L=L,Fit=Fit,balanced=balanced,verbose=verbose,...)
		} else {
			fit <- .buildNselect(X,Y,d,dMat=dMat,fold,method,L=L,Fit=Fit,balanced=balanced,verbose=verbose,bNb=bNb,forward=forward,maxVars=maxVars,...)
			}
			
	YHat <- Fit(Y,as.factor(fit$lab))[as.factor(fit$lab)]
	
	return(new("HIP",labels = as.factor(fit$lab),
	                 k = (1:ncol(X))[fit$k],
	                 alpha = fit$alpha,
	                 M = fit$medoids[order(levels(as.factor(fit$lab)))],
	                 YHat = c(YHat),#as.numeric(as.factor(c(YHat)))
	                 fold = fit$fold,
	                 hop = fit$hop,
	                 dMat = as.matrix(fit$dMat),
	                 Residuals = t(t(L(Y,c(YHat)[fit$lab]))),
	                 Risk = fit$err[1],
	                 call = call,
	                 args = argz))
	
	}
	
	
#################################
### Hopach Inference Partition

setClass("HIP",representation=representation(labels = "factor",
                                              k = "vector",
                                              alpha = "numeric",
                                              M = "numeric",
                                              YHat = "vector",
                                              fold = "numeric",
                                              hop = "list",
                                              dMat = "matrix",
                                              Residuals = "array",
                                              Risk = "numeric",
                                              call = "call",
                                              args = "list"))


##### lm Example...	
#Call:
#lm(formula = weight ~ group)
#
#Coefficients:
#(Intercept)     groupTrt  
#      5.032       -0.371  

setMethod("show",signature("HIP"),function(object){
	cat("Call:\n")
	print(object@call)
	cat("\nOptimal alpha: ")
	cat(object@alpha,"\n")
	cat("Number of Clusters: ")
	cat(length(levels(object@labels)),"\n")
	cat("Number of Variables: ")
	cat(length(object@k),"\n")
	})




#### WORK IN PROGRESS...
setMethod("plot", signature = "HIP", definition = function (x,whichOne=1L,main,...)  {
	
	show <- rep(FALSE, 3)
    show[whichOne] <- TRUE
	
	if(length(whichOne)>1)par(ask=TRUE)
	
	if(missing(main))main <- ""#bquote(paste(lambda,"=",.(x@lambda),sep=""))
	
	lab <- as.integer(x@labels)
	
	if(show[1L]){

		if(length(whichOne)>1) main <- "Clustered Y"
		
		if(is.numeric(x@args$Y)){
			barplot(x@args$Y[order(lab)],col=brewer.pal(n=max(3,length(unique(lab))),name="Blues")[sort(lab)],main=main,ylab=expression(Y[i(lab)]),xlab=expression(i(lab)))#BuGn
			}else{
				barplot(as.numeric(as.factor(x@args$Y[order(lab)])),col= brewer.pal(n=max(3,length(unique(lab))),name="Blues")[sort(lab)],main=main,ylab=expression(Y[i(lab)]),xlab=expression(i(lab)))
				}
				
		}	

	### Use pam's plotting trick....
	#main <- format(x@call)
    if(show[2L]){ 
    	sil <- silhouette(lab, dmatrix=as.matrix(x@dMat))
	 	
	 	if(length(whichOne)>1) main <- "Silhouette Plot"
		plot(sil,main=main)
		}
	
		
	if(show[3L])dplot(x@dMat, x@hop,ord="final", col=brewer.pal(n=11,name="Spectral"),showclusters=TRUE,labels=x@ labels,main=main)
							
	
	})
	
	
setMethod("[[", signature = "HIP", definition = function (x, i, j, ..., drop = TRUE) {
	if (!missing(j)) {stop("Wrong number of dimensions.")}
	if (!missing(i)) {
		return(switch(class(i),"character" = attributes(x)[[i]],
		                       "integer" = attributes(x)[[i]],
		                       "numeric" = attributes(x)[[i]],
		                       "logical" = attributes(x)[c(i,rep(FALSE,length(attributes(x))-length(i)))],
		                        stop("Subsetting object needs to be either a character, a numeric/integer or a logical.")
		                        ))
		}else{return(NULL)}
  
   })


setMethod("$", signature = "HIP", definition = function(x, name) {
	 x[[name]]
	 })
	 
setMethod("names", signature = "HIP", definition = function(x) {
	 slotNames(x)
	 })
	 
	 
setMethod("print",signature="HIP",function(x,...){
	cat("Call:\n")
	Kall <- lapply(1:length(x@call),function(j){
		if(j==1) return(x@call[[j]])
		if(length(x@args[[j-1]])<3){
			return(x@args[[j-1]])
			}else{
				return(x@call[[j]])
				}
			})
	names(Kall) <- names(x@call)
	Kall <- as.call(Kall)
	print(Kall)
	cat("\nOptimal alpha: ")
	cat(x@alpha,"\n")
	cat("Number of Clusters: ")
	cat(length(levels(x@labels)),"\n")
	cat("Number of Variables: ")
	cat(length(object@k),"\n")
	})


setMethod("summary",signature="HIP",function(object,...){
	
	print(object)
	cat("\nResiduals:\n")
	print(summary(c(object@Residuals)))
	
	})
	