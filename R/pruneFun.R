.findAlphas <- function(labels,Y,L,Fit){
	
	maxLev <- length(strsplit(labels[1],"")[[1]])
	
	groups <- lapply(1:maxLev,function(K){
		#browser()
		#cat(K," ")
		labs <- as.factor(substring(as.character(labels),1,K))
		fit <- Fit(Y,labs)
		YHat <- fit[labs]
		loss <- L(Y,c(YHat))
		loss<-tapply(loss,labs,mean)
		
		if(K<maxLev){
			child<-lapply(levels(labs),function(lab){
				lll<-unique(substring(as.character(labels),1,K+1))
				lll[substring(lll,1,K)==lab]
				})
			names(child)<-names(loss)<-as.character(levels(labs))
			}else{
				child<-NULL
				}
		
		list(loss=loss,children=child)
		
		})
		
	alphas <- unique(sort(unlist(lapply(maxLev:2,function(K){	
		SIC <- sort(unlist(lapply(names(groups[[K-1]]$loss),function(group){
			SIC1 <- 2*(1+groups[[K-1]]$loss[group])
			SIC2 <- 2*(length(groups[[K-1]]$children[group])+unlist(lapply(groups[[K-1]]$children[group],function(child)sum(groups[[K]]$loss[child]))))
			return(SIC1/SIC2)
			})))
		#cbind(AIC[-length(AIC)],AIC[-1])
		return(SIC)
		}))))
	
	return(alphas)
	
	}

.pruneTree <- function(labels,Y,L,Fit,alpha=1){
	
	maxLev <- length(strsplit(labels[1],"")[[1]])
	
	
	groups <- lapply(1:maxLev,function(K){
		#browser()
		#cat(K," ")
		labs <- as.factor(substring(as.character(labels),1,K))
		fit <- Fit(Y,labs)
		YHat <- fit[labs]
		loss <- L(Y,c(YHat))
		loss<-tapply(loss,labs,mean)
		
		if(K<maxLev){
			child<-lapply(levels(labs),function(lab){
				lll<-unique(substring(as.character(labels),1,K+1))
				lll[substring(lll,1,K)==lab]
				})
			names(child)<-names(loss)<-as.character(levels(labs))
			}else{
				child<-NULL
				}
		
		list(loss=loss,children=child)
		
		})
	
	doNotSplit <- lapply(maxLev:2,function(K){
		
		names(groups[[K-1]]$loss)[unlist(lapply(names(groups[[K-1]]$loss),function(group){
			
			SIC1 <- 2*(1+groups[[K-1]]$loss[group])
			SIC2 <- 2*(length(groups[[K-1]]$children[group])+unlist(lapply(groups[[K-1]]$children[group],function(child)sum(groups[[K]]$loss[child]))))
			return(SIC1<=alpha*SIC2)
			}))]
		})
		
	lab <- labels
	for(K in 1:(maxLev-1)){
		if(length(doNotSplit[[K]])==0) break
		for(grp in doNotSplit[[K]]){
			lab[substring(as.character(lab),1,maxLev-K)==grp]<-grp
			#lab[lab%in%groups[[maxLev-K]]$children[grp]] <- grp
			}
		}
		
	return(lab)
	
	
	}
#unlist(lapply(1:(maxLev-1),function(K)all(names(groups[[maxLev-K]]$children)%in%doNotSplit[[K]])))