L0 <- function(Y,YHat){
	1*(Y!=YHat)
	}

L2 <- function(Y,YHat){
	(Y-YHat)^2
	}

majVote <- function(Y,lab){

	tab <- table(lab, Y)
	
	return(colnames(tab)[apply(as.matrix(tab),1,which.max)])
		
	}

meanFit <- function(Y,lab){
	tapply(Y,lab,mean)
	}

medianFit <- function(Y,lab){
	tapply(Y,lab,median)
	}

maxFit <- function(Y,lab){
	tapply(Y,lab,max)
	}
