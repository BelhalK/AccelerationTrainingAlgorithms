data <- function(res)
{
  warfarin.saemix <- res$y1
	warfarin.saemix["amount"] <- 0
	treat <- res$treatment
	treat["y1"] <- 0
	treat <- treat[c(1,2,4,3)]

	j <- 1
	l<-c()
	for (i in 1:nrow(warfarin.saemix)) {
	    
	    if(t(warfarin.saemix["id"])[i]==t(treat["id"])[j]){
	        print(rownames(warfarin.saemix[i,]))
	        l <- rbind(l,rownames(warfarin.saemix[i,]))
	        j<-j+1
	      } 
	}

	warfarin.saemix <- rbind(treat[1,], warfarin.saemix)
	j <- 2
	for (i in l[-1]){
	  # print(typeof(as.numeric(i)))
	  warfarin.saemix <- rbind(warfarin.saemix[1:(as.numeric(i)-1),], treat[j,], warfarin.saemix[(as.numeric(i)+1):nrow(warfarin.saemix),])
	  j <- j +1
	}

	rownames(warfarin.saemix) <- 1:nrow(warfarin.saemix)
	return(warfarin.saemix)
}


warfarin.saemix <- data(res)