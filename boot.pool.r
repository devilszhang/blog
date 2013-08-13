
#Dadong Zhang
#http://devilzhang.blogspot.com/
#University of Virginia
	
	
##############A function to pool results from multiple glm
#####fitall: an object from multiple glm
	boot.pool<-function(fitall){
		call<-match.call()
		m<-length(fitall)
		k<-nrow(fitall[[1]]$coefficients)
		names<-rownames(fitall[[1]]$coefficients)
		qhat<-matrix(NA, nrow=m, ncol=k, dimnames=list(1:m, names))
		u<-array(NA, dim=c(m, k, k), dimnames=list(1:m, names, names))	
		dfcom<-NA
		for(i in 1:m){
			qhat[i, ]<-fitall[[i]]$coef[,"Estimate"]
			u[i,,]<-as.matrix(vcov(fitall[[i]]))		
			dfcom[i]<-df.residual(fitall[[i]])
						}
		#add a function to calculate pooled degree of freedom 
		#method=Barnard, J. and Rubin, D.B. 1999, biometrika, 86, 948-955
		df<-function (m, lambda, dfcom, method){
    					lambda[lambda < 1e-04] <- 1e-04
   						dfold <- (m - 1)/lambda^2
         				dfobs <- (dfcom + 1)/(dfcom + 3) * dfcom * (1 - lambda)
               			df <- dfold * dfobs/(dfold + dfobs)
				        	 return(df)
				        if (method != "smallsample") 
        				df <- dfold
                		 }
		dfcom<-dfcom[1]
		qbar<-apply(qhat, 2, mean)
		ubar<-apply(u, c(2,3), mean)
		e<-qhat-matrix(qbar, nrow=m, ncol=6, byrow=T)
		b<-(t(e)%*%e)/(m-1)
		t<-ubar+(1+1/m)*b
		r<-(1+1/m)*diag(b/ubar)
		lambda<-(1+1/m)*diag(b/t)
		df<-df(m, lambda, dfcom, "smallsample")	
		fmi<-(r+2/(df+3))/(r+1)
		names(r)<-names(df)<-names(fmi)<-names(lambda)<-names
		fit<-list(call=call, qhat=qhat, u=u, qbar=qbar, ubar=ubar, b=b, t=t, r=r, dfcom=dfcom, df=df, fmi=fmi, lambda=lambda)
		oldClass(fit)<-c("mipo", oldClass(fitall))
		return(fit)
					}
					
		
