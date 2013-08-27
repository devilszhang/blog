#Dadong Zhang
#dadong.zhang@virginia.edu
#CPHG, University of Virginia
#(434)982-0197

#################################Functions####################################################
####function for the imputation of the time1 based on the binomial distribution of the observed rate.
#x: a vector or data frame with NAs with column format ( id time1---timex)
#m: number of replicates
time1<-function(x, m){ 
							rep1<-function(x)
							pt1<-sum(x[, 2], na.rm=T)/length(which(!is.na(x[, 2])))
							t1ber<-rbinom(nrow(x), 1, pt1)
							t1imp<-ifelse(!is.na(x[, 2]), x[, 2], t1ber)
							return(t1imp)						
					}


####function for the each of the rest columns
#inacol: a list of the column names before timex without NAs, e.g. inacol<-c("t1imp", "t2imp").
#nacol: the column name with NA to be imputed, e.g. "day11"
#level: the number of times (columns) ahead of the time to be imputed. For example level c(0,1)=2, level c(0, 1, 2)=3.
timex<-function(x, inacol, nacol, m)
		{			
			if (length(inacol)<=1){x$ps=x[, inacol]}
			else {x$ps<-rowSums(x[, inacol])}
					level<-max(x$ps+1)
					rep1t2<-function(x, inacol, nacol) 
					{ 
							x<-x[order(x$ps), ]
							ms<-ns<-xs<-p1<-p2<-p<-NULL
							x$mi<-c()
								for (i in 1:level){
								ms[i]<-length(which(is.na(x[, nacol]) & x$ps==(i-1)))
								ns[i]<-length(which(x$ps==i-1))-ms[i]
								xs[i]<-length(which(x[, nacol]==1 & x$ps==(i-1)))	
													}
								for (i in 1:level){
								if (ms[i]>0 & xs[i]>0 & xs[i] !=ns[i]){
						    	p1[i]<-xs[i]-xs[i]/ns[i]
    		 				    p2[i]<-ns[i]-xs[i]-(ns[i]-xs[i])/ns[i]
    		 	    						  		 	    }
   								else
    								if (ms[i]>0 & xs[i]!=0){
    								r=0.1
      								p1[i]<-xs[i]+r
      								p2[i]<-ns[i]-xs[i]+r}   
    									 else { p1[i]<--9				
      									   p2[i]<--9}
			          								}
								for (i in 1:level){
								if (p1[i]!=-9 & p2[i]!=-9) {
								p[i]<-rbeta(1, p1[i], p2[i])}
								else {p[i]=-9}
													}
								for (i in 1:level) {
								if (p[i]!=-9){	bera<-rbinom(nrow(x), 1, p[i])}
								else {berb<-x[, nacol]
										bera<-rep(0, nrow(x))}
													}
								mss<-c()
								for (i in 1: level){
								mss<-c(mss, rep(length(which(is.na(x[, nacol]) & x$ps==i-1)), length(which(x$ps==i-1))))}
								x$mss<-mss
						ber<-ifelse(x$mss>0, bera, berb)
						x$oneimp<-ifelse(!is.na(x[, nacol]), x[, nacol], ber)
						return(x$oneimp)
					}
				repm<-function(m, x, inacol, nacol) replicate(m, rep1t2(x, inacol, nacol))
				rep.imp<-repm(m, x, inacol, nacol)
				x$prop<-rowSums(rep.imp)/m
				x$comimp<-ifelse(x$prop>0.5, 1, 0)
				x$imp<-ifelse(!is.na(x[, nacol]), x[, nacol], x$comimp)
		return(x$imp)
		}		

	
	
	
