	#Dadong Zhang
	#http://devilzhang.blogspot.com/
	#University of Virginia
	
	#######BOOTSTRAP IMPUTATION####
	##check the dataset before to apply the boot_imp function:
	##(1) code the variable to impute as: 0=miss, 1=not miss;
	##(2) group the miss variable according to the probs=c(0, 0.2, 0.4, 0.6, 0.8, 1) according to propensity scores;
	##(3) sort the data by group and the variable to impute;
	
	###function:
	#data: vector or data.frame
	#group, the groups according to propensity score
	#imp, the variable to impute
	#m, number of replicates
	boot.imp=function(data, group, imp, m){
		once=function(data, group, imp){
		g1=data[data[, group]==1, imp]
			grp1=g1[!is.na(g1)]
			n0g1=length(g1[is.na(g1)])
			impg1=sample(sample(grp1, length(grp1), replace=T), n0g1, replace=T)
		g2=data[data[, group]==2, imp]
			grp2=g2[!is.na(g2)]
			n0g2=length(g2[is.na(g2)])
			impg2=sample(sample(grp2, length(grp2), replace=T), n0g2, replace=T)
		g3=data[data[, group]==3, imp]
			grp3=g3[!is.na(g3)]
			n0g3=length(g3[is.na(g3)])
			impg3=sample(sample(grp3, length(grp3), replace=T), n0g3, replace=T)
		g4=data[data[, group]==4, imp]
			grp4=g4[!is.na(g4)]
			n0g4=length(g4[is.na(g4)])
			impg4=sample(sample(grp4, length(grp4), replace=T), n0g4, replace=T)
		g5=data[data[, group]==5, imp]
			grp5=g5[!is.na(g5)]
			n0g5=length(g5[is.na(g5)])
			impg5=sample(sample(grp5, length(grp5), replace=T), n0g5, replace=T)
		imp=c(impg1, grp1, impg2, grp2, impg3, grp3, impg4, grp4, impg5, grp5)
		   }
	bootm=function(data, group, imp, m) replicate(m, once(data, group, imp))	  
	impx=bootm(data, group, imp, m)
	return(impx)
		}
