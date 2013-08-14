se<-function(x) sd(x)/sqrt(r)
up<-function(x) mean(x)+qnorm(0.975)*se(x)
lo<-function(x) mean(x)-qnorm(0.975)*se(x)
#function for reference size# type1=alpha; type2=beta; pi1=specified rate of event in trt1; delta=minimum clinical effect size
tra.size<-function(type1, type2, pi1, pi2,delta)
    {
    t.size<-2*(qnorm(type1/2)+qnorm(type2))^2*(pi1*(1-pi1)+pi2*(1-pi2))/delta^2
    return(ceiling(t.size))
    }
#function for Gould size # opi: specified over all events across two trts.
Gould.size<-function(type1, type2, opi, delta)
    {
    ssize<-(4*(qnorm(type1/2)+qnorm(type2))^2)*(opi*(1-opi))/(delta^2)
    return (ceiling(ssize))
    }
#function for updating sample size after interim analysis;obspi: observed overall events rate. truepi: specifed rate; n: planned sample size
int.size<-function(obspi,truepi,n)
    {
    intn<-(obspi*(1-obspi))/(truepi*(1-truepi))*n
    return (ceiling(intn))
    }
#function for Shih size
Shih.size<-function(n,r, type1,type2,tpi1,tpi2,p,interim){
    hlam1<-(rbinom(r,as.integer(n*interim*p),tpi1)+rbinom(1,as.integer(n*interim*(1-p)),tpi2))/(n*interim)
    hlam2<-(rbinom(r,as.integer(n*interim*p),tpi2)+rbinom(1,as.integer(n*interim*(1-p)),tpi1))/(n*interim)
    epi1<-abs((p*hlam1-(1-p)*hlam2)/(2*p-1))
    epi2<-abs((p*hlam2-(1-p)*hlam1)/(2*p-1))
    eovp<-(epi1+epi2)/2
    size<-as.integer(int.size(eovp,(tpi1+tpi2)/2,n))
    return(ceiling(size))
}
#Function for Chisq test
chitest<-function(n1,n2,n11,n21)
    {
    stat<-(n1+n2)*(n11*(n2-n21)-(n1-n11)*n21)^2/(n1*n2*(n11+n21)*((n1+n2)-(n11+n21)))
    return(stat)
    }

#using T derivation
Tfunc<-function(Nr,Ng, opi=NULL, pi1=NULL, delta, pow)
	{
		Tg<-abs(delta)*sqrt(Ng/(4*opi*(1-opi)))-qnorm(1-pow)
		Tr<-abs(delta)*sqrt(Nr/(2*(pi1*(1-pi1)+(pi1+delta)*(1-pi1-delta))))-qnorm(1-pow)
		return(list(Tref=Tr, Tgould=Tg))
	}
Alphaerr<-function(estt, trut)
	{
	erro<-2*dnorm(estt)*(trut-estt)-estt*dnorm(estt)*(estt-trut)^2
	return(list(differror=erro))
	}

SAP<-function(r, type1, type2, interim, tpi1, delta,p)
{     
      tpi2<-tpi1+delta 
      #reference method: 
      			Ts<-tra.size(type1,type2,tpi1,tpi1+delta, delta)  # planned reference sample size
        		Tn1<-as.integer(Ts*interim)   #subjuects in each group at stage interim
          		Tnint1<-rbinom(r,Tn1,tpi1)    # number of events at interim from trt1
          		Tnint2<-rbinom(r,Tn1,tpi2) # number of events at interim from trt2 (under alternative hypothesis)
            	Tnint2n<-rbinom(r,Tn1,tpi1)    # number of events at interim from trt2 (under null hypothesis)
              	T.hpi<-(Tnint1+Tnint2)/as.integer(Ts*2*interim) # observed overall events rate
                chpi1<-abs(((2*T.hpi)-delta)/2) # observed event rate in trt1;
                chpi2<-((2*T.hpi)+delta)/2     # observed event rate in trt2;
      			Thsb<-tra.size(type1, type2,chpi1,chpi2,delta) #re-estimated sample size at interim
      #Gould methods:
      			Gs<-Gould.size(type1,type2,tpi1+0.5*delta,delta) # planned Gould sample size
       			Gn1<-as.integer(Gs*interim)
          		GSd<-Gn1-Tn1  
          		Gnint1<-Tnint1+rbinom(r,GSd,tpi1)    # Gould method: number of events at interim from trt1
           		Gnint2<-Tnint2+rbinom(r,GSd,tpi2) # Gould method:  number of events at interim from trt2 (under alternative hypothesis)
            	Gnint2n<-Tnint2n+rbinom(r,GSd,tpi1) # Gould method:  number of events at interim from trt2 (under null hypothesis)
            	G.hpi1<-(Gnint1)/as.integer(Gs*interim)
            	G.hpi2<-(Gnint2)/as.integer(Gs*interim)
              	G.hpi<-(Gnint1+Gnint2)/as.integer(Gs*2*interim) # Gould method:  observed overall events rate
      			Ghsb<-Gould.size(type1, type2,G.hpi,delta)   #re-estimated sample size at interim
      #Shih method based on the TS size
      			eventTA1<-rbinom(r,as.integer(Ts*interim*p),tpi1)
        		eventTA2<-rbinom(r,as.integer(Ts*interim*(1-p)),tpi2)
          		eventTAn2<-rbinom(r,as.integer(Ts*interim*(1-p)),tpi1)
          		thetaT1<-(eventTA1+eventTA2)/(Ts*interim) #stata1
      			eventTB1<-rbinom(r,as.integer(Ts*interim*p),tpi2)
        		eventTB2<-rbinom(r,as.integer(Ts*interim*(1-p)),tpi1)
          		eventTBn1<-rbinom(r,as.integer(Ts*interim*p),tpi1)
          		thetaT2<-(eventTB1+eventTB2)/(Ts*interim)  #stata2
      			Tepi1<-abs((p*thetaT1-(1-p)*thetaT2)/(2*p-1))  #pi1_hat
        		Tepi2<-abs((p*thetaT2-(1-p)*thetaT1)/(2*p-1)) #pi2_hat
      			Teovp<-(Tepi1+Tepi2)/2 #pi_hat
     			Sthsb<-as.integer(int.size(Teovp,(tpi1+tpi2)/2,Ts)) 
      #Shih method based on the GS size
      			eventGA1<-rbinom(r,as.integer(Gs*interim*p),tpi1)
        		eventGA2<-rbinom(r,as.integer(Gs*interim*(1-p)),tpi2)
          		eventGAn2<-rbinom(r,as.integer(Gs*interim*(1-p)),tpi1)
          		thetaG1<-(eventGA1+eventGA2)/(Gs*interim) #stata1
      			eventGB1<-rbinom(r,as.integer(Gs*interim*p),tpi2)
        		eventGBn1<-rbinom(r,as.integer(Gs*interim*p),tpi1)
          		eventGB2<-rbinom(r,as.integer(Gs*interim*(1-p)),tpi1)
            	thetaG2<-(eventGB1+eventGB2)/(Gs*interim)  #stata2
      			Gepi1<-abs((p*thetaG1-(1-p)*thetaG2)/(2*p-1))  #pi1_hat
        		Gepi2<-abs((p*thetaG2-(1-p)*thetaG1)/(2*p-1)) #pi2_hat
      			Geovp<-(Gepi1+Gepi2)/2 #pi_hat
      			Sghsb<-as.integer(int.size(Geovp,(tpi1+tpi2)/2,Gs)) 
      #Under and within proportion
       			#T.p.under<-sum(ifelse(Thsb< Ts, 1,0))/r
        		#G.p.under<-sum(ifelse(Ghsb< Gs, 1,0))/r
         		#Sht.p.under<-sum(ifelse(Sthsb< Ts, 1,0))/r
         		#Shg.p.under<-sum(ifelse(Sghsb< Gs, 1,0))/r
      			#T.p.in<-sum(ifelse((0.95*Ts<Thsb)&(Thsb<1.05*Ts),1,0))/r
          		#G.p.in<-sum(ifelse((0.95*Gs<Ghsb)&(Ghsb<1.05*Gs),1,0))/r
         		#Sht.p.in<-sum(ifelse((0.95*Ts<Sthsb)&(Sthsb< 1.05*Ts),1,0))/r
          		#Shg.p.in<-sum(ifelse((0.95*Gs<Sghsb)&(Sghsb< 1.05*Gs),1,0))/r
      			#put restriction on re-est SIZE >> planning size
      			#Ths<-ifelse(Thsb>=Ts, Thsb, Ts)   #Put restriction of N_hat>=N planning
      			#Sths<-ifelse(Sthsb>=Ts, Sthsb, Ts) 
      			#Ghs<-ifelse(Ghsb>=Gs, Ghsb, Gs)   #Put restriction of N_hat>=N planning
      			#Sghs<-ifelse(Sghsb>=Gs, Sghsb, Gs)
					Ths<-Thsb
      				Sths<-Sthsb 
      				Ghs<-Ghsb 
      				Sghs<-Sghsb
      #at the end of thr trial (finish all Ths or Ghs)
      #counting for alpha/power based on r replicates
      ##################################################################
      #n1n2:subjuects in each group at stage2
      				Tn2<-as.integer(Ths/2-Tn1)
        			Tn2<-ifelse(Tn2<=0,1,Tn2) #has to recruit at least on subject
      				Gn2<-as.integer(Ghs/2-Gn1)
        			Gn2<-ifelse(Gn2<=0,1,Gn2) #has to recruit at least on subject
      				STn2<-as.integer(Sths/2-Tn1)
        			STn2<-ifelse(STn2<=0,1,STn2) #has to recruit at least on subject
      				SGn2<-as.integer(Sghs/2-Gn1)
        			SGn2<-ifelse(SGn2<=0,1,SGn2) #has to recruit at least on subject
      ##################################################################
      				Tp1<-Tp2<-Tp2n<-Gads1<-Gads2<-Gads2n<-NULL
      				Tbin<-Gchi<-Talp<-Galp<-NULL
      				STp1<-STp2<-STp2n<-SGads1<-SGads2<-SGads2n<-NULL
      				STbin<-SGchi<-STalp<-SGalp<-NULL
      				Tn=Tn1+Tn2
      				Gn=Gn1+Gn2
        for (i in 1:r)
        	{ #reference method
         			Tp1[i]<-(rbinom(1,Tn2[i],tpi1)+Tnint1[i])/as.integer(Ths[i]*0.5) #final observed event rate in trt1
            		Tp2[i]<-(rbinom(1,Tn2[i],tpi2)+Tnint2[i])/as.integer(Ths[i]*0.5) #final observed event rate in trt2
              		Tp2n[i]<-(rbinom(1,Tn2[i],tpi1)+Tnint2n[i])/as.integer(Ths[i]*0.5) #final observed event rate in trt2(null hypothesis)
                	Tbin[i]<-abs((Tp2[i]-Tp1[i]))/sqrt((Tp1[i]*(1-Tp1[i]))/Tn[i]+(Tp2[i]*(1-Tp2[i]))/Tn[i]) #statistics for power
          			Talp[i]<-abs((Tp2n[i]-Tp1[i]))/sqrt((Tp1[i]*(1-Tp1[i]))/Tn[i]+(Tp2n[i]*(1-Tp2n[i]))/Tn[i])  #statistics for alpha
          	#Gould method
          			Gads1[i]<-rbinom(1,Gn2[i],tpi1)+Gnint1[i]#final number of events in trt1 (n11)
            		Gads2[i]<-rbinom(1,Gn2[i],tpi2)+Gnint2[i]#final number of events in trt2 (n21)
              		Gads2n[i]<-rbinom(1,Gn2[i],tpi1)+Gnint2n[i]#final number of events in trt2 (n21) null hypothesis
                	Gchi[i]<-chitest(Gn[i],Gn[i],Gads1[i],Gads2[i]) #statistics for power
                  	Galp[i]<-chitest(Gn[i],Gn[i],Gads1[i],Gads2n[i])#statistics for alpha
          	#shih method on Ts
          			STp1[i]<-(rbinom(1,STn2,tpi1)+eventTA1[i]+eventTB2[i])/as.integer(Sths[i]*0.5) #final observed event rate in trt1
            		STp2[i]<-(rbinom(1,STn2,tpi2)+eventTA2[i]+eventTB1[i])/as.integer(Sths[i]*0.5) #final observed event rate in trt2
              		STp2n[i]<-(rbinom(1,STn2,tpi1)+eventTAn2[i]+eventTBn1[i])/as.integer(Sths[i]*0.5) #final observed event rate in trt2(null hypothesis)
                	STbin[i]<-abs((STp2[i]-STp1[i]))/sqrt((STp1[i]*(1-STp1[i]))/as.integer(Sths[i]*0.5)+(STp2[i]*(1-STp2[i]))/as.integer(Sths[i]*0.5)) # power
          			STalp[i]<-abs((STp2n[i]-STp1[i]))/sqrt((STp1[i]*(1-STp1[i]))/as.integer(Sths[i]*0.5)+(STp2n[i]*(1-STp2n[i]))/as.integer(Sths[i]*0.5))  #alpha
          	#shih method on Gs
          			SGads1[i]<-(rbinom(1,SGn2[i],tpi1)+eventGA1[i]+eventGB2[i]) #final observed event rate in trt1
            		SGads2[i]<-(rbinom(1,SGn2[i],tpi2)+eventGA2[i]+eventGB1[i]) #final observed event rate in trt2
              		SGads2n[i]<-(rbinom(1,SGn2[i],tpi1)+eventGAn2[i]+eventGBn1[i]) #final observed event rate in trt2(null hypothesis)
                	SGchi[i]<-chitest(as.integer(Sghs[i]*0.5),as.integer(Sghs[i]*0.5),SGads1[i],SGads2[i]) #statistics for power
                  	SGalp[i]<-chitest(as.integer(Sghs[i]*0.5),as.integer(Sghs[i]*0.5),SGads1[i],SGads2n[i])#statistics for alpha    
          	}  
            #complete
      				Tbin<-Tbin[complete.cases(Tbin)]  #remove the NaNs
        			Gchi<-Gchi[complete.cases(Gchi)]  #remove the NaNs
          			Talp<-Talp[complete.cases(Talp)]  #remove the NaNs
            		Galp<-Galp[complete.cases(Galp)]  #remove the NaNs
              		STbin<-STbin[complete.cases(STbin)]   
                	STalp<-STalp[complete.cases(STalp)]
                 	SGchi<-SGchi[complete.cases(SGchi)] 
                    SGalp<-SGalp[complete.cases(SGalp)]
      		#Power
      				Tpower<-sum (ifelse(Tbin>qnorm(1-type1/2),1,0))/length(Tbin) #count number of rejects (power of reference method)
        			Gpower<-sum (ifelse(Gchi>qchisq(1-type1,1),1,0))/length(Gchi)  #count number of rejects  (power of Gould method)
      				STpower<-sum (ifelse(STbin>qnorm(1-type1/2),1,0))/length(STbin)
					SGpower<-sum(ifelse(SGchi>qchisq(1-type1,1),1,0))/length(SGchi)
      		#Alpha
      				Talpha<-sum (ifelse(Talp>qnorm(1-type1/2),1,0))/length(Talp)   #count number of rejects   (alpha of reference method)
        			Galpha<-sum (ifelse(Galp>qchisq(1-type1,1),1,0))/length(Galp)  #count number of rejects    (power of Gould method)
       				STalpha<-sum (ifelse(STalp>qnorm(1-type1/2),1,0))/length(STalp)
        	    	SGalpha<-sum(ifelse(SGalp>qchisq(1-type1,1),1,0))/length(SGalp) 
      		#CI		
      				Tr=mean(Ths)+c(-1,1)*qnorm(type1/2)*se(Ths) #cofidence interval of re-estimated reference sample size
      				Gr=mean(Ghs)+c(-1,1)*qnorm(type1/2)*se(Ghs) #cofidence interval of re-estimated Gould sample size
outcome<-rbind(
c(Ts, mean(Ths), se(Ths), mean(chpi1), se(chpi1), mean(chpi2), se(chpi2), mean(T.hpi), se(T.hpi), Talpha, Tpower),
c(Gs, mean(Ghs), se(Ghs), mean(G.hpi1), se(G.hpi1), mean(G.hpi2), se(G.hpi2), mean(G.hpi), se(G.hpi), Galpha, Gpower),
c(Ts, mean(Sths), se(Sths), mean(Tepi1), se(Tepi1), mean(Tepi2), se(Tepi2), mean(Teovp), se(Teovp), STalpha, STpower),
c(Gs, mean(Sghs), se(Sghs), mean(Gepi1), se(Gepi1), mean(Gepi2), se(Gepi2), mean(Geovp), se(Geovp), SGalpha, SGpower)
)
colnames(outcome)<-c("Planned(N)", "Smulated(N)", "Se(N)", "Simulated(p1)", "se(p1)", "Simulated(p2)", "se(p2)", "Simulated(p)", "se(p)", "Simulated(alpha)", "Smulated(Power)")
rownames(outcome)<-c("Chow", "Gould", "Shih-Chow", "Shih-Gould")
return(outcome)
}





