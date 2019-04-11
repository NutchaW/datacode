dir <- "/Users/Xi/Desktop/seminar/rdata/"
dir <-"~/Google Drive/seminar/rdata/"
source(paste0(dir,"SIRS.R"))
# year=05
# city="NewYorkNY"
# prediction("NewYorkNY",05)
Peakweekprediction <- prediction("NewYorkNY",05)
# write.csv(Peakweekprediction, file=paste0(dir,'Peakweekprediction1.csv'))


prediction <- function(city,year) {
  num_ens=100 #number of ensembles
  num_times=40 #number of weeks in each season
  num_pred=1 #number of prediction realizations
  scale=1 #scaling parameter
  
  # load absolute humidity
  AbsHumidity <- read.csv(paste0(dir,"AbsHumidity.csv"))
  AH = AbsHumidity[,city]
  AH =c(AH,AH)
  
  # load ILI+, prepare observations
  ILIplus <- read.csv(paste0(dir,"ILIplus.csv"))
  ILIp = ILIplus[,city]
  ILIp = ILIp*scale
  obstime = ILIplus[,"Date"]
  obstimeday = as.numeric(as.Date(obstime, format = "%m/%d/%y"))
  s=paste0("10/01/",year) # start training from Oct 1st
  temp = which(obstimeday > as.numeric(as.Date(s, format = "%m/%d/%y")))
  obs=ILIp[temp[1]:(temp[1]+num_times-1)] #observations
  startdate=obstime[temp[1]]
  s=paste0("01/01/",year) # the first day of the year
  ts=as.numeric(as.Date(startdate, format = "%m/%d/%y"))-as.numeric(as.Date(s, format = "%m/%d/%y"))+1 # start date number in the year
  
  # load the initial conditions for ensemble members
  phi <- read.csv(paste0(dir,"params_statespace_RK_ext_seed_phi.csv"),header=FALSE)
  susceps = matrix(as.matrix(phi[,(5:31)]),nrow=100000*27,ncol=1)
  infects = matrix(as.matrix(phi[,(32:58)]),nrow=100000*27,ncol=1)
  params = phi[,(1:4)]
  # ???????????? warning('off','all') 
  # output prediction results of EAKFC
  Tpkpred = matrix(0,num_pred,num_times)
  
  # mapping operator
  HH = diag(7)
  H = HH[3,]
  
  N=500000 # total population
  seed=0.1
  discrete=0
  dt=1
  lambda=1.02 # inflation factor
  tic <- Sys.time()
  for (p in 1:num_pred) {
    rnd=cbind(ceiling(1e5*runif(100)),ceiling(27*1e5*runif(100)),ceiling(27*1e5*runif(100)))
    x=matrix(0,7,num_ens)
    x[1,]=susceps[rnd[,2],1] # S
    x[2,]=infects[rnd[,3],1] # I
    x[4,]=params[rnd[,1],3] # R0max
    x[5,]=params[rnd[,1],4] # R0min
    x[6,]=params[rnd[,1],1]*365 # L
    x[7,]=params[rnd[,1],2]*365 # D
    ### Run far enough so model spreads
    tmstep=7
    x=SIRS(x,ts,dt,tmstep,N,AH,seed,discrete)
    ### Begin looping through observations
    xprior=replicate(num_times,matrix(NaN,7,num_ens))
    xpost=xprior
    for (tt in 1:num_times) {
      # Get the variance of the ensemble
      ave=obs[max(1,tt-1)]
      ave=ave+obs[max(1,tt-2)]
      ave=ave+obs[max(1,tt-3)]
      ave=ave/3
      obs_var = 10^5+ave^2/5
      # inflation of x before assimilation
      x=as.matrix(apply(x,1,mean))%*%rep(1,num_ens)+lambda*(x-as.matrix(apply(x,1,mean))%*%rep(1,num_ens))
      prior_var = var(as.vector(H%*%x))
      post_var = prior_var*obs_var/(prior_var+obs_var)
      if (prior_var==0){
        post_var=0
      } 
    prior_mean = mean(H%*%x)
    post_mean = post_var*(prior_mean/prior_var + obs[tt]/obs_var)
    # Compute alpha and adjust distribution to conform to posterior moments
    alpha = (obs_var/(obs_var+prior_var))^0.5
    dy = post_mean + alpha*((H%*%x)-prior_mean)-H%*%x
    #  Loop over each state variable
    rr=matrix(0,1,dim(x)[1])
    for (j in 1:dim(x)[1]) {
      A=cov(x[j,],as.vector(H%*%x))
      rr[j]=A/prior_var
    }
    dx=t(rr)%*%dy
    #  Get the new ensemble and save prior and posterior
    xprior[,,tt]=x
    xnew = x + dx
    #  Corrections to DA produced aphysicalities
    xnew=checkDA(xnew)
    xpost[,,tt]=xnew
    #  Integrate forward one time step
    tcurrent = ts+tmstep*tt
    x=SIRS(xnew,tcurrent,dt,tmstep,N,AH,seed,discrete)
    # EAKF update is finished
    # Forecast and evaluation
    xpred=replicate(num_times+1,matrix(0,7,num_ens))
    # assign xprd before tt as posterior
    for (t in 1:tt) {
      xpred[,,t]=xpost[,,t]
      xpred[3,,t]=obs[t]*rep(1,num_ens) # assign as true observations
    }
    # forecast with error correction
    #################################
    # Start error correction
    # correct R0max
    # Integrate backward for one tmstep
    tcurrent = ts+tmstep*tt
    xtemp=SIRS(xpost[,,tt],tcurrent,-dt,-tmstep,N,AH,seed,discrete)
    # breeding error
    magnitude=0.4
    num_err=20
    tcurrent = ts+tmstep*(tt-1)
    xcnt = xpost[,,max(1,tt-1)]
    # compensate S error
    xcnt[1,]=xtemp[1,]
    for (k in 1:num_ens) {
      xbred = as.matrix(xcnt[,k])%*%rep(1,num_err+1)
      # impose perturbation on R0max, first column store the
      # unperturbed trajectory
      xbred[4,(2:(num_err+1))]=xbred[4,1]*(1+rnorm(num_err)*magnitude)
      # run for one step
      xbred=Re(xbred)
      xbred=SIRS(xbred,tcurrent,dt,tmstep,N,AH,seed,discrete)
      # inflation of xbred
      xbred=as.matrix(apply(xbred,1,mean))%*%rep(1,(num_err+1))+lambda*(xbred-as.matrix(apply(xbred,1,mean))%*%rep(1,(num_err+1)))
      # calculate error
      R0maxerr=xbred[4,2:(num_err+1)]-xbred[4,1]
      obserr=xbred[3,2:(num_err+1)]-xbred[3,1]
      # normalize
      xscale=max(abs(obserr))
      yscale=max(abs(R0maxerr))
      if (xscale!=0) {
        obserr=obserr/xscale
        R0maxerr=R0maxerr/yscale
        # checker of phase transition
        # avoid phase transision
        # R0maxerr=R0maxerr[obserr>(min(obserr)+0.02)]
        R0maxerr=R0maxerr[Mod(obserr)>min(Mod(obserr[which(Mod(obserr)==min(Mod(obserr)))]+0.02))]
        # obserr=obserr[obserr>(min(obserr)+0.02)]
        obserr=obserr[Mod(obserr)>min(Mod(obserr[which(Mod(obserr)==min(Mod(obserr)))]+0.02))]
        # fit
        R0maxerr1=Re(R0maxerr); obserr1=Re(obserr)
        P=lm(R0maxerr1~poly(obserr1,3,raw=T))
        # P1=lm(Re(R0maxerr)~poly(Re(obserr),3,raw=T))
        # P2=lm(Im(R0maxerr)~poly(Im(obserr),3,raw=T))
        # infer error in R0max
        pcoef=as.numeric(P$coefficients)
        pcoef[is.na(pcoef)]=0
        deltaobs=xbred[3,1]-xprior[3,k,tt]-dy[k]
        deltaR0max=(pcoef[4]*(deltaobs/xscale)^3+pcoef[3]*(deltaobs/xscale)^2+pcoef[2]*(deltaobs/xscale)+pcoef[1])*yscale
      } else {
        R0maxerr=R0maxerr/yscale
        # checker of phase transition
        # avoid phase transision
        R0maxerr=c()
        obserr=c()
        # fit
        deltaobs=xbred[3,1]-xprior[3,k,tt]-dy[k]
        deltaR0max=0
      }

      if (Mod(deltaR0max)>(0.25*xprior[4,k,tt])) {
        # deltaR0max=sign(deltaR0max)*0.25*xprior[4,k,tt]
        deltaR0max=deltaR0max/Mod(deltaR0max)*0.25*xprior[4,k,tt]
      }
      xpred[4,k,tt]=xprior[4,k,tt]-deltaR0max 
    }
    # correct S
    magnitude=0.4
    num_err=20
    tcurrent = ts+tmstep*(tt-1)
    xcnt = xpost[,,max(1,tt-1)]
    xcnt[4,]=xpred[4,,tt] # compensate correction of R0max
    for (k in 1:num_ens) {
      xbred=as.matrix(xcnt[,k])%*%rep(1,num_err+1)
      # impose perturbation on S, first column store the
      # unperturbed trajectory
      xbred[1,2:(num_err+1)]=xbred[1,1]*(1+rnorm(num_err)*magnitude)
      # run for one step
      xbred=Re(xbred)
      xbred=SIRS(xbred,tcurrent,dt,tmstep,N,AH,seed,discrete)
      # inflation of xbred
      xbred=as.matrix(apply(xbred,1,mean))%*%rep(1,(num_err+1))+lambda*(xbred-as.matrix(apply(xbred,1,mean))%*%rep(1,(num_err+1)))
      # calculate error
      Serr=xbred[1,2:(num_err+1)]-xbred[1,1]
      obserr=xbred[3,2:(num_err+1)]-xbred[3,1]
      # normalize
      xscale=max(abs(obserr))
      yscale=max(abs(Serr))
      if (xscale!=0) {
        obserr=obserr/xscale
        Serr=Serr/yscale
        # checker of phase transition
        # avoid phase transision
        obserr=Re(obserr);Serr=Re(Serr)
        Serr=Serr[obserr>min(obserr)+0.02]
        # Serr=Serr[Mod(obserr)>min(Mod(obserr[which(Mod(obserr)==min(Mod(obserr)))]+0.02))]
        obserr=obserr[obserr>min(obserr)+0.02]
        # obserr=obserr[Mod(obserr)>min(Mod(obserr[which(Mod(obserr)==min(Mod(obserr)))]+0.02))]
        # fit
        if (length(obserr)>=3) {
          d_val=3} else {d_val=length(obserr)}
        P=lm(Serr~poly(obserr,d_val,raw=T))
        # infer error in S
        pcoef=as.numeric(P$coefficients)
        pcoef[is.na(pcoef)]=0
        if (d_val<3) {pcoef=c(pcoef,rep(0,(3-d_val)))}
        deltaobs=xbred[3,1]-xprior[3,k,tt]-dy[k]
        deltaS=(pcoef[4]*(deltaobs/xscale)^3+pcoef[3]*(deltaobs/xscale)^2+pcoef[2]*(deltaobs/xscale)+pcoef[1])*yscale
      } else {
        obserr=0
        Serr=Serr/yscale
        # checker of phase transition
        # avoid phase transision
        Serr=c()
        obserr=c()
        # fit
        deltaobs=xbred[3,1]-xprior[3,k,tt]-dy[k]
        deltaS=0
      }
      # if (abs(deltaS)>0.25*xprior[1,k,tt]) {
      #   deltaS=sign(deltaS)*0.25*xprior[1,k,tt]
      # }
      if (Mod(deltaS)>0.25*xprior[1,k,tt]) {
        deltaS=deltaS/Mod(deltaS)*0.25*xprior[1,k,tt]
      }
      xpred[1,k,tt]=xprior[1,k,tt]-deltaS   
    }
    xpred[,,tt]=checkDA(Re(xpred[,,tt]))
    # Forecast with error correction
    tpred=ts+tmstep*tt
    for (t in tt:num_times) {
      xpred[,,t+1]=SIRS(Re(xpred[,,t]),tpred,dt,tmstep,N,AH,seed,discrete)
      tpred=tpred+tmstep
    }
    obspred=aperm(xpred,c(2,3,1))
    obspred=Re(obspred[,,3]) # num_ens*num_times
    # evaluation
    # peak week
    pwens=rep(NaN,num_ens)
    for (i in 1:num_ens) {
      temp=which(obspred[i,]==max(obspred[i,]))
      pwens[i]=temp[1]
    }
    pwp=floor(mean(pwens))
    Tpkpred[p,tt]=pwp
  }
    
  }
  Tpkpred <- as.numeric(Tpkpred)
  Peakweekprediction <- cbind(obs,Tpkpred)
  return (Peakweekprediction)
}
toc <- Sys.time()

checkDA = function(x) {
  N=500000
  ug=max(x[1,]) # Corrects if S>N
  if (ug>N) {
    x[1,x[1,]>N]=N-1
  }
  ug=max(x[2,]) # Corrects if I>N
  if (ug>N) {
    x[2,x[2,]>N]=median(x[2,])
  }
  ug=min(min(x)) # Corrects if any state or parameter nudges negative
  if (ug<=0) {
    for (i in 1:7) {
      x[i,x[i,]<0]=mean(x[i,])
    }
  }
  ug=min(x[6,]) # Corrects if L < 200 days
  if (ug<200) {
    x[6,x[6,]<200]=median(x[6,])
  }
  ug=min(x[7,]) # Corrects if D < .5
  if (ug<=.5) {
    x[7,x[7,]<0.5]=median(x[7,])
  }
  ug=min(x[4,]-x[5,]) # Corrects if R0mx <= R0mn
  if (ug<=0) {
    x[4,x[4,]<x[5,]]=x[5,x[4,]<x[5,]]+0.01
  }
  return (x)
}
