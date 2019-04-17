library(ForecastFramework)
library(R6)
library(forecast)

source(paste0(dir,"SIRS.R"))

# load data
dir <- paste0(getwd(), "/")
AbsHumidity <- read.csv(paste0(dir,"AbsHumidity.csv"))
ILIplus <- read.csv(paste0(dir,"ILIplus.csv"))
phi <- read.csv(paste0(dir,"params_statespace_RK_ext_seed_phi.csv"),header=FALSE)

### fit
fitData = SIRS_EAKF_model_fit(ILIplus, AbsHumidity, phi, year=05, city="NewYorkNY")

### 3-week ahead
forecasts <- SIRS_EAKF_model_forecast(AbsHumidity, city="NewYorkNY", ts=275, fitData=fitData, steps = 3)

### functions

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


SIRS_EAKF_model_fit = function(data, AbsHumidity,phi,year, city) {
  
  num_ens=100 #number of ensembles
  num_times=40 #number of weeks in each season
  num_pred=1 #number of prediction realizations
  scale=1 #scaling parameter
  
  AH = AbsHumidity[,city]
  AH =c(AH,AH)
  ILIp = data[,city]
  ILIp = ILIp*scale
  obstime = data[,"Date"]
  obstimeday = as.numeric(as.Date(obstime, format = "%m/%d/%y"))
  s=paste0("10/01/",year) # start training from Oct 1st
  temp = which(obstimeday > as.numeric(as.Date(s, format = "%m/%d/%y")))
  obs=ILIp[temp[1]:(temp[1]+num_times-1)] #observations
  startdate=obstime[temp[1]]
  s=paste0("01/01/",year) # the first day of the year
  ts=as.numeric(as.Date(startdate, format = "%m/%d/%y"))-as.numeric(as.Date(s, format = "%m/%d/%y"))+1 # start date number in the year
  
  # load the initial conditions for ensemble members
  susceps = matrix(as.matrix(phi[,(5:31)]),nrow=100000*27,ncol=1)
  infects = matrix(as.matrix(phi[,(32:58)]),nrow=100000*27,ncol=1)
  params = phi[,(1:4)]
  
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
    }
    fitData=x
  }
  return(fitData)
}

### forecast
SIRS_EAKF_model_forecast = function(AbsHumidity, city, ts, fitData, steps) {
  
  AH = AbsHumidity[,city]
  AH =c(AH,AH)
  
  num_ens=100
  num_times=steps
  num_pred=1
  scale=1
  
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
  
  for (p in 1:num_pred) {
    x=fitData
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
      # xpred=replicate(num_times+1,matrix(0,7,num_ens))
      xpred=replicate(num_times,matrix(0,7,num_ens))
      # assign xprd before tt as posterior
      for (t in 1:tt) {
        xpred[,,t]=xpost[,,t]
        xpred[3,,t]=obs[t]*rep(1,num_ens) # assign as true observations
      }
      
      obspred=aperm(xpred,c(2,3,1))
      obspred=Re(obspred[,,2]) # num_ens*num_times
      # evaluation
      # incidence
      forecasts <- round(colMeans(obspred))
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
  return(forecasts)
}


