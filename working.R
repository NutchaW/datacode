library(ForecastFramework)
library(R6)
library(forecast)
source_github('https://raw.githubusercontent.com/reichlab/forecast-framework-demos/master/models/ContestModel.R')

# load data
AbsHumidity <- read.csv(paste0(dir,"AbsHumidity.csv"))
data <- read.csv(paste0(dir,"ILIplus.csv"))
phi <- read.csv(paste0(dir,"params_statespace_RK_ext_seed_phi.csv"),header=FALSE)

## source R6 SIRS_EAKF model
## create IncidenceMatrix object out of data
## use $fit method to fit all remaining location(s)
## instantiate a new SIRS_EAKF model using $new
## subset the data to only a few locations

### fit
SIRS_EAKF_model$fit(ILIplus, AbsHumidity, phi, year=05, city="NewYorkNY")

 ### forecast
steps <- 15 # forecast ahead `step` number of weeks
forecast_X <- SIRS_EAKF_model$forecast(fitData, steps = steps)

SIRS_EAKF_model <- R6Class(
  inherit = ContestModel,
  private = list(
    .data = NULL,        ## every model should have this
    .hudimity=NULL,
    .phi=NULL,
    .models = list(),    ## specific to models that are fit separately for each location
    .nsim = integer(0),     ## models that are simulating forecasts need this
    .period = integer(0), ## specific to SARIMA models
    .num_ens=integer(0), #number of ensembles
    .num_times=integer(0), #number of weeks in each season
    .num_pred=integer(0), #number of prediction realizations
    .scale=integer(0),
    .N=integer(0), # total population
    .seed=integer(0),
    .discrete=integer(0),
    .dt=integer(0),
    .lambda=integer(0) # inflation factor
  ),
  public = list(
    ## data will be MatrixData
    fit = function(data, AbsHumidity,phi,year, city) {
      if("fit" %in% private$.debug){browser()}
      ## stores data for easy access and checks to make sure it's the right class
      # private$.data <- IncidenceMatrix$new(data)
      # AbsHumidity <- read.csv(paste0(dir,"AbsHumidity.csv"))
      AH = AbsHumidity[,city]
      AH =c(AH,AH)
      
      # load ILI+, prepare observations
      # ILIplus <- read.csv(paste0(dir,"ILIplus.csv"))
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
      # phi <- read.csv(paste0(dir,"params_statespace_RK_ext_seed_phi.csv"),header=FALSE)
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
          # Forecast and evaluation
          xpred=replicate(num_times+1,matrix(0,7,num_ens))
          # assign xprd before tt as posterior
          for (t in 1:tt) {
            xpred[,,t]=xpost[,,t]
            xpred[3,,t]=obs[t]*rep(1,num_ens) # assign as true observations
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
        fitData=xnew
      }
    },
    forecast = function(fitData, steps) {
      ## include for debugging
      if("forecast" %in% private$.debug){browser()} 
      
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
        x=fitData
        ### Run far enough so model spreads
        tmstep=7
        x=SIRS(x,ts,dt,tmstep,N,AH,seed,discrete)
        ### Begin looping through observations
        xprior=replicate(num_times,matrix(NaN,7,num_ens))
        xpost=xprior
        for (tt in 1:steps) {
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
    },
    initialize = function(num_ens=100,num_times=40,N=500000,lambda=1.02,num_pred=1,dt=1,seed=0.1,scale=1) { 
      ## this code is run during SARIMAModel$new()
      ## need to store these arguments within the model object
      private$.nsim <- nsim
      private$.period <- period
      private$.num_ens=num_ens #number of ensembles
      private$.num_times=num_times #number of weeks in each season
      private$.num_pred=num_pred #number of prediction realizations
      private$.scale=scale
      private$.N=N # total population
      private$.seed=seed
      private$.discrete=discrete=0
      private$.dt=dt
      private$.lambda=lambda
    },
    predict = function(newdata) {
      stop("predict method has not been written.")
    }
  ),
  active = list(
    ## This list determines how you can access pieces of your model object
    data = function(value) {
      ## use this form when you want this parameter to be un-modifiable
      if(!missing(value))
        stop("Writing directly to the data is not allowed.")
      return(private$.data)
    },
    models = function(value) {
      ## use this form when you want this parameter to be un-modifiable
      if(!missing(value))
        stop("Writing directly to the models is not allowed.")
      return(private$.models)
    },
    nsim = function(value) {
      ## use this form when you want to be able to change this parameter
      private$defaultActive(type="private", ".nsim", val=value)
    },
    period = function(value) {
      ## use this form when you want this parameter to be un-modifiable
      if(!missing(value))
        stop("Writing directly to the model period is not allowed.")
      return(private$.period)
    }
  )
)
