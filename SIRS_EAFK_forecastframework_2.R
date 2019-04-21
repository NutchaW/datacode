library(ForecastFramework)
library(R6)
library(forecast)

source_github <- function(u) {
  library(RCurl)
  # read script lines from website and evaluate
  script <- getURL(u, ssl.verifypeer = FALSE)
  eval(parse(text = script),envir=.GlobalEnv)
}

source_github('https://raw.githubusercontent.com/reichlab/forecast-framework-demos/master/models/ContestModel.R')

SIRS_EAKF_model <- R6Class(
  inherit = ContestModel,
  
  private = list(
    .data = NULL, # ILIplus
    .initial = NULL, # initial states
    .num_ens = integer(0), # number of ensembles
    .num_times = integer(0), # number of weeks in each season
    .N = integer(0), # total population
    .scale = integer(0), # scale parameter of ILIplus
    .seed = integer(0), # seed in SIRS model
    .discrete = integer(0), # discrete (0/1) in SIRS model
    .dt = integer(0), # time unit
    .lambda = integer(0), # inflation factor
    .tmstep = integer(0), # number of runs, far enough so model spreads
    
    .obs = c(), # number of obsevations
    .ts = NULL, # start date of that year
    .AH = c(), # AbsHumidity
    
    .fitData = NULL # fitted state variables at the end of num_times
  ),
  
  public = list(
    fit = function(training_data) {
      if ("fit" %in% private$.debug) {
        browser()
      } # include for debugging
      
      # load data
      private$.data <- training_data
      AbsHumidity <- private$.data$metaData$humidity_dat
      
      x_fit <- replicate(nrow(training_data$mat),matrix(0,7,private$.num_ens))
      
      for (city_id in 1:nrow(training_data$mat)) {
      
      # load ILI+, prepare observations
      # city_id <- which(training_data$rnames == private$.city)
      ILIp <- private$.data$subset(rows = city_id, mutate = FALSE)
      ILIp <- as.numeric(ILIp$mat)
      private$.obs <- ILIp * private$.scale
      
      # load AbsHumidity
      AH <- AbsHumidity[, city_id]
      AH <- rep(AH, ceiling(length(ILIp)*7/365))
      private$.AH <- AH
      
      year <- substr(training_data$cnames, 2, 5)
      month <- substr(training_data$cnames, 20, 21)
      day <- substr(training_data$cnames, 23, 24)
      date <- paste0(year,"/",month,"/",day)
      private$.ts <- yday(date[1]) # start date number in the year

      # mapping operator
      HH <- diag(7)
      H <- HH[3, ]
      
      x <- private$.initial
      x <- SIRS(
        x, private$.ts, private$.dt, private$.tmstep,
        private$.N, private$.AH, private$.seed, private$.discrete
      )
      
      ### Begin looping through observations
      private$.num_times <- length(ILIp)
      xprior <- replicate(private$.num_times, matrix(NaN, 7, private$.num_ens))
      xpost <- xprior
      for (tt in 1:private$.num_times) {
        # Get the variance of the ensemble
        ave <- private$.obs[max(1, tt - 1)]
        ave <- ave + private$.obs[max(1, tt - 2)]
        ave <- ave + private$.obs[max(1, tt - 3)]
        ave <- ave / 3
        obs_var <- 10^5 + ave^2 / 5
        # inflation of x before assimilation
        x <- as.matrix(apply(x, 1, mean)) %*% rep(1, private$.num_ens) +
          private$.lambda * (x - as.matrix(apply(x, 1, mean)) %*% rep(1, private$.num_ens))
        prior_var <- var(as.vector(H %*% x))
        post_var <- prior_var * obs_var / (prior_var + obs_var)
        if (prior_var == 0) {
          post_var <- 0
        }
        prior_mean <- mean(H %*% x)
        post_mean <- post_var * (prior_mean / prior_var + private$.obs[tt] / obs_var)
        # Compute alpha and adjust distribution to conform to posterior moments
        alpha <- (obs_var / (obs_var + prior_var))^0.5
        dy <- post_mean + alpha * ((H %*% x) - prior_mean) - H %*% x
        #  Loop over each state variable
        rr <- matrix(0, 1, dim(x)[1])
        for (j in 1:dim(x)[1]) {
          A <- cov(x[j, ], as.vector(H %*% x))
          rr[j] <- A / prior_var
        }
        dx <- t(rr) %*% dy
        #  Get the new ensemble and save prior and posterior
        xprior[, , tt] <- x
        xnew <- x + dx
        #  Corrections to DA produced aphysicalities
        xnew <- checkDA(xnew)
        xpost[, , tt] <- xnew
        #  Integrate forward one time step
        tcurrent <- private$.ts + private$.tmstep * tt
        x <- SIRS(
          xnew, tcurrent, private$.dt, private$.tmstep,
          private$.N, private$.AH, private$.seed, private$.discrete
        )
        # EAKF update is finished
      }
       x_fit[,,city_id] <- x
      }
      private$.fitData <- x_fit
    },
    
    forecast = function(steps) {
      x_fit <- private$.fitData
      
      x_forecast <- replicate(nrow(training_data$mat),matrix(0,steps,private$.num_ens+1))
      
      for (city_id in 1:nrow(training_data$mat)) {

      forecasts <- matrix(0,steps,private$.num_ens+1)
      
      x <- x_fit[,,city_id]
      
      for (tt in 1:steps) {
        #  Integrate forward
        tcurrent <- private$.ts + private$.tmstep * (tt + private$.num_times)
        x <- SIRS(
          x, tcurrent, private$.dt, private$.tmstep,
          private$.N, private$.AH, private$.seed, private$.discrete
        )
        forecasts[tt,1:private$.num_ens] <- x[2,]
        forecasts[tt,private$.num_ens+1] <- rowMeans(x)[2]
      }
      x_forecast[,,city_id] <- forecasts
      }
      return(x_forecast)
    },
    
    initialize = function(num_ens = 100, N = 500000, lambda = 1.02, dt = 1, seed = 0.1,
                          scale = 1, discrete = 0, tmstep = 7, 
                          initial = initial) {
      private$.num_ens <- num_ens # number of ensembles
      private$.N <- N # total population
      private$.lambda <- lambda # inflation factor
      private$.dt <- dt # time unit
      private$.seed <- seed # seed in SIRS model
      private$.scale <- scale # scale parameter of ILIplus
      private$.discrete <- discrete # discrete (0/1) in SIRS model
      private$.tmstep <- tmstep # number of runs, far enough so model spreads
      private$.initial <- initial # initial states
    }
  )
)

dir <- paste0(getwd(), "/")
source(paste0(dir, "SIRS.R"))
source(paste0(dir, "checkDA.R"))


# load data
training_data <- readRDS(paste0(dir, "training_data.rds"))
initial <- as.matrix(read.csv(paste0(dir, "initials_original.csv"), header = TRUE))

# Create a new SIRS_EAKF model
SIRS.EAKF.model <- SIRS_EAKF_model$new(initial = initial)

### fit
SIRS.EAKF.model$fit(training_data)

### forecast
steps <- 6 # forecast ahead `step` number of weeks
SIRS.EAKF.model$forecast(steps = steps)
