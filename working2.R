library(ForecastFramework)
library(R6)
library(forecast)

dir <- paste0(getwd(), "/")
source(paste0(dir, "SIRS.R"))
source(paste0(dir, "checkDA.R"))

# source_github <- function(u) {
#   library(RCurl)
#   # read script lines from website and evaluate
#   script <- getURL(u, ssl.verifypeer = FALSE)
#   eval(parse(text = script),envir=.GlobalEnv)
# }
#

SIRS_EAKF_model <- R6Class(
  inherit = ContestModel,

  private = list(
    .data = NULL, # ILIplus
    .AbsHumidity = NULL, # humidity
    .phi = NULL, # initial states

    # .models = list(),    # model for each location

    # .num_pred=integer(0), # number of prediction realizations
    .city = NULL, # city whose ILIplus is predicted
    .year = integer(0), # year starts training from
    .num_ens = integer(0), # number of ensembles
    .num_times = integer(0), # number of weeks in each season
    .N = integer(0), # total population
    .scale = integer(0), # scale parameter of ILIplus
    .seed = integer(0), # seed in SIRS model
    .discrete = integer(0), # discrete (0/1) in SIRS model
    .dt = integer(0), # time unit
    .lambda = integer(0), # inflation factor
    .tmstep = integer(0), # number of runs, far enough so model spreads
    .Tpkpred = c(), # peak week
    .obs = c() # number of obsevations
  ),

  public = list(
    fit = function(data, AbsHumidity, phi) {
      if ("fit" %in% private$.debug) {
        browser()
      } # include for debugging

      # load data
      private$.data <- IncidenceMatrix$new(data)
      private$.AbsHumidity <- AbsHumidity
      private$.phi <- phi

      # load AbsHumidity
      AH <- private$.AbsHumidity[, private$.city]
      AH <- c(AH, AH)

      # load ILI+, prepare observations
      ILIp <- private$.data[, private$.city]
      ILIp <- ILIp * private$.scale
      obstime <- ILIplus[, "Date"]
      obstimeday <- as.numeric(as.Date(obstime, format = "%m/%d/%y"))
      s <- paste0("10/01/", year) # start training from Oct 1st
      temp <- which(obstimeday > as.numeric(as.Date(s, format = "%m/%d/%y")))
      obs <- ILIp[temp[1]:(temp[1] + num_times - 1)] # observations
      startdate <- obstime[temp[1]]
      s <- paste0("01/01/", private$.year) # the first day of the year
      ts <- as.numeric(as.Date(startdate, format = "%m/%d/%y")) -
        as.numeric(as.Date(s, format = "%m/%d/%y")) + 1 # start date number in the year

      # load the initial conditions for ensemble members
      susceps <- matrix(as.matrix(phi[, (5:31)]), nrow = 100000 * 27, ncol = 1)
      infects <- matrix(as.matrix(phi[, (32:58)]), nrow = 100000 * 27, ncol = 1)
      params <- phi[, (1:4)]

      # output prediction results of EAKFC
      Tpkpred <- rep(0, private$.num_times)

      # mapping operator
      HH <- diag(7)
      H <- HH[3, ]

      rnd <- cbind(
        ceiling(1e5 * runif(100)),
        ceiling(27 * 1e5 * runif(100)),
        ceiling(27 * 1e5 * runif(100))
      )
      x <- matrix(0, 7, private$.num_ens)
      x[1, ] <- susceps[rnd[, 2], 1] # S
      x[2, ] <- infects[rnd[, 3], 1] # I
      x[4, ] <- params[rnd[, 1], 3] # R0max
      x[5, ] <- params[rnd[, 1], 4] # R0min
      x[6, ] <- params[rnd[, 1], 1] * 365 # L
      x[7, ] <- params[rnd[, 1], 2] * 365 # D

      x <- SIRS(
        x, ts, private$.dt, private$.tmstep,
        private$.N, AH, private$.seed, private$.discrete
      )
      ### Begin looping through observations
      xprior <- replicate(private$.num_times, matrix(NaN, 7, private$.num_ens))
      xpost <- xprior
      for (tt in 1:private$.num_times) {
        # Get the variance of the ensemble
        ave <- obs[max(1, tt - 1)]
        ave <- ave + obs[max(1, tt - 2)]
        ave <- ave + obs[max(1, tt - 3)]
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
        post_mean <- post_var * (prior_mean / prior_var + obs[tt] / obs_var)
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
        tcurrent <- ts + private$.tmstep * tt
        x <- SIRS(
          xnew, tcurrent, private$.dt, private$.tmstep,
          private$.N, AH, private$.seed, private$.discrete
        )
        # EAKF update is finished

        # Forecast and evaluation
        xpred <- replicate(private$.num_times + 1, matrix(0, 7, private$.num_ens))
        # assign xpred before tt as posterior
        for (t in 1:tt) {
          xpred[, , t] <- xpost[, , t]
          xpred[3, , t] <- obs[t] * rep(1, private$.num_ens) # assign as true observations
        }
        ####################

        obspred <- aperm(xpred, c(2, 3, 1))
        obspred <- Re(obspred[, , 3]) # num_ens*num_times
        # evaluation
        # peak week
        pwens <- rep(NaN, num_ens)
        for (i in 1:num_ens) {
          temp <- which(obspred[i, ] == max(obspred[i, ]))
          pwens[i] <- temp[1]
        }
        pwp <- floor(mean(pwens))
        Tpkpred[tt] <- pwp
      }
      private$.Tpkpred <- Tpkpred
      private$.obs <- obs
    },

    forecast = function(steps) {
      Tpkpred <- as.numeric(private$.Tpkpred)
      Peakweekprediction <- cbind(private$.obs, Tpkpred)
      Peakweekprediction <- Peakweekprediction[1:steps, ]
      return(Peakweekprediction)
    },

    initialize = function(num_ens = 100, num_times = 40, N = 500000, lambda = 1.02, dt = 1, seed = 0.1,
                              scale = 1, discrete = 0, year = 05,
                              city = "NewYorkNY", tmstep = 7) {
      private$.num_ens <- num_ens # number of ensembles
      private$.num_times <- num_times # number of weeks in each season
      private$.N <- N # total population

      private$.lambda <- lambda # inflation factor
      private$.dt <- dt # time unit
      private$.seed <- seed # seed in SIRS model
      private$.scale <- scale # scale parameter of ILIplus
      private$.discrete <- discrete # discrete (0/1) in SIRS model
      private$.year <- year # year starts training from
      private$.city <- city # city whose ILIplus is predicted
      private$.tmstep <- tmstep # number of runs, far enough so model spreads
    }
  )
)

# load data
AbsHumidity <- read.csv(paste0(dir, "AbsHumidity.csv"))
data <- read.csv(paste0(dir, "ILIplus.csv"))
phi <- read.csv(paste0(dir, "params_statespace_RK_ext_seed_phi.csv"), header = FALSE)

# Create a new SIRS_EAKF model
SIRS_EAKF_model <- SIRS_EAKF_model$new()
#######ERROR###########

### fit
SIRS_EAKF_model$fit(data, AbsHumidity, phi)

### forecast
steps <- 15 # forecast ahead `step` number of weeks
forecast_X <- SIRS_EAKF_model$forecast(fitData, steps = steps)