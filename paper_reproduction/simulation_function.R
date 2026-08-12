
simulateInteractionRecords <- function(n_records_A,
                                       n_records_B,
                                       n_days,
                                       effectDurationDays,
                                       oddsRatio,
                                       doPlots = FALSE,
                                       family  = c("uniform", "von Mises", "von Mises Mixture"),
                                       speciesOffsetHours = 0,
                                       densityFunctionParameters
)
{
  
  require(CircStats)
  
  # check function input
  if(oddsRatio < 0) stop("oddsRatio must be be positive")
  if(oddsRatio < 1) stop("oddsRatio must be < 1 suggests attraction. This is not fully supported in this function.")
  
  family <- match.arg(family, choices =  c("uniform", "von Mises", "von Mises Mixture"))
  
  if(speciesOffsetHours > 24) stop("speciesOffsetHours can't be greater than 24")
  
  #  extract distrbution parameters of the given distribution function (uniform, von Mises, von Mises Mixture) from list
  #  and check if the distribution function arguments are correct
  
  if(family == "uniform")  {
    densityFunctionParameters <- list()
    mu    <- 0
    kappa <- 0
  } else {
    if(!hasArg(densityFunctionParameters))  stop("argument 'densityFunctionParameters' is missing")
    if(!is.list(densityFunctionParameters)) stop("'densityFunctionParameters' must be a list")
  }
  
  if(family == "von Mises") {
    stopifnot(all(c("mu", "kappa") %in% names(densityFunctionParameters)))
    
    if(densityFunctionParameters$mu < 0 | 
       densityFunctionParameters$mu > 2 * pi) stop("mu must be between 0 and 2*pi")
    
    mu    <- densityFunctionParameters$mu     # mu is mean direction of the von Mises distribution (in radians)
    kappa <- densityFunctionParameters$kappa  # kappa is concentration parameter of the von Mises distribution
  }
  
  if(family == "von Mises Mixture") {
    stopifnot(all(c("mu1", "mu2", "kappa1", "kappa2") %in% names(densityFunctionParameters)))
    
    if(densityFunctionParameters$mu1 < 0 | 
       densityFunctionParameters$mu1 > 2 * pi) stop("mu1 must be between 0 and 2*pi")
    if(densityFunctionParameters$mu2 < 0 | 
       densityFunctionParameters$mu2 > 2 * pi) stop("mu2 must be between 0 and 2*pi")
    
    mu1    <- densityFunctionParameters$mu1      # mean direction of the first von Mises distribution (in radians)
    mu2    <- densityFunctionParameters$mu2      # mean direction of the second von Mises distribution (in radians)
    kappa1 <- densityFunctionParameters$kappa1   # concentration parameter of the first von Mises distribution
    kappa2 <- densityFunctionParameters$kappa2   # concentration parameter of the second von Mises distribution
  }
  
  
  # generate a time sequence with 1-minute intervals for 1 day (00:00 - 23:59) in radial and clock time
  seq0 <- seq(0, (2 * pi), length.out = 1441)    # generate values from 0 to 2pi ((00:00 - 00:00 the next day, hence 1441 time steps) in radial time)
  seq0 <- seq0[-length(seq0)]                    # remove last value to crop that sequence to 00:00 - 23:59 in radial time
  
  species_offset_rad <- speciesOffsetHours / 24 * 2 * pi   # offset of mu (activity peak) for species B if speciesOffsetHours is defined
  
  # calculate probability densities of daily activity patterns of species A and B for the given probability distribution and its parameters
  # this is for each minute of the day, and only incorporates daily activity patterns, not yet spatiotemporal avoidance
  # if species activity peak offset = 0, density distributions are identical for both species
  
  # unimodal and uniform activiy patterns - von Mises Distributions (von Mises distribution is a uniform distribution if kappa = 0)
  if(family %in% c("uniform", "von Mises")){
    
    # for species A
    density_distribution_A <- dvm(theta = seq0,
                                  mu    = mu,
                                  kappa = kappa)
    
    if(speciesOffsetHours == 0){
      density_distribution_B <- density_distribution_A
    } else {
      if(mu + species_offset_rad > 2*pi) stop("mu + speciesOffsetHours is > 2 * pi")
      
      # for species B
      density_distribution_B <- dvm(theta = seq0,
                                    mu    = mu + species_offset_rad,
                                    kappa = kappa)
    }
  }
  
  # bimodal activiy pattern - mixture of 2 von Mises Distributions
  if(family == "von Mises Mixture"){
    # for species A
    density_distribution_A <- dmixedvm(theta  = seq0,
                                       mu1    = mu1,
                                       mu2    = mu2,
                                       kappa1 = kappa1,
                                       kappa2 = kappa2,
                                       p      = 0.5)
    
    if(speciesOffsetHours == 0){
      density_distribution_B <- density_distribution_A
    } else {
      if(mu1 + species_offset_rad > 2*pi) stop("mu1 + speciesOffsetHours is > 2 * pi")
      if(mu2 + species_offset_rad > 2*pi) stop("mu2 + speciesOffsetHours is > 2 * pi")
      
      # for species B
      density_distribution_B <- dmixedvm(theta  = seq0,
                                         mu1    = mu1 + species_offset_rad,
                                         mu2    = mu2 + species_offset_rad,
                                         kappa1 = kappa1,
                                         kappa2 = kappa2,
                                         p      = 0.5)
    }
  }
  
  
  
  # generate a sequence of time points with the specified number of days, each with 1440 minutes, in radial time.
  
  # create a sequence of dates for each minute of the study period (needed for the date/time functions)
  # the day of the origin is generic and does not influence results
  tz       <- "UTC"             # use a generic time zone for the date/time objects (main reason: no daylight saving time)
  date0    <- base::as.Date(rep(seq(0, n_days - 1, by = 1), each = length(seq0)), origin = "1970-01-01",  tz = tz)   # 0... n_days - 1 because origin = day 0
  n_events <- length(date0)     # the number of events (n_days * 1440 minutes)
  
  # calculate observation probabilites for each minute in the study period incorporating daily activity patterns
  # this is a relative measure, not an actual probability, and will be used as a probability weight for the random sample which realises the species records
  # takes into account daily activity patterns, but not yet avoidance
  
  observation_prob_A <- rep(1, times = n_events) * density_distribution_A
  observation_prob_B <- rep(1, times = n_events) * density_distribution_B
  
  # generate records of primary species A conditional on observation probability of A
  time_observations_A <- sort(sample(x       = seq(1, n_events),
                                     size    = n_records_A,
                                     prob    = observation_prob_A,
                                     replace = FALSE))
  
  # modify observation probability for subordinate species B dependent on records of species A and the strength of the avoidance
  
  # find all times (1-minute intervals) in which B is affected by A
  which_have_reduced_p <- unique(as.vector(sapply(time_observations_A, FUN = function(X){seq(from = X, to = X + (effectDurationDays * 1440))})))
  
  # remove those which are after the end of the study period
  if(max(which_have_reduced_p) > n_events){
    which_have_reduced_p <- which_have_reduced_p[-which(which_have_reduced_p > n_events)]
  }
  
  # for each time point affected by observation of A, calculate time since last record of A
  delta_time_since_last_record <- sapply(time_observations_A, FUN = function(x) {which_have_reduced_p - x})
  delta_time_since_last_record[delta_time_since_last_record < 0] <- NA                              # remove all values before first record
  distance_from_event <- apply(delta_time_since_last_record, MARGIN = 1, FUN = min, na.rm = TRUE)   # time steps since last observation of A
  
  # minutes until full recovery of detection probability (for next step)
  n_steps_until_recovery <- effectDurationDays * 1440
  
  # calculate observation probability of B taking into account linear recovery after observation of A
  observation_prob_B[which_have_reduced_p] <- observation_prob_B[which_have_reduced_p] / oddsRatio +
    (distance_from_event / n_steps_until_recovery) * (1 - (1 / oddsRatio)) * observation_prob_B[which_have_reduced_p]
  
  # generate records of subordinate species B conditional on observation probability of B and avoidance of A
  time_observations_B <- sort(sample(x       = seq(1, n_events),
                                     size    = n_records_B,
                                     replace = FALSE,
                                     prob    = observation_prob_B))
  
  # make output data frame with observations of species
  ObservationTime <- as.POSIXct(c(time_observations_A, time_observations_B) * 60,
                                origin="1970-01-01", tz=tz)
  
  outtable <- data.frame(Species          = rep(c("Species A", "Species B"), times = c(n_records_A, n_records_B)),
                         ObservationTime  = ObservationTime,
                         Date             = base::as.Date(ObservationTime, tz = tz),
                         Time             = format(ObservationTime, format = "%H:%M:%S"),
                         Time_sec         = as.numeric(format(ObservationTime, format = "%s")),
                         Time_rad         = ClocktimeToRadialTime(ObservationTime))
  
  
  # create plots of the observation probabilities of both species
  
  if(isTRUE(doPlots)){
    
    # reset graphics parameters on exiting function
    mfrow0 <- par()$mfrow
    on.exit(par(mfrow = mfrow0))
    
    # generate title and subtitles for plots
    main_title <- paste("family = ", family,
                        ";  odds ratio = ", oddsRatio, sep = "")
    subtitle   <- paste("n_records A = ", n_records_A, "; total length = ", n_days, " days", ";  effect_duration = ", effectDurationDays, " days", sep = "")
    subtitle2  <- paste("n_records B = ", n_records_B, "; speciesOffsetHours = ", speciesOffsetHours, sep = "")
    
    if(family == "von Mises") subtitle <- paste(subtitle, "; mu = ", round(mu, 2), " (= ", round(mu * 24/(2*pi)), " o'clock); kappa = ", kappa,
                                                "; odds ratio (B(no A) / B(A)) =", oddsRatio , sep = "")
    if(family == "von Mises Mixture")  subtitle <- paste(subtitle, "; mu1 = ", round(mu1, 2), " (= ", round(mu1 * 24/(2*pi)), " o'clock)",
                                                         "; mu2 = ", round(mu2, 2), " (= ", round(mu2 * 24/(2*pi)), " o'clock)",
                                                         "; kappa1 = ", kappa1, "; kappa2 = ", kappa2,
                                                         "; odds ratio (B(no A) / B(A)) =", oddsRatio, sep = "")
    
    # set graphical parameters
    col_abline  <- rgb(0, 0, 0, 0.2)     # light grey, transparent
    col_abline2 <- rgb(0, 0, 0, 0.1)     # light grey, even more transparent
    col_rug_A   <- "red"                 # rug colour for species A
    col_rug_B   <- "blue"                # rug colour for species B
    lwd_rug     <- 3                     # line width of the rugs
    par(mfrow = c(2,1))                  # plot layout: 2 rows, 1 column
    
    
    # define some plot elements shared between both plots
    x_axis_label_location       <- seq(720, length(observation_prob_A), by = 1440)                            # location for x axis label (noon)
    x_axis_midnight_locations <- seq(1,   length(observation_prob_A), by = 1440)                            # location for abline at midnight each day
    polygonToPlot_A <- data.frame(x = c(1, 1:length(observation_prob_A), length(observation_prob_A), 1),    # coordinates for a polygon to fill the area under the curve for species A
                                  y = c(0, observation_prob_A, 0, 0))
    polygonToPlot_B <- data.frame(x = c(1, 1:length(observation_prob_B), length(observation_prob_B), 1),    # coordinates for a polygon to fill the area under the curve for species B
                                  y = c(0, observation_prob_B, 0, 0))
    
    
    # create top plot: primary species A
    plot(observation_prob_A, type = "l", axes = F, ylim = c(0, max(observation_prob_A)),
         main = main_title, sub = subtitle, xlab = "", ylab = "probability weight (species A)")
    axis(1, at = x_axis_midnight_locations, labels = FALSE, tick = TRUE)             		# make ticks for days along x axis
    axis(1, at = x_axis_label_location, labels = paste("day", seq(1, (length(x_axis_label_location)))), tick = FALSE)   # day labels along x axis
    axis(2)                                                                          		# add y axis
    abline(v = x_axis_midnight_locations, col = col_abline)                           	# vertical lines separating days
    rug(time_observations_A, lwd = lwd_rug, col = col_rug_A); box()       			  		  # add rug to indicate observation of species A
    polygon(x = polygonToPlot_A$x, y = polygonToPlot_A$y, 								             	# plot the polygon for species A
            border = NA, col = col_abline)
    
    # create bottom plot: subordinate species B
    plot(observation_prob_B, type = "l", axes = F, ylim = c(0, max(observation_prob_B)), xlab = "", ylab = "probability weight (species B)", sub = subtitle2)
    axis(1, at = x_axis_midnight_locations, labels = FALSE, tick = TRUE)    				# make ticks for days along x axis
    axis(1, at = x_axis_label_location, labels = paste("day", seq(1, (length(x_axis_label_location)))), tick = FALSE)   # day labels along x axis
    axis(2)																					# add y axis
    abline(v = x_axis_midnight_locations, col = col_abline)									# vertical lines separating days
    rug(time_observations_A, lwd = lwd_rug, col = col_rug_A)								# add rug to indicate observation of species A
    rug(time_observations_B, lwd = lwd_rug, col = col_rug_B); box()					# add rug to indicate observation of species B
    polygon(x = polygonToPlot_A$x, y = polygonToPlot_A$y, 									# plot the polygon for species A
            border = col_abline, col = col_abline2)
    polygon(x = polygonToPlot_B$x, y = polygonToPlot_B$y, 									# plot the polygon for species A
            border = NA, col = col_abline)
    
  }
  # return data frame with record times
  return(outtable)
}


# Helper function to  convert clock time to radial time (0...2*pi)

ClocktimeToRadialTime <- function(Clocktime,
                                  timeformat = "%Y-%m-%d %H:%M:%S"
){
  DateTime2 <- strptime(as.character(Clocktime), format = timeformat, tz = "UTC")
  Time2     <- format(DateTime2, format = "%H:%M:%S", usetz = FALSE)
  Time.rad  <- (as.numeric(as.POSIXct(strptime(Time2, format = "%H:%M:%S", tz = "UTC"))) -
                  as.numeric(as.POSIXct(strptime("0", format = "%S", tz = "UTC")))) / 3600 * (pi/12)
  return(Time.rad)
}