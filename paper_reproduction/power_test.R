library(CircStats)
library(dplyr)

source("scripts/simulation_function.R")


run_one_test <- function(records) {
  
  # -------------------------------------------------------
  # 1. Sort observations chronologically
  # -------------------------------------------------------
  
  records_sorted <- records[
    order(records$ObservationTime),
  ]
  
  # -------------------------------------------------------
  # 2. Remove same-species detections within 60 minutes
  # -------------------------------------------------------
  
  records_filtered <- records_sorted[0, ]
  
  last_A <- as.POSIXct(NA)
  last_B <- as.POSIXct(NA)
  
  for (i in seq_len(nrow(records_sorted))) {
    
    current_species <- records_sorted$Species[i]
    current_time <- records_sorted$ObservationTime[i]
    
    if (current_species == "Species A") {
      
      if (
        is.na(last_A) ||
        as.numeric(
          difftime(
            current_time,
            last_A,
            units = "mins"
          )
        ) > 60
      ) {
        
        records_filtered <- rbind(
          records_filtered,
          records_sorted[i, ]
        )
        
        last_A <- current_time
      }
      
    } else if (current_species == "Species B") {
      
      if (
        is.na(last_B) ||
        as.numeric(
          difftime(
            current_time,
            last_B,
            units = "mins"
          )
        ) > 60
      ) {
        
        records_filtered <- rbind(
          records_filtered,
          records_sorted[i, ]
        )
        
        last_B <- current_time
      }
    }
  }
  
  # -------------------------------------------------------
  # 3. Calculate AB and BA intervals
  # -------------------------------------------------------
  
  intervals <- data.frame()
  
  for (i in 1:(nrow(records_filtered) - 1)) {
    
    current_species <- records_filtered$Species[i]
    next_species <- records_filtered$Species[i + 1]
    
    if (current_species != next_species) {
      
      current_time <- records_filtered$ObservationTime[i]
      next_time <- records_filtered$ObservationTime[i + 1]
      
      interval_type <- ifelse(
        current_species == "Species A",
        "AB",
        "BA"
      )
      
      interval_hours <- as.numeric(
        difftime(
          next_time,
          current_time,
          units = "hours"
        )
      )
      
      intervals <- rbind(
        intervals,
        data.frame(
          interval_type = interval_type,
          interval_hours = interval_hours
        )
      )
    }
  }
  
  # -------------------------------------------------------
  # 4. Need both AB and BA to fit the model
  # -------------------------------------------------------
  
  if (
    nrow(intervals) < 2 ||
    length(unique(intervals$interval_type)) < 2
  ) {
    return(NA)
  }
  
  # log(0) is undefined
  if (any(intervals$interval_hours <= 0)) {
    return(NA)
  }
  
  # BA = reference, matching the paper
  intervals$interval_type <- factor(
    intervals$interval_type,
    levels = c("BA", "AB")
  )
  
  # -------------------------------------------------------
  # 5. Fit log-transformed linear model
  # -------------------------------------------------------
  
  model <- lm(
    log(interval_hours) ~ interval_type,
    data = intervals
  )
  
  # p-value for AB vs BA
  p_value <- summary(model)$coefficients[
    "interval_typeAB",
    "Pr(>|t|)"
  ]
  
  return(p_value)
}

set.seed(42)

n_simulations <- 100

p_values <- numeric(n_simulations)

for (s in 1:n_simulations) {
  
  records_sim <- simulateInteractionRecords(
    n_records_A = 20,
    n_records_B = 20,
    n_days = 100,
    effectDurationDays = 1,
    oddsRatio = 10,
    doPlots = FALSE,
    family = "von Mises",
    speciesOffsetHours = 0,
    densityFunctionParameters = list(
      mu = pi,
      kappa = 2
    )
  )
  
  p_values[s] <- run_one_test(records_sim)
  
  cat("simulation", s, "of", n_simulations, "\n")
}

valid_p <- p_values[!is.na(p_values)]

power <- mean(valid_p < 0.05)

power


sum(valid_p < 0.05)
length(valid_p)


power_summary <- data.frame(
  n_simulations = n_simulations,
  n_valid = length(valid_p),
  n_significant = sum(valid_p < 0.05),
  power = power
)

power_summary

hist(
  valid_p,
  breaks = 20,
  main = "P-values across simulations",
  xlab = "P-value"
)

abline(
  v = 0.05,
  col = "red",
  lwd = 2
)

