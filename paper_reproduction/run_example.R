library(CircStats)

source("scripts/simulation_function.R")

set.seed(42)

records <- simulateInteractionRecords(
  n_records_A = 20,
  n_records_B = 20,
  n_days = 30,
  effectDurationDays = 1,
  oddsRatio = 10,
  doPlots = TRUE,
  family = "von Mises",
  speciesOffsetHours = 0,
  densityFunctionParameters = list(
    mu = pi,
    kappa = 2
  )
)

View(records)

head(records)

dim(records)

records_sorted <- records[order(records$ObservationTime), ]

View(records_sorted)

records_sorted$record_id <- seq_len(nrow(records_sorted))

records_sorted <- records_sorted[
  c(
    "record_id",
    "Species",
    "ObservationTime",
    "Date",
    "Time",
    "Time_sec",
    "Time_rad"
  )
]

View(records_sorted)

write.csv(
  records_sorted,
  "output/simulated_records.csv",
  row.names = FALSE
)


# ---------------------------------------------------------
# Remove repeated same-species detections within 60 minutes
# ---------------------------------------------------------

records_filtered <- records_sorted[0, ]

last_A <- as.POSIXct(NA)
last_B <- as.POSIXct(NA)

for (i in 1:nrow(records_sorted)) {
  
  current_species <- records_sorted$Species[i]
  current_time <- records_sorted$ObservationTime[i]
  
  if (current_species == "Species A") {
    
    if (is.na(last_A) ||
        as.numeric(difftime(current_time, last_A, units = "mins")) > 60) {
      
      records_filtered <- rbind(
        records_filtered,
        records_sorted[i, ]
      )
      
      last_A <- current_time
    }
    
  } else if (current_species == "Species B") {
    
    if (is.na(last_B) ||
        as.numeric(difftime(current_time, last_B, units = "mins")) > 60) {
      
      records_filtered <- rbind(
        records_filtered,
        records_sorted[i, ]
      )
      
      last_B <- current_time
    }
  }
}

records_filtered$record_id <- seq_len(nrow(records_filtered))


View(records_filtered)


# =========================================================
# Calculate AB and BA intervals
# =========================================================

intervals <- data.frame()

for (i in 1:(nrow(records_filtered) - 1)) {
  
  current_species <- records_filtered$Species[i]
  next_species    <- records_filtered$Species[i + 1]
  
  current_time <- records_filtered$ObservationTime[i]
  next_time    <- records_filtered$ObservationTime[i + 1]
  
  # only calculate an interval when species changes
  if (current_species != next_species) {
    
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
        start_species = current_species,
        start_time = current_time,
        end_species = next_species,
        end_time = next_time,
        interval_hours = interval_hours
      )
    )
  }
}

View(intervals)

aggregate(
  interval_hours ~ interval_type,
  data = intervals,
  FUN = mean
)


write.csv(
  intervals,
  "output/AB_BA_intervals.csv",
  row.names = FALSE
)


# =========================================================
# Linear model: compare log-transformed AB and BA intervals
# =========================================================

# Make BA the reference category, matching the paper
intervals$interval_type <- factor(
  intervals$interval_type,
  levels = c("BA", "AB")
)

model_log <- lm(
  log(interval_hours) ~ interval_type,
  data = intervals
)

summary(model_log)

boxplot(
  log(interval_hours) ~ interval_type,
  data = intervals,
  xlab = "Interval type",
  ylab = "Log time interval (hours)",
  main = "AB vs BA intervals"
)

boxplot(
  interval_hours ~ interval_type,
  data = intervals,
  xlab = "Interval type",
  ylab = "Time interval (hours)",
  main = "AB vs BA intervals"
)
