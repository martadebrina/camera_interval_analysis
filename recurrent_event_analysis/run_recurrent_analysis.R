# =========================================================
# RECURRENT EVENT ANALYSIS
# =========================================================


# =========================================================
# LOAD PACKAGES
# =========================================================

library(readr)
library(dplyr)
library(lubridate)
library(ctrecurrent)
library(pammtools)
library(ggplot2)


# =========================================================
# LOAD CONFIGURATION
# =========================================================

source("recurrent_event_analysis/config_recurrent.R")


# =========================================================
# 1. LOAD CLEANED CAMERA DATA
# =========================================================

cat("\nLoading cleaned camera data...\n")

camera_data <- read_csv(
  cleaned_file,
  show_col_types = FALSE
)

cat(
  "Rows:",
  nrow(camera_data),
  "\n"
)


# =========================================================
# 2. CHECK DATA
# =========================================================

print(names(camera_data))

print(head(camera_data))

print(
  camera_data %>%
    count(Species, sort = TRUE)
)


# =========================================================
# 3. TRANSFORM TO RECURRENT EVENT DATA
# =========================================================

cat("\nTransforming to recurrent event format...\n")

recurrent_data <- ct_to_recurrent(
  data = camera_data,
  primary = primary_species,
  secondary = secondary_species,
  datetime_var = "Timestamp",
  species_var = "Species",
  site_var = "Site Name",
  survey_duration = survey_duration
)

cat(
  "Recurrent-event rows:",
  nrow(recurrent_data),
  "\n"
)

print(
  recurrent_data %>%
    select(
      `Site Name`,
      survey_id,
      primary,
      secondary,
      t.start,
      t.stop,
      event,
      status,
      enum
    ) %>%
    head(20)
)


# =========================================================
# 4. SUMMARIZE RECURRENT EVENT DATA
# =========================================================

cat("\nRecurrent event summary:\n")

cat(
  "Number of surveys:",
  n_distinct(recurrent_data$survey_id),
  "\n"
)

cat(
  "Number of snow leopard events:",
  sum(recurrent_data$event),
  "\n"
)

survey_summary <- recurrent_data %>%
  group_by(survey_id, `Site Name`) %>%
  summarise(
    n_events = sum(event),
    survey_length_days = max(t.stop),
    .groups = "drop"
  )

print(
  survey_summary %>%
    count(n_events, name = "n_surveys")
)

# =========================================================
# 5. TRANSFORM TO PED FORMAT
# =========================================================

cat("\nTransforming recurrent events to PED format...\n")

recurrent_data <- recurrent_data %>%
  rename(
    Site = `Site Name`
  )

ped <- recurrent_data %>%
  as_ped(
    formula = Surv(t.start, t.stop, event) ~ Site,
    id = "survey_id",
    transition = "enum",
    timescale = "calendar"
  )

cat(
  "PED rows:",
  nrow(ped),
  "\n"
)

print(
  ped %>%
    select(
      survey_id,
      Site,
      tstart,
      tend,
      enum,
      ped_status
    ) %>%
    head(20)
)


# =========================================================
# 6. FIT BASELINE PAMM
# =========================================================

cat("\nFitting baseline PAMM...\n")

model_null <- pamm(
  formula = ped_status ~ s(tend),
  data = ped
)

print(
  summary(model_null)
)


# =========================================================
# 7. PLOT HAZARD OVER TIME
# =========================================================

cat("\nCreating hazard plot...\n")

hazard_data <- ped %>%
  make_newdata(
    tend = unique(tend)
  ) %>%
  add_hazard(model_null)

hazard_plot <- ggplot(
  hazard_data,
  aes(
    x = tend,
    y = hazard
  )
) +
  geom_line() +
  geom_ribbon(
    aes(
      ymin = ci_lower,
      ymax = ci_upper
    ),
    alpha = 0.3
  ) +
  labs(
    title = paste(
      secondary_species,
      "occurrence after",
      primary_species
    ),
    x = paste(
      "Time after",
      primary_species,
      "detection (days)"
    ),
    y = paste(
      secondary_species,
      "hazard"
    )
  ) +
  theme_minimal()

print(hazard_plot)


# =========================================================
# 8. SAVE RESULTS
# =========================================================

result_dir <- file.path(
  "results",
  "recurrent_events",
  paste0(
    primary_species,
    "_",
    secondary_species
  )
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  recurrent_data,
  file.path(
    result_dir,
    "recurrent_events.csv"
  )
)

write_csv(
  hazard_data,
  file.path(
    result_dir,
    "hazard_estimates.csv"
  )
)

ggsave(
  file.path(
    result_dir,
    "hazard_plot.png"
  ),
  hazard_plot,
  width = 8,
  height = 5,
  dpi = 300
)

cat(
  "\nResults saved to:",
  result_dir,
  "\n"
)