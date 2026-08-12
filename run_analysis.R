# =========================================================
# CAMERA-TRAP AB / BA INTERVAL ANALYSIS
# =========================================================


# =========================================================
# REQUIRED PACKAGES
# =========================================================

required_packages <- c(
  "readr",
  "dplyr",
  "lubridate",
  "stringr",
  "ggplot2"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them with:\n",
    "install.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "),
    "))"
  )
}

library(readr)
library(dplyr)
library(lubridate)
library(stringr)
library(ggplot2)


# =========================================================
# LOAD CONFIGURATION + FUNCTIONS
# =========================================================

source("config.R")
source("scripts/functions.R")


# =========================================================
# CHECK CONFIGURATION
# =========================================================

if (!file.exists(raw_file)) {
  stop(
    "Raw data file not found: ",
    raw_file
  )
}

if (species_A == species_B) {
  stop(
    "species_A and species_B must be different."
  )
}

if (collapse_seconds <= 0) {
  stop(
    "collapse_seconds must be greater than 0."
  )
}


# =========================================================
# 1. CLEAN + COMPRESS CAMERA DATA
# =========================================================

cat("\nCleaning camera data...\n")

cleaned <- clean_camera_data(
  raw_file = raw_file,
  collapse_seconds = collapse_seconds,
  species_aliases = species_aliases
)

cat(
  "Cleaned events:",
  nrow(cleaned),
  "\n"
)


# =========================================================
# CHECK FOCAL SPECIES
# =========================================================

available_species <- unique(cleaned$Species)

if (!species_A %in% available_species) {
  stop(
    "Species A not found in cleaned data: ",
    species_A
  )
}

if (!species_B %in% available_species) {
  stop(
    "Species B not found in cleaned data: ",
    species_B
  )
}

cat("\nIndependent focal-species events:\n")

cleaned %>%
  filter(
    Species %in% c(species_A, species_B)
  ) %>%
  count(Species) %>%
  print()



# =========================================================
# 2. CREATE OUTPUT FOLDERS
# =========================================================

dir.create(
  "data/processed",
  recursive = TRUE,
  showWarnings = FALSE
)


pair_name <- paste0(
  safe_name(species_A),
  "_",
  safe_name(species_B)
)


result_dir <- file.path(
  "results",
  pair_name
)


dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# =========================================================
# 3. SAVE CLEANED EVENT DATA
# =========================================================

cleaned_file <- file.path(
  "data/processed",
  paste0(
    tools::file_path_sans_ext(
      basename(raw_file)
    ),
    "_cleaned.csv"
  )
)


write_csv(
  cleaned,
  cleaned_file
)


cat(
  "Cleaned data saved to:",
  cleaned_file,
  "\n"
)



# =========================================================
# 4. BUILD AB / BA RELATIONSHIPS
# =========================================================

cat(
  "\nFocal pair:\n",
  "A =", species_A, "\n",
  "B =", species_B, "\n"
)


relationships <- build_relationships(
  cleaned_data = cleaned,
  species_A = species_A,
  species_B = species_B
)


cat(
  "\nNumber of relationships:",
  nrow(relationships),
  "\n"
)


print(
  relationships %>%
    count(direction)
)


# =========================================================
# 5. SAVE RELATIONSHIP TABLE
# =========================================================

write_csv(
  relationships,
  file.path(
    result_dir,
    "relationships.csv"
  )
)


# =========================================================
# 6. SUMMARY
# =========================================================

relationship_summary <-
  summarize_relationships(
    relationships
  )


print(relationship_summary)


write_csv(
  relationship_summary,
  file.path(
    result_dir,
    "summary.csv"
  )
)




# =========================================================
# 7. HISTOGRAM
# =========================================================

histogram <- ggplot(
  relationships,
  aes(
    x = time_difference_days
  )
) +
  
  geom_histogram(
    bins = 30
  ) +
  
  facet_wrap(
    ~ direction,
    ncol = 1
  ) +
  
  labs(
    
    title = paste(
      species_A,
      "and",
      species_B,
      "time intervals"
    ),
    
    subtitle = paste0(
      "AB = ",
      species_A,
      " → ",
      species_B,
      " | BA = ",
      species_B,
      " → ",
      species_A
    ),
    
    x = "Time difference (days)",
    
    y = "Number of intervals"
  ) +
  
  theme_minimal()


ggsave(
  file.path(
    result_dir,
    "histogram_days.png"
  ),
  histogram,
  width = 8,
  height = 7,
  dpi = 300
)



# =========================================================
# 8. BOXPLOT
# =========================================================

boxplot_graph <- ggplot(
  relationships,
  aes(
    x = direction,
    y = time_difference_days
  )
) +
  
  geom_boxplot() +
  
  scale_y_log10() +
  
  labs(
    
    title = paste(
      species_A,
      "and",
      species_B,
      "AB vs BA intervals"
    ),
    
    subtitle = paste0(
      "AB = ",
      species_A,
      " → ",
      species_B,
      " | BA = ",
      species_B,
      " → ",
      species_A
    ),
    
    x = "Direction",
    
    y = "Time difference (days, log scale)"
  ) +
  
  theme_minimal()


ggsave(
  file.path(
    result_dir,
    "boxplot_days.png"
  ),
  boxplot_graph,
  width = 7,
  height = 5,
  dpi = 300
)



# =========================================================
# DONE
# =========================================================

cat("\nAnalysis complete!\n")

cat(
  "Results saved in:",
  result_dir,
  "\n"
)



