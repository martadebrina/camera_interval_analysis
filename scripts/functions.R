# =========================================================
# CAMERA INTERVAL ANALYSIS FUNCTIONS
# =========================================================


# =========================================================
# 1. CLEAN + COMPRESS CAMERA DATA
# =========================================================

clean_camera_data <- function(
    raw_file,
    collapse_seconds = 60,
    species_aliases = NULL
) {
  
  # -------------------------------------------------------
  # Read raw CSV
  # -------------------------------------------------------
  
  raw <- readr::read_csv(
    raw_file,
    show_col_types = FALSE
  )
  
  
  # -------------------------------------------------------
  # Parse timestamp and clean species labels
  # -------------------------------------------------------
  
  data <- raw %>%
    dplyr::mutate(
      
      datetime = lubridate::ymd_hms(
        Timestamp,
        tz = "UTC",
        quiet = TRUE
      ),
      
      species_clean = stringr::str_trim(Species)
    )
  
  
  # -------------------------------------------------------
  # Apply species aliases
  # -------------------------------------------------------
  
  if (!is.null(species_aliases)) {
    
    for (old_label in names(species_aliases)) {
      
      new_label <- species_aliases[[old_label]]
      
      data$species_clean[
        data$species_clean == old_label
      ] <- new_label
    }
  }
  
  
  # -------------------------------------------------------
  # Remove unusable records
  # -------------------------------------------------------
  
  data <- data %>%
    dplyr::filter(
      !is.na(datetime),
      !is.na(`Site Name`),
      !is.na(species_clean),
      species_clean != "None"
    )
  
  
  # -------------------------------------------------------
  # Compress repeated same-species detections
  # -------------------------------------------------------
  
  compressed <- data %>%
    
    dplyr::arrange(
      `Site Name`,
      species_clean,
      datetime
    ) %>%
    
    dplyr::group_by(
      `Site Name`,
      species_clean
    ) %>%
    
    dplyr::mutate(
      
      gap_seconds = as.numeric(
        difftime(
          datetime,
          dplyr::lag(datetime),
          units = "secs"
        )
      ),
      
      new_event =
        is.na(gap_seconds) |
        gap_seconds > collapse_seconds,
      
      event_id = cumsum(new_event)
    ) %>%
    
    dplyr::group_by(
      `Site Name`,
      species_clean,
      event_id
    ) %>%
    
    dplyr::summarise(
      Timestamp = min(datetime),
      .groups = "drop"
    ) %>%
    
    dplyr::rename(
      Species = species_clean
    ) %>%
    
    dplyr::arrange(
      `Site Name`,
      Timestamp
    ) %>%
    
    dplyr::select(
      `Site Name`,
      Timestamp,
      Species
    )
  
  
  return(compressed)
}


# =========================================================
# 2. BUILD AB / BA RELATIONSHIPS
# =========================================================

build_relationships <- function(
    cleaned_data,
    species_A,
    species_B
) {
  
  # -------------------------------------------------------
  # Keep only focal species
  # -------------------------------------------------------
  
  focal <- cleaned_data %>%
    
    dplyr::filter(
      Species %in% c(
        species_A,
        species_B
      )
    ) %>%
    
    dplyr::mutate(
      species_code = dplyr::case_when(
        Species == species_A ~ "A",
        Species == species_B ~ "B"
      )
    ) %>%
    
    dplyr::arrange(
      `Site Name`,
      Timestamp
    )
  
  
  # -------------------------------------------------------
  # Calculate transitions within each camera
  # -------------------------------------------------------
  
  relationships <- focal %>%
    
    dplyr::group_by(`Site Name`) %>%
    
    dplyr::arrange(
      Timestamp,
      .by_group = TRUE
    ) %>%
    
    dplyr::mutate(
      previous_species = dplyr::lag(Species),
      previous_code = dplyr::lag(species_code),
      previous_time = dplyr::lag(Timestamp)
    ) %>%
    
    # Paper-style approach:
    # keep only transitions from one focal species
    # to the other.
    dplyr::filter(
      !is.na(previous_species),
      Species != previous_species
    ) %>%
    
    dplyr::mutate(
      
      direction = paste0(
        previous_code,
        species_code
      ),
      
      time_difference_minutes = as.numeric(
        difftime(
          Timestamp,
          previous_time,
          units = "mins"
        )
      ),
      
      time_difference_hours =
        time_difference_minutes / 60,
      
      time_difference_days =
        time_difference_hours / 24
    ) %>%
    
    dplyr::ungroup() %>%
    
    dplyr::transmute(
      
      site_name = `Site Name`,
      
      direction,
      
      start_species = previous_species,
      start_time = previous_time,
      
      end_species = Species,
      end_time = Timestamp,
      
      time_difference_minutes,
      time_difference_hours,
      time_difference_days
    )
  
  
  return(relationships)
}


# =========================================================
# 3. SUMMARIZE RELATIONSHIPS
# =========================================================

summarize_relationships <- function(
    relationships
) {
  
  summary_table <- relationships %>%
    
    dplyr::group_by(direction) %>%
    
    dplyr::summarise(
      
      n = dplyr::n(),
      
      mean_hours =
        mean(time_difference_hours),
      
      median_hours =
        median(time_difference_hours),
      
      mean_days =
        mean(time_difference_days),
      
      median_days =
        median(time_difference_days),
      
      min_days =
        min(time_difference_days),
      
      max_days =
        max(time_difference_days),
      
      .groups = "drop"
    )
  
  
  return(summary_table)
}


# =========================================================
# 4. CREATE SAFE FILE/FOLDER NAMES
# =========================================================

safe_name <- function(x) {
  
  x %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all(
      "[^a-z0-9]+",
      "_"
    ) %>%
    stringr::str_replace_all(
      "^_|_$",
      ""
    )
}