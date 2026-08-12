# =========================================================
# USER CONFIGURATION
# =========================================================
#
# This is the main file users should edit when running
# the analysis on a new dataset or species pair.
# =========================================================


# ---------------------------------------------------------
# INPUT DATA
# ---------------------------------------------------------

raw_file <- "data/raw/SPNP2017-18_Cycle1_TTE.csv"


# ---------------------------------------------------------
# DETECTION COMPRESSION
# ---------------------------------------------------------
#
# Consecutive records of the same species at the same
# camera within this many seconds are treated as one event.
#
# 60 seconds = 1 minute
# ---------------------------------------------------------

collapse_seconds <- 60


# ---------------------------------------------------------
# SPECIES ALIASES
# ---------------------------------------------------------
#
# Different labels that should be treated as the same
# species/category.
#
# Format:
# "original label" = "new standardized label"
# ---------------------------------------------------------

species_aliases <- c(
  "Vehicles/Humans/Livestock|human" = "human"
)


# ---------------------------------------------------------
# FOCAL SPECIES
# ---------------------------------------------------------
#
# A = first species
# B = second species
#
# AB = A followed by B
# BA = B followed by A
# ---------------------------------------------------------

species_A <- "human"
species_B <- "snow leopard"