# Camera Interval Analysis

This project calculates directional time intervals between two focal species detected at the same camera-trap site.

The workflow is based on the AB/BA interval approach used in Niedballa et al. for examining spatiotemporal relationships between species.

## How the analysis works

For two species:

* **A**
* **B**

the script calculates:

* **AB** = time from Species A to the next Species B detection
* **BA** = time from Species B to the next Species A detection

All intervals are calculated within the same camera site.

For example:

```text
A
A
A
B
B
A
```

The relationships are:

```text
last A → first B = AB

last B → next A = BA
```

Only transitions between the two focal species are used.

---

## Project structure

```text
camera_interval_analysis/
├── README.md
├── config.R
├── run_analysis.R
├── camera_interval_analysis.Rproj
│
├── data/
│   ├── raw/
│   └── processed/
│
├── scripts/
│   └── functions.R
│
├── results/
│
└── paper_reproduction/
```

### Main files

`config.R`
Contains the settings for the analysis, including the input file and focal species.

`run_analysis.R`
Runs the complete analysis.

`scripts/functions.R`
Contains the functions used for data cleaning, detection compression, and AB/BA calculations.

`data/raw/`
Contains the original camera-trap data.

`data/processed/`
Contains the cleaned event-level data.

`results/`
Contains the output for each species pair.

`paper_reproduction/`
Contains scripts used to reproduce and explore the Niedballa et al. simulation. These scripts are not required for the real-data analysis.

---

# Running the analysis

## 1. Open the project

Open:

```text
camera_interval_analysis.Rproj
```

in RStudio.

---

## 2. Install required packages

If needed, install:

```r
install.packages(c(
  "readr",
  "dplyr",
  "lubridate",
  "stringr",
  "ggplot2"
))
```

---

## 3. Put the raw CSV in `data/raw/`

The input CSV should contain at least:

```text
Site Name
Timestamp
Species
```

For example:

```text
data/raw/SPNP2017-18_Cycle1_TTE.csv
```

---

## 4. Edit `config.R`

Set the input file:

```r
raw_file <- "data/raw/SPNP2017-18_Cycle1_TTE.csv"
```

The current workflow compresses consecutive same-species detections within **1 minute**:

```r
collapse_seconds <- 60
```

Species aliases can also be standardized:

```r
species_aliases <- c(
  "Vehicles/Humans/Livestock|human" = "human"
)
```

Then choose the two focal species:

```r
species_A <- "human"
species_B <- "snow leopard"
```

In this example:

```text
AB = human → snow leopard
BA = snow leopard → human
```

To analyze another pair, only change `species_A` and `species_B`.

For example:

```r
species_A <- "snow leopard"
species_B <- "wolf"
```

---

## 5. Run the analysis

Open:

```text
run_analysis.R
```

and click **Source**.

The script will:

1. clean the camera-trap data
2. remove missing timestamps and `Species = "None"`
3. standardize configured species names
4. compress repeated same-species detections within 1 minute
5. calculate AB and BA intervals within each camera site
6. save summary tables and figures

---

# Outputs

The cleaned dataset is saved in:

```text
data/processed/
```

Results are saved separately for each species pair.

For example:

```text
results/human_snow_leopard/
```

contains:

```text
relationships.csv
summary.csv
histogram_days.png
boxplot_days.png
```

### `relationships.csv`

Contains each AB or BA transition and its time difference in:

```text
minutes
hours
days
```

### `summary.csv`

Contains summary statistics for AB and BA, including:

```text
number of intervals
mean
median
minimum
maximum
```

### Figures

`histogram_days.png` shows the distribution of time intervals.

`boxplot_days.png` compares AB and BA intervals.

---

# Quick start

For a new species pair:

1. Open `config.R`
2. Change `species_A` and `species_B`
3. Save the file
4. Run `run_analysis.R`

Example:

```r
species_A <- "snow leopard"
species_B <- "wolf"
```

Results will be saved automatically in:

```text
results/snow_leopard_wolf/
```

## Note

The AB/BA intervals describe temporal relationships between detections at the same camera site.

Differences between AB and BA should not automatically be interpreted as evidence of avoidance, competition, attraction, or causation.
