# Recurrent Event Analysis

This folder contains an exploratory recurrent-event / PAMM analysis using the camera-trap data.

The workflow follows the `ctrecurrent` method, which converts camera-trap detections into recurrent-event surveys and then into PED format for PAMM analysis.

## Current analysis

The current test uses:

```r
primary_species <- "human"
secondary_species <- "snow leopard"
survey_duration <- 30
```

The goal is to look at how snow leopard occurrence changes with time after a human detection.

## What I did

1. Loaded the cleaned camera-trap data from `data/processed/`.

2. Used `human` as the primary species and `snow leopard` as the secondary species.

3. Used `ct_to_recurrent()` to transform the camera detections into recurrent-event surveys.

   A human detection starts a survey, and snow leopard detections after that are treated as recurrent events.

4. Summarized the recurrent-event data.

   Current result:

```text
32 surveys
44 snow leopard events
```

5. Used `as_ped()` to transform the recurrent-event data into Piece-wise Exponential Data (PED), which is the format needed for PAMM.

```text
1359 PED rows
```

6. Fitted a baseline PAMM:

```r
ped_status ~ s(tend)
```

This tests whether snow leopard occurrence changes with time after a human detection.

7. Generated a hazard plot showing the estimated snow leopard occurrence hazard over the 30-day period after a human detection.

The current result shows a decreasing hazard over time, but this is still preliminary and should not yet be interpreted as attraction or avoidance.

## File structure

```text
camera_interval_analysis/
│
├── recurrent_event_analysis/
│   ├── README.md
│   ├── config_recurrent.R
│   └── run_recurrent_analysis.R
│
├── data/
│   ├── raw/
│   └── processed/
│
└── results/
```

### `config_recurrent.R`

Contains the settings for the recurrent-event analysis:

```r
cleaned_file
primary_species
secondary_species
survey_duration
```

Change these settings when testing a different species pair or survey duration.

### `run_recurrent_analysis.R`

Runs the analysis:

```text
load cleaned data
        ↓
recurrent-event transformation
        ↓
PED transformation
        ↓
baseline PAMM
        ↓
hazard plot
```

### `data/processed/`

Contains the cleaned camera-trap data used as input.

### `results/`

Will contain outputs generated from the recurrent-event analysis.

## Current status

The basic recurrent-event transformation and baseline PAMM are working.

The next step is to inspect the recurrent-event surveys and review the survey assumptions, especially:

* whether other species should end/censor a survey
* whether 30 days is an appropriate survey duration

After that, the model can be refined and tested with other species pairs.
