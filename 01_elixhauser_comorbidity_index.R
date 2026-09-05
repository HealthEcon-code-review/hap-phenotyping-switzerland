# Phenotyping Hospital-Acquired Pneumonia in Administrative Data
# 01 – Elixhauser comorbidity index
#
# Expected input:  dia_long.RData containing an object named `data_long`
# Created output: dia_eci.RData containing an object named `dia_eci`
#
# Set PNEUMONIA_DATA_DIR locally to the authorised data directory. For example,
# add the following line to a local .Renviron file that is not committed:
# PNEUMONIA_DATA_DIR=/path/to/authorised/data

suppressPackageStartupMessages({
  library(dplyr)
  library(comorbidity)
})
 
# Analysis settings ------------------------------------------------------------

last_study_year <- 2019L

apply_severity_hierarchy <- FALSE


# File paths ------------------------------------------------------------------

data_dir <- Sys.getenv("PNEUMONIA_DATA_DIR", unset = "")

if (!nzchar(data_dir)) {
  stop(
    "PNEUMONIA_DATA_DIR is not set. Point it to the authorised data ",
    "directory containing dia_long.RData."
  )
}

data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
input_file <- file.path(data_dir, "dia_long.RData")
output_file <- file.path(data_dir, "dia_eci.RData")

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}


# Load and validate the diagnosis data ----------------------------------------

input_environment <- new.env(parent = emptyenv())
loaded_objects <- load(input_file, envir = input_environment)

if (!"data_long" %in% loaded_objects) {
  stop("The input file must contain an object named `data_long`.")
}

data_long <- input_environment$data_long

required_columns <- c("patient_id", "id_jahr", "dia")
missing_columns <- setdiff(required_columns, names(data_long))

if (length(missing_columns)) {
  stop(
    "Missing required columns in `data_long`: ",
    paste(missing_columns, collapse = ", ")
  )
}

data_long <- data_long %>%
  mutate(
    patient_id = as.character(patient_id),
    id_jahr = as.integer(id_jahr),
    dia = as.character(dia)
  )


# Construct a unique patient-year identifier ----------------------------------

patient_years <- data_long %>%
  filter(id_jahr <= last_study_year) %>%
  distinct(patient_id, id_jahr) %>%
  arrange(patient_id, id_jahr) %>%
  mutate(patient_year_id = row_number())

diagnosis_input <- data_long %>%
  filter(
    id_jahr <= last_study_year,
    !is.na(dia),
    nzchar(trimws(dia))
  ) %>%
  inner_join(
    patient_years,
    by = c("patient_id", "id_jahr")
  ) %>%
  transmute(
    patient_year_id,
    diagnosis_code = dia
  ) %>%
  distinct()

if (!nrow(diagnosis_input)) {
  stop("No non-missing diagnosis codes were found in the study period.")
}


# Calculate the Swiss-weighted Elixhauser score -------------------------------

# ICD-10 mapping: Quan et al. (2005)
# Swiss weights: Sharma et al. (2021)

elixhauser_domains <- comorbidity::comorbidity(
  x = diagnosis_input,
  id = "patient_year_id",
  code = "diagnosis_code",
  map = "elixhauser_icd10_quan",
  assign0 = apply_severity_hierarchy,
  labelled = FALSE,
  tidy.codes = TRUE
)

elixhauser_scores <- tibble(
  patient_year_id = elixhauser_domains$patient_year_id,
  eci = as.integer(
    comorbidity::score(
      x = elixhauser_domains,
      weights = "swiss",
      assign0 = apply_severity_hierarchy
    )
  )
)

patient_year_scores <- patient_years %>%
  left_join(elixhauser_scores, by = "patient_year_id") %>%
  mutate(eci = coalesce(eci, 0L)) %>%
  select(patient_id, id_jahr, eci)


# Attach the score to the diagnosis-level data and save ------------------------

dia_eci <- data_long %>%
  left_join(
    patient_year_scores,
    by = c("patient_id", "id_jahr")
  ) %>%
  mutate(eci = coalesce(eci, 0L))

save(dia_eci, file = output_file)

