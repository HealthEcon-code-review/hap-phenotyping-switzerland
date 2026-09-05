# Phenotyping Hospital-Acquired Pneumonia in Administrative Data
# Phenotyping Hospital-Acquired Pneumonia in Administrative Data
# 02 - Cohort construction and data cleaning
#
# Expected inputs:
#   PNEUMONIA_DATA_DIR/dia_eci.RData containing an object named `dia_eci`
#   PNEUMONIA_RAW_DATA_DIR/min.RData containing an object named `min`
#
# Created outputs:
#   PNEUMONIA_DATA_DIR/df.RData containing `df`
#   PNEUMONIA_DATA_DIR/df_CV.RData containing `df_CV`
#   PNEUMONIA_DATA_DIR/cohort_checks.RData containing dynamically calculated
#   cohort denominators and annual ICD-10 HAP frequencies
#
# Set the directories locally, for example in an untracked .Renviron file:
#   PNEUMONIA_DATA_DIR=/path/to/derived/data
#   PNEUMONIA_RAW_DATA_DIR=/path/to/authorised/raw/data
#
# If PNEUMONIA_RAW_DATA_DIR is not set, PNEUMONIA_DATA_DIR is used for both.

suppressPackageStartupMessages({
  library(tidyverse)
})


# Analysis settings ------------------------------------------------------------

last_study_year <- 2019L


# File paths -------------------------------------------------------------------

data_dir <- Sys.getenv("PNEUMONIA_DATA_DIR", unset = "")
raw_data_dir <- Sys.getenv("PNEUMONIA_RAW_DATA_DIR", unset = data_dir)

if (!nzchar(data_dir)) {
  stop(
    "PNEUMONIA_DATA_DIR is not set. Point it to the directory containing ",
    "dia_eci.RData and used for derived output files."
  )
}

if (!nzchar(raw_data_dir)) {
  stop(
    "PNEUMONIA_RAW_DATA_DIR is not set. Point it to the authorised raw-data ",
    "directory containing min.RData."
  )
}

data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
raw_data_dir <- normalizePath(raw_data_dir, winslash = "/", mustWork = TRUE)

dia_eci_file <- file.path(data_dir, "dia_eci.RData")
minimum_data_file <- file.path(raw_data_dir, "min.RData")

df_file <- file.path(data_dir, "df.RData")
df_cv_file <- file.path(data_dir, "df_CV.RData")
cohort_checks_file <- file.path(data_dir, "cohort_checks.RData")


# Helper functions -------------------------------------------------------------

load_named_object <- function(file, object_name) {
  if (!file.exists(file)) {
    stop("Input file not found: ", file)
  }
  
  input_environment <- new.env(parent = emptyenv())
  loaded_objects <- load(file, envir = input_environment)
  
  if (!object_name %in% loaded_objects) {
    stop(
      "The file ", basename(file),
      " must contain an object named `", object_name, "`."
    )
  }
  
  input_environment[[object_name]]
}

require_columns <- function(x, required_columns, object_name) {
  missing_columns <- setdiff(required_columns, names(x))
  
  if (length(missing_columns)) {
    stop(
      "Missing required columns in `", object_name, "`: ",
      paste(missing_columns, collapse = ", ")
    )
  }
}


# Load diagnosis data with Elixhauser scores ----------------------------------

dia_eci <- load_named_object(dia_eci_file, "dia_eci")

require_columns(
  dia_eci,
  c("fall_id", "id_jahr", "patient_id", "dia", "rank", "eci"),
  "dia_eci"
)

dia_eci <- dia_eci %>%
  rename(
    case_id = fall_id,
    year = id_jahr
  ) %>%
  mutate(
    case_id = as.character(case_id),
    patient_id = as.character(patient_id),
    year = as.integer(year),
    dia = as.character(dia),
    rank = as.character(rank),
    eci = as.integer(eci)
  ) %>%
  filter(year <= last_study_year)


# Load and prepare the minimum dataset ----------------------------------------
minimum_data <- load_named_object(minimum_data_file, "min")

minimum_columns <- c(
  "fall_id",
  "id_jahr",
  "hospitaltyp",
  "LOS_sdrg",
  "tage_bis_hospitalisierung",
  "sequenz_num",
  "bur",
  "kanton_spital",
  "aufenthaltineinerIntensivstation",
  "patient_id",
  "Geschlecht",
  "alter_5",
  "monat_eintritt",
  "AufenthaltsortvordemEintritt",
  "Eintrittsart",
  "entscheidfuerAustritt",
  "aufenthaltnach_austritt",
  "dauer_kuenstlichen_beatmung",
  "tage_bis_tod",
  "list_home",
  "list _other",
  "list_netzwerk"
)

require_columns(minimum_data, minimum_columns, "min")

min2 <- minimum_data %>%
  filter(id_jahr <= last_study_year) %>%
  select(all_of(minimum_columns)) %>%
  rename(
    case_id = fall_id,
    year = id_jahr,
    hospital_type = hospitaltyp,
    length_of_stay = LOS_sdrg,
    days_to_hospitalisation = tage_bis_hospitalisierung,
    sequence_number = sequenz_num,
    hospital_id = bur,
    hospital_canton = kanton_spital,
    icu_stay = aufenthaltineinerIntensivstation,
    sex = Geschlecht,
    age = alter_5,
    admission_month = monat_eintritt,
    pre_admission_location = AufenthaltsortvordemEintritt,
    admission_type = Eintrittsart,
    discharge_decision = entscheidfuerAustritt,
    post_discharge_destination = aufenthaltnach_austritt,
    ventilation_hours = dauer_kuenstlichen_beatmung,
    days_to_death = tage_bis_tod,
    hospital_list_home = list_home,
    hospital_list_other = `list _other`,
    hospital_list_network = list_netzwerk
  ) %>%
  mutate(
    case_id = as.character(case_id),
    patient_id = as.character(patient_id),
    year = as.integer(year)
  )


# Eligible-admission denominator ----------------------------------------------

eligible_by_year <- min2 %>%
  filter(
    !as.character(age) %in% c("0-4", "5-9", "10-14", "15-19"),
    as.character(admission_type) != "3"
  ) %>%
  distinct(case_id, year) %>%
  count(year, name = "n_all_eligible")

n_all_eligible <- sum(eligible_by_year$n_all_eligible)


# Define pneumonia cases -------------------------------------------------------

pneumonia_exact_codes <- c(
  "A481", "J13", "J14", "J851", "J849"
)

pneumonia_prefix_pattern <- "^(J12|J15|J16|J17|J18|U690)"

is_pneumonia_code <- function(x) {
  !is.na(x) & (
    x %in% pneumonia_exact_codes |
      grepl(pneumonia_prefix_pattern, x)
  )
}

pneumonia_keys <- dia_eci %>%
  filter(is_pneumonia_code(dia)) %>%
  distinct(case_id, year)


min_pneu <- min2 %>%
  semi_join(pneumonia_keys, by = c("case_id", "year"))

duplicate_minimum_keys <- min_pneu %>%
  count(case_id, year) %>%
  filter(n > 1L)

if (nrow(duplicate_minimum_keys)) {
  stop("The minimum dataset contains duplicate case_id-year records.")
}

dia_pneu <- dia_eci %>%
  semi_join(pneumonia_keys, by = c("case_id", "year")) %>%
  mutate(
    is_pneumonia = as.integer(is_pneumonia_code(dia)),
    
    pneumonia_main_diagnosis = as.integer(
      is_pneumonia == 1L & !is.na(rank) & rank == "H"
    ),
    
    pneumonia_secondary_diagnosis = as.integer(
      is_pneumonia == 1L & !is.na(rank) & rank != "H"
    ),
    
    hap = as.integer(
      !is.na(dia) & grepl("^U690", dia)
    )
  ) %>%
  select(
    -any_of(c("Seitigkeitdia", "Tumoraktivitaetdia", "rank", "sequenz_num"))
  ) %>%
  distinct()

dia_pneu_case <- dia_pneu %>%
  group_by(case_id, patient_id, year) %>%
  summarise(
    pneumonia_main_diagnosis = max(pneumonia_main_diagnosis),
    pneumonia_secondary_diagnosis = max(pneumonia_secondary_diagnosis),
    hap = max(hap),
    pneumonia_case = max(is_pneumonia),
    .groups = "drop"
  )

pneu2 <- min_pneu %>%
  left_join(
    dia_pneu_case,
    by = c("case_id", "patient_id", "year")
  )

if (nrow(pneu2) != nrow(min_pneu)) {
  stop("Joining diagnosis flags changed the number of pneumonia admissions.")
}

if (anyNA(pneu2$pneumonia_case)) {
  stop("At least one selected pneumonia admission has no diagnosis summary.")
}


# Attach one Elixhauser score per patient-year --------------------------------

eci_py <- dia_eci %>%
  distinct(patient_id, year, eci)

duplicate_eci <- eci_py %>%
  count(patient_id, year) %>%
  filter(n > 1L)

if (nrow(duplicate_eci)) {
  stop("More than one distinct ECI value exists for at least one patient-year.")
}

n_before_eci_join <- nrow(pneu2)

pneu2 <- pneu2 %>%
  left_join(eci_py, by = c("patient_id", "year"))

if (nrow(pneu2) != n_before_eci_join) {
  stop("Joining ECI values changed the number of pneumonia admissions.")
}

if (anyNA(pneu2$eci)) {
  stop("At least one pneumonia admission has no patient-year ECI value.")
}


# Construct the adult pneumonia cohort ---------------------------------------

df <- pneu2 %>%
  mutate(
    gender = factor(
      as.character(sex),
      levels = c("1", "2"),
      labels = c("Male", "Female")
    ),
    admission_type = forcats::fct_collapse(
      factor(as.character(admission_type)),
      Emergency = "1",
      Scheduled = "2",
      Birth = "3",
      Transfer = c("4", "5"),
      Other_Unknown = c("8", "9")
    ),
    pre_admission_code = as.character(pre_admission_location),
    place_of_stay = case_when(
      pre_admission_code %in% c("1", "2") ~ "Home",
      pre_admission_code %in% c("3", "4") ~ "Nursing/Residential",
      pre_admission_code %in% c("6", "66") ~ "Acute",
      pre_admission_code %in% c("5", "55", "83", "84") ~ "Psych/Rehab",
      pre_admission_code %in% c("7", "8", "9") ~ "Other/Unknown",
      TRUE ~ pre_admission_code
    ),
    place_of_stay = relevel(factor(place_of_stay), ref = "Home"),
    list = pmax(
      hospital_list_home,
      hospital_list_other,
      hospital_list_network,
      na.rm = TRUE
    )
  ) %>%
  select(-pre_admission_code) %>%
  filter(
    !as.character(age) %in% c("0-4", "5-9", "10-14", "15-19"),
    admission_type != "Birth"
  ) %>%
  mutate(
    admission_type = forcats::fct_drop(admission_type)
  ) %>%
  droplevels()

if (anyDuplicated(df[c("case_id", "year")])) {
  stop("The final pneumonia cohort contains duplicate case_id-year records.")
}

save(df, file = df_file)


# Add annual hospital pneumonia volume ----------------------------------------

hospital_volume_df <- df %>%
  count(hospital_id, year, name = "volume") %>%
  group_by(year) %>%
  mutate(
    volume_quartile = ntile(volume, 4L),
    volume_group = factor(
      volume_quartile,
      levels = 1:4,
      labels = c("Low", "Medium-low", "Medium-high", "High")
    )
  ) %>%
  ungroup()

df_CV <- df %>%
  left_join(hospital_volume_df, by = c("hospital_id", "year")) %>%
  mutate(
    death_in_hospital = as.integer(
      as.character(discharge_decision) == "5"
    ),
    death_30 = as.integer(days_to_death <= 30),
    death_90 = as.integer(days_to_death <= 90),
    death_365 = as.integer(days_to_death <= 365),
    readm_30 = as.integer(days_to_hospitalisation <= 30),
    readm_90 = as.integer(days_to_hospitalisation <= 90),
    readm_365 = as.integer(days_to_hospitalisation <= 365)
  ) %>%
  select(
    -volume_quartile,
    -any_of(c(
      "hospital_list_home",
      "hospital_list_other",
      "hospital_list_network",
      "sequence_number"
    ))
  )

if (nrow(df_CV) != nrow(df)) {
  stop("Joining hospital volumes changed the number of pneumonia admissions.")
}

save(df_CV, file = df_cv_file)


# Reproducible cohort checks ---------------------------------------------------

pneumonia_icd_hap_by_year <- df %>%
  group_by(year) %>%
  summarise(
    n_pneumonia = n(),
    n_u69 = sum(hap == 1L, na.rm = TRUE),
    .groups = "drop"
  )

u69_all_admissions_by_year <- eligible_by_year %>%
  left_join(pneumonia_icd_hap_by_year, by = "year") %>%
  mutate(
    across(
      c(n_pneumonia, n_u69),
      ~ dplyr::coalesce(.x, 0L)
    ),
    pneumonia_per_1000 = round(
      1000 * n_pneumonia / n_all_eligible, 2
    ),
    u69_per_1000 = round(
      1000 * n_u69 / n_all_eligible, 2
    ),
    pct_u69_among_pneumonia = round(
      100 * n_u69 / n_pneumonia, 2
    )
  )

cohort_counts <- tibble(
  group = c(
    "All eligible admissions",
    "Pneumonia cohort",
    "ICD-10 HAP"
  ),
  n = c(
    n_all_eligible,
    nrow(df),
    sum(df$hap == 1L, na.rm = TRUE)
  )
) %>%
  mutate(percentage_of_all_eligible = 100 * n / n_all_eligible)

panel_a_check <- u69_all_admissions_by_year %>%
  select(
    year,
    eligible_admissions = n_all_eligible,
    pneumonia = n_pneumonia,
    pneumonia_per_1000,
    u69_pneumonia = n_u69,
    u69_per_1000,
    pct_u69_among_pneumonia
  )

panel_b_check <- u69_volume_check %>%
  select(
    year,
    hospitals_ge100,
    hospitals_ge100_with_u69,
    pct_ge100_with_u69
  ) %>%
  left_join(
    u69_rule_by_year %>%
      select(
        year,
        pct_rule_hap,
        pct_rule_added_without_u69
      ),
    by = "year"
  )

panel_b_final <- panel_b_check %>%
  transmute(
    Year = as.character(year),
    `Hospitals with ≥100 pneumonia hospitalisations, n` =
      hospitals_ge100,
    `Hospitals recording U69.0*, n (%)` = sprintf(
      "%d (%.1f%%)",
      hospitals_ge100_with_u69,
      pct_ge100_with_u69
    ),
    `Rule-based HAP, %` = round(pct_rule_hap, 1),
    `Additional rule-based HAP without U69.0*, %` =
      round(pct_rule_added_without_u69, 1)
  )

save(
  n_all_eligible,
  cohort_counts,
  u69_all_admissions_by_year,
  file = cohort_checks_file
)


adult_pneumonia_before_birth_exclusion <- pneu2 %>%
  filter(
    !as.character(age) %in%
      c("0-4", "5-9", "10-14", "15-19")
  )

cohort_flow_counts <- tibble(
  stage = c(
    "Diagnosis records",
    "Patient-information records",
    "Pneumonia diagnosis records",
    "Unique pneumonia hospitalisations",
    "Hospitalisations linked to patient data",
    "Excluded: age <20 years",
    "Excluded: admission type birth",
    "Final analytical cohort"
  ),
  n = c(
    nrow(dia_eci),
    nrow(min2),
    sum(is_pneumonia_code(dia_eci$dia)),
    nrow(pneumonia_keys),
    nrow(pneu2),
    sum(
      as.character(pneu2$age) %in%
        c("0-4", "5-9", "10-14", "15-19")
    ),
    sum(
      as.character(
        adult_pneumonia_before_birth_exclusion$admission_type
      ) == "3",
      na.rm = TRUE
    ),
    nrow(df)
  )
)

cohort_flow_counts
