# Phenotyping Hospital-Acquired Pneumonia in Administrative Data
# 05 - Outcome analysis, model tables, and visualisations
#
# Expected input:
#   PNEUMONIA_DATA_DIR/df_classified.RData containing an object named df
#
# Created outputs:
#   PNEUMONIA_OUTPUT_DIR/tables/Appendix_model_tables.docx
#   PNEUMONIA_OUTPUT_DIR/figures/Figure1_adjusted_hazard_ratios.pdf
#   PNEUMONIA_OUTPUT_DIR/figures/Figure1_adjusted_hazard_ratios.tiff
#   PNEUMONIA_OUTPUT_DIR/figures/Figure2_CAP_HAP_hospital_mortality.pdf
#   PNEUMONIA_OUTPUT_DIR/figures/Figure2_CAP_HAP_hospital_mortality.tiff
#   PNEUMONIA_OUTPUT_DIR/outcome_analysis_results.RData
#
# Set the directories locally, for example in an untracked .Renviron file:
#   PNEUMONIA_DATA_DIR=/path/to/derived/data
#   PNEUMONIA_OUTPUT_DIR=/path/to/analysis/output
#   PNEUMONIA_PLOT_FONT=Times New Roman
#
# If PNEUMONIA_OUTPUT_DIR is not set, a sibling directory named 03_Output is
# used. The portable default plot font is serif.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(splines)
  library(survival)
  library(lme4)
  library(glmmTMB)
  library(psych)
  library(flextable)
  library(officer)
  library(ggplot2)
  library(scales)
})

# File paths -------------------------------------------------------------------

data_dir <- Sys.getenv("PNEUMONIA_DATA_DIR", unset = "")

if (!nzchar(data_dir)) {
  stop(
    "PNEUMONIA_DATA_DIR is not set. Point it to the directory containing ",
    "df_classified.RData."
  )
}

data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)

output_dir <- Sys.getenv(
  "PNEUMONIA_OUTPUT_DIR",
  unset = file.path(dirname(data_dir), "03_Output")
)

tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
tables_dir <- normalizePath(tables_dir, winslash = "/", mustWork = TRUE)
figures_dir <- normalizePath(figures_dir, winslash = "/", mustWork = TRUE)

input_file <- file.path(data_dir, "df_classified.RData")

model_tables_file <- file.path(
  tables_dir,
  "Appendix_model_tables.docx"
)

figure1_pdf <- file.path(
  figures_dir,
  "Figure1_adjusted_hazard_ratios.pdf"
)

figure1_tiff <- file.path(
  figures_dir,
  "Figure1_adjusted_hazard_ratios.tiff"
)

figure2_pdf <- file.path(
  figures_dir,
  "Figure2_CAP_HAP_hospital_mortality.pdf"
)

figure2_tiff <- file.path(
  figures_dir,
  "Figure2_CAP_HAP_hospital_mortality.tiff"
)

analysis_results_file <- file.path(
  output_dir,
  "outcome_analysis_results.RData"
)

plot_font_family <- Sys.getenv(
  "PNEUMONIA_PLOT_FONT",
  unset = "serif"
)


# Load and validate the classified cohort --------------------------------------

load_named_object <- function(file, object_name) {
  if (!file.exists(file)) {
    stop("Input file not found: ", file)
  }

  input_environment <- new.env(parent = emptyenv())
  loaded_objects <- load(file, envir = input_environment)

  if (!object_name %in% loaded_objects) {
    stop(
      "The file ", basename(file),
      " must contain an object named ", object_name, "."
    )
  }

  input_environment[[object_name]]
}

require_columns <- function(x, required_columns, object_name) {
  missing_columns <- setdiff(required_columns, names(x))

  if (length(missing_columns)) {
    stop(
      "Missing required columns in ", object_name, ": ",
      paste(missing_columns, collapse = ", ")
    )
  }
}

assert_binary <- function(x, variable_name) {
  observed_values <- sort(unique(stats::na.omit(as.integer(x))))

  if (anyNA(x) || !all(observed_values %in% c(0L, 1L))) {
    stop(
      variable_name,
      " must contain only 0 and 1 without missing values."
    )
  }
}

df <- load_named_object(input_file, "df")

required_columns <- c(
  "case_id",
  "year",
  "hospital_id",
  "hospital_type",
  "hospital_canton",
  "volume_group",
  "list",
  "age",
  "gender",
  "eci",
  "admission_month",
  "death_in_hospital",
  "days_to_death",
  "days_to_hospitalisation",
  "hap",
  "HAP_rule",
  "HAP_PU80"
)

require_columns(df, required_columns, "df")

if (anyDuplicated(df[c("case_id", "year")])) {
  stop("The classified cohort contains duplicate case_id-year records.")
}

assert_binary(df$hap, "hap")
assert_binary(df$HAP_rule, "HAP_rule")
assert_binary(df$HAP_PU80, "HAP_PU80")
assert_binary(df$death_in_hospital, "death_in_hospital")


# Prepare analysis variables ---------------------------------------------------

age_lookup <- tibble(
  age = c(
    "20-24", "25-29", "30-34", "35-39", "40-44", "45-49",
    "50-54", "55-59", "60-64", "65-69", "70-74", "75-79",
    "80-84", "85-89", "90-94", "95+"
  ),
  age_num = c(
    22, 27, 32, 37, 42, 47, 52, 57,
    62, 67, 72, 77, 82, 87, 92, 97
  )
)

observed_age_bands <- unique(stats::na.omit(as.character(df$age)))
unmapped_age_bands <- setdiff(observed_age_bands, age_lookup$age)

if (length(unmapped_age_bands)) {
  stop(
    "Unmapped age bands: ",
    paste(sort(unmapped_age_bands), collapse = ", ")
  )
}

# Negative follow-up intervals are impossible for post-discharge outcomes.
# Retain their counts for reproducibility, then treat the invalid intervals as
# missing. Admissions themselves remain in the cohort.
follow_up_data_checks <- tibble(
  n_negative_death = sum(
    df$days_to_death < 0,
    na.rm = TRUE
  ),
  n_negative_hospitalisation = sum(
    df$days_to_hospitalisation < 0,
    na.rm = TRUE
  ),
  n_negative_hospitalisation_after_alive_discharge = sum(
    df$days_to_hospitalisation < 0 &
      df$death_in_hospital == 0L,
    na.rm = TRUE
  )
)

df <- df %>%
  select(-any_of("age_num")) %>%
  mutate(
    age = as.character(age),
    days_to_death = if_else(
      !is.na(days_to_death) &
        (days_to_death < 0 | days_to_death > 365),
      NA_real_,
      as.numeric(days_to_death)
    ),
    days_to_hospitalisation = if_else(
      !is.na(days_to_hospitalisation) &
        (
          days_to_hospitalisation < 0 |
            days_to_hospitalisation > 365
        ),
      NA_real_,
      as.numeric(days_to_hospitalisation)
    ),
    admission_month = factor(
      as.integer(as.character(admission_month)),
      levels = 1:12
    ),
    year = factor(as.character(year)),
    hospital_id = factor(as.character(hospital_id)),
    hospital_canton = factor(as.character(hospital_canton)),
    volume_group = factor(
      as.character(volume_group),
      levels = c("Low", "Medium-low", "Medium-high", "High")
    ),
    gender = relevel(
      factor(as.character(gender)),
      ref = "Male"
    ),
    list = as.numeric(as.character(list)),
    hap = as.integer(hap),
    HAP_rule = as.integer(HAP_rule),
    HAP_PU80 = as.integer(HAP_PU80),
    death_in_hospital = as.integer(death_in_hospital)
  ) %>%
  left_join(age_lookup, by = "age")

if ("AG" %in% levels(df$hospital_canton)) {
  df$hospital_canton <- relevel(df$hospital_canton, ref = "AG")
}

negative_follow_up <- c(
  df$days_to_death[!is.na(df$days_to_death)],
  df$days_to_hospitalisation[!is.na(df$days_to_hospitalisation)]
)

if (any(negative_follow_up < 0)) {
  stop("Follow-up times must not be negative.")
}

hospital_group_levels <- c(
  "University hospital",
  "Central care",
  "Maximum care",
  "Supra-regional basic care",
  "Regional basic care",
  "Psychiatry",
  "Rehabilitation",
  "Specialized clinic"
)

df <- df %>%
  mutate(
    hospital_group = case_when(
      hospital_type == "K111" ~ "University hospital",
      hospital_type == "K112" ~ "Central care",
      hospital_type == "K121" ~ "Maximum care",
      hospital_type == "K122" ~ "Supra-regional basic care",
      hospital_type == "K123" ~ "Regional basic care",
      hospital_type %in% c("K211", "K212") ~ "Psychiatry",
      hospital_type == "K221" ~ "Rehabilitation",
      hospital_type %in%
        c("K231", "K232", "K233", "K234", "K235") ~
        "Specialized clinic",
      TRUE ~ NA_character_
    ),
    hospital_group = factor(
      hospital_group,
      levels = hospital_group_levels
    )
  )

unmapped_hospital_types <- sort(unique(
  as.character(df$hospital_type[is.na(df$hospital_group)])
))

if (length(unmapped_hospital_types)) {
  stop(
    "Unmapped or missing hospital types: ",
    paste(unmapped_hospital_types, collapse = ", ")
  )
}

complete_covariates <- c(
  "hospital_id",
  "hospital_canton",
  "volume_group",
  "hospital_group",
  "age_num",
  "gender",
  "eci",
  "admission_month",
  "year"
)

# Missing hospital-list status is retained in the cohort and handled by
# complete-case analysis only in models that include this covariate.
model_covariate_missingness <- tibble(
  variable = c(complete_covariates, "list"),
  n_missing = vapply(
    df[c(complete_covariates, "list")],
    function(x) sum(is.na(x)),
    integer(1)
  )
) %>%
  mutate(
    percentage_missing = 100 * n_missing / nrow(df)
  )

missing_covariates <- vapply(
  df[complete_covariates],
  anyNA,
  logical(1)
)

if (any(missing_covariates)) {
  stop(
    "Missing values in required model covariates: ",
    paste(names(missing_covariates)[missing_covariates], collapse = ", ")
  )
}

if (!all(is.finite(df$age_num)) ||
    !all(is.finite(df$eci)) ||
    !all(is.finite(df$list[!is.na(df$list)]))) {
  stop("Non-missing age_num, eci, and list values must be finite.")
}



# ICD-10 HAP coding checks -----------------------------------------------------

u69_hospital_year <- df %>%
  group_by(year, hospital_id) %>%
  summarise(
    n_pneumonia = n(),
    n_u69 = sum(hap == 1L),
    pct_u69 = 100 * n_u69 / n_pneumonia,
    .groups = "drop"
  )

u69_hospital_use <- u69_hospital_year %>%
  group_by(year) %>%
  summarise(
    n_hospitals = n(),
    n_hospitals_with_u69 = sum(n_u69 > 0L),
    pct_hospitals_with_u69 =
      100 * n_hospitals_with_u69 / n_hospitals,
    .groups = "drop"
  )

u69_volume_check <- u69_hospital_year %>%
  group_by(year) %>%
  summarise(
    hospitals_ge50 = sum(n_pneumonia >= 50L),
    hospitals_ge50_with_u69 =
      sum(n_pneumonia >= 50L & n_u69 > 0L),
    pct_ge50_with_u69 = if_else(
      hospitals_ge50 > 0L,
      100 * hospitals_ge50_with_u69 / hospitals_ge50,
      NA_real_
    ),
    hospitals_ge100 = sum(n_pneumonia >= 100L),
    hospitals_ge100_with_u69 =
      sum(n_pneumonia >= 100L & n_u69 > 0L),
    pct_ge100_with_u69 = if_else(
      hospitals_ge100 > 0L,
      100 * hospitals_ge100_with_u69 / hospitals_ge100,
      NA_real_
    ),
    median_u69_pct_among_coders =
      median(pct_u69[n_u69 > 0L], na.rm = TRUE),
    .groups = "drop"
  )

n_years <- n_distinct(df$year)

stable_large_ids <- u69_hospital_year %>%
  group_by(hospital_id) %>%
  summarise(
    years_observed = n_distinct(year),
    minimum_volume = min(n_pneumonia),
    .groups = "drop"
  ) %>%
  filter(
    years_observed == n_years,
    minimum_volume >= 100L
  ) %>%
  pull(hospital_id)

u69_stable_large <- u69_hospital_year %>%
  filter(hospital_id %in% stable_large_ids) %>%
  group_by(year) %>%
  summarise(
    n_hospitals = n(),
    n_pneumonia = sum(n_pneumonia),
    n_u69 = sum(n_u69),
    pct_u69_weighted = 100 * n_u69 / n_pneumonia,
    median_hospital_pct = median(pct_u69),
    q1_hospital_pct = unname(quantile(pct_u69, 0.25)),
    q3_hospital_pct = unname(quantile(pct_u69, 0.75)),
    .groups = "drop"
  )

u69_rule_by_year <- df %>%
  group_by(year) %>%
  summarise(
    n_pneumonia = n(),
    pct_icd_hap = 100 * mean(hap == 1L),
    pct_rule_hap = 100 * mean(HAP_rule == 1L),
    pct_rule_added_without_u69 =
      100 * mean(HAP_rule == 1L & hap == 0L),
    .groups = "drop"
  )


# Spline basis and model formulas ---------------------------------------------

add_spline_basis <- function(data) {
  age_basis <- splines::ns(data$age_num, df = 5)
  eci_basis <- splines::ns(data$eci, df = 4)

  data %>%
    mutate(
      age_ns1 = age_basis[, 1],
      age_ns2 = age_basis[, 2],
      age_ns3 = age_basis[, 3],
      age_ns4 = age_basis[, 4],
      age_ns5 = age_basis[, 5],
      eci_ns1 = eci_basis[, 1],
      eci_ns2 = eci_basis[, 2],
      eci_ns3 = eci_basis[, 3],
      eci_ns4 = eci_basis[, 4]
    )
}

vars_patient <- c(
  "gender",
  "year",
  "admission_month"
)

vars_hospital <- c(
  "volume_group",
  "hospital_canton",
  "list"
)

vars_hospital_type <- c(
  "hospital_group",
  "hospital_canton",
  "list"
)

make_piecewise_formula <- function(
    response,
    hap_variable,
    adjustment_variables) {
  time_varying_terms <- c(
    paste0(hap_variable, ":strata(interval)"),
    paste0("age_ns", 1:5, ":strata(interval)"),
    paste0("eci_ns", 1:4, ":strata(interval)")
  )

  as.formula(
    paste(
      response,
      "~",
      paste(
        c(
          "strata(interval)",
          adjustment_variables,
          time_varying_terms
        ),
        collapse = " + "
      )
    )
  )
}


# Post-discharge mortality -----------------------------------------------------

df_surv_death <- df %>%
  filter(death_in_hospital == 0L) %>%
  mutate(
    time_death = if_else(
      is.na(days_to_death),
      365,
      pmin(days_to_death, 365)
    ),
    event_death = as.integer(
      !is.na(days_to_death) & days_to_death <= 365
    ),
    time_death = pmax(time_death, 1e-3)
  ) %>%
  add_spline_basis()

df_pw_death <- survival::survSplit(
  data = df_surv_death,
  cut = c(30, 90),
  end = "time_death",
  event = "event_death",
  start = "tstart",
  episode = "interval"
) %>%
  mutate(
    interval = factor(
      interval,
      levels = 1:3,
      labels = c("0-30", "31-90", "91-365")
    )
  )

f_ttd_m1 <- make_piecewise_formula(
  "Surv(tstart, time_death, event_death)",
  "HAP_rule",
  vars_patient
)

f_ttd_m2 <- make_piecewise_formula(
  "Surv(tstart, time_death, event_death)",
  "HAP_rule",
  c(vars_patient, vars_hospital)
)

ttd_m1 <- survival::coxph(
  f_ttd_m1,
  data = df_pw_death,
  cluster = hospital_id
)

ttd_m2 <- survival::coxph(
  f_ttd_m2,
  data = df_pw_death,
  cluster = hospital_id
)


# Post-discharge readmission ---------------------------------------------------

df_surv_readm <- df %>%
  filter(death_in_hospital == 0L) %>%
  mutate(
    t_readm = if_else(
      is.na(days_to_hospitalisation),
      Inf,
      days_to_hospitalisation
    ),
    t_death = if_else(
      is.na(days_to_death),
      Inf,
      days_to_death
    ),
    time_readm = pmin(t_readm, t_death, 365),
    event_readm = as.integer(
      t_readm <= t_death & t_readm <= 365
    ),
    time_readm = pmax(time_readm, 1e-3)
  ) %>%
  add_spline_basis()

df_pw_readm <- survival::survSplit(
  data = df_surv_readm,
  cut = c(30, 90),
  end = "time_readm",
  event = "event_readm",
  start = "tstart",
  episode = "interval"
) %>%
  mutate(
    interval = factor(
      interval,
      levels = 1:3,
      labels = c("0-30", "31-90", "91-365")
    )
  )

f_ttr_m1 <- make_piecewise_formula(
  "Surv(tstart, time_readm, event_readm)",
  "HAP_rule",
  vars_patient
)

f_ttr_m2 <- make_piecewise_formula(
  "Surv(tstart, time_readm, event_readm)",
  "HAP_rule",
  c(vars_patient, vars_hospital)
)

ttr_m1 <- survival::coxph(
  f_ttr_m1,
  data = df_pw_readm,
  cluster = hospital_id
)

ttr_m2 <- survival::coxph(
  f_ttr_m2,
  data = df_pw_readm,
  cluster = hospital_id
)


# In-hospital mortality --------------------------------------------------------

glm_ihm <- glm(
  death_in_hospital ~
    HAP_rule +
    gender +
    ns(age_num, df = 5) +
    ns(eci, df = 4) +
    year +
    admission_month,
  data = df,
  family = binomial()
)


# Hospital variation and profiling --------------------------------------------

calc_icc_mor <- function(fit) {
  hospital_sd <- attr(
    lme4::VarCorr(fit)$cond$hospital_id,
    "stddev"
  )

  hospital_variance <- as.numeric(hospital_sd)^2
  icc <- hospital_variance / (hospital_variance + pi^2 / 3)
  mor <- exp(0.6745 * sqrt(2 * hospital_variance))

  tibble(
    hospital_variance = hospital_variance,
    ICC = icc,
    MOR = mor
  )
}

df_var <- df %>%
  filter(death_in_hospital == 0L) %>%
  mutate(
    death_30 = as.integer(
      !is.na(days_to_death) & days_to_death <= 30
    ),
    readm_30 = as.integer(
      !is.na(days_to_hospitalisation) &
        days_to_hospitalisation <= 30 &
        (
          is.na(days_to_death) |
            days_to_hospitalisation <= days_to_death
        )
    )
  ) %>%
  add_spline_basis()

fit_var_mort_30 <- glmmTMB::glmmTMB(
  death_30 ~
    HAP_rule +
    gender +
    age_ns1 + age_ns2 + age_ns3 + age_ns4 + age_ns5 +
    eci_ns1 + eci_ns2 + eci_ns3 + eci_ns4 +
    year +
    admission_month +
    volume_group +
    hospital_canton +
    list +
    (1 | hospital_id),
  data = df_var,
  family = binomial()
)

fit_var_readm_30 <- glmmTMB::glmmTMB(
  readm_30 ~
    HAP_rule +
    gender +
    age_ns1 + age_ns2 + age_ns3 + age_ns4 + age_ns5 +
    eci_ns1 + eci_ns2 + eci_ns3 + eci_ns4 +
    year +
    admission_month +
    volume_group +
    hospital_canton +
    list +
    (1 | hospital_id),
  data = df_var,
  family = binomial()
)

icc_mor_table <- bind_rows(
  calc_icc_mor(fit_var_mort_30) %>%
    mutate(outcome = "30-day mortality"),
  calc_icc_mor(fit_var_readm_30) %>%
    mutate(outcome = "30-day readmission")
) %>%
  select(outcome, hospital_variance, ICC, MOR)

run_hospital_profiling_30d <- function(data, hap_value) {
  data_30 <- data %>%
    filter(HAP_rule == hap_value) %>%
    droplevels()

  if (!nrow(data_30)) {
    stop("No observations available for HAP_rule = ", hap_value, ".")
  }

  fit <- glmmTMB::glmmTMB(
    death_30 ~
      gender +
      age_ns1 + age_ns2 + age_ns3 + age_ns4 + age_ns5 +
      eci_ns1 + eci_ns2 + eci_ns3 + eci_ns4 +
      year +
      admission_month +
      (1 | hospital_id),
    data = data_30,
    family = binomial()
  )

  random_effects <- lme4::ranef(
    fit,
    condVar = TRUE
  )$cond$hospital_id

  post_mean <- random_effects[, "(Intercept)"]
  post_variance <- attr(random_effects, "condVar")
  post_se <- sqrt(as.numeric(post_variance[1, 1, ]))

  hospital_effects <- tibble(
    hospital_id = rownames(random_effects),
    post_mean = as.numeric(post_mean),
    post_se = post_se,
    q2.5 = post_mean - 1.96 * post_se,
    q97.5 = post_mean + 1.96 * post_se
  ) %>%
    mutate(
      class = case_when(
        q97.5 < 0 ~ "lower-mortality",
        q2.5 > 0 ~ "higher-mortality",
        TRUE ~ "average-mortality"
      ),
      class = factor(
        class,
        levels = c(
          "lower-mortality",
          "average-mortality",
          "higher-mortality"
        )
      ),
      or = exp(post_mean),
      or_low = exp(q2.5),
      or_high = exp(q97.5)
    )

  # Apply every hospital random intercept to the same national case mix within
  # the respective CAP or HAP group. This produces shrinkage-adjusted,
  # risk-standardised mortality rather than an average over each hospital's own
  # patients.
  fixed_linear_predictor <- predict(
    fit,
    newdata = data_30,
    type = "link",
    re.form = NA
  )

  standardised_risk <- vapply(
    hospital_effects$post_mean,
    function(hospital_intercept) {
      mean(plogis(fixed_linear_predictor + hospital_intercept))
    },
    numeric(1)
  )

  hospital_counts <- data_30 %>%
    mutate(hospital_id = as.character(hospital_id)) %>%
    count(hospital_id, name = "n_admissions")

  hospital_predictions <- hospital_effects %>%
    transmute(
      hospital_id,
      risk = standardised_risk,
      risk_pct = 100 * standardised_risk
    ) %>%
    left_join(hospital_counts, by = "hospital_id")

  list(
    hospital_effects = hospital_effects,
    hospital_predictions = hospital_predictions,
    hospital_profile = left_join(
      hospital_effects,
      hospital_predictions,
      by = "hospital_id"
    )
  )
}

profile_30d_cap <- run_hospital_profiling_30d(
  df_var,
  hap_value = 0L
)

profile_30d_hap <- run_hospital_profiling_30d(
  df_var,
  hap_value = 1L
)

rm(df_var)


# Sensitivity analyses ---------------------------------------------------------

f_ttd_pu_m2 <- make_piecewise_formula(
  "Surv(tstart, time_death, event_death)",
  "HAP_PU80",
  c(vars_patient, vars_hospital)
)

f_ttr_pu_m2 <- make_piecewise_formula(
  "Surv(tstart, time_readm, event_readm)",
  "HAP_PU80",
  c(vars_patient, vars_hospital)
)

ttd_m2_80 <- survival::coxph(
  f_ttd_pu_m2,
  data = df_pw_death,
  cluster = hospital_id
)

ttr_m2_80 <- survival::coxph(
  f_ttr_pu_m2,
  data = df_pw_readm,
  cluster = hospital_id
)

f_ttd_m2_type <- make_piecewise_formula(
  "Surv(tstart, time_death, event_death)",
  "HAP_rule",
  c(vars_patient, vars_hospital_type)
)

f_ttr_m2_type <- make_piecewise_formula(
  "Surv(tstart, time_readm, event_readm)",
  "HAP_rule",
  c(vars_patient, vars_hospital_type)
)

ttd_m2_type <- survival::coxph(
  f_ttd_m2_type,
  data = df_pw_death,
  cluster = hospital_id
)

ttr_m2_type <- survival::coxph(
  f_ttr_m2_type,
  data = df_pw_readm,
  cluster = hospital_id
)

ttd_frailty <- survival::coxph(
  update(
    f_ttd_m1,
    . ~ . + frailty(hospital_id, distribution = "gamma")
  ),
  data = df_pw_death
)

ttr_frailty <- survival::coxph(
  update(
    f_ttr_m1,
    . ~ . + frailty(hospital_id, distribution = "gamma")
  ),
  data = df_pw_readm
)

df_finegray <- df %>%
  filter(death_in_hospital == 0L) %>%
  mutate(
    t_readm = if_else(
      is.na(days_to_hospitalisation),
      Inf,
      days_to_hospitalisation
    ),
    t_death = if_else(
      is.na(days_to_death),
      Inf,
      days_to_death
    ),
    time = pmax(pmin(t_readm, t_death, 365), 1e-3),
    status = case_when(
      t_readm <= t_death & t_readm <= 365 ~ 1L,
      t_death < t_readm & t_death <= 365 ~ 2L,
      TRUE ~ 0L
    )
  ) %>%
  add_spline_basis()

df_finegray_long <- survival::finegray(
  Surv(time, status, type = "mstate") ~
    HAP_rule +
    gender +
    age_ns1 + age_ns2 + age_ns3 + age_ns4 + age_ns5 +
    eci_ns1 + eci_ns2 + eci_ns3 + eci_ns4 +
    year +
    admission_month +
    volume_group +
    hospital_canton +
    list +
    hospital_id,
  data = df_finegray,
  etype = 1
)

ttr_finegray <- survival::coxph(
  Surv(fgstart, fgstop, fgstatus) ~
    HAP_rule +
    gender +
    age_ns1 + age_ns2 + age_ns3 + age_ns4 + age_ns5 +
    eci_ns1 + eci_ns2 + eci_ns3 + eci_ns4 +
    year +
    admission_month +
    volume_group +
    hospital_canton +
    list,
  data = df_finegray_long,
  weights = fgwt,
  cluster = hospital_id
)

rm(df_finegray, df_finegray_long)


# Model diagnostics ------------------------------------------------------------

extract_concordance <- function(fit, model_name) {
  concordance_values <- summary(fit)$concordance

  tibble(
    model = model_name,
    concordance = unname(concordance_values[1]),
    standard_error = unname(concordance_values[2])
  )
}

concordance_table <- bind_rows(
  extract_concordance(ttd_m1, "Mortality model 1"),
  extract_concordance(ttd_m2, "Mortality model 2"),
  extract_concordance(ttr_m1, "Readmission model 1"),
  extract_concordance(ttr_m2, "Readmission model 2")
)

proportional_hazards_tests <- list(
  mortality_model_2 = survival::cox.zph(ttd_m2),
  readmission_model_2 = survival::cox.zph(ttr_m2)
)

rm(df_surv_death, df_surv_readm, df_pw_death, df_pw_readm)


# Hospital benchmarking summaries ---------------------------------------------

cap_predictions <- profile_30d_cap$hospital_predictions %>%
  transmute(
    hospital_id,
    cap_risk_pct = risk_pct,
    cap_n = n_admissions
  )

hap_predictions <- profile_30d_hap$hospital_predictions %>%
  transmute(
    hospital_id,
    hap_risk_pct = risk_pct,
    hap_n = n_admissions
  )

plot_df <- cap_predictions %>%
  inner_join(hap_predictions, by = "hospital_id") %>%
  mutate(total_n = cap_n + hap_n)

cap_profile <- profile_30d_cap$hospital_profile %>%
  select(
    hospital_id,
    cap_risk = risk_pct,
    cap_class = class
  )

hap_profile <- profile_30d_hap$hospital_profile %>%
  select(
    hospital_id,
    hap_risk = risk_pct,
    hap_class = class
  )

rank_df <- inner_join(
  cap_profile,
  hap_profile,
  by = "hospital_id"
)

if (nrow(rank_df) < 4L) {
  stop("At least four hospitals with both CAP and HAP are required.")
}

spearman_test <- cor.test(
  rank_df$cap_risk,
  rank_df$hap_risk,
  method = "spearman",
  exact = FALSE
)

rank_df <- rank_df %>%
  mutate(
    cap_quartile = ntile(cap_risk, 4L),
    hap_quartile = ntile(hap_risk, 4L),
    cap_best_quartile = cap_quartile == 1L,
    hap_best_quartile = hap_quartile == 1L,
    cap_worst_quartile = cap_quartile == 4L,
    hap_worst_quartile = hap_quartile == 4L,
    switcher = as.integer(cap_class != hap_class)
  )

quartile_agreement_table <- table(
  CAP_quartile = rank_df$cap_quartile,
  HAP_quartile = rank_df$hap_quartile
)

quartile_kappa <- psych::cohen.kappa(
  cbind(
    rank_df$cap_quartile,
    rank_df$hap_quartile
  )
)

kappa_confidence <- as.data.frame(quartile_kappa$confid)
weighted_kappa_row <- grep(
  "^weighted kappa$",
  rownames(kappa_confidence),
  ignore.case = TRUE
)

if (length(weighted_kappa_row) != 1L) {
  stop("Could not extract the weighted kappa confidence interval.")
}

weighted_kappa_values <- as.numeric(
  kappa_confidence[weighted_kappa_row, c("lower", "estimate", "upper")]
)

performance_class_table <- table(
  CAP_class = rank_df$cap_class,
  HAP_class = rank_df$hap_class
)

n_same_quartile <- sum(
  rank_df$cap_quartile == rank_df$hap_quartile
)

n_same_class <- sum(
  rank_df$cap_class == rank_df$hap_class
)

n_cap_best <- sum(rank_df$cap_best_quartile)
n_best_overlap <- sum(
  rank_df$cap_best_quartile & rank_df$hap_best_quartile
)

n_cap_worst <- sum(rank_df$cap_worst_quartile)
n_worst_overlap <- sum(
  rank_df$cap_worst_quartile & rank_df$hap_worst_quartile
)

profiling_summary <- tibble(
  metric = c(
    "Hospitals represented in both groups",
    "Spearman rank correlation",
    "Same mortality quartile",
    "Weighted kappa",
    "Same mortality-performance class",
    "CAP best-quartile hospitals also best for HAP",
    "CAP worst-quartile hospitals also worst for HAP"
  ),
  estimate = c(
    nrow(rank_df),
    unname(spearman_test$estimate),
    n_same_quartile / nrow(rank_df),
    weighted_kappa_values[2],
    n_same_class / nrow(rank_df),
    n_best_overlap / n_cap_best,
    n_worst_overlap / n_cap_worst
  ),
  lower = c(
    NA_real_,
    NA_real_,
    NA_real_,
    weighted_kappa_values[1],
    NA_real_,
    NA_real_,
    NA_real_
  ),
  upper = c(
    NA_real_,
    NA_real_,
    NA_real_,
    weighted_kappa_values[3],
    NA_real_,
    NA_real_,
    NA_real_
  ),
  numerator = c(
    nrow(rank_df),
    NA_real_,
    n_same_quartile,
    NA_real_,
    n_same_class,
    n_best_overlap,
    n_worst_overlap
  ),
  denominator = c(
    nrow(rank_df),
    NA_real_,
    nrow(rank_df),
    NA_real_,
    nrow(rank_df),
    n_cap_best,
    n_cap_worst
  )
)

hospital_characteristics <- df %>%
  mutate(hospital_id = as.character(hospital_id)) %>%
  group_by(hospital_id) %>%
  summarise(
    n_cases = n(),
    n_cap = sum(HAP_rule == 0L),
    n_hap = sum(HAP_rule == 1L),
    hap_share = mean(HAP_rule),
    mean_age = mean(age_num),
    mean_eci = mean(eci),
    volume_group = first(volume_group),
    list = first(list),
    .groups = "drop"
  )

switch_df <- rank_df %>%
  left_join(hospital_characteristics, by = "hospital_id")

switcher_characteristics <- switch_df %>%
  group_by(switcher) %>%
  summarise(
    n_hospitals = n(),
    mean_n_cases = mean(n_cases),
    mean_hap_share = mean(hap_share),
    median_hap_share = median(hap_share),
    mean_age = mean(mean_age),
    mean_eci = mean(mean_eci),
    .groups = "drop"
  )

switcher_volume_table <- table(
  Switcher = switch_df$switcher,
  Volume_group = switch_df$volume_group
)


# Model-table labels and formatting -------------------------------------------

month_labels <- c(
  "1" = "January",
  "2" = "February",
  "3" = "March",
  "4" = "April",
  "5" = "May",
  "6" = "June",
  "7" = "July",
  "8" = "August",
  "9" = "September",
  "10" = "October",
  "11" = "November",
  "12" = "December"
)

canton_labels <- c(
  "AG" = "Aargau",
  "AI" = "Appenzell Innerrhoden",
  "AR" = "Appenzell Ausserrhoden",
  "BE" = "Bern",
  "BL" = "Basel-Landschaft",
  "BS" = "Basel-Stadt",
  "FR" = "Fribourg",
  "GE" = "Geneva",
  "GL" = "Glarus",
  "GR" = "Graubünden",
  "JU" = "Jura",
  "LU" = "Lucerne",
  "NE" = "Neuchâtel",
  "NW" = "Nidwalden",
  "OW" = "Obwalden",
  "SG" = "St Gallen",
  "SH" = "Schaffhausen",
  "SO" = "Solothurn",
  "SZ" = "Schwyz",
  "TG" = "Thurgau",
  "TI" = "Ticino",
  "UR" = "Uri",
  "VD" = "Vaud",
  "VS" = "Valais",
  "ZG" = "Zug",
  "ZH" = "Zurich"
)

hospital_type_labels <- c(
  "centralcare" = "central care",
  "maximumcare" = "maximum care",
  "psychiatry" = "psychiatry",
  "regionalbasiccare" = "regional basic care",
  "rehabilitation" = "rehabilitation",
  "specializedclinic" = "specialized clinic",
  "supraregionalbasiccare" = "supra-regional basic care",
  "universityhospital" = "university hospital"
)

interval_labels <- c(
  "0-30" = "0–30 days",
  "31-90" = "31–90 days",
  "91-365" = "91–365 days"
)

extract_interval <- function(term) {
  interval <- stringr::str_extract(
    term,
    "0[-–]30|31[-–]90|91[-–]365"
  )

  if (is.na(interval)) {
    return(NA_character_)
  }

  interval <- stringr::str_replace_all(interval, "–", "-")
  unname(interval_labels[interval])
}

extract_spline_basis <- function(
    compact_term,
    variable = c("age", "eci")) {
  variable <- match.arg(variable)

  if (variable == "age") {
    patterns <- c(
      "(?:age_s|age_ns)([1-5])",
      "ns\\(age_num,df=5\\)([1-5])"
    )
  } else {
    patterns <- c(
      "(?:eci_s|eci_ns)([1-4])",
      "ns\\(eci,df=4\\)([1-4])"
    )
  }

  for (pattern in patterns) {
    matched_value <- stringr::str_match(
      compact_term,
      pattern
    )

    if (!is.na(matched_value[1, 2])) {
      return(matched_value[1, 2])
    }
  }

  NA_character_
}

label_model_term <- function(term) {
  clean_term <- stringr::str_replace_all(
    term,
    intToUtf8(96),
    ""
  )
  compact_term <- stringr::str_replace_all(clean_term, "\\s+", "")
  lower_compact <- stringr::str_to_lower(compact_term)
  interval <- extract_interval(clean_term)

  if (stringr::str_detect(lower_compact, "hap_rule")) {
    label <- "Rule-based HAP"
    if (!is.na(interval)) {
      label <- paste0(label, " (", interval, ")")
    }
    return(label)
  }

  if (stringr::str_detect(lower_compact, "hap_pu80")) {
    label <- "PU-learning HAP"
    if (!is.na(interval)) {
      label <- paste0(label, " (", interval, ")")
    }
    return(label)
  }

  age_basis <- extract_spline_basis(lower_compact, "age")
  if (!is.na(age_basis)) {
    label <- paste("Age spline basis", age_basis)
    if (!is.na(interval)) {
      label <- paste0(label, " (", interval, ")")
    }
    return(label)
  }

  eci_basis <- extract_spline_basis(lower_compact, "eci")
  if (!is.na(eci_basis)) {
    label <- paste("ECI spline basis", eci_basis)
    if (!is.na(interval)) {
      label <- paste0(label, " (", interval, ")")
    }
    return(label)
  }

  if (stringr::str_detect(
    lower_compact,
    "^gender(?:female|f|1)$"
  )) {
    return("Sex: female")
  }

  year_match <- stringr::str_match(
    lower_compact,
    "^year([12][0-9]{3})$"
  )

  if (!is.na(year_match[1, 2])) {
    return(paste0("Admission year: ", year_match[1, 2]))
  }

  if (stringr::str_detect(lower_compact, "^admission_month")) {
    month_value <- stringr::str_remove(
      lower_compact,
      "^admission_month"
    )

    if (month_value %in% names(month_labels)) {
      return(
        paste0(
          "Admission month: ",
          unname(month_labels[month_value])
        )
      )
    }
  }

  if (stringr::str_detect(lower_compact, "^volume_group")) {
    volume_value <- stringr::str_remove(
      lower_compact,
      "^volume_group"
    )
    volume_key <- stringr::str_replace_all(
      volume_value,
      "[^a-z]",
      ""
    )
    volume_map <- c(
      "low" = "low",
      "mediumlow" = "medium-low",
      "mediumhigh" = "medium-high",
      "high" = "high"
    )

    if (volume_key %in% names(volume_map)) {
      return(
        paste0(
          "Annual hospital volume: ",
          unname(volume_map[volume_key])
        )
      )
    }
  }

  if (stringr::str_detect(
    lower_compact,
    "^hospital_canton"
  )) {
    canton_code <- stringr::str_remove(
      compact_term,
      stringr::regex(
        "^hospital_canton",
        ignore_case = TRUE
      )
    )
    canton_code <- stringr::str_to_upper(canton_code)

    if (canton_code %in% names(canton_labels)) {
      return(
        paste0(
          "Hospital canton: ",
          canton_labels[canton_code],
          " (",
          canton_code,
          ")"
        )
      )
    }
  }

  if (stringr::str_detect(
    lower_compact,
    "^hospital_group"
  )) {
    hospital_type <- stringr::str_remove(
      lower_compact,
      "^hospital_group"
    )
    hospital_type <- stringr::str_replace_all(
      hospital_type,
      "[^a-z]",
      ""
    )

    if (hospital_type %in% names(hospital_type_labels)) {
      return(
        paste0(
          "Hospital type: ",
          unname(hospital_type_labels[hospital_type])
        )
      )
    }
  }

  if (stringr::str_detect(
    lower_compact,
    "^list(?:1|yes|true|listed|hospital_listed)?$"
  )) {
    return("Cantonal hospital-list status: listed")
  }

  NA_character_
}

extract_model_estimates <- function(fit, conf_level = 0.95) {
  beta <- stats::coef(fit)
  covariance <- stats::vcov(fit)

  if (is.null(names(beta)) ||
      is.null(rownames(covariance))) {
    stop("The fitted model does not provide named coefficients.")
  }

  coefficient_names <- intersect(
    names(beta),
    rownames(covariance)
  )
  beta <- beta[coefficient_names]
  covariance <- covariance[
    coefficient_names,
    coefficient_names,
    drop = FALSE
  ]
  standard_error <- sqrt(diag(covariance))

  keep <- is.finite(beta) &
    is.finite(standard_error) &
    standard_error > 0 &
    !stringr::str_detect(
      coefficient_names,
      "^\\(Intercept\\)$"
    ) &
    !stringr::str_detect(
      coefficient_names,
      "^frailty\\("
    )

  beta <- beta[keep]
  standard_error <- standard_error[keep]
  coefficient_names <- coefficient_names[keep]
  critical_value <- stats::qnorm(
    1 - (1 - conf_level) / 2
  )
  z_value <- beta / standard_error

  tibble(
    term = coefficient_names,
    estimate = exp(beta),
    conf_low = exp(
      beta - critical_value * standard_error
    ),
    conf_high = exp(
      beta + critical_value * standard_error
    ),
    p_value = 2 * stats::pnorm(-abs(z_value))
  )
}

use_decimal_point <- function(x) {
  stringr::str_replace_all(
    x,
    stringr::fixed("."),
    "."
  )
}

format_effect_ci <- function(
    estimate,
    conf_low,
    conf_high) {
  formatted <- sprintf(
    "%.2f (%.2f–%.2f)",
    estimate,
    conf_low,
    conf_high
  )
  use_decimal_point(formatted)
}

format_p <- function(p_value) {
  vapply(
    p_value,
    function(p) {
      if (is.na(p)) {
        return("")
      }
      if (p < 0.0001) {
        return("<0.0001")
      }

      formatted <- format(
        signif(p, digits = 2),
        scientific = FALSE,
        trim = TRUE
      )
      use_decimal_point(formatted)
    },
    character(1)
  )
}

make_model_table <- function(
    fit,
    effect_type = c(
      "OR",
      "HR",
      "subdistribution HR"
    ),
    strict_labels = TRUE) {
  effect_type <- match.arg(effect_type)
  estimates <- extract_model_estimates(fit)
  labels <- vapply(
    estimates$term,
    label_model_term,
    character(1)
  )

  if (strict_labels && anyNA(labels)) {
    unknown_terms <- estimates$term[is.na(labels)]
    stop(
      "No publication label is defined for: ",
      paste(unknown_terms, collapse = ", "),
      ". Add the term to label_model_term before exporting."
    )
  }

  labels[is.na(labels)] <-
    estimates$term[is.na(labels)]

  effect_header <- switch(
    effect_type,
    "OR" = "Adjusted OR (95% CI)",
    "HR" = "Adjusted HR (95% CI)",
    "subdistribution HR" =
      "Adjusted subdistribution HR (95% CI)"
  )

  output <- tibble(
    Variable = labels,
    effect = format_effect_ci(
      estimates$estimate,
      estimates$conf_low,
      estimates$conf_high
    ),
    p_value = format_p(estimates$p_value)
  )

  names(output) <- c(
    "Variable",
    effect_header,
    "p value"
  )
  output
}

as_flextable <- function(table_data) {
  effect_column <- names(table_data)[2]

  flextable::flextable(table_data) %>%
    flextable::theme_booktabs() %>%
    flextable::font(
      fontname = "Times New Roman",
      part = "all"
    ) %>%
    flextable::fontsize(
      size = 8.5,
      part = "all"
    ) %>%
    flextable::bold(part = "header") %>%
    flextable::align(
      j = "Variable",
      align = "left",
      part = "all"
    ) %>%
    flextable::align(
      j = c(effect_column, "p value"),
      align = "center",
      part = "all"
    ) %>%
    flextable::valign(
      valign = "center",
      part = "all"
    ) %>%
    flextable::padding(
      padding = 2.5,
      part = "all"
    ) %>%
    flextable::width(
      j = "Variable",
      width = 3.65
    ) %>%
    flextable::width(
      j = effect_column,
      width = 2.05
    ) %>%
    flextable::width(
      j = "p value",
      width = 0.75
    ) %>%
    flextable::set_table_properties(
      layout = "fixed",
      opts_word = list(
        split = TRUE,
        keep_with_next = FALSE,
        repeat_headers = TRUE
      )
    )
}

model_table_specifications <- tibble(
  table_no = c(
    "S8", "S9", "S10", "S11", "S12", "S13",
    "S14", "S17", "S18", "S19", "S20", "S21"
  ),
  model_name = c(
    "glm_ihm",
    "ttd_m1",
    "ttd_m2",
    "ttd_m2_80",
    "ttr_m1",
    "ttr_m2",
    "ttr_m2_80",
    "ttd_frailty",
    "ttr_frailty",
    "ttd_m2_type",
    "ttr_m2_type",
    "ttr_finegray"
  ),
  effect_type = c(
    "OR",
    rep("HR", 10),
    "subdistribution HR"
  ),
  caption = c(
    "Adjusted model for in-hospital mortality",
    "Adjusted model 1 for post-discharge mortality",
    "Adjusted model 2 for post-discharge mortality",
    "PU-learning model for post-discharge mortality",
    "Adjusted model 1 for readmission",
    "Adjusted model 2 for readmission",
    "PU-learning model for readmission",
    paste0(
      "Shared-frailty sensitivity analysis for ",
      "post-discharge mortality"
    ),
    "Shared-frailty sensitivity analysis for readmission",
    paste0(
      "Hospital-type sensitivity analysis for ",
      "post-discharge mortality"
    ),
    "Hospital-type sensitivity analysis for readmission",
    paste0(
      "Fine–Gray competing-risk sensitivity analysis ",
      "for readmission"
    )
  ),
  fit = list(
    glm_ihm,
    ttd_m1,
    ttd_m2,
    ttd_m2_80,
    ttr_m1,
    ttr_m2,
    ttr_m2_80,
    ttd_frailty,
    ttr_frailty,
    ttd_m2_type,
    ttr_m2_type,
    ttr_finegray
  )
)

model_table_note <- function(table_no) {
  if (table_no == "S8") {
    return(
      paste(
        "OR=odds ratio; CI=confidence interval;",
        "ECI=Elixhauser Comorbidity Index;",
        "HAP=hospital-acquired pneumonia."
      )
    )
  }
  
  if (table_no %in% c("S17", "S18")) {
    return(
      paste(
        "HR=hazard ratio; CI=confidence interval;",
        "ECI=Elixhauser Comorbidity Index;",
        "HAP=hospital-acquired pneumonia.",
        paste0(
          "The model includes a hospital-level gamma ",
          "shared-frailty term."
        )
      )
    )
  }
  
  if (table_no == "S21") {
    return(
      paste(
        "Subdistribution HR=subdistribution hazard ratio;",
        "CI=confidence interval;",
        "ECI=Elixhauser Comorbidity Index;",
        "HAP=hospital-acquired pneumonia.",
        "Death was treated as a competing event;",
        "standard errors were clustered by hospital."
      )
    )
  }
  
  paste(
    "HR=hazard ratio; CI=confidence interval;",
    "ECI=Elixhauser Comorbidity Index;",
    "HAP=hospital-acquired pneumonia.",
    "Standard errors were clustered by hospital."
  )
}

model_table_results <- purrr::pmap(
  model_table_specifications,
  function(
      table_no,
      model_name,
      effect_type,
      caption,
      fit) {
    table_data <- make_model_table(
      fit,
      effect_type = effect_type,
      strict_labels = TRUE
    )

    list(
      table_no = table_no,
      model_name = model_name,
      caption = caption,
      note = model_table_note(table_no),
      data = table_data,
      flextable = as_flextable(table_data)
    )
  }
)

model_table_document <- officer::read_docx()

for (i in seq_along(model_table_results)) {
  current_table <- model_table_results[[i]]
  caption_text <- paste0(
    "Table ",
    current_table$table_no,
    ". ",
    current_table$caption
  )

  model_table_document <- officer::body_add_fpar(
    model_table_document,
    officer::fpar(
      officer::ftext(
        caption_text,
        officer::fp_text(
          font.family = "Times New Roman",
          font.size = 10,
          bold = TRUE
        )
      )
    )
  )

  model_table_document <- flextable::body_add_flextable(
    model_table_document,
    current_table$flextable
  )

  model_table_document <- officer::body_add_fpar(
    model_table_document,
    officer::fpar(
      officer::ftext(
        current_table$note,
        officer::fp_text(
          font.family = "Times New Roman",
          font.size = 8
        )
      )
    )
  )

  if (i < length(model_table_results)) {
    model_table_document <- officer::body_add_break(
      model_table_document
    )
  }
}

print(
  model_table_document,
  target = model_tables_file
)

model_table_data <- purrr::map(
  model_table_results,
  "data"
)

names(model_table_data) <-
  model_table_specifications$table_no

rank_df %>%
  summarise(
    n_hospitals = n(),
    n_switchers = sum(switcher == 1L),
    pct_switchers = 100 * mean(switcher == 1L)
  )

# Figure 1: Interval-specific adjusted hazard ratios ---------------------------

extract_hap_interval_effects <- function(
    fit,
    outcome_label) {
  effects <- extract_model_estimates(fit) %>%
    filter(stringr::str_detect(term, "HAP_rule")) %>%
    mutate(
      interval = vapply(
        term,
        extract_interval,
        character(1)
      )
    ) %>%
    filter(!is.na(interval)) %>%
    transmute(
      outcome = outcome_label,
      interval,
      hr = estimate,
      lo = conf_low,
      hi = conf_high
    )

  if (nrow(effects) != 3L ||
      n_distinct(effects$interval) != 3L) {
    stop(
      "Expected three interval-specific HAP effects for ",
      outcome_label,
      "."
    )
  }

  effects
}

hr_plot_df <- bind_rows(
  extract_hap_interval_effects(
    ttd_m2,
    "Mortality"
  ),
  extract_hap_interval_effects(
    ttr_m2,
    "Readmission"
  )
) %>%
  mutate(
    outcome = factor(
      outcome,
      levels = c("Mortality", "Readmission")
    ),
    interval = factor(
      interval,
      levels = c(
        "91–365 days",
        "31–90 days",
        "0–30 days"
      )
    ),
    label = paste0(
      use_decimal_point(sprintf("%.2f", hr)),
      " (",
      use_decimal_point(sprintf("%.2f", lo)),
      "–",
      use_decimal_point(sprintf("%.2f", hi)),
      ")"
    )
  ) %>%
  group_by(outcome) %>%
  mutate(label_x = max(hi) + 0.04) %>%
  ungroup()

hr_x_min <- min(0.95, min(hr_plot_df$lo) - 0.03)
hr_x_max <- max(hr_plot_df$label_x) + 0.28

figure1_plot <- ggplot(
  hr_plot_df,
  aes(x = hr, y = interval)
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey50"
  ) +
  geom_errorbar(
    aes(xmin = lo, xmax = hi),
    orientation = "y",
    width = 0.12,
    linewidth = 0.45,
    colour = "grey45"
  ) +
  geom_point(
    size = 2.3,
    colour = "black"
  ) +
  geom_text(
    aes(
      x = label_x,
      label = label
    ),
    hjust = 0,
    size = 3.5,
    family = plot_font_family
  ) +
  facet_wrap(
    ~ outcome,
    ncol = 2,
    labeller = as_labeller(
      c(
        Mortality = "a  Mortality",
        Readmission = "b  Readmission"
      )
    )
  ) +
  scale_x_continuous(
    limits = c(hr_x_min, hr_x_max),
    breaks = scales::breaks_pretty(n = 4),
    labels = scales::label_number(
      accuracy = 0.1,
      decimal.mark = "."
    )
  ) +
  labs(
    x = "Adjusted hazard ratio (HAP relative to CAP)",
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(
    base_family = plot_font_family,
    base_size = 10
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    axis.text.y = element_text(face = "plain"),
    axis.title.x = element_text(
      margin = margin(t = 8)
    ),
    plot.margin = margin(10, 60, 10, 10),
    legend.position = "none"
  )

ggplot2::ggsave(
  filename = figure1_pdf,
  plot = figure1_plot,
  width = 180,
  height = 75,
  units = "mm",
  device = grDevices::cairo_pdf
)

ggplot2::ggsave(
  filename = figure1_tiff,
  plot = figure1_plot,
  width = 180,
  height = 75,
  units = "mm",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)


# Figure 2: CAP versus HAP hospital mortality ----------------------------------

max_mortality_limit <- max(
  plot_df$cap_risk_pct,
  plot_df$hap_risk_pct,
  na.rm = TRUE
)

max_mortality_limit <- ceiling(
  max_mortality_limit * 2
) / 2

figure2_plot <- ggplot(
  plot_df,
  aes(
    x = cap_risk_pct,
    y = hap_risk_pct
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    colour = "grey55"
  ) +
  geom_point(
    aes(size = total_n),
    shape = 21,
    fill = "grey35",
    colour = "black",
    stroke = 0.25,
    alpha = 0.65
  ) +
  scale_size_continuous(
    range = c(1, 4),
    name = "Pneumonia\nadmissions",
    labels = scales::comma
  ) +
  scale_x_continuous(
    limits = c(0, max_mortality_limit),
    expand = expansion(
      mult = c(0.02, 0.05)
    )
  ) +
  scale_y_continuous(
    limits = c(0, max_mortality_limit),
    expand = expansion(
      mult = c(0.02, 0.05)
    )
  ) +
  labs(
    x = "CAP risk-standardised 30-day mortality (%)",
    y = "HAP risk-standardised 30-day mortality (%)"
  ) +
  coord_equal(clip = "off") +
  theme_minimal(
    base_family = plot_font_family,
    base_size = 10
  ) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(
      linewidth = 0.2,
      colour = "grey90"
    ),
    axis.title.y = element_text(
      margin = margin(r = 8)
    ),
    axis.title.x = element_text(
      margin = margin(t = 8)
    ),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(6, 12, 6, 6)
  )

ggplot2::ggsave(
  filename = figure2_pdf,
  plot = figure2_plot,
  width = 107,
  height = 85,
  units = "mm",
  device = grDevices::cairo_pdf
)

ggplot2::ggsave(
  filename = figure2_tiff,
  plot = figure2_plot,
  width = 107,
  height = 85,
  units = "mm",
  dpi = 300,
  device = "tiff",
  compression = "lzw"
)


# Save compact, dynamically calculated analysis results -----------------------

model_formulas <- list(
  in_hospital_mortality = formula(glm_ihm),
  mortality_model_1 = f_ttd_m1,
  mortality_model_2 = f_ttd_m2,
  readmission_model_1 = f_ttr_m1,
  readmission_model_2 = f_ttr_m2,
  mortality_pu80 = f_ttd_pu_m2,
  readmission_pu80 = f_ttr_pu_m2,
  mortality_hospital_type = f_ttd_m2_type,
  readmission_hospital_type = f_ttr_m2_type
)

analysis_session_info <- sessionInfo()

save(
  u69_hospital_year,
  u69_hospital_use,
  u69_volume_check,
  u69_stable_large,
  u69_rule_by_year,
  icc_mor_table,
  concordance_table,
  proportional_hazards_tests,
  plot_df,
  rank_df,
  spearman_test,
  quartile_agreement_table,
  quartile_kappa,
  performance_class_table,
  profiling_summary,
  switch_df,
  switcher_characteristics,
  switcher_volume_table,
  model_table_data,
  hr_plot_df,
  model_formulas,
  analysis_session_info,
  file = analysis_results_file
)





