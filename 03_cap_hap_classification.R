# Phenotyping Hospital-Acquired Pneumonia in Administrative Data
# 03 - CAP/HAP phenotyping
#
# Expected input:
#   PNEUMONIA_DATA_DIR/df_CV.RData containing an object named `df_CV`
#
# Created outputs:
#   PNEUMONIA_DATA_DIR/df_classified.RData containing `df`
#   PNEUMONIA_DATA_DIR/classification_checks.RData containing dynamically
#   calculated rule-based and PU-learning checks
#   PNEUMONIA_OUTPUT_DIR/Table2_reclassification_phenotyping.docx
#   PNEUMONIA_OUTPUT_DIR/pu_feature_importance.docx
#
# Set the directories locally, for example in an untracked .Renviron file:
#   PNEUMONIA_DATA_DIR=/path/to/derived/data
#   PNEUMONIA_OUTPUT_DIR=/path/to/analysis/output
#
# If PNEUMONIA_OUTPUT_DIR is not set, a sibling directory named `03_Output`
# is used.

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(xgboost)
  library(irr)
  library(flextable)
})


# Analysis settings ------------------------------------------------------------

random_seed <- 123L
n_folds <- 5L
n_rounds <- 2000L
recall_targets <- c(PU70 = 0.70, PU80 = 0.80, PU90 = 0.90)


# File paths -------------------------------------------------------------------

data_dir <- Sys.getenv("PNEUMONIA_DATA_DIR", unset = "")

if (!nzchar(data_dir)) {
  stop(
    "PNEUMONIA_DATA_DIR is not set. Point it to the directory containing ",
    "df_CV.RData."
  )
}

data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)

output_dir <- Sys.getenv(
  "PNEUMONIA_OUTPUT_DIR",
  unset = file.path(dirname(data_dir), "03_Output")
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

input_file <- file.path(data_dir, "df_CV.RData")
classified_file <- file.path(data_dir, "df_classified.RData")
checks_file <- file.path(data_dir, "classification_checks.RData")

reclassification_table_file <- file.path(
  output_dir,
  "Table2_reclassification_phenotyping.docx"
)

feature_importance_file <- file.path(
  output_dir,
  "pu_feature_importance.docx"
)


# Load and validate the analysis cohort ---------------------------------------

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

df_CV <- load_named_object(input_file, "df_CV")

require_columns(
  df_CV,
  c(
    "case_id",
    "year",
    "hospital_id",
    "hap",
    "admission_type",
    "place_of_stay",
    "pneumonia_main_diagnosis",
    "pneumonia_secondary_diagnosis",
    "age",
    "gender",
    "eci",
    "length_of_stay",
    "icu_stay",
    "ventilation_hours"
  ),
  "df_CV"
)

if (anyDuplicated(df_CV[c("case_id", "year")])) {
  stop("The analysis cohort contains duplicate case_id-year records.")
}


# Rule-based HAP definition ----------------------------------------------------

rule_source_levels <- c(
  "ICD-10 HAP",
  "Scheduled admission",
  "Transfer and secondary diagnosis",
  "Prior acute stay and secondary diagnosis",
  "Prior psychiatric/rehabilitation stay",
  "CAP"
)

df <- df_CV %>%
  mutate(
    HAP_rule_source = case_when(
      hap == 1L ~ "ICD-10 HAP",
      admission_type == "Scheduled" ~ "Scheduled admission",
      admission_type == "Transfer" &
        pneumonia_secondary_diagnosis == 1L ~ "Transfer and secondary diagnosis",
      place_of_stay == "Acute" &
        pneumonia_secondary_diagnosis == 1L ~
        "Prior acute stay and secondary diagnosis",
      place_of_stay == "Psych/Rehab" ~
        "Prior psychiatric/rehabilitation stay",
      TRUE ~ "CAP"
    ),
    HAP_rule_source = factor(HAP_rule_source, levels = rule_source_levels),
    HAP_rule = as.integer(HAP_rule_source != "CAP")
  )

rule_counts <- df %>%
  count(HAP_rule_source, HAP_rule, name = "n") %>%
  mutate(share = n / sum(n))


# Positive-unlabelled classifier ----------------------------------------------

df_pu <- df %>%
  mutate(
    hospital_id = as.character(hospital_id),
    y = as.integer(hap == 1L),
    year = factor(year)
  )

if (anyNA(df_pu$hospital_id)) {
  stop("Hospital identifiers must not be missing for hospital-grouped CV.")
}

categorical_variables <- c(
  "age",
  "admission_type",
  "place_of_stay",
  "gender",
  "year"
)

numeric_variables <- c(
  "pneumonia_main_diagnosis",
  "eci"
)

pu_formula <- as.formula(
  paste0(
    "~ -1 + ",
    paste(c(categorical_variables, numeric_variables), collapse = " + ")
  )
)

X <- Matrix::sparse.model.matrix(pu_formula, data = df_pu)
y <- df_pu$y

if (nrow(X) != nrow(df_pu)) {
  stop("The PU model matrix does not have one row per admission.")
}

set.seed(random_seed)

hospital_ids <- sort(unique(df_pu$hospital_id))
fold_id <- sample(rep(seq_len(n_folds), length.out = length(hospital_ids)))

fold_map <- tibble(
  hospital_id = hospital_ids,
  fold = fold_id
)

df_pu <- df_pu %>%
  left_join(fold_map, by = "hospital_id")

if (anyNA(df_pu$fold)) {
  stop("At least one admission was not assigned to a CV fold.")
}

fold_counts <- df_pu %>%
  count(fold, y, name = "n")

positive_counts_by_fold <- df_pu %>%
  group_by(fold) %>%
  summarise(n_positive = sum(y == 1L), .groups = "drop")

if (nrow(positive_counts_by_fold) != n_folds ||
    any(positive_counts_by_fold$n_positive == 0L)) {
  stop("Every CV fold must contain at least one ICD-10 HAP-positive admission.")
}

base_parameters <- list(
  objective = "binary:logistic",
  eval_metric = "logloss",
  eta = 0.05,
  max_depth = 5L,
  min_child_weight = 1,
  subsample = 0.8,
  colsample_bytree = 0.8,
  seed = random_seed
)

oof_prediction <- rep(NA_real_, nrow(df_pu))

for (k in seq_len(n_folds)) {
  training_rows <- which(df_pu$fold != k)
  test_rows <- which(df_pu$fold == k)
  
  training_labels <- y[training_rows]
  
  if (!length(test_rows) ||
      !any(training_labels == 0L) ||
      !any(training_labels == 1L)) {
    stop("Each CV split must contain test observations and both training classes.")
  }
  
  scale_positive_weight <-
    sum(training_labels == 0L) / sum(training_labels == 1L)
  
  fold_parameters <- c(
    base_parameters,
    list(scale_pos_weight = scale_positive_weight)
  )
  
  training_matrix <- xgboost::xgb.DMatrix(
    data = X[training_rows, ],
    label = training_labels
  )
  
  test_matrix <- xgboost::xgb.DMatrix(data = X[test_rows, ])
  
  fold_model <- xgboost::xgb.train(
    params = fold_parameters,
    data = training_matrix,
    nrounds = n_rounds,
    verbose = 0
  )
  
  oof_prediction[test_rows] <- predict(fold_model, test_matrix)
}

if (anyNA(oof_prediction)) {
  stop("Out-of-fold predictions are missing for at least one admission.")
}

df_pu$hap_pu_score <- oof_prediction


# Recall-based PU thresholds ---------------------------------------------------

positive_scores <- df_pu$hap_pu_score[df_pu$y == 1L]

pu_cutoffs <- vapply(
  recall_targets,
  function(target) {
    as.numeric(
      quantile(
        positive_scores,
        probs = 1 - target,
        na.rm = TRUE,
        names = FALSE
      )
    )
  },
  numeric(1)
)

df_pu <- df_pu %>%
  mutate(
    HAP_PU70 = as.integer(hap_pu_score >= pu_cutoffs[["PU70"]]),
    HAP_PU80 = as.integer(hap_pu_score >= pu_cutoffs[["PU80"]]),
    HAP_PU90 = as.integer(hap_pu_score >= pu_cutoffs[["PU90"]])
  )

pu_thresholds <- tibble(
  definition = names(recall_targets),
  recall_target = unname(recall_targets),
  cutoff = unname(pu_cutoffs),
  recall_achieved = c(
    mean(df_pu$HAP_PU70[df_pu$y == 1L] == 1L),
    mean(df_pu$HAP_PU80[df_pu$y == 1L] == 1L),
    mean(df_pu$HAP_PU90[df_pu$y == 1L] == 1L)
  ),
  n_hap = c(
    sum(df_pu$HAP_PU70 == 1L),
    sum(df_pu$HAP_PU80 == 1L),
    sum(df_pu$HAP_PU90 == 1L)
  )
) %>%
  mutate(share = n_hap / nrow(df_pu))

# Agreement and reclassification ---------------------------------------------

calculate_kappa <- function(first, second, comparison) {
  fit <- irr::kappa2(
    data.frame(first = first, second = second),
    weight = "unweighted"
  )
  
  tibble(
    comparison = comparison,
    kappa = unname(fit$value),
    p_value = unname(fit$p.value)
  )
}

agreement_results <- bind_rows(
  calculate_kappa(df_pu$y, df_pu$HAP_rule, "ICD-10 HAP vs rule-based HAP"),
  calculate_kappa(df_pu$y, df_pu$HAP_PU80, "ICD-10 HAP vs PU80 HAP"),
  calculate_kappa(
    df_pu$HAP_rule,
    df_pu$HAP_PU80,
    "Rule-based HAP vs PU80 HAP"
  )
)

classification_summary <- tibble(
  definition = c(
    "ICD-10 HAP",
    "Rule-based HAP",
    "PU70 HAP",
    "PU80 HAP",
    "PU90 HAP"
  ),
  n_hap = c(
    sum(df_pu$y == 1L),
    sum(df_pu$HAP_rule == 1L),
    sum(df_pu$HAP_PU70 == 1L),
    sum(df_pu$HAP_PU80 == 1L),
    sum(df_pu$HAP_PU90 == 1L)
  )
) %>%
  mutate(share = n_hap / nrow(df_pu))

rule_pu80_table <- table(
  Rule_based = df_pu$HAP_rule,
  PU80 = df_pu$HAP_PU80
)

df <- df_pu %>%
  mutate(
    class_group = case_when(
      HAP_rule == 0L & HAP_PU80 == 0L ~ "CAP_both",
      HAP_rule == 1L & HAP_PU80 == 1L ~ "HAP_both",
      HAP_rule == 0L & HAP_PU80 == 1L ~ "PU_only_HAP",
      HAP_rule == 1L & HAP_PU80 == 0L ~ "Rule_only_HAP"
    ),
    class_group = factor(
      class_group,
      levels = c(
        "CAP_both",
        "Rule_only_HAP",
        "PU_only_HAP",
        "HAP_both"
      )
    ),
    icu_admission = as.integer(dplyr::coalesce(icu_stay, 0) > 0),
    mechanical_ventilation = as.integer(
      dplyr::coalesce(ventilation_hours, 0) > 0
    )
  ) %>%
  select(-y, -fold)

if (anyNA(df$class_group)) {
  stop("At least one admission could not be assigned to a reclassification group.")
}


# Reclassification table -------------------------------------------------------

fmt <- function(x, digits = 2L) {
  x <- round(x, digits)
  x[abs(x) < 10^(-digits)] <- 0
  output <- formatC(x, format = "f", digits = digits)
  output <- sub("^-", "−", output)
  gsub("\\.", ".", output)
}

fmt_pct <- function(x, digits = 1L) {
  paste0(fmt(100 * x, digits), "%")
}

reclassification_table <- df %>%
  group_by(class_group) %>%
  summarise(
    n = n(),
    eci = median(eci, na.rm = TRUE),
    length_of_stay = median(length_of_stay, na.rm = TRUE),
    icu = mean(icu_admission, na.rm = TRUE),
    ventilation = mean(mechanical_ventilation, na.rm = TRUE),
    secondary_diagnosis = mean(
      pneumonia_secondary_diagnosis == 1L,
      na.rm = TRUE
    ),
    emergency_admission = mean(
      admission_type == "Emergency",
      na.rm = TRUE
    ),
    home_before_admission = mean(place_of_stay == "Home", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    class_group = factor(
      class_group,
      levels = c(
        "CAP_both",
        "Rule_only_HAP",
        "PU_only_HAP",
        "HAP_both"
      ),
      labels = c(
        "CAP by both approaches",
        "Rule-based HAP only",
        "PU learning HAP only",
        "HAP by both approaches"
      )
    )
  ) %>%
  arrange(class_group) %>%
  transmute(
    `Phenotyping group` = class_group,
    Hospitalisations = formatC(n, format = "d", big.mark = ","),
    `ECI, median` = fmt(eci, 0L),
    `Length of stay, median days` = fmt(length_of_stay, 0L),
    `ICU admission` = fmt_pct(icu, 1L),
    `Mechanical ventilation` = fmt_pct(ventilation, 1L),
    `Pneumonia as secondary diagnosis` =
      fmt_pct(secondary_diagnosis, 1L),
    `Emergency admission` = fmt_pct(emergency_admission, 1L),
    `Home before admission` = fmt_pct(home_before_admission, 1L)
  )

reclassification_flextable <-
  flextable::flextable(reclassification_table) %>%
  flextable::set_header_labels(`Phenotyping group` = "") %>%
  flextable::align(part = "header", align = "center") %>%
  flextable::align(j = 1L, align = "left", part = "all") %>%
  flextable::align(j = 2:9, align = "right", part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::autofit() %>%
  flextable::add_footer_lines(
    paste0(
      "Values are n, median, or %. CAP = community-acquired pneumonia; ",
      "HAP = hospital-acquired pneumonia; PU = positive–unlabelled; ",
      "ECI = Elixhauser Comorbidity Index; ICU = intensive care unit."
    )
  ) %>%
  flextable::add_footer_lines(
    paste0(
      "Groups are defined by agreement or disagreement between the primary ",
      "rule-based algorithm and the PU-learning HAP definition using the ",
      "80% recall threshold."
    )
  )

flextable::save_as_docx(
  "Reclassification between phenotyping approaches" =
    reclassification_flextable,
  path = reclassification_table_file
)

# Supplementary PU feature importance -----------------------------------------

all_data_matrix <- xgboost::xgb.DMatrix(data = X, label = y)
scale_positive_weight_all <- sum(y == 0L) / sum(y == 1L)

all_data_parameters <- c(
  base_parameters,
  list(scale_pos_weight = scale_positive_weight_all)
)

all_data_model <- xgboost::xgb.train(
  params = all_data_parameters,
  data = all_data_matrix,
  nrounds = n_rounds,
  verbose = 0
)

feature_importance <- xgboost::xgb.importance(
  feature_names = colnames(X),
  model = all_data_model
)

term_labels <- c(
  "age20-24" = "20–24 years",
  "age25-29" = "25–29 years",
  "age30-34" = "30–34 years",
  "age35-39" = "35–39 years",
  "age40-44" = "40–44 years",
  "age45-49" = "45–49 years",
  "age50-54" = "50–54 years",
  "age55-59" = "55–59 years",
  "age60-64" = "60–64 years",
  "age65-69" = "65–69 years",
  "age70-74" = "70–74 years",
  "age75-79" = "75–79 years",
  "age80-84" = "80–84 years",
  "age85-89" = "85–89 years",
  "age90-94" = "90–94 years",
  "age95+" = "95+ years",
  "genderFemale" = "Female",
  "admission_typeScheduled" = "Scheduled admission",
  "admission_typeTransfer" = "Transfer admission",
  "admission_typeEmergency" = "Emergency admission",
  "admission_typeOther_Unknown" = "Other/unknown admission",
  "place_of_stayHome" = "Home",
  "place_of_stayNursing/Residential" = "Nursing/residential care",
  "place_of_stayAcute" = "Acute care hospital",
  "place_of_stayPsych/Rehab" = "Psychiatry/rehabilitation",
  "place_of_stayOther/Unknown" = "Place of prior stay unknown",
  "eci" = "Elixhauser comorbidity index (ECI)",
  "pneumonia_main_diagnosis" = "Pneumonia as main diagnosis",
  "year2010" = "2010",
  "year2011" = "2011",
  "year2012" = "2012",
  "year2013" = "2013",
  "year2014" = "2014",
  "year2015" = "2015",
  "year2016" = "2016",
  "year2017" = "2017",
  "year2018" = "2018",
  "year2019" = "2019"
)

importance_table <- feature_importance %>%
  mutate(
    Variable = if_else(
      Feature %in% names(term_labels),
      unname(term_labels[Feature]),
      Feature
    ),
    Gain_share = Gain / sum(Gain)
  ) %>%
  arrange(desc(Gain_share)) %>%
  transmute(
    Variable,
    `Share of total model gain` = fmt_pct(Gain_share, 1L),
    Gain = fmt(Gain, 3L),
    Cover = fmt(Cover, 3L),
    Frequency = fmt(Frequency, 3L)
  )

importance_flextable <- flextable::flextable(importance_table) %>%
  flextable::set_header_labels(Variable = "") %>%
  flextable::align(part = "header", align = "center") %>%
  flextable::align(j = 1L, align = "left", part = "all") %>%
  flextable::align(j = 2:5, align = "right", part = "all") %>%
  flextable::bold(part = "header") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::autofit() %>%
  flextable::add_footer_lines(
    paste0(
      "Feature importance is gain-based and does not imply causal effects. ",
      "Predictors were derived from administrative discharge data. ",
      "ECI = Elixhauser Comorbidity Index."
    )
  )

flextable::save_as_docx(
  "PU feature importance (XGBoost)" = importance_flextable,
  path = feature_importance_file
)


# Save classified data and reproducibility checks -----------------------------

pu_model_settings <- list(
  random_seed = random_seed,
  n_folds = n_folds,
  n_rounds = n_rounds,
  base_parameters = base_parameters,
  recall_targets = recall_targets
)


save(
  rule_counts,
  fold_counts,
  pu_thresholds,
  agreement_results,
  classification_summary,
  rule_pu80_table,
  reclassification_table,
  importance_table,
  pu_model_settings,
  file = checks_file
)

rule_counts
classification_summary
pu_thresholds
agreement_results

rule_detail_counts <- df %>%
  mutate(
    rule_step = case_when(
      hap == 1L ~ "ICD-10 HAP",
      
      admission_type == "Scheduled" ~
        "Scheduled admission",
      
      admission_type == "Transfer" &
        pneumonia_secondary_diagnosis == 1L ~
        "Transfer and secondary pneumonia",
      
      place_of_stay == "Acute" &
        pneumonia_secondary_diagnosis == 1L ~
        "Prior acute stay and secondary pneumonia",
      
      place_of_stay == "Psych/Rehab" ~
        "Prior Psych/Rehab stay",
      
      TRUE ~ "CAP"
    )
  ) %>%
  count(rule_step, name = "n") %>%
  mutate(
    percentage = round(100 * n / nrow(df), 1)
  )

rule_detail_counts

