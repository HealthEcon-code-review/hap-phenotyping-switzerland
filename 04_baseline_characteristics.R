# Phenotyping Hospital-Acquired Pneumonia in Administrative Data
# 04 - Baseline characteristics and construct-validity tables
#
# Expected input:
#   PNEUMONIA_DATA_DIR/df_classified.RData containing an object named `df`
#
# Created outputs:
#   PNEUMONIA_OUTPUT_DIR/Table1_patient_characteristics_rule_based.docx
#   PNEUMONIA_OUTPUT_DIR/TableS1_alternative_HAP_definitions.docx
#   PNEUMONIA_OUTPUT_DIR/TableS2_rule_based_components.docx
#
# Set the directories locally, for example in an untracked .Renviron file:
#   PNEUMONIA_DATA_DIR=/path/to/derived/data
#   PNEUMONIA_OUTPUT_DIR=/path/to/analysis/output
#
# If PNEUMONIA_OUTPUT_DIR is not set, a sibling directory named `03_Output`
# is used.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(flextable)
  library(officer)
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

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

input_file <- file.path(data_dir, "df_classified.RData")

table1_file <- file.path(
  output_dir,
  "Table1_patient_characteristics_rule_based.docx"
)

table_s1_file <- file.path(
  output_dir,
  "TableS1_alternative_HAP_definitions.docx"
)

table_s2_file <- file.path(
  output_dir,
  "TableS2_rule_based_components.docx"
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

assert_binary <- function(x, variable_name) {
  observed_values <- sort(unique(stats::na.omit(x)))
  
  if (anyNA(x) || !all(observed_values %in% c(0, 1))) {
    stop("`", variable_name, "` must contain only 0 and 1 without missing values.")
  }
}

df <- load_named_object(input_file, "df")

require_columns(
  df,
  c(
    "case_id",
    "year",
    "hospital_id",
    "age",
    "gender",
    "eci",
    "pneumonia_secondary_diagnosis",
    "admission_type",
    "place_of_stay",
    "icu_stay",
    "ventilation_hours",
    "length_of_stay",
    "hap",
    "HAP_rule",
    "HAP_rule_source",
    "HAP_PU80"
  ),
  "df"
)

if (anyDuplicated(df[c("case_id", "year")])) {
  stop("The classified cohort contains duplicate case_id-year records.")
}

assert_binary(df$hap, "hap")
assert_binary(df$HAP_rule, "HAP_rule")
assert_binary(df$HAP_PU80, "HAP_PU80")

if (anyNA(df$HAP_rule_source)) {
  stop("`HAP_rule_source` contains missing values.")
}

rule_indicator_from_source <- as.integer(
  as.character(df$HAP_rule_source) != "CAP"
)

if (!identical(rule_indicator_from_source, as.integer(df$HAP_rule))) {
  stop("`HAP_rule` and `HAP_rule_source` are inconsistent.")
}


# Age-band midpoints -----------------------------------------------------------

# Age is available only in five-year bands. Midpoints are used to approximate
# the cohort median and IQR; the open-ended 95+ group is represented by 97.
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

df_characteristics <- df %>%
  mutate(
    age = as.character(age),
    icu_any = dplyr::coalesce(icu_stay, 0) > 0,
    ventilation_any = dplyr::coalesce(ventilation_hours, 0) > 0
  ) %>%
  left_join(age_lookup, by = "age")


# Formatting functions --------------------------------------------------------

fmt <- function(x, digits = 1L) {
  x <- round(x, digits)
  x[abs(x) < 10^(-digits)] <- 0
  output <- formatC(x, format = "f", digits = digits)
  output <- sub("^-", "−", output)
  gsub("\\.", ".", output)
}

fmt_n <- function(x) {
  format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_n_pct <- function(n, denominator, digits = 1L) {
  paste0(
    fmt_n(n),
    " (",
    fmt(100 * n / denominator, digits),
    "%)"
  )
}

fmt_median_iqr <- function(median_value, q1, q3, digits = 1L) {
  paste0(
    fmt(median_value, digits),
    " (",
    fmt(q1, digits),
    "–",
    fmt(q3, digits),
    ")"
  )
}

format_characteristic_table <- function(summary_data, variable_levels) {
  summary_data %>%
    pivot_longer(
      -group_label,
      names_to = "Variable",
      values_to = "Value"
    ) %>%
    pivot_wider(names_from = group_label, values_from = Value) %>%
    mutate(
      Variable = factor(Variable, levels = variable_levels)
    ) %>%
    arrange(Variable) %>%
    mutate(Variable = as.character(Variable))
}

save_manuscript_table <- function(table, title, footer, path) {
  output_flextable <- flextable::flextable(table) %>%
    flextable::theme_vanilla() %>%
    flextable::set_header_labels(Variable = "") %>%
    flextable::bold(part = "header") %>%
    flextable::autofit() %>%
    flextable::add_footer_lines(footer)
  
  output_document <- officer::read_docx() %>%
    officer::body_add_par(title, style = "heading 1") %>%
    flextable::body_add_flextable(output_flextable, align = "center")
  
  print(output_document, target = path)
}


# Table 1: Primary rule-based classification ----------------------------------

table1_summary <- bind_rows(
  df_characteristics %>%
    mutate(
      group = if_else(
        HAP_rule == 1L,
        "Rule-based HAP",
        "Rule-based CAP"
      )
    ),
  df_characteristics %>%
    mutate(group = "All")
) %>%
  group_by(group) %>%
  summarise(
    n_cases = n(),
    age_median = median(age_num, na.rm = TRUE),
    age_q1 = unname(quantile(age_num, 0.25, na.rm = TRUE)),
    age_q3 = unname(quantile(age_num, 0.75, na.rm = TRUE)),
    female_n = sum(gender == "Female", na.rm = TRUE),
    eci_median = median(eci, na.rm = TRUE),
    eci_q1 = unname(quantile(eci, 0.25, na.rm = TRUE)),
    eci_q3 = unname(quantile(eci, 0.75, na.rm = TRUE)),
    secondary_pneumonia_n = sum(
      pneumonia_secondary_diagnosis == 1L,
      na.rm = TRUE
    ),
    emergency_n = sum(admission_type == "Emergency", na.rm = TRUE),
    home_n = sum(place_of_stay == "Home", na.rm = TRUE),
    icu_n = sum(icu_any, na.rm = TRUE),
    ventilation_n = sum(ventilation_any, na.rm = TRUE),
    los_median = median(length_of_stay, na.rm = TRUE),
    los_q1 = unname(quantile(length_of_stay, 0.25, na.rm = TRUE)),
    los_q3 = unname(quantile(length_of_stay, 0.75, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c("All", "Rule-based CAP", "Rule-based HAP")
    ),
    group_label = paste0(as.character(group), " (n=", fmt_n(n_cases), ")"),
    Age = fmt_median_iqr(age_median, age_q1, age_q3),
    Female = fmt_n_pct(female_n, n_cases),
    ECI = fmt_median_iqr(eci_median, eci_q1, eci_q3),
    `Pneumonia as secondary diagnosis` = fmt_n_pct(
      secondary_pneumonia_n,
      n_cases
    ),
    `Emergency admission` = fmt_n_pct(emergency_n, n_cases),
    `Home before admission` = fmt_n_pct(home_n, n_cases),
    `ICU admission` = fmt_n_pct(icu_n, n_cases),
    `Mechanical ventilation` = fmt_n_pct(ventilation_n, n_cases),
    `Length of stay, days` = fmt_median_iqr(
      los_median,
      los_q1,
      los_q3
    )
  ) %>%
  arrange(group) %>%
  select(
    group_label,
    Age,
    Female,
    ECI,
    `Pneumonia as secondary diagnosis`,
    `Emergency admission`,
    `Home before admission`,
    `ICU admission`,
    `Mechanical ventilation`,
    `Length of stay, days`
  )

table1 <- format_characteristic_table(
  table1_summary,
  variable_levels = c(
    "Age",
    "Female",
    "ECI",
    "Pneumonia as secondary diagnosis",
    "Emergency admission",
    "Home before admission",
    "ICU admission",
    "Mechanical ventilation",
    "Length of stay, days"
  )
)

save_manuscript_table(
  table = table1,
  title = paste0(
    "Table 1: Patient characteristics by primary rule-based pneumonia ",
    "classification"
  ),
  footer = paste0(
    "Values are median (IQR) for continuous variables and n (%) for ",
    "categorical variables. Age was available in five-year bands; medians ",
    "and IQRs were approximated using category midpoints. ",
    "ECI = Elixhauser Comorbidity Index; ICU = intensive care unit."
  ),
  path = table1_file
)


# Supplementary Table S1: Alternative HAP definitions -------------------------

# These groups intentionally overlap: each column describes one classification
# definition rather than mutually exclusive patient groups.
classification_summary <- bind_rows(
  df_characteristics %>%
    filter(hap == 1L) %>%
    mutate(group = "ICD-10 HAP"),
  df_characteristics %>%
    filter(HAP_PU80 == 1L) %>%
    mutate(group = "PU-learning HAP"),
  df_characteristics %>%
    filter(HAP_rule == 1L) %>%
    mutate(group = "Rule-based HAP"),
  df_characteristics %>%
    filter(HAP_rule == 0L) %>%
    mutate(group = "Rule-based CAP")
) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "ICD-10 HAP",
        "PU-learning HAP",
        "Rule-based HAP",
        "Rule-based CAP"
      )
    )
  ) %>%
  group_by(group) %>%
  summarise(
    n_cases = n(),
    eci_median = median(eci, na.rm = TRUE),
    eci_q1 = unname(quantile(eci, 0.25, na.rm = TRUE)),
    eci_q3 = unname(quantile(eci, 0.75, na.rm = TRUE)),
    secondary_n = sum(
      pneumonia_secondary_diagnosis == 1L,
      na.rm = TRUE
    ),
    emergency_n = sum(admission_type == "Emergency", na.rm = TRUE),
    icu_n = sum(icu_any, na.rm = TRUE),
    ventilation_n = sum(ventilation_any, na.rm = TRUE),
    los_median = median(length_of_stay, na.rm = TRUE),
    los_q1 = unname(quantile(length_of_stay, 0.25, na.rm = TRUE)),
    los_q3 = unname(quantile(length_of_stay, 0.75, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    group_label = paste0(as.character(group), " (n=", fmt_n(n_cases), ")"),
    ECI = fmt_median_iqr(eci_median, eci_q1, eci_q3),
    `Pneumonia as secondary diagnosis` = fmt_n_pct(secondary_n, n_cases),
    `Emergency admission` = fmt_n_pct(emergency_n, n_cases),
    `ICU admission` = fmt_n_pct(icu_n, n_cases),
    `Mechanical ventilation` = fmt_n_pct(ventilation_n, n_cases),
    `Length of stay, days` = fmt_median_iqr(
      los_median,
      los_q1,
      los_q3
    )
  ) %>%
  arrange(group) %>%
  select(
    group_label,
    ECI,
    `Pneumonia as secondary diagnosis`,
    `Emergency admission`,
    `ICU admission`,
    `Mechanical ventilation`,
    `Length of stay, days`
  )

classification_table <- format_characteristic_table(
  classification_summary,
  variable_levels = c(
    "ECI",
    "Pneumonia as secondary diagnosis",
    "Emergency admission",
    "ICU admission",
    "Mechanical ventilation",
    "Length of stay, days"
  )
)

save_manuscript_table(
  table = classification_table,
  title = paste0(
    "Supplementary Table S1: Clinical plausibility of alternative HAP ",
    "definitions"
  ),
  footer = paste0(
    "Values are median (IQR) for continuous variables and n (%) for ",
    "categorical variables. ECI = Elixhauser Comorbidity Index; ",
    "HAP = hospital-acquired pneumonia; ICU = intensive care unit; ",
    "PU-learning HAP uses the threshold corresponding to 80% recall of ",
    "ICD-10 HAP cases."
  ),
  path = table_s1_file
)


# Supplementary Table S2: Rule-based classification components ----------------

# HAP_rule_source was created by the ordered decision rule in script 03. Using
# that variable here prevents the component definitions from diverging from the
# actual rule-based classification.
rule_component_labels <- c(
  "ICD-10 HAP" = "ICD-10 HAP",
  "Scheduled admission" = "Scheduled admission",
  "Transfer and secondary diagnosis" =
    "Transfer admission and secondary pneumonia",
  "Prior acute stay and secondary diagnosis" =
    "Prior acute stay and secondary pneumonia",
  "Prior psychiatric/rehabilitation stay" =
    "Prior psychiatric or rehabilitation stay",
  "CAP" = "Rule-based CAP"
)

unexpected_rule_sources <- setdiff(
  unique(as.character(df_characteristics$HAP_rule_source)),
  names(rule_component_labels)
)

if (length(unexpected_rule_sources)) {
  stop(
    "Unexpected values in `HAP_rule_source`: ",
    paste(sort(unexpected_rule_sources), collapse = ", ")
  )
}

rule_component_summary <- df_characteristics %>%
  mutate(
    group = unname(
      rule_component_labels[as.character(HAP_rule_source)]
    ),
    group = factor(group, levels = unname(rule_component_labels))
  ) %>%
  group_by(group) %>%
  summarise(
    n_cases = n(),
    eci_median = median(eci, na.rm = TRUE),
    eci_q1 = unname(quantile(eci, 0.25, na.rm = TRUE)),
    eci_q3 = unname(quantile(eci, 0.75, na.rm = TRUE)),
    secondary_n = sum(
      pneumonia_secondary_diagnosis == 1L,
      na.rm = TRUE
    ),
    emergency_n = sum(admission_type == "Emergency", na.rm = TRUE),
    icu_n = sum(icu_any, na.rm = TRUE),
    ventilation_n = sum(ventilation_any, na.rm = TRUE),
    los_median = median(length_of_stay, na.rm = TRUE),
    los_q1 = unname(quantile(length_of_stay, 0.25, na.rm = TRUE)),
    los_q3 = unname(quantile(length_of_stay, 0.75, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    group_label = paste0(as.character(group), " (n=", fmt_n(n_cases), ")"),
    ECI = fmt_median_iqr(eci_median, eci_q1, eci_q3),
    `Pneumonia as secondary diagnosis` = fmt_n_pct(secondary_n, n_cases),
    `Emergency admission` = fmt_n_pct(emergency_n, n_cases),
    `ICU admission` = fmt_n_pct(icu_n, n_cases),
    `Mechanical ventilation` = fmt_n_pct(ventilation_n, n_cases),
    `Length of stay, days` = fmt_median_iqr(
      los_median,
      los_q1,
      los_q3
    )
  ) %>%
  arrange(group) %>%
  select(
    group_label,
    ECI,
    `Pneumonia as secondary diagnosis`,
    `Emergency admission`,
    `ICU admission`,
    `Mechanical ventilation`,
    `Length of stay, days`
  )

rule_component_table <- format_characteristic_table(
  rule_component_summary,
  variable_levels = c(
    "ECI",
    "Pneumonia as secondary diagnosis",
    "Emergency admission",
    "ICU admission",
    "Mechanical ventilation",
    "Length of stay, days"
  )
)

save_manuscript_table(
  table = rule_component_table,
  title = paste0(
    "Supplementary Table S2: Characteristics of rule-based classification ",
    "components"
  ),
  footer = paste0(
    "Values are median (IQR) for continuous variables and n (%) for ",
    "categorical variables. Components are mutually exclusive and follow ",
    "the ordered rule-based decision path. ECI = Elixhauser Comorbidity ",
    "Index; HAP = hospital-acquired pneumonia; ICU = intensive care unit."
  ),
  path = table_s2_file
)

