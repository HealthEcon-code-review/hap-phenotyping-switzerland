# Code Repository – Phenotyping Hospital-Acquired Pneumonia in Administrative Data

This repository contains the R code used for cohort construction, pneumonia phenotyping, statistical analyses, and figure generation in the manuscript **“Phenotyping hospital-acquired pneumonia in administrative data: implications for outcome measurement and hospital performance assessment in a Swiss nationwide cohort study”**, submitted for peer review to *The Lancet Regional Health – Europe*.

The code includes:

- Construction of the nationwide cohort of adult pneumonia hospitalisations
- Implementation and comparison of three approaches to identifying hospital-acquired pneumonia: ICD-based classification, a rule-based algorithm, and positive–unlabelled learning
- Analyses of in-hospital mortality, post-discharge mortality, and hospital readmission
- Hospital-level outcome profiling
- Sensitivity analyses and generation of manuscript tables and figures

## Usage

All analyses were conducted in R version 4.5.2 (R Foundation for Statistical Computing, Vienna, Austria). Scripts are numbered according to their intended execution order. Because the underlying administrative data are not publicly available, the code cannot be executed without authorised access to the original dataset.

## Structure

- 01_elixhauser_comorbidity_index.R – Derivation of the Elixhauser Comorbidity Index from ICD-10-GM diagnoses using the Quan mapping and Swiss-specific weights.
- 02_data_cleaning.R – Data preparation, cohort construction, variable derivation, and application of the study eligibility criteria.
- 03_cap_hap_classification.R – Implementation and comparison of the ICD-based, rule-based, and positive–unlabelled learning HAP phenotypes.
- 04_baseline_characteristics.R – Descriptive analyses and generation of baseline characteristic tables.
- 05_outcome_analysis_tables_visualisations.R – Outcome modelling, hospital-level analyses, sensitivity analyses, and generation of manuscript tables and figures.
- README.md – Repository overview, usage information, and data-availability statement.

## Data availability

The analyses use administrative hospital data provided by the Swiss Federal Statistical Office under a data-use agreement. These data are not publicly available, and access requires approval from the Swiss Federal Statistical Office. The authors are not permitted to distribute the individual-level data. No patient-level data are included in this repository.

> ⚠️ This repository contains the peer-review version of the analysis code. A citable, versioned release will be archived on Zenodo upon acceptance.
