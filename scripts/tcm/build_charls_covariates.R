#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

argval <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

input_dta <- argval(
  "--harmonized",
  "data/raw/charls/downloads/harmonized/regular_waves/H_CHARLS_D_Data.dta"
)
linked_tsv <- argval(
  "--linked",
  "local CHARLS linked table"
)
context_tsv <- argval(
  "--context",
  "results/tcm/province_year_contextual_covariates.tsv"
)
covariates_tsv <- argval(
  "--covariates-output",
  "results/tcm/charls_covariates_harmonized.tsv"
)
analysis_tsv <- argval(
  "--analysis-output",
  "results/tcm/person_wave_tcm_core_density_analysis.tsv"
)
qc_tsv <- argval(
  "--qc-output",
  "results/tcm/charls_covariate_missingness.tsv"
)
standardization_tsv <- argval(
  "--standardization-output",
  "results/tcm/province_year_standardization_constants.tsv"
)

if (!file.exists(input_dta)) stop(sprintf("Missing Harmonized CHARLS file: %s", input_dta))
if (!file.exists(linked_tsv)) stop(sprintf("Missing linked TCM table: %s", linked_tsv))
if (!file.exists(context_tsv)) stop(sprintf("Missing contextual table: %s", context_tsv))

dir.create(dirname(covariates_tsv), recursive = TRUE, showWarnings = FALSE)

value <- function(df, var) {
  if (!var %in% names(df)) return(rep(NA_real_, nrow(df)))
  suppressWarnings(as.numeric(df[[var]]))
}

safe_log1p_positive <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x) & x > 0
  out[ok] <- log1p(x[ok])
  out
}

yes01 <- function(x) {
  ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, ifelse(x == 0, 0L, NA_integer_)))
}

income_var <- function(wave) paste0("hh", wave, "itot")
wealth_var <- function(wave) {
  if (wave %in% 1:2) paste0("hh", wave, "atotb") else paste0("h", wave, "atotb")
}

h <- read_dta(input_dta) %>%
  mutate(ID = as.character(ID))

chronic_suffixes <- c(
  "hibpe", "diabe", "cancre", "lunge", "hearte", "stroke",
  "arthre", "dyslipe", "livere", "kidneye", "digeste", "memrye"
)
chronic_names <- c(
  "hypertension", "diabetes", "cancer", "chronic_lung", "heart", "stroke",
  "arthritis", "dyslipidemia", "liver", "kidney", "digestive", "memory"
)

covariates <- bind_rows(lapply(1:4, function(wave) {
  year <- c(2011L, 2013L, 2015L, 2018L)[[wave]]

  chronic_vars <- paste0("r", wave, chronic_suffixes)
  chronic_mat <- sapply(chronic_vars, function(v) yes01(value(h, v)))
  if (is.null(dim(chronic_mat))) chronic_mat <- matrix(chronic_mat, ncol = length(chronic_vars))
  chronic_dx <- as_tibble(chronic_mat, .name_repair = "minimal")
  names(chronic_dx) <- paste0("chronic_dx_", chronic_names)

  age <- value(h, paste0("r", wave, "agey"))
  gender <- value(h, "ragender")
  education <- value(h, "raeducl")
  marital <- value(h, paste0("r", wave, "mstat"))
  hukou <- value(h, paste0("r", wave, "hukou"))
  rural_community <- value(h, paste0("h", wave, "rural"))
  rural_hukou <- value(h, paste0("r", wave, "rural2"))
  self_rated_health <- value(h, paste0("r", wave, "shlta"))
  adl_count <- value(h, paste0("r", wave, "adla_c"))
  iadl_count <- value(h, paste0("r", wave, "iadla"))
  public_insurance <- value(h, paste0("r", wave, "higov"))
  household_income <- value(h, income_var(wave))
  household_wealth <- value(h, wealth_var(wave))
  respondent_weight <- value(h, paste0("r", wave, "wtresp"))
  respondent_weight_adjusted <- value(h, paste0("r", wave, "wtrespb"))
  general_outpatient_visit <- value(h, paste0("r", wave, "doctor1m"))
  general_hospital_stay <- value(h, paste0("r", wave, "hosp1y"))

  chronic_count <- rowSums(chronic_mat, na.rm = TRUE)
  chronic_observed <- rowSums(!is.na(chronic_mat))
  chronic_count <- ifelse(chronic_observed == 0, NA_real_, chronic_count)

  bind_cols(tibble(
    harmonized_charls_id = h$ID,
    wave = wave,
    year = year,
    age_harmonized = age,
    female = ifelse(gender == 2, 1L, ifelse(gender == 1, 0L, NA_integer_)),
    education_level = education,
    education_group = case_when(
      education == 0 ~ "none",
      education == 1 ~ "less_than_lower_secondary",
      education == 2 ~ "upper_secondary_or_vocational",
      education == 3 ~ "tertiary",
      TRUE ~ NA_character_
    ),
    married_or_partnered = case_when(
      marital %in% c(1, 3) ~ 1L,
      marital %in% c(4, 5, 7, 8) ~ 0L,
      TRUE ~ NA_integer_
    ),
    agricultural_hukou = case_when(
      hukou == 1 ~ 1L,
      hukou %in% c(2, 3, 4) ~ 0L,
      TRUE ~ NA_integer_
    ),
    rural_hukou = case_when(
      rural_hukou == 1 ~ 1L,
      rural_hukou == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),
    rural_community = case_when(
      rural_community == 1 ~ 1L,
      rural_community == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),
    self_rated_health_alt = self_rated_health,
    poor_self_rated_health = case_when(
      self_rated_health %in% c(4, 5) ~ 1L,
      self_rated_health %in% c(1, 2, 3) ~ 0L,
      TRUE ~ NA_integer_
    ),
    adl_count = adl_count,
    any_adl_limitation = case_when(
      adl_count > 0 ~ 1L,
      adl_count == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),
    iadl_count = iadl_count,
    any_iadl_limitation = case_when(
      iadl_count > 0 ~ 1L,
      iadl_count == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),
    chronic_count = chronic_count,
    multi_chronic = case_when(
      chronic_count >= 2 ~ 1L,
      chronic_count < 2 ~ 0L,
      TRUE ~ NA_integer_
    ),
    public_insurance = case_when(
      public_insurance == 1 ~ 1L,
      public_insurance == 0 ~ 0L,
      TRUE ~ NA_integer_
    ),
    household_income = household_income,
    household_wealth = household_wealth,
    log_household_income = safe_log1p_positive(household_income),
    log_household_wealth = safe_log1p_positive(household_wealth),
    respondent_weight = respondent_weight,
    respondent_weight_adjusted = respondent_weight_adjusted,
    general_outpatient_visit_last_month = yes01(general_outpatient_visit),
    general_hospital_stay_last_year = yes01(general_hospital_stay),
    covariate_source = "Harmonized CHARLS Version D"
  ), chronic_dx)
}))

covariates <- covariates %>%
  group_by(harmonized_charls_id, wave, year) %>%
  summarize(across(everything(), ~ dplyr::first(.x)), .groups = "drop") %>%
  group_by(year) %>%
  mutate(
    household_income_quartile = case_when(
      is.na(log_household_income) ~ "missing",
      ntile(log_household_income, 4) == 1 ~ "Q1_lowest",
      ntile(log_household_income, 4) == 2 ~ "Q2",
      ntile(log_household_income, 4) == 3 ~ "Q3",
      ntile(log_household_income, 4) == 4 ~ "Q4_highest",
      TRUE ~ "missing"
    ),
    household_wealth_quartile = case_when(
      is.na(log_household_wealth) ~ "missing",
      ntile(log_household_wealth, 4) == 1 ~ "Q1_lowest",
      ntile(log_household_wealth, 4) == 2 ~ "Q2",
      ntile(log_household_wealth, 4) == 3 ~ "Q3",
      ntile(log_household_wealth, 4) == 4 ~ "Q4_highest",
      TRUE ~ "missing"
    )
  ) %>%
  ungroup()

linked <- read_tsv(linked_tsv, show_col_types = FALSE) %>%
  mutate(
    ID = as.character(ID),
    wave = as.integer(wave),
    year = as.integer(year),
    harmonized_charls_id = case_when(
      year == 2011L & nchar(ID) == 11L ~ paste0(substr(ID, 1, 9), "0", substr(ID, 10, 11)),
      TRUE ~ ID
    )
  )

context <- read_tsv(context_tsv, show_col_types = FALSE) %>%
  transmute(
    province_supply_key = province,
    year = as.integer(year),
    value_total_health_institution_beds,
    value_hospital_beds,
    value_comprehensive_hospital_beds,
    value_per_10000_population_total_health_institution_beds,
    value_per_10000_population_hospital_beds,
    value_per_10000_population_comprehensive_hospital_beds,
    gdp_per_capita_current_yuan,
    log_gdp_per_capita_current_yuan,
    urban_population_percent,
    population_age_65_plus_percent,
    context_z_tcm_hospital_beds = z_value_per_10000_population_tcm_hospital_beds,
    context_z_tcm_physicians = z_value_per_10000_population_tcm_practicing_assistant_physicians,
    context_z_comprehensive_hospital_beds = z_value_per_10000_population_comprehensive_hospital_beds,
    context_z_all_hospital_beds = z_value_per_10000_population_hospital_beds,
    context_z_log_gdp_per_capita = z_log_gdp_per_capita_current_yuan,
    context_z_urban_population = z_urban_population_percent,
    context_z_population_age_65_plus = z_population_age_65_plus_percent
  )

analysis <- linked %>%
  left_join(covariates, by = c("harmonized_charls_id", "wave", "year"), suffix = c("", "_cov")) %>%
  left_join(context, by = c("province_supply_key", "year")) %>%
  mutate(
    person_id = harmonized_charls_id,
    age_model = coalesce(age_harmonized, age),
    main_model_age60 = main_core_density_linked_panel &
      age_60plus == 1 &
      !is.na(primary_condition_tcm_any)
  )

py_standardization_vars <- tribble(
  ~source_variable, ~standardized_variable, ~label,
  "value_per_10000_population_tcm_hospital_beds", "z_py_tcm_beds_per_10000", "TCM hospital beds per 10,000 population",
  "value_per_10000_population_tcm_practicing_assistant_physicians", "z_py_tcm_physicians_per_10000", "TCM physicians per 10,000 population",
  "value_per_10000_population_comprehensive_hospital_beds", "z_py_comprehensive_hospital_beds", "Comprehensive-hospital beds per 10,000 population",
  "value_per_10000_population_hospital_beds", "z_py_hospital_beds", "Hospital beds per 10,000 population",
  "log_gdp_per_capita_current_yuan", "z_py_log_gdp", "Log GDP per capita",
  "urban_population_percent", "z_py_urbanization", "Urban population share",
  "population_age_65_plus_percent", "z_py_age_65_plus", "Population aged 65 years or older"
)

standardization_constants <- bind_rows(lapply(seq_len(nrow(py_standardization_vars)), function(i) {
  source_var <- py_standardization_vars$source_variable[[i]]
  values <- analysis %>%
    filter(main_model_age60) %>%
    distinct(province_supply_key, year, value = .data[[source_var]]) %>%
    filter(!is.na(value))
  tibble(
    source_variable = source_var,
    standardized_variable = py_standardization_vars$standardized_variable[[i]],
    label = py_standardization_vars$label[[i]],
    standardization_population = "unique linked CHARLS province-year units among age-eligible observations with a nonmissing primary TCM outcome",
    n_province_years = nrow(values),
    mean = mean(values$value),
    sd = sd(values$value)
  )
}))

for (i in seq_len(nrow(standardization_constants))) {
  source_var <- standardization_constants$source_variable[[i]]
  standardized_var <- standardization_constants$standardized_variable[[i]]
  center <- standardization_constants$mean[[i]]
  spread <- standardization_constants$sd[[i]]
  analysis[[standardized_var]] <- (analysis[[source_var]] - center) / spread
}

qc_vars <- c(
  "age_model", "female", "education_group", "married_or_partnered",
  "rural_community", "agricultural_hukou", "public_insurance",
  "chronic_count", "any_adl_limitation", "poor_self_rated_health",
  "household_income_quartile", "household_wealth_quartile",
  "log_household_income", "log_household_wealth",
  "general_outpatient_visit_last_month", "general_hospital_stay_last_year"
)

qc <- bind_rows(lapply(qc_vars, function(v) {
  analysis %>%
    filter(main_core_density_linked_panel) %>%
    group_by(year) %>%
    summarize(
      variable = v,
      n = n(),
      missing = sum(is.na(.data[[v]])),
      missing_pct = round(100 * missing / n, 2),
      .groups = "drop"
    )
}))

write_tsv(covariates, covariates_tsv)
write_tsv(analysis, analysis_tsv)
write_tsv(qc, qc_tsv)
write_tsv(standardization_constants, standardization_tsv)

cat(sprintf("Wrote %s\n", covariates_tsv))
cat(sprintf("Wrote %s\n", analysis_tsv))
cat(sprintf("Wrote %s\n", qc_tsv))
cat(sprintf("Wrote %s\n", standardization_tsv))
print(analysis %>% count(year, main_model_age60, name = "n_person_waves"))
print(qc %>% filter(variable %in% c("female", "education_group", "chronic_count", "log_household_income")))
