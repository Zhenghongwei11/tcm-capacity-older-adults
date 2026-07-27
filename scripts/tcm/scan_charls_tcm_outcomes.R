#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(readr)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

argval <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

input_dir <- argval("--input-dir", "data/raw/charls")
crosswalk_path <- argval("--crosswalk", "data/processed/charls/community_city_crosswalk.tsv")
out_dir <- argval("--out-dir", "results/tcm")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(input_dir)) {
  stop(sprintf("Input directory not found: %s", input_dir))
}

path_first_existing <- function(paths) {
  found <- paths[file.exists(paths)]
  if (length(found) == 0) return(NA_character_)
  found[[1]]
}

module_path <- function(wave_dir, candidates) {
  path_first_existing(file.path(
    input_dir,
    c(
      wave_dir,
      file.path("downloads", wave_dir, "all_modules")
    ),
    rep(candidates, each = 2)
  ))
}

read_dta_safe <- function(path) {
  x <- try(read_dta(path), silent = TRUE)
  if (!inherits(x, "try-error")) return(x)
  x <- try(read_dta(path, encoding = "UTF-8"), silent = TRUE)
  if (!inherits(x, "try-error")) return(x)
  read_dta(path, encoding = "GB18030")
}

var_label <- function(df, var) {
  if (!var %in% names(df)) return("")
  lab <- attr(df[[var]], "label")
  if (is.null(lab)) "" else as.character(lab)
}

selected_multichoice <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(v) & v != 0, 1L, 0L)
}

yes01 <- function(x) {
  v <- suppressWarnings(as.numeric(x))
  ifelse(v == 1, 1L, ifelse(v == 0, 0L, NA_integer_))
}

make_age_group <- function(age) {
  case_when(
    is.na(age) ~ NA_character_,
    age < 45 ~ "<45",
    age < 60 ~ "45-59",
    age < 70 ~ "60-69",
    age < 80 ~ "70-79",
    TRUE ~ "80+"
  )
}

files <- tibble::tribble(
  ~wave, ~year, ~module, ~path,
  1L, 2011L, "health_status", module_path("wave1_2011", c("health_status_and_functioning.dta", "Health_Status_and_Functioning.dta")),
  1L, 2011L, "health_care", module_path("wave1_2011", c("health_care_and_insurance.dta", "Health_Care_and_Insurance.dta")),
  2L, 2013L, "health_status", module_path("wave2_2013", c("health_status_and_functioning.dta", "Health_Status_and_Functioning.dta")),
  2L, 2013L, "health_care", module_path("wave2_2013", c("health_care_and_insurance.dta", "Health_Care_and_Insurance.dta")),
  3L, 2015L, "health_status", module_path("wave3_2015", c("health_status_and_functioning.dta", "Health_Status_and_Functioning.dta")),
  3L, 2015L, "health_care", module_path("wave3_2015", c("health_care_and_insurance.dta", "Health_Care_and_Insurance.dta")),
  4L, 2018L, "health_status", module_path("wave4_2018", c("health_status_and_functioning.dta", "Health_Status_and_Functioning.dta")),
  4L, 2018L, "health_care", module_path("wave4_2018", c("health_care_and_insurance.dta", "Health_Care_and_Insurance.dta"))
)

demo_files <- tibble::tribble(
  ~wave, ~year, ~path,
  1L, 2011L, module_path("wave1_2011", c("demographic_background.dta", "Demographic_Background.dta")),
  2L, 2013L, module_path("wave2_2013", c("demographic_background.dta", "Demographic_Background.dta")),
  3L, 2015L, module_path("wave3_2015", c("demographic_background.dta", "Demographic_Background.dta")),
  4L, 2018L, module_path("wave4_2018", c("demographic_background.dta", "Demographic_Background.dta"))
)

harm_path <- path_first_existing(file.path(
  input_dir,
  c(
    file.path("harmonized", "regular_waves"),
    file.path("downloads", "harmonized", "regular_waves")
  ),
  rep(c("H_CHARLS_D_Data.dta", "H_CHARLS_D_data.dta"), each = 2)
))
if (is.na(harm_path)) {
  stop("Harmonized CHARLS file not found under input_dir/harmonized/regular_waves")
}

missing_files <- files %>% filter(is.na(path))
if (nrow(missing_files) > 0) {
  stop(sprintf(
    "Missing required CHARLS module files: %s",
    paste(sprintf("wave %s %s", missing_files$wave, missing_files$module), collapse = "; ")
  ))
}

read_demographics <- function(wave_i, year_i) {
  demo_path <- demo_files %>% filter(wave == wave_i) %>% pull(path)
  if (is.na(demo_path)) {
    return(tibble(ID = character(), wave = integer(), year = integer(), age_raw = numeric()))
  }
  demo <- read_dta_safe(demo_path) %>% mutate(ID = as.character(ID))
  age_raw <- rep(NA_real_, nrow(demo))
  if ("ba002_1" %in% names(demo)) {
    birth_year <- suppressWarnings(as.numeric(demo$ba002_1))
    age_raw <- ifelse(!is.na(birth_year) & birth_year > 1900 & birth_year <= year_i, year_i - birth_year, NA_real_)
  } else if ("ba004" %in% names(demo)) {
    age_raw <- suppressWarnings(as.numeric(demo$ba004))
  } else if ("ba004_w3_1" %in% names(demo)) {
    birth_year <- suppressWarnings(as.numeric(demo$ba004_w3_1))
    age_raw <- ifelse(!is.na(birth_year) & birth_year > 1900 & birth_year <= year_i, year_i - birth_year, NA_real_)
  }
  tibble(ID = demo$ID, wave = wave_i, year = year_i, age_raw = age_raw) %>%
    group_by(ID, wave, year) %>%
    summarize(age_raw = dplyr::first(age_raw[!is.na(age_raw)]), .groups = "drop")
}

disease_map <- tibble::tribble(
  ~condition_id, ~condition_label, ~w1, ~w2, ~w3, ~w4,
  "dyslipidemia", "Dyslipidemia", "da010_2_s1", "da010_2_s1", "da010_2_s1", "da010_w4_2__s1",
  "chronic_lung", "Chronic lung disease", "da010_5_s1", "da010_5_s1", "da010_5_s1", "da010_w4_5__s1",
  "liver", "Liver disease", "da010_6_s1", "da010_6_s1", "da010_6_s1", "da010_w4_6__s1",
  "heart", "Heart disease", "da010_7_s1", "da010_7_s1", "da010_7_s1", "da010_w4_7__s1",
  "kidney", "Kidney disease", "da010_9_s1", "da010_9_s1", "da010_9_s1", "da010_w4_9__s1",
  "digestive", "Stomach or digestive disease", "da010_10_s1", "da010_10_s1", "da010_10_s1", "da010_w4_10__s1",
  "memory", "Memory-related disease", "da010_12_s1", "da010_12_s1", "da010_12_s1", "da010_w4_12__s1",
  "arthritis", "Arthritis or rheumatism", "da010_13_s1", "da010_13_s1", "da010_13_s1", "da010_w4_13__s1",
  "hypertension", "Hypertension", "da011s1", "da011s1", "da011s1", "da011_w4_s1",
  "diabetes", "Diabetes", "da014s1", "da014s1", "da014s1", "da014_w4_s1",
  "cancer", "Cancer", "da018s1", "da018s1", "da018s1", "da018_w4_s1",
  "stroke_chinese_medicine", "Stroke treated by Chinese medicine", "da019s1", "da019s1", "da019s1", "da019_w4_s1",
  "stroke_acupuncture_moxibustion", "Stroke treated by acupuncture or moxibustion", "da019s4", "da019s4", "da019s4", "da019_w4_s4"
)

raw_healthcare_map <- tibble::tribble(
  ~outcome, ~outcome_tier, ~wave, ~year, ~var,
  "raw_tcm_hospital_visit", "secondary_broader_raw", 1L, 2011L, "ed004s3",
  "raw_traditional_outpatient_treatment", "secondary_broader_raw", 1L, 2011L, "ed021s7",
  "raw_traditional_inpatient_treatment", "secondary_broader_raw", 1L, 2011L, "ee021s8",
  "raw_self_treatment_herbs_medicines", "secondary_broader_raw", 1L, 2011L, "ef001s3",
  "raw_tcm_hospital_visit", "secondary_broader_raw", 2L, 2013L, "ed004s3",
  "raw_traditional_outpatient_treatment", "secondary_broader_raw", 2L, 2013L, "ed021s7",
  "raw_self_treatment_herbs_medicines", "secondary_broader_raw", 2L, 2013L, "ef001s3",
  "raw_tcm_hospital_visit", "secondary_broader_raw", 3L, 2015L, "ed004s3",
  "raw_tcm_hospital_visit_count", "secondary_broader_raw", 3L, 2015L, "ed005_3_",
  "raw_traditional_outpatient_treatment", "secondary_broader_raw", 3L, 2015L, "ed021s7",
  "raw_traditional_inpatient_treatment", "secondary_broader_raw", 3L, 2015L, "ee021s8",
  "raw_self_treatment_herbs_medicines", "secondary_broader_raw", 3L, 2015L, "ef001s3",
  "raw_tcm_hospital_visit", "context_incomplete_2018", 4L, 2018L, "ed004_w4_s3"
)

harmonized_map <- tibble::tribble(
  ~outcome, ~outcome_tier, ~wave, ~year, ~visit_var, ~count_var,
  "strict_tcm_hospital_visit", "strict_secondary", 1L, 2011L, "r1trdmed1m", "r1trdmdtim1m",
  "strict_tcm_hospital_visit", "strict_secondary", 2L, 2013L, "r2trdmed1m", "r2trdmdtim1m",
  "strict_tcm_hospital_visit", "strict_secondary", 3L, 2015L, "r3trdmed1m", "r3trdmdtim1m",
  "strict_tcm_hospital_visit", "strict_secondary", 4L, 2018L, "r4trdmed1m", "r4trdmdtim1m"
)

harm <- read_dta(harm_path)
harm <- harm %>% mutate(ID = as.character(ID))

harm_long <- bind_rows(lapply(seq_len(nrow(harmonized_map)), function(i) {
  row <- harmonized_map[i, ]
  tibble(
    ID = harm$ID,
    wave = row$wave,
    year = row$year,
    age = suppressWarnings(as.numeric(harm[[paste0("r", row$wave, "agey")]])),
    strict_tcm_hospital_visit = yes01(harm[[row$visit_var]]),
    strict_tcm_hospital_visit_count = suppressWarnings(as.numeric(harm[[row$count_var]]))
  )
}))

crosswalk <- NULL
if (file.exists(crosswalk_path)) {
  crosswalk <- read_tsv(crosswalk_path, show_col_types = FALSE) %>%
    mutate(communityID = as.character(communityID)) %>%
    select(any_of(c("communityID", "province", "city", "province_short", "city_short", "urban_nbs")))
}

harmonization_rows <- list()
outcome_rows <- list()

for (wave_i in 1:4) {
  year_i <- c(`1` = 2011L, `2` = 2013L, `3` = 2015L, `4` = 2018L)[[as.character(wave_i)]]
  hs_path <- files %>% filter(wave == wave_i, module == "health_status") %>% pull(path)
  hc_path <- files %>% filter(wave == wave_i, module == "health_care") %>% pull(path)
  hs <- read_dta_safe(hs_path) %>%
    mutate(ID = as.character(ID), communityID = as.character(communityID))
  hc <- read_dta_safe(hc_path) %>%
    mutate(ID = as.character(ID), communityID = as.character(communityID))

  dvars <- disease_map[[paste0("w", wave_i)]]
  disease_outcomes <- tibble(
    ID = hs$ID,
    householdID = if ("householdID" %in% names(hs)) as.character(hs$householdID) else NA_character_,
    communityID = hs$communityID,
    wave = wave_i,
    year = year_i
  )
  for (j in seq_along(dvars)) {
    v <- dvars[[j]]
    condition_id <- disease_map$condition_id[[j]]
    present <- v %in% names(hs)
    disease_outcomes[[paste0("condition_tcm_", condition_id)]] <- if (present) selected_multichoice(hs[[v]]) else NA_integer_
    harmonization_rows[[length(harmonization_rows) + 1]] <- tibble(
      wave = wave_i,
      year = year_i,
      source_file = hs_path,
      module = "health_status",
      raw_variable = v,
      label = if (present) var_label(hs, v) else "",
      selected_value_coding = "multi-select item; non-missing/non-zero means selected",
      denominator_rule = "All respondents with health-status module record; condition-specific denominator to be checked using diagnosis variables",
      missingness_rule = "For multi-select items, missing usually means not selected or not applicable; composite treats not selected as 0 within module records",
      outcome = paste0("condition_tcm_", condition_id),
      outcome_tier = "primary_candidate"
    )
  }

  condition_cols <- names(disease_outcomes)[startsWith(names(disease_outcomes), "condition_tcm_")]
  disease_outcomes <- disease_outcomes %>%
    mutate(primary_condition_tcm_any = as.integer(rowSums(across(all_of(condition_cols)), na.rm = TRUE) > 0))

  hc_map_i <- raw_healthcare_map %>% filter(wave == wave_i)
  hc_out <- tibble(
    ID = hc$ID,
    householdID = if ("householdID" %in% names(hc)) as.character(hc$householdID) else NA_character_,
    communityID = hc$communityID,
    wave = wave_i,
    year = year_i
  )
  for (i in seq_len(nrow(hc_map_i))) {
    v <- hc_map_i$var[[i]]
    out <- hc_map_i$outcome[[i]]
    present <- v %in% names(hc)
    if (present && grepl("count", out)) {
      hc_out[[out]] <- suppressWarnings(as.numeric(hc[[v]]))
    } else {
      hc_out[[out]] <- if (present) selected_multichoice(hc[[v]]) else NA_integer_
    }
    harmonization_rows[[length(harmonization_rows) + 1]] <- tibble(
      wave = wave_i,
      year = year_i,
      source_file = hc_path,
      module = "health_care",
      raw_variable = v,
      label = if (present) var_label(hc, v) else "",
      selected_value_coding = ifelse(grepl("count", out), "numeric visit count where available", "multi-select item; non-missing/non-zero means selected"),
      denominator_rule = ifelse(year_i == 2018L, "2018 health-care module has incomplete comparable broader TCM components", "All respondents with health-care module record"),
      missingness_rule = "For multi-select items, missing usually means not selected or not applicable; counts retain numeric missingness",
      outcome = out,
      outcome_tier = hc_map_i$outcome_tier[[i]]
    )
  }

  hc_component_cols <- intersect(
    c("raw_tcm_hospital_visit", "raw_traditional_outpatient_treatment", "raw_traditional_inpatient_treatment", "raw_self_treatment_herbs_medicines"),
    names(hc_out)
  )
  if (year_i %in% c(2011L, 2013L, 2015L)) {
    hc_out <- hc_out %>%
      mutate(broader_tcm_use_2011_2015 = as.integer(rowSums(across(all_of(hc_component_cols)), na.rm = TRUE) > 0))
  } else {
    hc_out <- hc_out %>% mutate(broader_tcm_use_2011_2015 = NA_integer_)
  }

  person_wave <- disease_outcomes %>%
    full_join(hc_out, by = c("ID", "householdID", "communityID", "wave", "year")) %>%
    left_join(harm_long %>% filter(wave == wave_i), by = c("ID", "wave", "year")) %>%
    left_join(read_demographics(wave_i, year_i), by = c("ID", "wave", "year")) %>%
    mutate(
      age = coalesce(age, age_raw),
      strict_tcm_hospital_visit = coalesce(strict_tcm_hospital_visit, raw_tcm_hospital_visit)
    ) %>%
    select(-age_raw)

  if (!is.null(crosswalk)) {
    person_wave <- person_wave %>% left_join(crosswalk, by = "communityID")
  }

  outcome_rows[[length(outcome_rows) + 1]] <- person_wave
}

harmonization <- bind_rows(harmonization_rows)
for (i in seq_len(nrow(harmonized_map))) {
  row <- harmonized_map[i, ]
  harmonization <- bind_rows(
    harmonization,
    tibble(
      wave = row$wave,
      year = row$year,
      source_file = harm_path,
      module = "harmonized_charls",
      raw_variable = row$visit_var,
      label = var_label(harm, row$visit_var),
      selected_value_coding = "0/1 harmonized indicator",
      denominator_rule = "Respondents with valid Harmonized CHARLS health-care response for wave",
      missingness_rule = "Keep missing as missing",
      outcome = "strict_tcm_hospital_visit",
      outcome_tier = "strict_secondary"
    ),
    tibble(
      wave = row$wave,
      year = row$year,
      source_file = harm_path,
      module = "harmonized_charls",
      raw_variable = row$count_var,
      label = var_label(harm, row$count_var),
      selected_value_coding = "numeric count",
      denominator_rule = "Respondents with valid Harmonized CHARLS health-care response for wave",
      missingness_rule = "Keep missing as missing",
      outcome = "strict_tcm_hospital_visit_count",
      outcome_tier = "strict_secondary"
    )
  )
}

person_wave <- bind_rows(outcome_rows) %>%
  mutate(
    age_group = make_age_group(age),
    age_60plus = ifelse(!is.na(age) & age >= 60, 1L, ifelse(!is.na(age), 0L, NA_integer_))
  )

outcome_vars <- c(
  "primary_condition_tcm_any",
  "strict_tcm_hospital_visit",
  "strict_tcm_hospital_visit_count",
  "broader_tcm_use_2011_2015",
  "raw_tcm_hospital_visit",
  "raw_traditional_outpatient_treatment",
  "raw_traditional_inpatient_treatment",
  "raw_self_treatment_herbs_medicines"
)
outcome_vars <- intersect(outcome_vars, names(person_wave))

summarize_events <- function(data, group_vars) {
  data %>%
    select(all_of(group_vars), all_of(outcome_vars)) %>%
    pivot_longer(cols = all_of(outcome_vars), names_to = "outcome", values_to = "value") %>%
    group_by(across(all_of(group_vars)), outcome) %>%
    summarize(
      n_rows = n(),
      n_nonmissing = sum(!is.na(value)),
      n_events = sum(value > 0, na.rm = TRUE),
      event_rate = ifelse(n_nonmissing > 0, n_events / n_nonmissing, NA_real_),
      n_missing = sum(is.na(value)),
      missing_rate = n_missing / n_rows,
      .groups = "drop"
    )
}

harmonization <- harmonization %>%
  arrange(year, module, outcome, raw_variable)

front_cols <- intersect(
  c("ID", "householdID", "communityID", "province", "city", "province_short", "city_short", "wave", "year", "age", "age_group", "age_60plus"),
  names(person_wave)
)
person_wave <- person_wave %>%
  arrange(year, ID) %>%
  relocate(all_of(front_cols))

write_tsv(harmonization, file.path(out_dir, "tcm_outcome_harmonization.tsv"))
write_tsv(person_wave, file.path(out_dir, "person_wave_tcm_outcomes.tsv"))
write_tsv(summarize_events(person_wave, c("wave", "year")), file.path(out_dir, "tcm_event_counts_by_wave.tsv"))
write_tsv(summarize_events(person_wave, c("age_group")), file.path(out_dir, "tcm_event_counts_by_age_group.tsv"))
if ("province" %in% names(person_wave)) {
  write_tsv(summarize_events(person_wave, c("province", "wave", "year")), file.path(out_dir, "tcm_event_counts_by_province.tsv"))
}
write_tsv(
  summarize_events(person_wave, c("wave", "year")) %>% select(wave, year, outcome, n_rows, n_missing, missing_rate),
  file.path(out_dir, "tcm_missingness_by_wave.tsv")
)

cat(sprintf("Wrote %s\n", file.path(out_dir, "tcm_outcome_harmonization.tsv")))
cat(sprintf("Wrote %s\n", file.path(out_dir, "person_wave_tcm_outcomes.tsv")))
cat(sprintf("Wrote event-count summaries under %s\n", out_dir))
