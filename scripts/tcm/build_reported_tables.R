#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

out_dir <- "tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits, big.mark = ","))
}

fmt_int <- function(x) {
  ifelse(is.na(x), "", formatC(round(x), format = "d", big.mark = ","))
}

fmt_ci <- function(effect, lo, hi, digits = 2) {
  paste0(fmt_num(effect, digits), " (", fmt_num(lo, digits), " to ", fmt_num(hi, digits), ")")
}

fmt_p <- function(x) {
  ifelse(is.na(x), "", ifelse(x < 0.001, "<0.001", fmt_num(x, 3)))
}

fmt_p4 <- function(x) {
  ifelse(is.na(x), "", ifelse(x < 0.0001, "<0.0001", fmt_num(x, 4)))
}

fmt_p_interaction <- function(x) {
  ifelse(is.na(x), "", ifelse(x < 0.01, fmt_p4(x), fmt_p(x)))
}

label_component <- function(x) {
  recode(
    x,
    dyslipidemia = "Dyslipidemia",
    chronic_lung = "Chronic lung disease",
    liver = "Liver disease",
    heart = "Heart disease",
    kidney = "Kidney disease",
    digestive = "Digestive disease",
    memory = "Memory-related disease",
    arthritis = "Arthritis or rheumatism",
    hypertension = "Hypertension",
    diabetes = "Diabetes",
    cancer = "Cancer",
    stroke_chinese_medicine = "Stroke, Chinese medicine",
    stroke_acupuncture_moxibustion = "Stroke, acupuncture or moxibustion",
    .default = str_to_sentence(str_replace_all(x, "_", " "))
  )
}

md_table <- function(df) {
  names(df) <- str_replace_all(names(df), "_", " ")
  header <- paste(names(df), collapse = " | ")
  sep <- paste(rep("---", ncol(df)), collapse = " | ")
  body <- apply(df, 1, function(row) paste(row, collapse = " | "))
  paste(c(paste0("| ", header, " |"), paste0("| ", sep, " |"), paste0("| ", body, " |")), collapse = "\n")
}

analysis <- read_tsv("results/tcm/person_wave_tcm_core_density_analysis.tsv", show_col_types = FALSE)
supply <- read_tsv("results/tcm/province_year_tcm_core_density.tsv", show_col_types = FALSE)
main <- read_tsv("results/tcm/tcm_supply_main_models.tsv", show_col_types = FALSE)
strict <- read_tsv("results/tcm/tcm_supply_strict_secondary_models.tsv", show_col_types = FALSE)
robust <- read_tsv("results/tcm/tcm_supply_robustness.tsv", show_col_types = FALSE)
interactions <- read_tsv("results/tcm/tcm_supply_equity_interactions.tsv", show_col_types = FALSE)
components <- read_tsv("results/tcm/tcm_supply_primary_component_models.tsv", show_col_types = FALSE)
cluster <- read_tsv("results/tcm/tcm_supply_cluster_sensitivity.tsv", show_col_types = FALSE)
contextual <- read_tsv("results/tcm/tcm_supply_contextual_models.tsv", show_col_types = FALSE)
contextual_diagnostics <- read_tsv("results/tcm/tcm_supply_contextual_diagnostics.tsv", show_col_types = FALSE)
falsification <- read_tsv("results/tcm/tcm_supply_falsification_models.tsv", show_col_types = FALSE)
bed_physician <- read_tsv("results/tcm/tcm_supply_bed_physician_models.tsv", show_col_types = FALSE)
weighted <- read_tsv("results/tcm/tcm_supply_weighted_attrition_models.tsv", show_col_types = FALSE)
identification <- read_tsv("results/tcm/tcm_supply_longitudinal_identification_models.tsv", show_col_types = FALSE)
composite_stability <- read_tsv("results/tcm/tcm_supply_composite_stability.tsv", show_col_types = FALSE)
flow <- read_tsv("results/tcm/charls_participant_flow.tsv", show_col_types = FALSE)
retention <- read_tsv("results/tcm/charls_retention_by_wave.tsv", show_col_types = FALSE)

age60 <- analysis %>% filter(main_model_age60)

table1 <- age60 %>%
  group_by(year) %>%
  summarise(
    observations = n(),
    respondents = n_distinct(person_id),
    `age, mean` = mean(age_model, na.rm = TRUE),
    `female, %` = 100 * mean(female == 1, na.rm = TRUE),
    `rural residence, %` = 100 * mean(rural_community == 1, na.rm = TRUE),
    `less than lower secondary education, %` = 100 * mean(education_group == "less_than_lower_secondary", na.rm = TRUE),
    `public insurance, %` = 100 * mean(public_insurance == 1, na.rm = TRUE),
    `chronic conditions, mean` = mean(chronic_count, na.rm = TRUE),
    `multimorbidity, %` = 100 * mean(multi_chronic == 1, na.rm = TRUE),
    `ADL limitation, %` = 100 * mean(any_adl_limitation == 1, na.rm = TRUE),
    `primary TCM use, %` = 100 * mean(primary_condition_tcm_any == 1, na.rm = TRUE),
    `TCM hospital visit, %` = 100 * mean(strict_tcm_hospital_visit == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    year = as.character(year),
    observations = fmt_int(observations),
    respondents = fmt_int(respondents),
    across(-c(year, observations, respondents), ~ fmt_num(.x, 1))
  )

overall <- age60 %>%
  summarise(
    year = "Overall",
    observations = fmt_int(n()),
    respondents = fmt_int(n_distinct(person_id)),
    `age, mean` = fmt_num(mean(age_model, na.rm = TRUE), 1),
    `female, %` = fmt_num(100 * mean(female == 1, na.rm = TRUE), 1),
    `rural residence, %` = fmt_num(100 * mean(rural_community == 1, na.rm = TRUE), 1),
    `less than lower secondary education, %` = fmt_num(100 * mean(education_group == "less_than_lower_secondary", na.rm = TRUE), 1),
    `public insurance, %` = fmt_num(100 * mean(public_insurance == 1, na.rm = TRUE), 1),
    `chronic conditions, mean` = fmt_num(mean(chronic_count, na.rm = TRUE), 1),
    `multimorbidity, %` = fmt_num(100 * mean(multi_chronic == 1, na.rm = TRUE), 1),
    `ADL limitation, %` = fmt_num(100 * mean(any_adl_limitation == 1, na.rm = TRUE), 1),
    `primary TCM use, %` = fmt_num(100 * mean(primary_condition_tcm_any == 1, na.rm = TRUE), 1),
    `TCM hospital visit, %` = fmt_num(100 * mean(strict_tcm_hospital_visit == 1, na.rm = TRUE), 1)
  )

table1 <- bind_rows(table1, overall)

table2 <- supply %>%
  group_by(year) %>%
  summarise(
    provinces = n(),
    `TCM hospital beds` = sum(value_tcm_hospital_beds, na.rm = TRUE),
    `beds per 10,000, mean` = mean(value_per_10000_population_tcm_hospital_beds, na.rm = TRUE),
    `beds per 10,000, range` = paste0(
      fmt_num(min(value_per_10000_population_tcm_hospital_beds, na.rm = TRUE), 2),
      " to ",
      fmt_num(max(value_per_10000_population_tcm_hospital_beds, na.rm = TRUE), 2)
    ),
    `TCM physicians` = sum(value_tcm_practicing_assistant_physicians, na.rm = TRUE),
    `physicians per 10,000, mean` = mean(value_per_10000_population_tcm_practicing_assistant_physicians, na.rm = TRUE),
    `physicians per 10,000, range` = paste0(
      fmt_num(min(value_per_10000_population_tcm_practicing_assistant_physicians, na.rm = TRUE), 2),
      " to ",
      fmt_num(max(value_per_10000_population_tcm_practicing_assistant_physicians, na.rm = TRUE), 2)
    ),
    .groups = "drop"
  ) %>%
  mutate(
    year = as.character(year),
    provinces = fmt_int(provinces),
    `TCM hospital beds` = fmt_int(`TCM hospital beds`),
    `beds per 10,000, mean` = fmt_num(`beds per 10,000, mean`, 2),
    `TCM physicians` = fmt_int(`TCM physicians`),
    `physicians per 10,000, mean` = fmt_num(`physicians per 10,000, mean`, 2)
  )

model_rows <- bind_rows(
  main %>% filter(model == "M2_covariate_adjusted_lpm"),
  strict %>%
    filter(model == "S2_covariate_adjusted_lpm") %>%
    filter(str_detect(supply_indicator, "beds|physicians")),
  robust %>%
    filter(check_type %in% c("age_threshold", "covariate_set", "lag_structure", "alternative_outcome")) %>%
    filter(str_detect(supply_indicator, "bed")) %>%
    filter(outcome != "TCM hospital visit in the past month")
)

table3 <- model_rows %>%
  transmute(
    outcome = outcome,
    exposure = supply_indicator,
    model = case_when(
      model == "M2_covariate_adjusted_lpm" ~ "Primary, fully adjusted",
      model == "S2_covariate_adjusted_lpm" ~ "Strict visit, fully adjusted",
      check_type == "age_threshold" ~ "Adults aged 45 years or older",
      check_type == "covariate_set" ~ "Without household income",
      check_type == "lag_structure" ~ "Lagged supply",
      check_type == "alternative_outcome" ~ "Alternative outcome",
      TRUE ~ model
    ),
    estimate = fmt_ci(effect, ci_lower, ci_upper, 2),
    `clustered P value` = fmt_p(pvalue),
    `wild-cluster P value` = fmt_p(wild_cluster_pvalue),
    observations = fmt_int(n_person_waves),
    provinces = fmt_int(n_clusters)
  )

table4 <- interactions %>%
  transmute(
    stratum = case_when(
      stratum == "rural_stratum" ~ "Rural residence",
      stratum == "education_stratum" ~ "Lower education",
      stratum == "income_stratum" & str_detect(contrast, "lower_income") ~ "Lower household income",
      stratum == "income_stratum" & str_detect(contrast, "income_missing") ~ "Income missing",
      stratum == "multichronic_stratum" ~ "Multimorbidity",
      stratum == "adl_stratum" ~ "ADL limitation",
      TRUE ~ stratum
    ),
    reference = case_when(
      reference_level == "urban" ~ "Urban residence",
      reference_level == "higher_education" ~ "Higher education",
      reference_level == "higher_income" ~ "Higher household income",
      reference_level == "no_multimorbidity" ~ "No multimorbidity",
      reference_level == "no_adl_limitation" ~ "No ADL limitation",
      TRUE ~ reference_level
    ),
    `interaction estimate` = fmt_ci(interaction_effect, ci_lower, ci_upper, 2),
    `unadjusted P value` = fmt_p_interaction(pvalue),
    `Holm-adjusted P value` = fmt_p_interaction(pvalue_holm),
    observations = fmt_int(n_person_waves),
    provinces = fmt_int(n_clusters)
  )

supp_component <- components %>%
  transmute(
    component = label_component(str_remove(outcome, "^condition_tcm_")),
    events = fmt_int(n_events),
    `event prevalence, %` = fmt_num(event_percent, 2),
    estimate = fmt_ci(effect, ci_lower, ci_upper, 2),
    `P value` = fmt_p(pvalue)
  )

supp_cluster <- cluster %>%
  filter(check_type != "leave_one_province_out") %>%
  transmute(
    check = case_when(
      check_type == "main_cluster_robust" ~ "Main province-clustered estimate",
      check_type == "leave_one_province_out_summary" ~ "Leave-one-province-out range",
      check_type == "province_cluster_bootstrap" ~ "Province resampling",
      TRUE ~ check_type
    ),
    detail = detail,
    estimate = case_when(
      check_type == "leave_one_province_out_summary" ~ fmt_num(effect, 2),
      TRUE ~ fmt_ci(effect, ci_lower, ci_upper, 2)
    ),
    `interval or range` = case_when(
      check_type == "leave_one_province_out_summary" ~ paste0(fmt_num(ci_lower, 2), " to ", fmt_num(ci_upper, 2), " (min-max)"),
      TRUE ~ paste0(fmt_num(ci_lower, 2), " to ", fmt_num(ci_upper, 2), " (95% CI)")
    ),
    observations = fmt_int(n_person_waves),
    provinces = fmt_int(n_clusters)
  )

supp_cluster$detail <- recode(
  supp_cluster$detail,
  `499 weighted cluster-bootstrap replicates` = "499 weighted province-level resamples",
  .default = supp_cluster$detail
)

contextual_bed <- contextual %>%
  filter(term == "z_tcm_beds") %>%
  left_join(
    contextual_diagnostics %>%
      filter(term == "z_tcm_beds") %>%
      select(model, two_way_fe_vif),
    by = "model"
  )

supp_contextual <- contextual_bed %>%
  filter(model %in% c("C0", "C1", "C2", "C3", "C4")) %>%
  transmute(
    model = recode(
      model,
      C0 = "Main individual-covariate model",
      C1 = "Add comprehensive-hospital bed density",
      C2 = "Add provincial GDP per capita and urbanization",
      C3 = "Comprehensive-hospital beds, GDP per capita, and urbanization",
      C4 = "Add province-level population aged 65 years or older"
    ),
    `TCM bed estimate` = fmt_ci(effect, ci_lower, ci_upper, 2),
    `conventional P value` = fmt_p(cluster_pvalue),
    `wild-cluster P value` = fmt_p(wild_cluster_pvalue),
    `TCM-bed VIF` = fmt_num(two_way_fe_vif, 2)
  )

supp_falsification <- bind_rows(
  falsification %>%
    filter(term == "z_future_tcm_beds") %>%
    transmute(
      check = "Future exposure for current TCM use",
      exposure = "Next-wave TCM hospital beds",
      estimate = fmt_ci(effect, ci_lower, ci_upper, 2),
      `P value` = fmt_p(cluster_pvalue),
      observations = fmt_int(n_person_waves)
    ),
  falsification %>%
    filter(term == "z_tcm_beds", model %in% c("F2_general_outpatient", "F3_general_inpatient")) %>%
    transmute(
      check = if_else(model == "F2_general_outpatient", "General outpatient use", "General hospitalization"),
      exposure = "Concurrent TCM hospital beds",
      estimate = fmt_ci(effect, ci_lower, ci_upper, 2),
      `P value` = fmt_p(cluster_pvalue),
      observations = fmt_int(n_person_waves)
    ),
  bed_physician %>%
    filter(model %in% c("J2", "J4"), str_detect(term, "tcm_beds|tcm_physicians")) %>%
    transmute(
      check = if_else(model == "J2", "Joint concurrent resource model", "Joint lagged resource model"),
      exposure = case_when(
        term == "z_tcm_beds" ~ "TCM hospital beds",
        term == "z_tcm_physicians" ~ "TCM physicians",
        term == "z_lag_tcm_beds" ~ "Lagged TCM hospital beds",
        term == "z_lag_tcm_physicians" ~ "Lagged TCM physicians"
      ),
      estimate = fmt_ci(effect, ci_lower, ci_upper, 2),
      `P value` = fmt_p(cluster_pvalue),
      observations = fmt_int(n_person_waves)
    )
)

supp_weight_identification <- bind_rows(
  weighted %>%
    transmute(
      model = recode(
        model,
        W0_unweighted_reproduction = "Unweighted primary model",
        W1_adjusted_cross_sectional_weight = "Adjusted respondent-weight model",
        W2_trimmed_adjusted_cross_sectional_weight = "Adjusted respondent weight, 1st-99th percentile truncation"
      ),
      estimate = fmt_ci(effect, ci_lower, ci_upper, 2),
      `P value` = fmt_p(pvalue),
      respondents = fmt_int(n_persons),
      observations = fmt_int(n_person_waves),
      `effective observations` = fmt_int(effective_person_waves)
    ),
  identification %>%
    filter(model != "I0_primary_reproduction") %>%
    transmute(
      model = model_label,
      estimate = fmt_ci(effect, ci_lower, ci_upper, 2),
      `P value` = fmt_p(pvalue),
      respondents = fmt_int(n_persons),
      observations = fmt_int(n_person_waves),
      `effective observations` = ""
    )
)

supp_composite_stability <- composite_stability %>%
  transmute(
    check = str_replace_all(scenario, "_", " ") %>% str_to_sentence(),
    `outcome prevalence, %` = fmt_num(outcome_prevalence_percent, 2),
    estimate = fmt_ci(effect, ci_lower, ci_upper, 2),
    `P value` = fmt_p(pvalue),
    observations = fmt_int(n_person_waves)
  )

supp_flow <- flow %>%
  transmute(stage, respondents = fmt_int(n_persons), observations = fmt_int(n_person_waves))

supp_retention <- retention %>%
  transmute(
    interval = paste(index_year, "to", next_year),
    `index respondents` = fmt_int(n_index),
    `observed next wave` = fmt_int(n_observed_next_wave),
    `retention, %` = fmt_num(retention_percent, 1)
  )

tables_md <- c(
  "# Manuscript Tables",
  "",
  "## Table 1. Characteristics of CHARLS respondents aged 60 years or older in the linked descriptive sample",
  "",
  md_table(table1),
  "",
  "Note: Values are person-wave summaries unless labelled as respondents. Percentages are unweighted. ADL denotes activities of daily living; TCM denotes traditional Chinese medicine. Fully adjusted models used the subset with complete model covariates.",
  "",
  "## Table 2. Provincial TCM hospital bed and physician supply, 2011-2018",
  "",
  md_table(table2),
  "",
  "Note: Provincial supply indicators were taken from official National Administration of Traditional Chinese Medicine statistical extracts. Physician counts refer to TCM physicians.",
  "",
  "## Table 3. Main and sensitivity associations between provincial TCM supply and realized TCM use",
  "",
  md_table(table3),
  "",
  "Note: Estimates are percentage-point differences per 1-SD higher province-year supply. Conventional P values use province-clustered standard errors. The wild-cluster P value for the prespecified primary bed-density model uses a null-imposed bootstrap with Webb six-point weights and 9,999 repetitions. Models include province and survey-year fixed effects. Fully adjusted models include age, age squared, sex, education, marital status, rural residence, agricultural hukou, public insurance, chronic disease count, ADL limitation, poor self-rated health, and household income quartile unless otherwise stated.",
  "",
  "## Table 4. Formal interaction contrasts for the association between TCM bed density and primary TCM treatment use",
  "",
  md_table(table4),
  "",
  "Note: Interaction estimates are additional percentage-point differences per 1-SD higher TCM hospital bed density compared with the reference group. Holm adjustment treats the six reported coefficients from five prespecified equity domains as one family. All models include province and survey-year fixed effects and the main covariate set, omitting the stratifying covariate where appropriate.",
  "",
  "## Supplementary Table 1. Component checks for the primary TCM treatment composite",
  "",
  md_table(supp_component),
  "",
  "Note: Estimates are percentage-point differences per 1-SD higher provincial TCM hospital bed density. Components with few events should be interpreted as descriptive checks of the composite outcome.",
  "",
  "## Supplementary Table 2. Small-cluster and province-influence sensitivity checks",
  "",
  md_table(supp_cluster),
  "",
  "Note: Estimates are percentage-point differences per 1-SD higher provincial TCM hospital bed density for the primary outcome. Parenthetical intervals for the main and province-resampling rows are 95% confidence intervals. The leave-one-province-out row reports the median estimate and the minimum-to-maximum range across 28 omitted-province models. Province resampling used 499 weighted province-level resamples.",
  "",
  "## Supplementary Table 3. General health-system and province-context adjustment",
  "",
  md_table(supp_contextual),
  "",
  "Note: Estimates are percentage-point differences per 1-SD higher provincial TCM hospital bed density. Wild-cluster P values use a null-imposed bootstrap with Webb six-point weights and 9,999 repetitions. Comprehensive-hospital beds are a co-exposure, not a placebo.",
  "",
  "## Supplementary Table 4. Temporal, utilization-specificity, and joint-resource checks",
  "",
  md_table(supp_falsification),
  "",
  "Note: General outpatient and hospitalization outcomes can include TCM care and are not strict non-TCM negative controls. Joint concurrent and lagged resource models include comprehensive-hospital bed density, log provincial GDP per capita, and urbanization.",
  "",
  "## Supplementary Table 5. Weighting and longitudinal identification sensitivity analyses",
  "",
  md_table(supp_weight_identification),
  "",
  "Note: Weighted models use wave-specific individual cross-sectional weights with household and individual non-response adjustment, normalized within wave. They are estimand-specific sensitivities rather than replacements for the unweighted association model. Individual fixed-effects and province-trend models use province-clustered uncertainty.",
  "",
  "## Supplementary Table 6. Leave-component-out stability of the primary composite outcome",
  "",
  md_table(supp_composite_stability),
  "",
  "Note: Each model reconstructs the primary outcome after omitting the named component or clinically coherent component family and retains the main fixed-effects and covariate specification.",
  "",
  "## Supplementary Table 7. Analytic participant flow",
  "",
  md_table(supp_flow),
  "",
  "## Supplementary Table 8. Wave-to-wave observation among age-eligible respondents",
  "",
  md_table(supp_retention),
  "",
  "Note: Subsequent observation combines continued eligibility and survey participation; it should not be interpreted as pure non-response because death and other reasons for non-observation are not separated in this table."
)

write_lines(tables_md, file.path(out_dir, "tables.md"))
cat("Wrote tables/tables.md\n")
