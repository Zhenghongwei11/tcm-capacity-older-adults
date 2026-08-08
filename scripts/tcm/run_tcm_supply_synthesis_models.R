#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clubSandwich)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

argval <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

analysis_tsv <- argval(
  "--analysis",
  "results/tcm/person_wave_tcm_core_density_analysis.tsv"
)
small_cluster_out <- argval(
  "--small-cluster-output",
  "results/tcm/tcm_supply_small_cluster_inference.tsv"
)
specificity_out <- argval(
  "--specificity-output",
  "results/tcm/tcm_supply_specificity_synthesis.tsv"
)

main_tsv <- argval("--main", "results/tcm/tcm_supply_main_models.tsv")
contextual_tsv <- argval("--contextual", "results/tcm/tcm_supply_contextual_models.tsv")
strict_tsv <- argval("--strict", "results/tcm/tcm_supply_strict_secondary_models.tsv")
falsification_tsv <- argval("--falsification", "results/tcm/tcm_supply_falsification_models.tsv")
bed_physician_tsv <- argval("--bed-physician", "results/tcm/tcm_supply_bed_physician_models.tsv")

required <- c(analysis_tsv, main_tsv, contextual_tsv, strict_tsv, falsification_tsv, bed_physician_tsv)
missing <- required[!file.exists(required)]
if (length(missing) > 0) stop(sprintf("Missing required input(s): %s", paste(missing, collapse = ", ")))

dir.create(dirname(small_cluster_out), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(specificity_out), recursive = TRUE, showWarnings = FALSE)

analysis <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  filter(main_model_age60) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    outcome_primary = as.numeric(primary_condition_tcm_any),
    outcome_strict_tcm_visit = as.numeric(strict_tcm_hospital_visit),
    outcome_general_outpatient = as.numeric(general_outpatient_visit_last_month),
    outcome_general_hospitalization = as.numeric(general_hospital_stay_last_year),
    z_tcm_beds = z_py_tcm_beds_per_10000,
    z_tcm_physicians = z_py_tcm_physicians_per_10000,
    z_comprehensive_beds = z_py_comprehensive_hospital_beds,
    z_log_gdp = z_py_log_gdp,
    z_urbanization = z_py_urbanization
  )

individual_covariates <- c(
  "age_model", "I(age_model^2)", "female", "education_group",
  "married_or_partnered", "rural_community", "agricultural_hukou",
  "public_insurance", "chronic_count", "any_adl_limitation",
  "poor_self_rated_health", "household_income_quartile"
)

cluster_vcov <- function(model, cluster) {
  cluster <- droplevels(as.factor(cluster))
  x <- model.matrix(model)
  u <- residuals(model)
  xtx_inv <- qr.solve(crossprod(x))
  meat <- matrix(0, ncol(x), ncol(x))
  for (g in levels(cluster)) {
    idx <- cluster == g
    score_g <- crossprod(x[idx, , drop = FALSE], u[idx])
    meat <- meat + tcrossprod(score_g)
  }
  n <- nrow(x)
  k <- ncol(x)
  g <- nlevels(cluster)
  correction <- (g / (g - 1)) * ((n - 1) / (n - k))
  out <- correction * xtx_inv %*% meat %*% xtx_inv
  dimnames(out) <- list(colnames(x), colnames(x))
  out
}

fit_terms <- function(model_id, model_label, outcome_var, exposures, contextual_covariates = character()) {
  rhs <- c(exposures, contextual_covariates, individual_covariates, "province_fe", "wave_fe")
  formula <- as.formula(paste(outcome_var, "~", paste(rhs, collapse = " + ")))
  needed <- all.vars(formula)
  d <- analysis %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
  fit <- lm(formula, data = d)
  conventional_vcov <- cluster_vcov(fit, d$province_supply_key)
  cr2 <- coef_test(
    fit,
    vcov = "CR2",
    cluster = d$province_supply_key,
    test = "Satterthwaite"
  )
  conventional_df <- n_distinct(d$province_supply_key) - 1

  bind_rows(lapply(exposures, function(term) {
    beta <- unname(coef(fit)[[term]])
    se <- sqrt(conventional_vcov[term, term])
    ci <- beta + c(-1, 1) * qt(0.975, df = conventional_df) * se
    cr2_row <- cr2[term, ]
    tibble(
      model = model_id,
      model_label = model_label,
      outcome = recode(
        outcome_var,
        outcome_primary = "Disease-specific TCM treatment use",
        outcome_strict_tcm_visit = "TCM hospital visit in the past month",
        outcome_general_outpatient = "Any outpatient or doctor visit",
        outcome_general_hospitalization = "Any hospitalization"
      ),
      exposure = recode(
        term,
        z_tcm_beds = "TCM hospital beds per 10,000 population",
        z_tcm_physicians = "TCM physicians per 10,000 population",
        z_comprehensive_beds = "Comprehensive-hospital beds per 10,000 population"
      ),
      effect_type = "percentage-point difference per 1-SD higher province-year measure",
      effect = 100 * beta,
      ci_lower = 100 * ci[[1]],
      ci_upper = 100 * ci[[2]],
      conventional_pvalue = 2 * pt(abs(beta / se), df = conventional_df, lower.tail = FALSE),
      cr2_se = 100 * as.numeric(cr2_row$SE),
      cr2_df = as.numeric(cr2_row$df),
      cr2_pvalue = as.numeric(cr2_row$p_Satt),
      small_cluster_method = "CR2 cluster-robust variance with Satterthwaite degrees of freedom",
      cluster_level = "province",
      n_clusters = n_distinct(d$province_supply_key),
      n_persons = n_distinct(d$person_id),
      n_person_waves = nrow(d),
      fixed_effects = "province and survey year",
      covariate_set = paste(c(individual_covariates, contextual_covariates), collapse = "; ")
    )
  }))
}

small_cluster <- bind_rows(
  fit_terms(
    "S1_primary",
    "Primary individual-covariate model",
    "outcome_primary",
    "z_tcm_beds"
  ),
  fit_terms(
    "S2_tcm_and_comprehensive_beds",
    "Primary model with comprehensive-hospital bed density",
    "outcome_primary",
    c("z_tcm_beds", "z_comprehensive_beds")
  ),
  fit_terms(
    "S3_tcm_beds_physicians_and_context",
    "Joint TCM resource model with health-system and socioeconomic context",
    "outcome_primary",
    c("z_tcm_beds", "z_tcm_physicians", "z_comprehensive_beds"),
    c("z_log_gdp", "z_urbanization")
  )
)

main <- read_tsv(main_tsv, show_col_types = FALSE)
contextual <- read_tsv(contextual_tsv, show_col_types = FALSE)

wild_lookup <- bind_rows(
  main %>%
    filter(model == "M2_covariate_adjusted_lpm", str_detect(supply_indicator, "beds")) %>%
    transmute(
      model = "S1_primary",
      exposure = "TCM hospital beds per 10,000 population",
      wild_cluster_pvalue
    ),
  contextual %>%
    filter(model == "C1", term == "z_tcm_beds") %>%
    transmute(
      model = "S2_tcm_and_comprehensive_beds",
      exposure = "TCM hospital beds per 10,000 population",
      wild_cluster_pvalue
    )
)

small_cluster <- small_cluster %>%
  left_join(wild_lookup, by = c("model", "exposure")) %>%
  mutate(across(
    c(effect, ci_lower, ci_upper, conventional_pvalue, cr2_se, cr2_df, cr2_pvalue, wild_cluster_pvalue),
    ~ round(.x, 5)
  ))

write_tsv(small_cluster, small_cluster_out)

strict <- read_tsv(strict_tsv, show_col_types = FALSE)
falsification <- read_tsv(falsification_tsv, show_col_types = FALSE)
bed_physician <- read_tsv(bed_physician_tsv, show_col_types = FALSE)

base_primary <- main %>%
  filter(model == "M2_covariate_adjusted_lpm") %>%
  transmute(
    analysis_domain = "Primary outcome",
    outcome = "Disease-specific TCM treatment use",
    exposure = supply_indicator,
    model = "Individual covariates, province fixed effects, and survey-year fixed effects",
    estimate = effect,
    ci_lower,
    ci_upper,
    pvalue,
    small_cluster_pvalue = wild_cluster_pvalue,
    interpretation_note = if_else(
      str_detect(supply_indicator, "beds"),
      "Institutional TCM capacity was positively associated with realized TCM treatment use, but small-cluster inference was less precise.",
      "Physician headcount did not show a parallel positive association."
    )
  )

strict_tcm <- strict %>%
  filter(model == "S2_covariate_adjusted_lpm", str_detect(supply_indicator, "beds|physicians")) %>%
  transmute(
    analysis_domain = "Strict TCM institution use",
    outcome = "TCM hospital visit in the past month",
    exposure = supply_indicator,
    model = "Individual covariates, province fixed effects, and survey-year fixed effects",
    estimate = effect,
    ci_lower,
    ci_upper,
    pvalue,
    small_cluster_pvalue = NA_real_,
    interpretation_note = "The strict institution-use endpoint was rare and showed smaller, imprecise associations."
  )

context_primary <- contextual %>%
  filter(model == "C1") %>%
  transmute(
    analysis_domain = "Health-system context",
    outcome = "Disease-specific TCM treatment use",
    exposure = term_label,
    model = "TCM and comprehensive-hospital bed densities entered jointly",
    estimate = effect,
    ci_lower,
    ci_upper,
    pvalue = cluster_pvalue,
    small_cluster_pvalue = wild_cluster_pvalue,
    interpretation_note = if_else(
      term == "z_tcm_beds",
      "The TCM bed association was not explained by comprehensive-hospital bed density in the same model.",
      "Comprehensive-hospital bed density did not show a positive association with TCM treatment use."
    )
  )

general_use <- falsification %>%
  filter(model %in% c("F2_general_outpatient", "F3_general_inpatient"), term %in% c("z_tcm_beds", "z_comprehensive_beds")) %>%
  transmute(
    analysis_domain = "Broader health-care use",
    outcome = outcome,
    exposure = term_label,
    model = "General utilization model with health-system and socioeconomic context",
    estimate = effect,
    ci_lower,
    ci_upper,
    pvalue = cluster_pvalue,
    small_cluster_pvalue = NA_real_,
    interpretation_note = "Broader utilization outcomes help judge whether TCM capacity is also marking wider care-seeking or service expansion."
  )

joint_resources <- bed_physician %>%
  filter(model == "J2", term %in% c("z_tcm_beds", "z_tcm_physicians", "z_comprehensive_beds")) %>%
  transmute(
    analysis_domain = "Joint resource model",
    outcome = "Disease-specific TCM treatment use",
    exposure = term_label,
    model = "TCM bed, TCM physician, comprehensive-hospital bed, GDP, and urbanization entered jointly",
    estimate = effect,
    ci_lower,
    ci_upper,
    pvalue = cluster_pvalue,
    small_cluster_pvalue = NA_real_,
    interpretation_note = case_when(
      term == "z_tcm_beds" ~ "Bed density retained the positive association in the joint resource model.",
      term == "z_tcm_physicians" ~ "Physician density did not translate into higher realized use in the same period.",
      TRUE ~ "Comprehensive-hospital bed density did not account for the TCM bed association."
    )
  )

specificity <- bind_rows(base_primary, strict_tcm, context_primary, general_use, joint_resources) %>%
  mutate(
    effect_type = "percentage-point difference per 1-SD higher province-year measure",
    n_clusters = 28L,
    across(c(estimate, ci_lower, ci_upper, pvalue, small_cluster_pvalue), ~ round(.x, 5))
  ) %>%
  select(
    analysis_domain, outcome, exposure, model, effect_type, estimate, ci_lower,
    ci_upper, pvalue, small_cluster_pvalue, n_clusters, interpretation_note
  )

write_tsv(specificity, specificity_out)

cat(sprintf("Wrote %s\n", small_cluster_out))
print(small_cluster)
cat(sprintf("Wrote %s\n", specificity_out))
print(specificity)
