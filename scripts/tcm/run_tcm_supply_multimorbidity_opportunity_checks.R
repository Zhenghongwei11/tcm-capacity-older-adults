#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clubSandwich)
  library(dplyr)
  library(readr)
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
out_tsv <- argval(
  "--output",
  "results/tcm/tcm_supply_multimorbidity_opportunity_checks.tsv"
)

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

cluster_vcov <- function(model, cluster) {
  cluster <- droplevels(as.factor(cluster))
  estimable <- !is.na(coef(model))
  x <- model.matrix(model)[, estimable, drop = FALSE]
  u <- residuals(model)
  x_qr <- qr(x)
  if (x_qr$rank < ncol(x)) stop("Model matrix remains rank deficient after removing aliased terms")
  r_inv <- backsolve(qr.R(x_qr), diag(ncol(x)))
  bread_pivoted <- tcrossprod(r_inv)
  xtx_inv <- matrix(0, ncol(x), ncol(x))
  xtx_inv[x_qr$pivot, x_qr$pivot] <- bread_pivoted
  meat <- matrix(0, ncol(x), ncol(x))
  for (g in levels(cluster)) {
    idx <- cluster == g
    score_g <- crossprod(x[idx, , drop = FALSE], u[idx])
    meat <- meat + tcrossprod(score_g)
  }
  n <- nrow(x)
  k <- ncol(x)
  groups <- nlevels(cluster)
  correction <- (groups / (groups - 1)) * ((n - 1) / (n - k))
  out <- correction * xtx_inv %*% meat %*% xtx_inv
  dimnames(out) <- list(colnames(x), colnames(x))
  out
}

analysis <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  filter(main_model_age60) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    z_tcm_beds = z_py_tcm_beds_per_10000,
    multimorbidity = as.numeric(multi_chronic)
  )

base_covariates <- c(
  "age_model", "I(age_model^2)", "female", "education_group",
  "married_or_partnered", "rural_community", "agricultural_hukou",
  "public_insurance", "any_adl_limitation", "poor_self_rated_health",
  "household_income_quartile"
)

extract_term <- function(data, formula, term, analysis_id, outcome_label,
                         denominator, interaction_scale, condition = NA_character_) {
  needed <- all.vars(formula)
  d <- data %>% filter(if_all(all_of(needed), ~ !is.na(.x))) %>% droplevels()
  fit <- lm(formula, data = d)
  if (!(term %in% names(coef(fit)))) return(tibble())

  vc <- cluster_vcov(fit, d$province_supply_key)
  df <- n_distinct(d$province_supply_key) - 1
  beta <- unname(coef(fit)[[term]])
  se <- sqrt(vc[term, term])
  ci <- beta + c(-1, 1) * qt(0.975, df = df) * se
  pvalue <- 2 * pt(abs(beta / se), df = df, lower.tail = FALSE)

  cr2 <- tryCatch(
    coef_test(fit, vcov = "CR2", cluster = d$province_supply_key, test = "Satterthwaite"),
    error = function(e) NULL
  )
  cr2_se <- cr2_df <- cr2_p <- NA_real_
  if (!is.null(cr2) && term %in% rownames(cr2)) {
    cr2_se <- as.numeric(cr2[term, "SE"])
    cr2_df <- as.numeric(cr2[term, "df_Satt"])
    cr2_p <- as.numeric(cr2[term, "p_Satt"])
  }
  cr2_ci <- if (is.finite(cr2_se) && is.finite(cr2_df)) {
    beta + c(-1, 1) * qt(0.975, df = cr2_df) * cr2_se
  } else {
    c(NA_real_, NA_real_)
  }

  tibble(
    analysis_id = analysis_id,
    condition = condition,
    outcome = outcome_label,
    denominator = denominator,
    interaction_scale = interaction_scale,
    interaction_effect = 100 * beta,
    ci_lower = 100 * ci[[1]],
    ci_upper = 100 * ci[[2]],
    conventional_pvalue = pvalue,
    cr2_ci_lower = 100 * cr2_ci[[1]],
    cr2_ci_upper = 100 * cr2_ci[[2]],
    cr2_pvalue = cr2_p,
    cr2_df = cr2_df,
    n_persons = n_distinct(d$person_id),
    n_person_waves = nrow(d),
    n_events = sum(model.response(model.frame(fit)) == 1),
    n_clusters = n_distinct(d$province_supply_key),
    fixed_effects = "province and survey year",
    covariate_note = "Individual covariates; diagnosis-eligible models additionally adjust for chronic-condition count"
  )
}

continuous_formula <- as.formula(paste(
  "primary_condition_tcm_any ~ z_tcm_beds * chronic_count +",
  paste(c(base_covariates, "province_fe", "wave_fe"), collapse = " + ")
))

binary_adjusted_formula <- as.formula(paste(
  "primary_condition_tcm_any ~ z_tcm_beds * multimorbidity + chronic_count +",
  paste(c(base_covariates, "province_fe", "wave_fe"), collapse = " + ")
))

results <- bind_rows(
  extract_term(
    analysis, continuous_formula, "z_tcm_beds:chronic_count",
    "composite_continuous_burden", "Any condition-specific TCM treatment use",
    "All age-eligible observations", "Additional bed-density association per additional diagnosed condition"
  ),
  extract_term(
    analysis, binary_adjusted_formula, "z_tcm_beds:multimorbidity",
    "composite_multimorbidity_adjusted_for_count", "Any condition-specific TCM treatment use",
    "All age-eligible observations", "Additional bed-density association for multimorbidity after adjustment for chronic-condition count"
  )
)

condition_map <- tribble(
  ~condition, ~diagnosis, ~outcome, ~label,
  "hypertension", "chronic_dx_hypertension", "condition_tcm_hypertension", "TCM treatment for hypertension",
  "diabetes", "chronic_dx_diabetes", "condition_tcm_diabetes", "TCM treatment for diabetes",
  "cancer", "chronic_dx_cancer", "condition_tcm_cancer", "TCM treatment for cancer",
  "chronic_lung", "chronic_dx_chronic_lung", "condition_tcm_chronic_lung", "TCM treatment for chronic lung disease",
  "heart", "chronic_dx_heart", "condition_tcm_heart", "TCM treatment for heart disease",
  "arthritis", "chronic_dx_arthritis", "condition_tcm_arthritis", "TCM treatment for arthritis",
  "dyslipidemia", "chronic_dx_dyslipidemia", "condition_tcm_dyslipidemia", "TCM treatment for dyslipidemia",
  "liver", "chronic_dx_liver", "condition_tcm_liver", "TCM treatment for liver disease",
  "kidney", "chronic_dx_kidney", "condition_tcm_kidney", "TCM treatment for kidney disease",
  "digestive", "chronic_dx_digestive", "condition_tcm_digestive", "TCM treatment for digestive disease",
  "memory", "chronic_dx_memory", "condition_tcm_memory", "TCM treatment for memory-related disease"
)

condition_results <- bind_rows(lapply(seq_len(nrow(condition_map)), function(i) {
  diagnosis <- condition_map$diagnosis[[i]]
  outcome <- condition_map$outcome[[i]]
  d <- analysis %>%
    filter(.data[[diagnosis]] == 1) %>%
    mutate(condition_outcome = as.numeric(.data[[outcome]]))
  if (sum(d$condition_outcome == 1, na.rm = TRUE) < 100) return(tibble())
  formula <- as.formula(paste(
    "condition_outcome ~ z_tcm_beds * multimorbidity + chronic_count +",
    paste(c(base_covariates, "province_fe", "wave_fe"), collapse = " + ")
  ))
  extract_term(
    d, formula, "z_tcm_beds:multimorbidity",
    paste0("diagnosis_eligible_", condition_map$condition[[i]]),
    condition_map$label[[i]], paste0("Respondents diagnosed with ", condition_map$condition[[i]]),
    "Additional bed-density association for multimorbidity within a fixed diagnosis-eligible denominator",
    condition_map$condition[[i]]
  )
})) %>%
  mutate(
    conventional_pvalue_holm = p.adjust(conventional_pvalue, method = "holm"),
    cr2_pvalue_holm = p.adjust(cr2_pvalue, method = "holm")
  )

results <- bind_rows(
  results %>% mutate(conventional_pvalue_holm = NA_real_, cr2_pvalue_holm = NA_real_),
  condition_results
) %>%
  mutate(across(
    c(interaction_effect, ci_lower, ci_upper, conventional_pvalue,
      cr2_ci_lower, cr2_ci_upper, cr2_pvalue, cr2_df,
      conventional_pvalue_holm, cr2_pvalue_holm),
    ~ round(.x, 5)
  ))

write_tsv(results, out_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
print(results %>% select(analysis_id, interaction_effect, ci_lower, ci_upper, conventional_pvalue, cr2_pvalue))
