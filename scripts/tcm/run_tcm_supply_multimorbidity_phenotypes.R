#!/usr/bin/env Rscript

suppressPackageStartupMessages({
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
  "results/tcm/tcm_supply_multimorbidity_phenotype_interactions.tsv"
)

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

required_dx <- c(
  "chronic_dx_hypertension", "chronic_dx_diabetes", "chronic_dx_dyslipidemia",
  "chronic_dx_heart", "chronic_dx_stroke", "chronic_dx_arthritis",
  "chronic_dx_chronic_lung", "chronic_dx_digestive", "chronic_dx_liver",
  "chronic_dx_kidney", "chronic_dx_memory"
)

raw <- read_tsv(analysis_tsv, show_col_types = FALSE)
missing_dx <- setdiff(required_dx, names(raw))
if (length(missing_dx) > 0) {
  stop(sprintf(
    "The analysis table is missing chronic diagnosis variables required for phenotype analysis: %s. Re-run scripts/tcm/build_charls_covariates.R first.",
    paste(missing_dx, collapse = ", ")
  ))
}

count_present <- function(df) {
  rowSums(df == 1, na.rm = TRUE)
}

count_observed <- function(df) {
  rowSums(!is.na(df))
}

analysis <- raw %>%
  filter(main_model_age60) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    z_tcm_beds_per_10000 = z_py_tcm_beds_per_10000,
    outcome_primary = as.numeric(primary_condition_tcm_any)
  )

cardio_vars <- c(
  "chronic_dx_hypertension", "chronic_dx_diabetes", "chronic_dx_dyslipidemia",
  "chronic_dx_heart", "chronic_dx_stroke"
)
visceral_vars <- c("chronic_dx_digestive", "chronic_dx_liver", "chronic_dx_kidney")
neuro_vars <- c("chronic_dx_memory", "chronic_dx_stroke")

analysis <- analysis %>%
  mutate(
    cardio_count = count_present(pick(all_of(cardio_vars))),
    cardio_observed = count_observed(pick(all_of(cardio_vars))),
    visceral_count = count_present(pick(all_of(visceral_vars))),
    visceral_observed = count_observed(pick(all_of(visceral_vars))),
    neuro_count = count_present(pick(all_of(neuro_vars))),
    neuro_observed = count_observed(pick(all_of(neuro_vars))),
    phenotype_cardiometabolic = case_when(
      cardio_observed == 0 ~ NA_integer_,
      cardio_count >= 2 ~ 1L,
      TRUE ~ 0L
    ),
    phenotype_arthritis_multimorbidity = case_when(
      is.na(chronic_dx_arthritis) | is.na(multi_chronic) ~ NA_integer_,
      chronic_dx_arthritis == 1 & multi_chronic == 1 ~ 1L,
      TRUE ~ 0L
    ),
    phenotype_respiratory_multimorbidity = case_when(
      is.na(chronic_dx_chronic_lung) | is.na(multi_chronic) ~ NA_integer_,
      chronic_dx_chronic_lung == 1 & multi_chronic == 1 ~ 1L,
      TRUE ~ 0L
    ),
    phenotype_digestive_liver_kidney = case_when(
      visceral_observed == 0 ~ NA_integer_,
      visceral_count >= 1 & multi_chronic == 1 ~ 1L,
      TRUE ~ 0L
    ),
    phenotype_neurocognitive_stroke = case_when(
      neuro_observed == 0 ~ NA_integer_,
      neuro_count >= 1 & multi_chronic == 1 ~ 1L,
      TRUE ~ 0L
    )
  )

cluster_vcov <- function(model, cluster) {
  ok <- !is.na(cluster)
  cluster <- as.factor(cluster[ok])
  x <- model.matrix(model)[ok, , drop = FALSE]
  u <- residuals(model)[ok]
  xtx_inv <- solve(crossprod(x))
  meat <- matrix(0, ncol(x), ncol(x))
  for (g in levels(cluster)) {
    idx <- cluster == g
    xg <- x[idx, , drop = FALSE]
    ug <- u[idx]
    sg <- crossprod(xg, ug)
    meat <- meat + tcrossprod(sg)
  }
  n <- nrow(x)
  k <- ncol(x)
  g <- nlevels(cluster)
  correction <- (g / (g - 1)) * ((n - 1) / (n - k))
  correction * xtx_inv %*% meat %*% xtx_inv
}

format_estimate <- function(effect, lo, hi) {
  sprintf("%.2f (%.2f to %.2f)", effect, lo, hi)
}

fit_phenotype <- function(data, phenotype_var, label, definition) {
  data <- data %>%
    mutate(
      phenotype = factor(
        case_when(
          .data[[phenotype_var]] == 1 ~ "Present",
          .data[[phenotype_var]] == 0 ~ "Absent",
          TRUE ~ NA_character_
        ),
        levels = c("Absent", "Present")
      )
    )

  covariates <- c(
    "age_model", "I(age_model^2)", "female", "education_group",
    "married_or_partnered", "rural_community", "agricultural_hukou",
    "public_insurance", "any_adl_limitation", "poor_self_rated_health",
    "household_income_quartile"
  )
  rhs <- c("z_tcm_beds_per_10000 * phenotype", covariates, "province_fe", "wave_fe")
  f <- as.formula(paste("outcome_primary ~", paste(rhs, collapse = " + ")))
  needed <- all.vars(f)
  d <- data %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
  if (n_distinct(d$phenotype) < 2) return(tibble())

  fit <- lm(f, data = d)
  vc <- cluster_vcov(fit, d$province_supply_key)
  coefs <- coef(fit)
  z_name <- "z_tcm_beds_per_10000"
  int_name <- "z_tcm_beds_per_10000:phenotypePresent"
  df <- n_distinct(d$province_supply_key) - 1

  absent_effect <- unname(coefs[[z_name]])
  absent_se <- sqrt(vc[z_name, z_name])
  absent_ci <- absent_effect + c(-1, 1) * qt(0.975, df = df) * absent_se

  interaction_effect <- unname(coefs[[int_name]])
  interaction_se <- sqrt(vc[int_name, int_name])
  interaction_ci <- interaction_effect + c(-1, 1) * qt(0.975, df = df) * interaction_se
  interaction_p <- 2 * pt(abs(interaction_effect / interaction_se), df = df, lower.tail = FALSE)

  contrast <- rep(0, length(coefs))
  names(contrast) <- names(coefs)
  contrast[z_name] <- 1
  contrast[int_name] <- 1
  present_effect <- sum(contrast * coefs)
  present_se <- sqrt(as.numeric(t(contrast) %*% vc %*% contrast))
  present_ci <- present_effect + c(-1, 1) * qt(0.975, df = df) * present_se
  present_p <- 2 * pt(abs(present_effect / present_se), df = df, lower.tail = FALSE)

  tibble(
    phenotype = label,
    definition = definition,
    phenotype_prevalence_percent = 100 * mean(d$phenotype == "Present"),
    tcm_use_prevalence_if_present_percent = 100 * mean(d$outcome_primary[d$phenotype == "Present"] == 1),
    tcm_use_prevalence_if_absent_percent = 100 * mean(d$outcome_primary[d$phenotype == "Absent"] == 1),
    bed_association_if_absent = format_estimate(100 * absent_effect, 100 * absent_ci[[1]], 100 * absent_ci[[2]]),
    bed_association_if_present = format_estimate(100 * present_effect, 100 * present_ci[[1]], 100 * present_ci[[2]]),
    interaction_estimate = format_estimate(100 * interaction_effect, 100 * interaction_ci[[1]], 100 * interaction_ci[[2]]),
    interaction_pvalue = interaction_p,
    bed_association_present_pvalue = present_p,
    observations = nrow(d),
    respondents = n_distinct(d$person_id),
    provinces = n_distinct(d$province_supply_key),
    fixed_effects = "province and survey year",
    covariate_set = paste(covariates, collapse = "; ")
  )
}

results <- bind_rows(
  fit_phenotype(
    analysis,
    "phenotype_cardiometabolic",
    "Cardiometabolic or cerebrovascular multimorbidity",
    "At least 2 diagnosed conditions among hypertension, diabetes, dyslipidemia, heart disease, and stroke"
  ),
  fit_phenotype(
    analysis,
    "phenotype_arthritis_multimorbidity",
    "Arthritis with multimorbidity",
    "Arthritis or rheumatism plus at least 1 other chronic condition"
  ),
  fit_phenotype(
    analysis,
    "phenotype_respiratory_multimorbidity",
    "Chronic lung disease with multimorbidity",
    "Chronic lung disease plus at least 1 other chronic condition"
  ),
  fit_phenotype(
    analysis,
    "phenotype_digestive_liver_kidney",
    "Digestive, liver, or kidney disease with multimorbidity",
    "Digestive disease, liver disease, or kidney disease among respondents with at least 2 chronic conditions"
  ),
  fit_phenotype(
    analysis,
    "phenotype_neurocognitive_stroke",
    "Memory disease or stroke with multimorbidity",
    "Memory-related disease or stroke among respondents with at least 2 chronic conditions"
  )
) %>%
  mutate(
    interaction_holm_pvalue = p.adjust(interaction_pvalue, method = "holm"),
    across(
      c(
        phenotype_prevalence_percent,
        tcm_use_prevalence_if_present_percent,
        tcm_use_prevalence_if_absent_percent,
        interaction_pvalue,
        bed_association_present_pvalue,
        interaction_holm_pvalue
      ),
      ~ round(.x, 4)
    )
  )

write_tsv(results, out_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
print(results)
