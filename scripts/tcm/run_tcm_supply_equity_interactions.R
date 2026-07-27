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
  "local CHARLS analytic table"
)
out_tsv <- argval(
  "--output",
  "results/tcm/tcm_supply_equity_interactions.tsv"
)

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

analysis <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  filter(main_model_age60) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    z_tcm_beds_per_10000 = as.numeric(scale(value_per_10000_population_tcm_hospital_beds)),
    outcome_primary = as.numeric(primary_condition_tcm_any),
    rural_stratum = factor(case_when(
      rural_community == 0 ~ "urban",
      rural_community == 1 ~ "rural",
      TRUE ~ NA_character_
    ), levels = c("urban", "rural")),
    education_stratum = factor(case_when(
      education_group %in% c("upper_secondary_or_vocational", "tertiary") ~ "higher_education",
      education_group %in% c("none", "less_than_lower_secondary") ~ "lower_education",
      TRUE ~ NA_character_
    ), levels = c("higher_education", "lower_education")),
    income_stratum = factor(case_when(
      household_income_quartile %in% c("Q3", "Q4_highest") ~ "higher_income",
      household_income_quartile %in% c("Q1_lowest", "Q2") ~ "lower_income",
      household_income_quartile == "missing" ~ "income_missing",
      TRUE ~ NA_character_
    ), levels = c("higher_income", "lower_income", "income_missing")),
    multichronic_stratum = factor(case_when(
      multi_chronic == 0 ~ "no_multimorbidity",
      multi_chronic == 1 ~ "multimorbidity",
      TRUE ~ NA_character_
    ), levels = c("no_multimorbidity", "multimorbidity")),
    adl_stratum = factor(case_when(
      any_adl_limitation == 0 ~ "no_adl_limitation",
      any_adl_limitation == 1 ~ "adl_limitation",
      TRUE ~ NA_character_
    ), levels = c("no_adl_limitation", "adl_limitation"))
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

fit_interaction <- function(data, stratum_var, omit_covariates = character()) {
  base_covariates <- c(
    "age_model", "I(age_model^2)", "female", "education_group",
    "married_or_partnered", "rural_community", "agricultural_hukou",
    "public_insurance", "chronic_count", "any_adl_limitation",
    "poor_self_rated_health", "household_income_quartile"
  )
  covariates <- setdiff(base_covariates, omit_covariates)
  rhs <- c(paste0("z_tcm_beds_per_10000 * ", stratum_var), covariates, "province_fe", "wave_fe")
  f <- as.formula(paste("outcome_primary ~", paste(rhs, collapse = " + ")))
  needed <- all.vars(f)
  d <- data %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
  fit <- lm(f, data = d)
  vc <- cluster_vcov(fit, d$province_supply_key)
  coefs <- coef(fit)
  se <- sqrt(diag(vc))
  interaction_terms <- names(coefs)[grepl("^z_tcm_beds_per_10000:", names(coefs))]
  if (length(interaction_terms) == 0) {
    return(tibble())
  }
  bind_rows(lapply(interaction_terms, function(term) {
    effect <- unname(coefs[[term]])
    se_term <- unname(se[[term]])
    df <- length(unique(d$province_supply_key)) - 1
    ci <- effect + c(-1, 1) * qt(0.975, df = df) * se_term
    p_value <- 2 * pt(abs(effect / se_term), df = df, lower.tail = FALSE)
    tibble(
      stratum = stratum_var,
      contrast = sub("^z_tcm_beds_per_10000:", "", term),
      reference_level = levels(d[[stratum_var]])[[1]],
      interaction_term = term,
      effect_type = "additional percentage-point difference per 1 SD higher bed density versus reference level",
      interaction_effect = 100 * effect,
      ci_lower = 100 * ci[[1]],
      ci_upper = 100 * ci[[2]],
      pvalue = p_value,
      n_persons = n_distinct(d$person_id),
      n_person_waves = nrow(d),
      n_clusters = n_distinct(d$province_supply_key),
      fixed_effects = "province and survey year",
      covariate_set = paste(covariates, collapse = "; "),
      source_table = analysis_tsv
    )
  }))
}

results <- bind_rows(
  fit_interaction(analysis, "rural_stratum", "rural_community"),
  fit_interaction(analysis, "education_stratum", "education_group"),
  fit_interaction(analysis, "income_stratum", "household_income_quartile"),
  fit_interaction(analysis, "multichronic_stratum", "chronic_count"),
  fit_interaction(analysis, "adl_stratum", "any_adl_limitation")
) %>%
  mutate(
    pvalue_holm = p.adjust(pvalue, method = "holm"),
    pvalue_fdr_bh = p.adjust(pvalue, method = "BH"),
    multiplicity_family = "five prespecified equity interaction domains",
    across(
      c(interaction_effect, ci_lower, ci_upper, pvalue, pvalue_holm, pvalue_fdr_bh),
      ~ round(.x, 5)
    )
  )

write_tsv(results, out_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
print(results)
