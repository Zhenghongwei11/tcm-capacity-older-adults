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
  "results/tcm/tcm_supply_equity_heterogeneity.tsv"
)

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

east <- c("北京市", "天津市", "河北省", "上海市", "江苏省", "浙江省", "福建省", "山东省", "广东省", "海南省")
central <- c("山西省", "安徽省", "江西省", "河南省", "湖北省", "湖南省")
northeast <- c("辽宁省", "吉林省", "黑龙江省")
west <- c("内蒙古自治区", "广西壮族自治区", "重庆市", "四川省", "贵州省", "云南省", "西藏自治区", "陕西省", "甘肃省", "青海省", "宁夏回族自治区", "新疆维吾尔自治区")

analysis <- read_tsv(analysis_tsv, show_col_types = FALSE) %>%
  filter(main_model_age60) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    z_tcm_physicians_per_10000 = z_py_tcm_physicians_per_10000,
    z_tcm_beds_per_10000 = z_py_tcm_beds_per_10000,
    outcome_primary = as.numeric(primary_condition_tcm_any),
    equity_rural_community = case_when(
      rural_community == 1 ~ "rural",
      rural_community == 0 ~ "urban",
      TRUE ~ NA_character_
    ),
    equity_education = case_when(
      education_group %in% c("none", "less_than_lower_secondary") ~ "lower_education",
      education_group %in% c("upper_secondary_or_vocational", "tertiary") ~ "higher_education",
      TRUE ~ NA_character_
    ),
    equity_income = case_when(
      household_income_quartile %in% c("Q1_lowest", "Q2") ~ "lower_income",
      household_income_quartile %in% c("Q3", "Q4_highest") ~ "higher_income",
      household_income_quartile == "missing" ~ "income_missing",
      TRUE ~ NA_character_
    ),
    equity_multichronic = case_when(
      multi_chronic == 1 ~ "multimorbidity",
      multi_chronic == 0 ~ "no_multimorbidity",
      TRUE ~ NA_character_
    ),
    equity_adl = case_when(
      any_adl_limitation == 1 ~ "adl_limitation",
      any_adl_limitation == 0 ~ "no_adl_limitation",
      TRUE ~ NA_character_
    ),
    equity_region = case_when(
      province_supply_key %in% east ~ "east",
      province_supply_key %in% central ~ "central",
      province_supply_key %in% northeast ~ "northeast",
      province_supply_key %in% west ~ "west",
      TRUE ~ NA_character_
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
  if (g <= 1 || n <= k) return(matrix(NA_real_, ncol(x), ncol(x)))
  correction <- (g / (g - 1)) * ((n - 1) / (n - k))
  correction * xtx_inv %*% meat %*% xtx_inv
}

fit_stratum <- function(data, stratum_var, stratum_level, supply_var) {
  base_covariates <- c(
    "age_model", "I(age_model^2)", "female", "education_group",
    "married_or_partnered", "rural_community", "agricultural_hukou",
    "public_insurance", "chronic_count", "any_adl_limitation",
    "poor_self_rated_health", "household_income_quartile"
  )
  omit <- switch(
    stratum_var,
    equity_rural_community = "rural_community",
    equity_education = "education_group",
    equity_income = "household_income_quartile",
    equity_multichronic = "chronic_count",
    equity_adl = "any_adl_limitation",
    equity_region = character(),
    character()
  )
  covariates <- setdiff(base_covariates, omit)
  rhs <- c(supply_var, covariates, "province_fe", "wave_fe")
  f <- as.formula(paste("outcome_primary ~", paste(rhs, collapse = " + ")))
  needed <- all.vars(f)
  d <- data %>%
    filter(.data[[stratum_var]] == stratum_level) %>%
    filter(if_all(all_of(needed), ~ !is.na(.x)))

  events <- sum(d$outcome_primary == 1, na.rm = TRUE)
  if (nrow(d) < 500 || events < 50 || n_distinct(d$province_supply_key) < 5) {
    return(tibble(
      stratum = stratum_var,
      stratum_level = stratum_level,
      supply_indicator = supply_var,
      model = "E1_covariate_adjusted_lpm",
      effect_type = "percentage-point difference per 1 SD higher province-year supply",
      effect = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      pvalue = NA_real_,
      n_persons = n_distinct(d$person_id),
      n_person_waves = nrow(d),
      events = events,
      cluster_level = "province",
      n_clusters = n_distinct(d$province_supply_key),
      status = "not_estimated_sparse_or_too_few_clusters"
    ))
  }

  fit <- lm(f, data = d)
  vc <- cluster_vcov(fit, d$province_supply_key)
  coefs <- coef(fit)
  se <- sqrt(diag(vc))
  effect <- unname(coefs[[supply_var]])
  se_term <- unname(se[[supply_var]])
  df <- length(unique(d$province_supply_key)) - 1
  ci <- effect + c(-1, 1) * qt(0.975, df = df) * se_term
  p_value <- 2 * pt(abs(effect / se_term), df = df, lower.tail = FALSE)

  tibble(
    stratum = stratum_var,
    stratum_level = stratum_level,
    supply_indicator = case_when(
      supply_var == "z_tcm_physicians_per_10000" ~ "TCM physicians per 10,000 population",
      supply_var == "z_tcm_beds_per_10000" ~ "TCM hospital beds per 10,000 population",
      TRUE ~ supply_var
    ),
    model = "E1_covariate_adjusted_lpm",
    effect_type = "percentage-point difference per 1 SD higher province-year supply",
    effect = 100 * effect,
    ci_lower = 100 * ci[[1]],
    ci_upper = 100 * ci[[2]],
    pvalue = p_value,
    n_persons = n_distinct(d$person_id),
    n_person_waves = nrow(d),
    events = events,
    cluster_level = "province",
    n_clusters = n_distinct(d$province_supply_key),
    status = "estimated"
  )
}

strata <- c(
  "equity_rural_community", "equity_education", "equity_income",
  "equity_multichronic", "equity_adl", "equity_region"
)
exposures <- c("z_tcm_physicians_per_10000", "z_tcm_beds_per_10000")

results <- bind_rows(lapply(strata, function(s) {
  levels_s <- sort(unique(analysis[[s]][!is.na(analysis[[s]])]))
  bind_rows(lapply(levels_s, function(level_s) {
    bind_rows(lapply(exposures, function(exposure) fit_stratum(analysis, s, level_s, exposure)))
  }))
})) %>%
  mutate(across(c(effect, ci_lower, ci_upper, pvalue), ~ round(.x, 5)))

write_tsv(results, out_tsv)
cat(sprintf("Wrote %s\n", out_tsv))
print(results)
