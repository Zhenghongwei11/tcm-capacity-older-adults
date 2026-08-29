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
  "results/tcm/tcm_supply_capital_labor_models.tsv"
)
diagnostics_tsv <- argval(
  "--diagnostics-output",
  "results/tcm/tcm_supply_capital_labor_diagnostics.tsv"
)

if (!file.exists(analysis_tsv)) stop(sprintf("Missing analysis table: %s", analysis_tsv))
dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

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
  groups <- nlevels(cluster)
  correction <- (groups / (groups - 1)) * ((n - 1) / (n - k))
  out <- correction * xtx_inv %*% meat %*% xtx_inv
  dimnames(out) <- list(colnames(x), colnames(x))
  out
}

analysis_raw <- read_tsv(analysis_tsv, show_col_types = FALSE)

standardization_panel <- analysis_raw %>%
  filter(main_model_age60) %>%
  distinct(
    province_supply_key,
    year,
    value_tcm_hospital_beds,
    value_tcm_practicing_assistant_physicians
  ) %>%
  mutate(
    tcm_beds_per_physician = value_tcm_hospital_beds /
      value_tcm_practicing_assistant_physicians
  )

ratio_mean <- mean(standardization_panel$tcm_beds_per_physician, na.rm = TRUE)
ratio_sd <- sd(standardization_panel$tcm_beds_per_physician, na.rm = TRUE)
if (!is.finite(ratio_sd) || ratio_sd <= 0) stop("Bed-to-physician ratio has no usable variation")

standardization_panel <- standardization_panel %>%
  mutate(
    z_tcm_beds_per_physician = (tcm_beds_per_physician - ratio_mean) / ratio_sd
  )

analysis <- analysis_raw %>%
  filter(main_model_age60) %>%
  left_join(
    standardization_panel %>%
      select(province_supply_key, year, tcm_beds_per_physician, z_tcm_beds_per_physician),
    by = c("province_supply_key", "year")
  ) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    outcome_primary = as.numeric(primary_condition_tcm_any),
    z_tcm_beds = z_py_tcm_beds_per_10000,
    z_tcm_physicians = z_py_tcm_physicians_per_10000,
    z_resource_sum = (z_py_tcm_beds_per_10000 + z_py_tcm_physicians_per_10000) / 2,
    z_resource_difference = (z_py_tcm_beds_per_10000 - z_py_tcm_physicians_per_10000) / 2,
    bed_physician_interaction = z_tcm_beds * z_tcm_physicians
  )

individual_covariates <- c(
  "age_model", "I(age_model^2)", "female", "education_group",
  "married_or_partnered", "rural_community", "agricultural_hukou",
  "public_insurance", "chronic_count", "any_adl_limitation",
  "poor_self_rated_health", "household_income_quartile"
)

model_specs <- list(
  CL1 = c("z_tcm_beds_per_physician"),
  CL2 = c("z_tcm_beds", "z_tcm_physicians"),
  CL3 = c("z_tcm_beds", "z_tcm_physicians", "z_tcm_beds_per_physician"),
  CL4 = c("z_tcm_beds", "z_tcm_physicians", "bed_physician_interaction")
)

model_labels <- c(
  CL1 = "Bed-to-physician ratio only",
  CL2 = "Bed density and physician density",
  CL3 = "Bed density, physician density, and bed-to-physician ratio",
  CL4 = "Bed density, physician density, and their interaction"
)

term_labels <- c(
  z_tcm_beds_per_physician = "TCM hospital beds per TCM practicing or assistant physician",
  z_tcm_beds = "TCM hospital beds per 10,000 population",
  z_tcm_physicians = "TCM practicing or assistant physicians per 10,000 population",
  bed_physician_interaction = "Interaction between TCM bed density and physician density"
)

term_effect_types <- c(
  z_tcm_beds_per_physician = "percentage-point difference per 1-SD higher province-year bed-to-physician ratio",
  z_tcm_beds = "percentage-point difference per 1-SD higher province-year bed density",
  z_tcm_physicians = "percentage-point difference per 1-SD higher province-year physician density",
  bed_physician_interaction = "percentage-point interaction between province-year bed and physician density"
)

lookup_chr <- function(x, key) {
  value <- unname(x[key])
  if (length(value) == 0 || is.na(value)) return(NA_character_)
  as.character(value[[1]])
}

scalar_dbl <- function(x) {
  value <- suppressWarnings(as.numeric(x))
  if (length(value) == 0) return(NA_real_)
  value[[1]]
}

fit_spec <- function(model_id, terms) {
  rhs <- c(terms, individual_covariates, "province_fe", "wave_fe")
  formula <- as.formula(paste("outcome_primary ~", paste(rhs, collapse = " + ")))
  needed <- all.vars(formula)
  d <- analysis %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
  fit <- lm(formula, data = d)
  conventional_vcov <- cluster_vcov(fit, d$province_supply_key)
  conventional_df <- n_distinct(d$province_supply_key) - 1
  cr2 <- tryCatch(
    coef_test(fit, vcov = "CR2", cluster = d$province_supply_key, test = "Satterthwaite"),
    error = function(e) NULL
  )

  rows <- list()
  model_terms <- names(coef(fit))
  for (term in terms) {
    if (!(term %in% model_terms)) {
      rows <- append(rows, list(tibble(
        model = model_id,
        model_label = lookup_chr(model_labels, model_id),
        exposure = lookup_chr(term_labels, term),
        effect_type = lookup_chr(term_effect_types, term),
        effect = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        conventional_pvalue = NA_real_,
        cr2_ci_lower = NA_real_,
        cr2_ci_upper = NA_real_,
        cr2_pvalue = NA_real_,
        cr2_df = NA_real_,
        n_persons = n_distinct(d$person_id),
        n_person_waves = nrow(d),
        n_clusters = n_distinct(d$province_supply_key),
        n_province_years = n_distinct(paste(d$province_supply_key, d$year)),
        interpretation_role = "Exploratory mechanism context",
        estimation_note = "Term was not estimable in this model"
      )))
      next
    }

    beta <- scalar_dbl(coef(fit)[term])
    se <- scalar_dbl(sqrt(conventional_vcov[term, term]))
    ci <- beta + c(-1, 1) * qt(0.975, df = conventional_df) * se
    pvalue <- if (is.na(beta) || is.na(se) || se <= 0) {
      NA_real_
    } else {
      2 * pt(abs(beta / se), df = conventional_df, lower.tail = FALSE)
    }
    cr2_p <- NA_real_
    cr2_df <- NA_real_
    cr2_ci <- c(NA_real_, NA_real_)
    if (!is.null(cr2) && term %in% rownames(cr2)) {
      if ("p_Satt" %in% colnames(cr2)) cr2_p <- scalar_dbl(cr2[term, "p_Satt", drop = TRUE])
      if ("df_Satt" %in% colnames(cr2)) cr2_df <- scalar_dbl(cr2[term, "df_Satt", drop = TRUE])
      if ("SE" %in% colnames(cr2) && is.finite(cr2_df)) {
        cr2_se <- scalar_dbl(cr2[term, "SE", drop = TRUE])
        cr2_ci <- beta + c(-1, 1) * qt(0.975, df = cr2_df) * cr2_se
      }
    }
    rows <- append(rows, list(tibble(
      model = model_id,
      model_label = lookup_chr(model_labels, model_id),
      exposure = lookup_chr(term_labels, term),
      effect_type = lookup_chr(term_effect_types, term),
      effect = 100 * beta,
      ci_lower = 100 * ci[[1]],
      ci_upper = 100 * ci[[2]],
      conventional_pvalue = pvalue,
      cr2_ci_lower = 100 * cr2_ci[[1]],
      cr2_ci_upper = 100 * cr2_ci[[2]],
      cr2_pvalue = cr2_p,
      cr2_df = cr2_df,
      n_persons = n_distinct(d$person_id),
      n_person_waves = nrow(d),
      n_clusters = n_distinct(d$province_supply_key),
      n_province_years = n_distinct(paste(d$province_supply_key, d$year)),
      interpretation_role = "Exploratory mechanism context",
      estimation_note = "Province and survey-year fixed effects with individual covariates"
    )))
  }
  bind_rows(rows)
}

fit_bed_physician_contrast <- function() {
  rhs <- c("z_resource_sum", "z_resource_difference", individual_covariates, "province_fe", "wave_fe")
  formula <- as.formula(paste("outcome_primary ~", paste(rhs, collapse = " + ")))
  needed <- all.vars(formula)
  d <- analysis %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
  fit <- lm(formula, data = d)
  term <- "z_resource_difference"
  conventional_vcov <- cluster_vcov(fit, d$province_supply_key)
  conventional_df <- n_distinct(d$province_supply_key) - 1
  beta <- scalar_dbl(coef(fit)[term])
  se <- scalar_dbl(sqrt(conventional_vcov[term, term]))
  ci <- beta + c(-1, 1) * qt(0.975, df = conventional_df) * se
  pvalue <- 2 * pt(abs(beta / se), df = conventional_df, lower.tail = FALSE)

  cr2 <- tryCatch(
    coef_test(fit, vcov = "CR2", cluster = d$province_supply_key, test = "Satterthwaite"),
    error = function(e) NULL
  )
  cr2_p <- cr2_df <- NA_real_
  cr2_ci <- c(NA_real_, NA_real_)
  if (!is.null(cr2) && term %in% rownames(cr2)) {
    cr2_se <- scalar_dbl(cr2[term, "SE", drop = TRUE])
    cr2_df <- scalar_dbl(cr2[term, "df_Satt", drop = TRUE])
    cr2_p <- scalar_dbl(cr2[term, "p_Satt", drop = TRUE])
    cr2_ci <- beta + c(-1, 1) * qt(0.975, df = cr2_df) * cr2_se
  }

  tibble(
    model = "CL2_CONTRAST",
    model_label = "Bed density and physician density coefficient contrast",
    exposure = "Difference between TCM bed-density and physician-density associations",
    effect_type = "percentage-point difference between coefficients per 1-SD higher province-year resource density",
    effect = 100 * beta,
    ci_lower = 100 * ci[[1]],
    ci_upper = 100 * ci[[2]],
    conventional_pvalue = pvalue,
    cr2_ci_lower = 100 * cr2_ci[[1]],
    cr2_ci_upper = 100 * cr2_ci[[2]],
    cr2_pvalue = cr2_p,
    cr2_df = cr2_df,
    n_persons = n_distinct(d$person_id),
    n_person_waves = nrow(d),
    n_clusters = n_distinct(d$province_supply_key),
    n_province_years = n_distinct(paste(d$province_supply_key, d$year)),
    interpretation_role = "Formal exploratory coefficient contrast",
    estimation_note = "Reparameterized joint model; coefficient equals bed-density coefficient minus physician-density coefficient"
  )
}

residualize_exposures <- function(panel, terms) {
  out <- lapply(terms, function(term) {
    residuals(lm(
      as.formula(paste(term, "~ factor(province_supply_key) + factor(year)")),
      data = panel
    ))
  })
  names(out) <- terms
  as.data.frame(out)
}

diagnose_spec <- function(model_id, terms) {
  panel <- analysis %>%
    distinct(province_supply_key, year, across(all_of(terms))) %>%
    filter(if_all(all_of(terms), ~ !is.na(.x)))

  x <- residualize_exposures(panel, terms)
  x <- x[, vapply(x, function(col) sd(col, na.rm = TRUE) > 0, logical(1)), drop = FALSE]

  if (ncol(x) == 0) {
    return(tibble(
      model = model_id,
      model_label = lookup_chr(model_labels, model_id),
      exposure = "No estimable exposure variation",
      two_way_fe_vif = NA_real_,
      maximum_absolute_pairwise_correlation = NA_real_,
      design_condition_index = NA_real_,
      n_province_years = nrow(panel),
      n_provinces = n_distinct(panel$province_supply_key),
      diagnostic_note = "No residual exposure variation after province and survey-year fixed effects"
    ))
  }

  x_scaled <- scale(as.matrix(x))
  correlation <- if (ncol(x_scaled) == 1) matrix(1, 1, 1) else cor(x_scaled)
  eigenvalues <- eigen(crossprod(x_scaled), symmetric = TRUE, only.values = TRUE)$values
  positive <- eigenvalues[eigenvalues > max(eigenvalues) * .Machine$double.eps]
  condition_index <- sqrt(max(positive) / min(positive))

  bind_rows(lapply(seq_len(ncol(x_scaled)), function(j) {
    term <- colnames(x_scaled)[[j]]
    other <- setdiff(seq_len(ncol(x_scaled)), j)
    vif <- if (length(other) == 0) {
      1
    } else {
      1 / (1 - summary(lm(x_scaled[, j] ~ x_scaled[, other, drop = FALSE]))$r.squared)
    }
    tibble(
      model = model_id,
      model_label = lookup_chr(model_labels, model_id),
      exposure = lookup_chr(term_labels, term),
      two_way_fe_vif = vif,
      maximum_absolute_pairwise_correlation = if (length(other) == 0) 0 else max(abs(correlation[j, other])),
      design_condition_index = condition_index,
      n_province_years = nrow(panel),
      n_provinces = n_distinct(panel$province_supply_key),
      diagnostic_note = "Province-year exposure residuals after province and survey-year fixed effects"
    )
  }))
}

result_list <- list()
for (model_id in names(model_specs)) {
  model_result <- fit_spec(model_id, unname(model_specs[[model_id]]))
  result_list <- append(result_list, list(model_result))
}
results <- dplyr::bind_rows(result_list)
results <- bind_rows(results, fit_bed_physician_contrast())
results <- results %>%
  mutate(across(
    c(effect, ci_lower, ci_upper, conventional_pvalue, cr2_ci_lower,
      cr2_ci_upper, cr2_pvalue, cr2_df),
    ~ round(.x, 5)
  ))

if (nrow(results) == 0) stop("Capital-labor model results are empty")

diagnostics <- bind_rows(lapply(names(model_specs), function(model_id) {
  diagnose_spec(model_id, unname(model_specs[[model_id]]))
})) %>%
  mutate(across(
    c(two_way_fe_vif, maximum_absolute_pairwise_correlation, design_condition_index),
    ~ round(.x, 5)
  ))

write_tsv(results, out_tsv)
write_tsv(diagnostics, diagnostics_tsv)

cat(sprintf("Wrote %s\n", out_tsv))
cat(sprintf("Wrote %s\n", diagnostics_tsv))
cat(sprintf(
  "Bed-to-physician ratio standardization: mean %.5f, SD %.5f across %d linked province-years\n",
  ratio_mean, ratio_sd, nrow(standardization_panel)
))
print(results %>% select(model, exposure, effect, ci_lower, ci_upper, conventional_pvalue, cr2_pvalue))
print(diagnostics %>% select(model, exposure, two_way_fe_vif, maximum_absolute_pairwise_correlation, design_condition_index))
