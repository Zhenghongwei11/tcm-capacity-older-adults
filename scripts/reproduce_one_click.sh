#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw/tcm_supply/official_excels data/raw/province_covariates results/tcm plots/publication tables

echo "Checking expected CHARLS directory structure."
if [[ ! -d data/raw/charls/downloads ]]; then
  echo "CHARLS files are not present. See data/README.md for access and placement instructions." >&2
  exit 1
fi

Rscript scripts/tcm/scan_charls_tcm_outcomes.R
Rscript scripts/tcm/download_natcm_core_density_sources.R
Rscript scripts/tcm/extract_natcm_core_density_panel.R
Rscript scripts/tcm/link_charls_tcm_core_density.R

Rscript scripts/tcm/extract_health_system_bed_panel.R
Rscript scripts/tcm/extract_nbs_province_covariates.R
Rscript scripts/tcm/build_province_year_contextual_panel.R
Rscript scripts/tcm/build_charls_covariates.R

Rscript scripts/tcm/run_tcm_supply_main_models.R
Rscript scripts/tcm/run_tcm_supply_strict_models.R
Rscript scripts/tcm/run_tcm_supply_broader_models.R
Rscript scripts/tcm/run_tcm_supply_equity_models.R
Rscript scripts/tcm/run_tcm_supply_equity_interactions.R
Rscript scripts/tcm/run_tcm_supply_robustness.R
Rscript scripts/tcm/run_tcm_supply_component_and_cluster_checks.R
Rscript scripts/tcm/run_tcm_supply_contextual_models.R
Rscript scripts/tcm/run_tcm_supply_falsification_models.R
Rscript scripts/tcm/run_tcm_supply_bed_physician_models.R
Rscript scripts/tcm/audit_charls_longitudinal_sample.R
Rscript scripts/tcm/run_tcm_supply_weighted_attrition_models.R
Rscript scripts/tcm/run_tcm_supply_longitudinal_identification_models.R
Rscript scripts/tcm/run_tcm_supply_composite_stability.R

Rscript scripts/tcm/build_reported_tables.R
Rscript scripts/tcm/make_publication_figures.R

echo "Rebuild complete."
