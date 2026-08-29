#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
})

out_path <- "tables/tables.md"

fmt_num <- function(x, digits = 2) ifelse(is.na(x), "", formatC(x, format = "f", digits = digits, big.mark = ","))
fmt_int <- function(x) ifelse(is.na(x), "", formatC(round(x), format = "d", big.mark = ","))
fmt_ci <- function(x, lo, hi, digits = 2) paste0(fmt_num(x, digits), " (", fmt_num(lo, digits), " to ", fmt_num(hi, digits), ")")
fmt_p <- function(x) ifelse(is.na(x), "", ifelse(x < 0.001, "<0.001", fmt_num(x, 3)))

md_table <- function(df) {
  header <- paste(names(df), collapse = " | ")
  separator <- paste(rep("---", ncol(df)), collapse = " | ")
  body <- apply(df, 1, function(row) paste(row, collapse = " | "))
  paste(c(paste0("| ", header, " |"), paste0("| ", separator, " |"), paste0("| ", body, " |")), collapse = "\n")
}

section <- function(title, table, note) c(paste0("## ", title), "", md_table(table), "", paste0("Note: ", note), "")

analysis <- read_tsv("results/tcm/person_wave_tcm_core_density_analysis.tsv", show_col_types = FALSE)
supply <- read_tsv("results/tcm/province_year_tcm_core_density.tsv", show_col_types = FALSE)
main <- read_tsv("results/tcm/tcm_supply_main_models.tsv", show_col_types = FALSE)
strict <- read_tsv("results/tcm/tcm_supply_strict_secondary_models.tsv", show_col_types = FALSE)
cluster <- read_tsv("results/tcm/tcm_supply_cluster_sensitivity.tsv", show_col_types = FALSE)
small <- read_tsv("results/tcm/tcm_supply_small_cluster_inference.tsv", show_col_types = FALSE)
contextual <- read_tsv("results/tcm/tcm_supply_contextual_models.tsv", show_col_types = FALSE)
joint <- read_tsv("results/tcm/tcm_supply_bed_physician_models.tsv", show_col_types = FALSE)
capital <- read_tsv("results/tcm/tcm_supply_capital_labor_models.tsv", show_col_types = FALSE)
opportunity <- read_tsv("results/tcm/tcm_supply_multimorbidity_opportunity_checks.tsv", show_col_types = FALSE)
interactions <- read_tsv("results/tcm/tcm_supply_equity_interactions.tsv", show_col_types = FALSE)
heterogeneity <- read_tsv("results/tcm/tcm_supply_equity_heterogeneity.tsv", show_col_types = FALSE)
components <- read_tsv("results/tcm/tcm_supply_primary_component_models.tsv", show_col_types = FALSE)
specificity <- read_tsv("results/tcm/tcm_supply_specificity_synthesis.tsv", show_col_types = FALSE)
weighted <- read_tsv("results/tcm/tcm_supply_weighted_attrition_models.tsv", show_col_types = FALSE)
identification <- read_tsv("results/tcm/tcm_supply_longitudinal_identification_models.tsv", show_col_types = FALSE)
stability <- read_tsv("results/tcm/tcm_supply_composite_stability.tsv", show_col_types = FALSE)
flow <- read_tsv("results/tcm/charls_participant_flow.tsv", show_col_types = FALSE)
retention <- read_tsv("results/tcm/charls_retention_by_wave.tsv", show_col_types = FALSE)
age_qc <- read_tsv("results/tcm/charls_age_harmonization_qc.tsv", show_col_types = FALSE)

age60 <- analysis %>% filter(main_model_age60)
table1 <- age60 %>%
  group_by(Year = year) %>%
  summarise(
    Observations = n(), Respondents = n_distinct(person_id),
    `Age, mean` = mean(age_model, na.rm = TRUE),
    `Female, %` = 100 * mean(female == 1, na.rm = TRUE),
    `Rural residence, %` = 100 * mean(rural_community == 1, na.rm = TRUE),
    `Public insurance, %` = 100 * mean(public_insurance == 1, na.rm = TRUE),
    `Chronic conditions, mean` = mean(chronic_count, na.rm = TRUE),
    `Multimorbidity, %` = 100 * mean(multi_chronic == 1, na.rm = TRUE),
    `ADL limitation, %` = 100 * mean(any_adl_limitation == 1, na.rm = TRUE),
    `Disease-specific TCM use, %` = 100 * mean(primary_condition_tcm_any == 1, na.rm = TRUE),
    `TCM hospital visit, %` = 100 * mean(strict_tcm_hospital_visit == 1, na.rm = TRUE), .groups = "drop"
  ) %>%
  mutate(
    Year = as.character(Year), Observations = fmt_int(Observations), Respondents = fmt_int(Respondents),
    across(-c(Year, Observations, Respondents), ~ fmt_num(.x, 1))
  )

overall <- age60 %>% summarise(
  Year = "Overall", Observations = fmt_int(n()), Respondents = fmt_int(n_distinct(person_id)),
  `Age, mean` = fmt_num(mean(age_model, na.rm = TRUE), 1),
  `Female, %` = fmt_num(100 * mean(female == 1, na.rm = TRUE), 1),
  `Rural residence, %` = fmt_num(100 * mean(rural_community == 1, na.rm = TRUE), 1),
  `Public insurance, %` = fmt_num(100 * mean(public_insurance == 1, na.rm = TRUE), 1),
  `Chronic conditions, mean` = fmt_num(mean(chronic_count, na.rm = TRUE), 1),
  `Multimorbidity, %` = fmt_num(100 * mean(multi_chronic == 1, na.rm = TRUE), 1),
  `ADL limitation, %` = fmt_num(100 * mean(any_adl_limitation == 1, na.rm = TRUE), 1),
  `Disease-specific TCM use, %` = fmt_num(100 * mean(primary_condition_tcm_any == 1, na.rm = TRUE), 1),
  `TCM hospital visit, %` = fmt_num(100 * mean(strict_tcm_hospital_visit == 1, na.rm = TRUE), 1)
)
table1 <- bind_rows(table1, overall)

pick <- function(df, ...) df %>% filter(...) %>% slice(1)
main_bed <- pick(main, model == "M2_covariate_adjusted_lpm", str_detect(supply_indicator, "hospital beds"))
main_phys <- pick(main, model == "M2_covariate_adjusted_lpm", str_detect(supply_indicator, "physicians"))
ctx_bed <- pick(contextual, model == "C1", term == "z_tcm_beds")
ctx_comp <- pick(contextual, model == "C1", term == "z_comprehensive_beds")
strict_bed <- pick(strict, model == "S2_covariate_adjusted_lpm", str_detect(supply_indicator, "hospital beds"))
j2 <- joint %>% filter(model == "J2")
contrast <- capital %>% filter(model == "CL2_CONTRAST") %>% slice(1)
multi_original <- interactions %>% filter(stratum == "multichronic_stratum") %>% slice(1)
multi_adjusted <- opportunity %>% filter(analysis_id == "composite_multimorbidity_adjusted_for_count") %>% slice(1)
outpatient <- specificity %>% filter(`evidence_domain` == "Broader health-care use", str_detect(outcome, "doctor"), str_detect(exposure, "TCM")) %>% slice(1)

table2 <- bind_rows(
  tibble(Analysis="Primary outcome", Outcome="Disease-specific TCM treatment use", Measure="TCM hospital beds", Estimate=fmt_ci(main_bed$effect,main_bed$ci_lower,main_bed$ci_upper), `Conventional P value`=fmt_p(main_bed$pvalue), `Small-cluster inference`=paste0("CR2 ",fmt_p(small$cr2_pvalue[small$model=="JHF1_primary"]),"; wild-cluster ",fmt_p(main_bed$wild_cluster_pvalue))),
  tibble(Analysis="Primary outcome", Outcome="Disease-specific TCM treatment use", Measure="TCM physicians", Estimate=fmt_ci(main_phys$effect,main_phys$ci_lower,main_phys$ci_upper), `Conventional P value`=fmt_p(main_phys$pvalue), `Small-cluster inference`="Not estimated"),
  tibble(Analysis="Health-system context", Outcome="Disease-specific TCM treatment use", Measure="TCM beds adjusted for comprehensive-hospital beds", Estimate=fmt_ci(ctx_bed$effect,ctx_bed$ci_lower,ctx_bed$ci_upper), `Conventional P value`=fmt_p(ctx_bed$cluster_pvalue), `Small-cluster inference`=paste0("CR2 ",fmt_p(ctx_bed$cr2_pvalue),"; wild-cluster ",fmt_p(ctx_bed$wild_cluster_pvalue))),
  tibble(Analysis="Health-system context", Outcome="Disease-specific TCM treatment use", Measure="Comprehensive-hospital beds adjusted for TCM beds", Estimate=fmt_ci(ctx_comp$effect,ctx_comp$ci_lower,ctx_comp$ci_upper), `Conventional P value`=fmt_p(ctx_comp$cluster_pvalue), `Small-cluster inference`=paste0("CR2 ",fmt_p(ctx_comp$cr2_pvalue))),
  tibble(Analysis="Strict institution use", Outcome="TCM hospital visit in past month", Measure="TCM hospital beds", Estimate=fmt_ci(strict_bed$effect,strict_bed$ci_lower,strict_bed$ci_upper), `Conventional P value`=fmt_p(strict_bed$pvalue), `Small-cluster inference`="Not estimated"),
  tibble(Analysis="Broader health-care use", Outcome="Any doctor or outpatient visit", Measure="TCM hospital beds", Estimate=fmt_ci(outpatient$estimate,outpatient$ci_lower,outpatient$ci_upper), `Conventional P value`=fmt_p(outpatient$pvalue), `Small-cluster inference`="Not estimated"),
  j2 %>% filter(term %in% c("z_tcm_beds","z_tcm_physicians","z_comprehensive_beds")) %>% transmute(Analysis="Joint resource model", Outcome="Disease-specific TCM treatment use", Measure=recode(term,z_tcm_beds="TCM hospital beds",z_tcm_physicians="TCM physicians",z_comprehensive_beds="Comprehensive-hospital beds"), Estimate=fmt_ci(effect,ci_lower,ci_upper), `Conventional P value`=fmt_p(cluster_pvalue), `Small-cluster inference`=paste0("CR2 ",fmt_p(cr2_pvalue))),
  tibble(Analysis="Coefficient contrast", Outcome="Disease-specific TCM treatment use", Measure="TCM bed minus physician coefficient", Estimate=fmt_ci(contrast$effect,contrast$ci_lower,contrast$ci_upper), `Conventional P value`=fmt_p(contrast$conventional_pvalue), `Small-cluster inference`=paste0("CR2 ",fmt_ci(contrast$effect,contrast$cr2_ci_lower,contrast$cr2_ci_upper),"; P = ",fmt_p(contrast$cr2_pvalue))),
  tibble(Analysis="Need-related interaction", Outcome="Disease-specific TCM treatment use", Measure="Multimorbidity, original specification", Estimate=fmt_ci(multi_original$interaction_effect,multi_original$ci_lower,multi_original$ci_upper), `Conventional P value`=fmt_p(multi_original$pvalue), `Small-cluster inference`=paste0("Holm-adjusted P = ",fmt_p(multi_original$pvalue_holm))),
  tibble(Analysis="Need-related interaction", Outcome="Disease-specific TCM treatment use", Measure="Multimorbidity, adjusted for condition count", Estimate=fmt_ci(multi_adjusted$interaction_effect,multi_adjusted$ci_lower,multi_adjusted$ci_upper), `Conventional P value`=fmt_p(multi_adjusted$conventional_pvalue), `Small-cluster inference`=paste0("CR2 P = ",fmt_p(multi_adjusted$cr2_pvalue)))
)

supply_table <- supply %>% group_by(Year=year) %>% summarise(
  Provinces=n(), `TCM hospital beds`=sum(value_tcm_hospital_beds),
  `Beds per 10,000, mean`=mean(value_per_10000_population_tcm_hospital_beds),
  `Beds per 10,000, range`=paste0(fmt_num(min(value_per_10000_population_tcm_hospital_beds))," to ",fmt_num(max(value_per_10000_population_tcm_hospital_beds))),
  `TCM physicians`=sum(value_tcm_practicing_assistant_physicians),
  `Physicians per 10,000, mean`=mean(value_per_10000_population_tcm_practicing_assistant_physicians),
  `Physicians per 10,000, range`=paste0(fmt_num(min(value_per_10000_population_tcm_practicing_assistant_physicians))," to ",fmt_num(max(value_per_10000_population_tcm_practicing_assistant_physicians))), .groups="drop") %>%
  mutate(Year=as.character(Year),Provinces=fmt_int(Provinces),`TCM hospital beds`=fmt_int(`TCM hospital beds`),`TCM physicians`=fmt_int(`TCM physicians`),`Beds per 10,000, mean`=fmt_num(`Beds per 10,000, mean`),`Physicians per 10,000, mean`=fmt_num(`Physicians per 10,000, mean`))

flow_table <- bind_rows(
  flow %>% transmute(Section="Participant flow", `Stage or year`=stage, Respondents=fmt_int(n_persons), `Observations or count`=fmt_int(n_person_waves), Detail=""),
  retention %>% transmute(Section="Wave-to-wave observation", `Stage or year`=paste(index_year,"to",next_year), Respondents=fmt_int(n_index), `Observations or count`=fmt_int(n_observed_next_wave), Detail=paste0(fmt_num(retention_percent,1),"%")),
  age_qc %>% transmute(Section="Age-source quality control", `Stage or year`=as.character(year), Respondents="", `Observations or count`=fmt_int(n_age60_eligibility_discordant), Detail=paste0("discordant eligibility; ",fmt_int(n_raw_age_fallback)," raw-age fallbacks"))
)

cluster_main <- cluster %>% filter(check_type=="main_cluster_robust") %>% slice(1)
cluster_loo <- cluster %>% filter(check_type=="leave_one_province_out_summary") %>% slice(1)
cluster_boot <- cluster %>% filter(check_type=="province_cluster_bootstrap") %>% slice(1)
small_main <- small %>% filter(model=="JHF1_primary") %>% slice(1)
cluster_table <- tribble(
  ~`Sensitivity analysis`,~Estimate,~`Interval or detail`,~`P value`,
  "Main province-clustered estimate",fmt_num(cluster_main$effect),paste0(fmt_num(cluster_main$ci_lower)," to ",fmt_num(cluster_main$ci_upper)," (95% CI)"),fmt_p(cluster_main$pvalue),
  "CR2 small-cluster inference",fmt_num(small_main$effect),paste0(fmt_num(small_main$cr2_se)," SE; ",fmt_num(small_main$cr2_df)," df"),fmt_p(small_main$cr2_pvalue),
  "Wild-cluster inference",fmt_num(main_bed$effect),"9,999 Webb-weight repetitions",fmt_p(main_bed$wild_cluster_pvalue),
  "Leave-one-province-out",fmt_num(cluster_loo$effect),paste0(fmt_num(cluster_loo$ci_lower)," to ",fmt_num(cluster_loo$ci_upper)," (min-max)"),"",
  "Province resampling",fmt_num(cluster_boot$effect),paste0(fmt_num(cluster_boot$ci_lower)," to ",fmt_num(cluster_boot$ci_upper)," (95% interval)"),""
)

resource_table <- bind_rows(
  j2 %>% filter(term %in% c("z_tcm_beds","z_tcm_physicians","z_comprehensive_beds")) %>% transmute(Model="Joint contextual resource model",Measure=recode(term,z_tcm_beds="TCM hospital beds",z_tcm_physicians="TCM physicians",z_comprehensive_beds="Comprehensive-hospital beds"),`Conventional estimate`=fmt_ci(effect,ci_lower,ci_upper),`Conventional P value`=fmt_p(cluster_pvalue),`CR2 estimate or P value`=paste0(fmt_ci(effect,cr2_ci_lower,cr2_ci_upper),"; P = ",fmt_p(cr2_pvalue))),
  capital %>% filter(model %in% c("CL1","CL2_CONTRAST")) %>% transmute(Model=model_label,Measure=exposure,`Conventional estimate`=fmt_ci(effect,ci_lower,ci_upper),`Conventional P value`=fmt_p(conventional_pvalue),`CR2 estimate or P value`=paste0(fmt_ci(effect,cr2_ci_lower,cr2_ci_upper),"; P = ",fmt_p(cr2_pvalue)))
)

condition_labels <- c(
  hypertension = "Hypertension", diabetes = "Diabetes", cancer = "Cancer",
  chronic_lung = "Chronic lung disease", heart = "Heart disease",
  arthritis = "Arthritis", dyslipidemia = "Dyslipidemia", liver = "Liver disease",
  kidney = "Kidney disease", digestive = "Digestive disease", memory = "Memory disorder"
)

opportunity_table <- opportunity %>% transmute(
  Analysis=case_when(
    analysis_id == "composite_continuous_burden" ~ "Continuous chronic-condition burden",
    analysis_id == "composite_multimorbidity_adjusted_for_count" ~ "Multimorbidity adjusted for condition count",
    TRUE ~ recode(condition, !!!condition_labels)
  ),
  Denominator=case_when(
    denominator == "All age-eligible observations" ~ denominator,
    !is.na(condition) ~ paste0("Respondents diagnosed with ", str_to_lower(recode(condition, !!!condition_labels))),
    TRUE ~ str_replace_all(denominator,"_"," ")
  ), `Interaction estimate`=fmt_ci(interaction_effect,ci_lower,ci_upper),
  `Conventional P value`=fmt_p(conventional_pvalue), `CR2 estimate`=fmt_ci(interaction_effect,cr2_ci_lower,cr2_ci_upper),
  `CR2 P value`=fmt_p(cr2_pvalue), `Holm-adjusted CR2 P value`=fmt_p(cr2_pvalue_holm), Events=fmt_int(n_events)
)

interaction_table <- interactions %>% transmute(
  Group=recode(stratum,rural_stratum="Rural residence",education_stratum="Lower education",income_stratum=if_else(str_detect(contrast,"lower_income"),"Lower income","Income missing"),multichronic_stratum="Multimorbidity",adl_stratum="ADL limitation"),
  `Interaction estimate`=fmt_ci(interaction_effect,ci_lower,ci_upper),`P value`=fmt_p(pvalue),`Holm-adjusted P value`=fmt_p(pvalue_holm),Provinces=fmt_int(n_clusters)
)
region_table <- heterogeneity %>% filter(stratum=="equity_region",str_detect(supply_indicator,"beds")) %>% transmute(
  Group=paste0(str_to_title(stratum_level)," region"),`Interaction estimate`=if_else(status=="estimated",fmt_ci(effect,ci_lower,ci_upper),"Not estimated"),`P value`=fmt_p(pvalue),`Holm-adjusted P value`="",Provinces=fmt_int(n_clusters)
)
interaction_table <- bind_rows(interaction_table,region_table)

component_table <- components %>% transmute(
  Component=case_when(
    outcome == "condition_tcm_stroke_chinese_medicine" ~ "Stroke, Chinese medicine",
    outcome == "condition_tcm_stroke_acupuncture_moxibustion" ~ "Stroke, acupuncture or moxibustion",
    TRUE ~ str_to_sentence(str_replace_all(str_remove(outcome,"^condition_tcm_"),"_"," "))
  ), Events=fmt_int(n_events),`Event prevalence, %`=fmt_num(event_percent),Estimate=fmt_ci(effect,ci_lower,ci_upper),`P value`=fmt_p(pvalue))
specificity_table <- specificity %>% transmute(Analysis=evidence_domain,Outcome=outcome,Exposure=exposure,Estimate=fmt_ci(estimate,ci_lower,ci_upper),`P value`=fmt_p(pvalue),`Small-cluster P value`=fmt_p(small_cluster_pvalue))
weighted_table <- bind_rows(
  weighted %>% transmute(Model=case_when(model=="W0_unweighted_reproduction"~"Unweighted primary model",model=="W1_adjusted_cross_sectional_weight"~"Adjusted respondent-weight model",TRUE~"Trimmed respondent-weight model"),Estimate=fmt_ci(effect,ci_lower,ci_upper),`P value`=fmt_p(pvalue),Respondents=fmt_int(n_persons),Observations=fmt_int(n_person_waves)),
  identification %>% filter(model!="I0_primary_reproduction") %>% transmute(Model=model_label,Estimate=fmt_ci(effect,ci_lower,ci_upper),`P value`=fmt_p(pvalue),Respondents=fmt_int(n_persons),Observations=fmt_int(n_person_waves))
)
stability_table <- stability %>% transmute(`Composite definition`=str_to_sentence(str_replace_all(scenario,"_"," ")),`Outcome prevalence, %`=fmt_num(outcome_prevalence_percent),Estimate=fmt_ci(effect,ci_lower,ci_upper),`P value`=fmt_p(pvalue),Observations=fmt_int(n_person_waves))

lines <- c("# Manuscript Tables", "",
  section("Table 1. Characteristics of CHARLS Respondents Aged 60 Years or Older", table1, "Values are unweighted person-wave summaries. The fully adjusted sample comprised 32,634 observations from 13,083 respondents. ADL indicates activities of daily living; TCM, traditional Chinese medicine."),
  section("Table 2. TCM Capacity, Resource Context, and Treatment Use Among Older Adults", table2, "Estimates are percentage-point differences per 1-SD higher province-year measure. CR2 and wild-cluster results address inference with 28 province clusters."),
  section("Supplementary Table 1. Participant Flow, Retention, and Age-Source Quality Control", flow_table, "Age eligibility and adjustment use Harmonized CHARLS age, with raw-wave age used only when harmonized age is unavailable."),
  section("Supplementary Table 2. Provincial TCM Hospital Bed and Physician Supply, 2011-2018", supply_table, "Data are from official National Administration of Traditional Chinese Medicine statistical extracts."),
  section("Supplementary Table 3. Small-Cluster and Province-Influence Sensitivity", cluster_table, "The leave-one-province-out row reports a median and min-max range, not a confidence interval."),
  section("Supplementary Table 4. Joint Resource and Capital-Labor Structure Models", resource_table, "Coefficients in joint models are mutually adjusted for correlated province-year resource measures."),
  section("Supplementary Table 5. Outcome and Resource Specificity Synthesis", specificity_table, "General outpatient and hospital outcomes are broader utilization checks, not strict negative controls."),
  section("Supplementary Table 6. Weighting and Longitudinal Identification", weighted_table, "Weighted models target wave-specific representation; individual and province-trend models use different identifying variation."),
  section("Supplementary Table 7. Component Checks for the Primary Composite", component_table, "Components with few events are descriptive checks of the composite outcome."),
  section("Supplementary Table 8. Sensitivity Analysis Excluding Outcome Components", stability_table, "Each model reconstructs the primary outcome after omitting the named component or component family."),
  section("Supplementary Table 9. Need-, Access-, and Region-Related Estimates", interaction_table, "The first six rows are formal interaction coefficients from five domains. Region rows are exploratory stratified estimates; the northeast model was not estimable with only three province clusters."),
  section("Supplementary Table 10. Multimorbidity Outcome-Opportunity Checks", opportunity_table, "Diagnosis-eligible models use a fixed condition denominator and adjust for chronic-condition count. Holm adjustment applies across the diagnosis-specific interaction family."))

while (length(lines) > 0 && tail(lines, 1) == "") lines <- head(lines, -1)

write_lines(lines, out_path)
cat(sprintf("Wrote %s\n", out_path))
