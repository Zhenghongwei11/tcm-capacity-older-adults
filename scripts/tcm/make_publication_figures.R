#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggalluvial)
  library(ggplot2)
  library(ggrepel)
  library(ggridges)
  library(grid)
  library(patchwork)
  library(readr)
  library(scales)
  library(stringr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

argval <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

out_dir <- argval("--out-dir", "plots/publication")
source_dir <- argval("--source-dir", "results/tcm/figure_source_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

pal <- c(
  ink = "#24282B",
  muted = "#6B747A",
  grid = "#E4E8E8",
  neutral = "#BCC5C7",
  neutral_light = "#E9ECEB",
  supply = "#2F7D78",
  supply_dark = "#205D59",
  supply_light = "#C8DEDB",
  use = "#C9663B",
  use_dark = "#8F4228",
  use_light = "#F0D2C4",
  reference = "#708F9E",
  reference_light = "#D9E3E7"
)

trajectory_colors <- c(
  "Never" = "#D7DCDA",
  "Intermittent" = "#9CB6B5",
  "Incident and sustained" = "#D98A5B",
  "Discontinued" = "#7694A3",
  "Persistent" = "#9B452B",
  "No use" = "#F1F3F2",
  "TCM use" = "#E7B59C"
)

theme_paper <- function(base_size = 6.7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = pal[["ink"]]),
      axis.ticks = element_line(linewidth = 0.3, colour = pal[["ink"]]),
      axis.ticks.length = unit(1.4, "mm"),
      axis.title = element_text(size = base_size, colour = pal[["ink"]]),
      axis.text = element_text(size = base_size - 0.5, colour = pal[["muted"]]),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.8, colour = pal[["ink"]]),
      plot.title = element_text(size = base_size + 0.7, face = "bold", colour = pal[["ink"]], margin = margin(b = 2.5)),
      plot.subtitle = element_text(size = base_size - 0.4, colour = pal[["muted"]], margin = margin(b = 4)),
      plot.tag.position = "topleft",
      plot.tag = element_text(size = 8, face = "bold", colour = pal[["ink"]], margin = margin(r = 4, b = 2)),
      panel.grid = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", colour = pal[["ink"]]),
      plot.margin = margin(4, 5, 4, 5)
    )
}

save_pub <- function(plot, stem, width_mm = 183, height_mm = 115, dpi = 600) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svg_path <- file.path(out_dir, paste0(stem, ".svg"))
  pdf_path <- file.path(out_dir, paste0(stem, ".pdf"))
  tiff_path <- file.path(out_dir, paste0(stem, ".tiff"))
  png_path <- file.path(out_dir, paste0(stem, ".png"))

  svglite::svglite(svg_path, width = width_in, height = height_in)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in, family = "Arial")
  print(plot)
  dev.off()

  ragg::agg_tiff(tiff_path, width = width_in, height = height_in, units = "in", res = dpi, compression = "lzw")
  print(plot)
  dev.off()

  ragg::agg_png(png_path, width = width_in, height = height_in, units = "in", res = 300, background = "white")
  print(plot)
  dev.off()
}

cluster_vcov <- function(model, cluster) {
  ok <- !is.na(cluster)
  cluster <- as.factor(cluster[ok])
  x <- model.matrix(model)[ok, , drop = FALSE]
  u <- residuals(model)[ok]
  xtx_inv <- solve(crossprod(x))
  meat <- matrix(0, ncol(x), ncol(x))
  for (g in levels(cluster)) {
    idx <- cluster == g
    score <- crossprod(x[idx, , drop = FALSE], u[idx])
    meat <- meat + tcrossprod(score)
  }
  n <- nrow(x)
  k <- ncol(x)
  groups <- nlevels(cluster)
  correction <- (groups / (groups - 1)) * ((n - 1) / (n - k))
  correction * xtx_inv %*% meat %*% xtx_inv
}

province_names <- c(
  "北京市" = "Beijing", "天津市" = "Tianjin", "河北省" = "Hebei",
  "山西省" = "Shanxi", "内蒙古自治区" = "Inner Mongolia", "辽宁省" = "Liaoning",
  "吉林省" = "Jilin", "黑龙江省" = "Heilongjiang", "上海市" = "Shanghai",
  "江苏省" = "Jiangsu", "浙江省" = "Zhejiang", "安徽省" = "Anhui",
  "福建省" = "Fujian", "江西省" = "Jiangxi", "山东省" = "Shandong",
  "河南省" = "Henan", "湖北省" = "Hubei", "湖南省" = "Hunan",
  "广东省" = "Guangdong", "广西壮族自治区" = "Guangxi", "海南省" = "Hainan",
  "重庆市" = "Chongqing", "四川省" = "Sichuan", "贵州省" = "Guizhou",
  "云南省" = "Yunnan", "西藏自治区" = "Tibet", "陕西省" = "Shaanxi",
  "甘肃省" = "Gansu", "青海省" = "Qinghai", "宁夏回族自治区" = "Ningxia",
  "新疆维吾尔自治区" = "Xinjiang"
)

supply <- read_tsv("results/tcm/province_year_tcm_core_density.tsv", show_col_types = FALSE)
panel_raw <- read_tsv("results/tcm/person_wave_tcm_core_density_analysis.tsv", show_col_types = FALSE)
main_models <- read_tsv("results/tcm/tcm_supply_main_models.tsv", show_col_types = FALSE)
robustness <- read_tsv("results/tcm/tcm_supply_robustness.tsv", show_col_types = FALSE)
cluster_checks <- read_tsv("results/tcm/tcm_supply_cluster_sensitivity.tsv", show_col_types = FALSE)
interactions <- read_tsv("results/tcm/tcm_supply_equity_interactions.tsv", show_col_types = FALSE)
weighted_models <- read_tsv("results/tcm/tcm_supply_weighted_attrition_models.tsv", show_col_types = FALSE)
identification_models <- read_tsv("results/tcm/tcm_supply_longitudinal_identification_models.tsv", show_col_types = FALSE)

# Figure 1: individual use streams and provincial capacity trajectories.
age60_use <- panel_raw %>%
  filter(main_model_age60, !is.na(primary_condition_tcm_any)) %>%
  distinct(harmonized_charls_id, year, .keep_all = TRUE)

balanced_ids <- age60_use %>%
  count(harmonized_charls_id, name = "waves") %>%
  filter(waves == 4)

sequence_counts <- age60_use %>%
  semi_join(balanced_ids, by = "harmonized_charls_id") %>%
  arrange(harmonized_charls_id, year) %>%
  group_by(harmonized_charls_id) %>%
  summarise(
    sequence = paste0(primary_condition_tcm_any, collapse = ""),
    uses = sum(primary_condition_tcm_any),
    .groups = "drop"
  ) %>%
  mutate(
    trajectory = case_when(
      uses == 0 ~ "Never",
      uses == 4 ~ "Persistent",
      sequence %in% c("0001", "0011", "0111") ~ "Incident and sustained",
      sequence %in% c("1000", "1100", "1110") ~ "Discontinued",
      TRUE ~ "Intermittent"
    ),
    trajectory = factor(
      trajectory,
      levels = c("Never", "Intermittent", "Incident and sustained", "Discontinued", "Persistent")
    )
  ) %>%
  count(sequence, trajectory, name = "n") %>%
  mutate(percent = 100 * n / sum(n))

sequence_long <- sequence_counts %>%
  crossing(wave_index = 1:4) %>%
  mutate(
    year = factor(c(2011, 2013, 2015, 2018)[wave_index], levels = c(2011, 2013, 2015, 2018)),
    state = if_else(substr(sequence, wave_index, wave_index) == "1", "TCM use", "No use"),
    state = factor(state, levels = c("No use", "TCM use"))
  )

write_tsv(
  sequence_counts %>% select(sequence, trajectory, n, percent),
  file.path(source_dir, "fig1_individual_use_trajectories.tsv")
)

fig1a <- ggplot(
  sequence_long,
  aes(x = year, stratum = state, alluvium = sequence, y = n, fill = trajectory)
) +
  geom_alluvium(width = 0.18, alpha = 0.82, colour = NA, knot.pos = 0.42) +
  geom_stratum(aes(fill = state), width = 0.18, colour = "white", linewidth = 0.3) +
  stat_stratum(aes(label = after_stat(stratum)), geom = "text", size = 2.0, colour = pal[["ink"]], lineheight = 0.9) +
  scale_fill_manual(
    values = trajectory_colors,
    breaks = levels(sequence_counts$trajectory),
    drop = FALSE
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
  labs(
    title = "Individual Chinese medicine use trajectories",
    subtitle = sprintf("Balanced four-wave panel, age 60+ (n = %s)", comma(nrow(balanced_ids))),
    x = NULL,
    y = NULL
  ) +
  theme_paper(6.5) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(3.5, "mm"),
    legend.spacing.x = unit(1.5, "mm")
  ) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(alpha = 1)))

supply_plot <- supply %>%
  transmute(
    province,
    province_label = recode(province, !!!province_names),
    year,
    beds = value_per_10000_population_tcm_hospital_beds,
    physicians = value_per_10000_population_tcm_practicing_assistant_physicians
  ) %>%
  arrange(province, year)

endpoints <- supply_plot %>%
  filter(year == 2018) %>%
  mutate(
    endpoint_score = as.numeric(scale(beds)) + as.numeric(scale(physicians)),
    endpoint_rank = rank(endpoint_score, ties.method = "first"),
    label_endpoint = endpoint_rank <= 2 | endpoint_rank >= n() - 1
  )

mean_path <- supply_plot %>%
  group_by(year) %>%
  summarise(beds = mean(beds), physicians = mean(physicians), .groups = "drop")

write_tsv(supply_plot, file.path(source_dir, "fig1_provincial_capacity_trajectories.tsv"))

fig1b <- ggplot(supply_plot, aes(beds, physicians, group = province)) +
  geom_path(
    colour = pal[["neutral"]], linewidth = 0.36, alpha = 0.58,
    arrow = arrow(type = "closed", length = unit(1.4, "mm"))
  ) +
  geom_point(
    data = filter(supply_plot, year %in% c(2013, 2015)),
    size = 0.85, shape = 21, fill = "white", colour = pal[["neutral"]], stroke = 0.28
  ) +
  geom_point(
    data = filter(supply_plot, year == 2011),
    size = 1.6, shape = 21, fill = "white", colour = pal[["supply_dark"]], stroke = 0.42
  ) +
  geom_point(
    data = filter(supply_plot, year == 2018),
    size = 1.8, shape = 21, fill = pal[["supply"]], colour = "white", stroke = 0.35
  ) +
  geom_path(
    data = mean_path,
    aes(beds, physicians, group = 1),
    colour = pal[["supply_dark"]], linewidth = 0.95,
    arrow = arrow(type = "closed", length = unit(2.0, "mm"))
  ) +
  geom_point(
    data = mean_path, aes(beds, physicians), inherit.aes = FALSE,
    shape = 21, fill = pal[["supply_light"]], colour = pal[["supply_dark"]], size = 1.8, stroke = 0.4
  ) +
  geom_text_repel(
    data = filter(endpoints, label_endpoint),
    aes(label = province_label),
    size = 1.9, colour = pal[["ink"]], min.segment.length = 0,
    segment.colour = pal[["neutral"]], segment.size = 0.25,
    box.padding = 0.25, point.padding = 0.12, max.overlaps = Inf, seed = 20260720
  ) +
  annotate("text", x = mean_path$beds[1] - 0.12, y = mean_path$physicians[1] + 0.18, label = "2011", hjust = 1, size = 1.9, colour = pal[["supply_dark"]]) +
  annotate("text", x = mean_path$beds[4] + 0.12, y = mean_path$physicians[4] - 0.20, label = "2018 mean", hjust = 0, size = 1.9, colour = pal[["supply_dark"]]) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(
    title = "Provincial capacity trajectories",
    subtitle = "All 31 provinces, 2011-2018",
    x = "TCM hospital beds per 10,000",
    y = "TCM physicians per 10,000"
  ) +
  theme_paper(6.5)

fig1 <- fig1a + fig1b +
  plot_layout(widths = c(1.45, 1), guides = "keep") +
  plot_annotation(tag_levels = "A")

save_pub(fig1, "fig1_two_longitudinal_systems", height_mm = 112)

# Figure 2: exact partial association from the primary model and province influence.
analysis <- panel_raw %>%
  filter(main_model_age60) %>%
  mutate(
    province_fe = factor(province_supply_key),
    wave_fe = factor(year),
    education_group = factor(education_group),
    household_income_quartile = factor(household_income_quartile),
    z_beds = z_py_tcm_beds_per_10000,
    outcome = as.numeric(primary_condition_tcm_any)
  )

controls <- c(
  "age_model", "I(age_model^2)", "female", "education_group",
  "married_or_partnered", "rural_community", "agricultural_hukou",
  "public_insurance", "chronic_count", "any_adl_limitation",
  "poor_self_rated_health", "household_income_quartile", "province_fe", "wave_fe"
)

main_formula <- as.formula(paste("outcome ~ z_beds +", paste(controls, collapse = " + ")))
needed <- all.vars(main_formula)
model_data <- analysis %>% filter(if_all(all_of(needed), ~ !is.na(.x)))
main_fit <- lm(main_formula, data = model_data)
main_vcov <- cluster_vcov(main_fit, model_data$province_supply_key)
main_beta <- unname(coef(main_fit)[["z_beds"]])
main_se <- sqrt(main_vcov["z_beds", "z_beds"])
main_df <- n_distinct(model_data$province_supply_key) - 1
main_ci <- main_beta + c(-1, 1) * qt(0.975, df = main_df) * main_se

y_control <- lm(as.formula(paste("outcome ~", paste(controls, collapse = " + "))), data = model_data)
x_control <- lm(as.formula(paste("z_beds ~", paste(controls, collapse = " + "))), data = model_data)

partial <- model_data %>%
  transmute(
    province = province_supply_key,
    x_residual = residuals(x_control),
    y_residual_pp = 100 * residuals(y_control)
  ) %>%
  mutate(bin = ntile(x_residual, 10))

clustered_bin <- function(df) {
  n_obs <- nrow(df)
  mean_y <- mean(df$y_residual_pp)
  scores <- df %>%
    mutate(centered = y_residual_pp - mean_y) %>%
    group_by(province) %>%
    summarise(score = sum(centered), .groups = "drop")
  groups <- nrow(scores)
  se <- sqrt((groups / (groups - 1)) * sum(scores$score^2) / n_obs^2)
  tibble(
    x = mean(df$x_residual), y = mean_y, n = n_obs, provinces = groups,
    se = se, lower = mean_y - qt(0.975, groups - 1) * se,
    upper = mean_y + qt(0.975, groups - 1) * se
  )
}

partial_bins <- partial %>%
  group_by(bin) %>%
  group_modify(~ clustered_bin(.x)) %>%
  ungroup()

fan <- tibble(x = seq(min(partial_bins$x) - 0.08, max(partial_bins$x) + 0.08, length.out = 200)) %>%
  mutate(
    fit = 100 * main_beta * x,
    edge_1 = 100 * main_ci[[1]] * x,
    edge_2 = 100 * main_ci[[2]] * x,
    lower = pmin(edge_1, edge_2),
    upper = pmax(edge_1, edge_2)
  )

write_tsv(partial_bins, file.path(source_dir, "fig2_partial_association_bins.tsv"))

fig2a <- ggplot() +
  geom_ribbon(data = fan, aes(x, ymin = lower, ymax = upper), fill = pal[["use_light"]], alpha = 0.78) +
  geom_line(data = fan, aes(x, fit), colour = pal[["use_dark"]], linewidth = 0.75) +
  geom_errorbar(
    data = partial_bins,
    aes(x = x, ymin = lower, ymax = upper),
    width = 0, linewidth = 0.38, colour = pal[["reference"]]
  ) +
  geom_point(
    data = partial_bins,
    aes(x, y, size = n), shape = 21, fill = "white",
    colour = pal[["reference"]], stroke = 0.55
  ) +
  geom_hline(yintercept = 0, colour = pal[["muted"]], linewidth = 0.28) +
  geom_vline(xintercept = 0, colour = pal[["grid"]], linewidth = 0.28) +
  annotate(
    "label", x = Inf, y = Inf,
    label = sprintf("+%.2f pp per SD\n95%% CI %.2f to %.2f", 100 * main_beta, 100 * main_ci[[1]], 100 * main_ci[[2]]),
    hjust = 1.04, vjust = 1.12, size = 2.05,
    fill = alpha("white", 0.88), colour = pal[["ink"]], linewidth = 0,
    label.padding = unit(1.4, "mm")
  ) +
  scale_size_continuous(range = c(1.7, 2.8), guide = "none") +
  labs(
    title = "Adjusted within-design association",
    subtitle = "Equal-count bins of model partial residuals; intervals clustered by province",
    x = "Adjusted TCM bed-density residual (province-year SD units)",
    y = "Adjusted difference in treatment use (percentage points)"
  ) +
  theme_paper(6.6)

loo <- cluster_checks %>%
  filter(check_type == "leave_one_province_out") %>%
  transmute(
    omitted_province = detail,
    omitted_label = recode(detail, !!!province_names),
    effect
  ) %>%
  arrange(effect) %>%
  mutate(label_extreme = row_number() %in% c(1, n()))

main_cluster <- cluster_checks %>% filter(check_type == "main_cluster_robust") %>% slice(1)
write_tsv(loo, file.path(source_dir, "fig2_leave_one_province_out.tsv"))

fig2b <- ggplot(loo, aes(effect)) +
  geom_dotplot(
    binaxis = "x", stackdir = "center", binwidth = 0.075,
    dotsize = 0.72, stackratio = 1.02,
    fill = pal[["supply"]], colour = "white", stroke = 0.3
  ) +
  annotate(
    "segment", x = main_cluster$ci_lower, xend = main_cluster$ci_upper,
    y = 0.72, yend = 0.72, linewidth = 0.65, colour = pal[["use_dark"]]
  ) +
  annotate(
    "point", x = main_cluster$effect, y = 0.72,
    shape = 21, size = 2.25, fill = pal[["use"]], colour = "white", stroke = 0.35
  ) +
  geom_text(
    data = filter(loo, label_extreme),
    aes(x = effect, y = -0.22, label = omitted_label),
    inherit.aes = FALSE, size = 1.9, colour = pal[["muted"]]
  ) +
  annotate("text", x = main_cluster$effect, y = 0.88, label = "Full sample", size = 1.9, colour = pal[["use_dark"]]) +
  scale_x_continuous(breaks = c(0.5, 1.5, 2.5, 3.5)) +
  scale_y_continuous(breaks = NULL) +
  coord_cartesian(xlim = c(0.05, 3.75), ylim = c(-0.42, 1.02), clip = "off") +
  labs(
    title = "Province influence",
    subtitle = "28 leave-one-province-out estimates",
    x = "Percentage-point difference per SD",
    y = NULL
  ) +
  theme_paper(6.6) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    panel.grid = element_blank()
  )

fig2 <- fig2a + fig2b +
  plot_layout(widths = c(1.75, 1)) +
  plot_annotation(tag_levels = "A")

save_pub(fig2, "fig2_supply_to_realized_use", height_mm = 96)

# Figure 3: formal interaction contrasts shown as response-path differences.
contrast_lookup <- tribble(
  ~stratum, ~contrast, ~panel_title, ~reference_text,
  "multichronic_stratum", "multichronic_stratummultimorbidity", "Multimorbidity", "versus no multimorbidity",
  "rural_stratum", "rural_stratumrural", "Rural residence", "versus urban residence",
  "education_stratum", "education_stratumlower_education", "Lower education", "versus higher education",
  "income_stratum", "income_stratumlower_income", "Lower income", "versus higher income",
  "adl_stratum", "adl_stratumadl_limitation", "ADL limitation", "versus no ADL limitation"
)

contrast_data <- interactions %>%
  inner_join(contrast_lookup, by = c("stratum", "contrast")) %>%
  mutate(
    panel_title = factor(panel_title, levels = contrast_lookup$panel_title),
    p_label = if_else(pvalue < 0.001, "P interaction < 0.001", sprintf("P interaction = %.3f", pvalue))
  )

contrast_paths <- contrast_data %>%
  select(panel_title, reference_text, interaction_effect, ci_lower, ci_upper, p_label) %>%
  crossing(x = seq(-1, 1, length.out = 101)) %>%
  mutate(
    fit = interaction_effect * x,
    edge_1 = ci_lower * x,
    edge_2 = ci_upper * x,
    lower = pmin(edge_1, edge_2),
    upper = pmax(edge_1, edge_2)
  )

write_tsv(
  contrast_data %>% select(panel_title, reference_text, interaction_effect, ci_lower, ci_upper, pvalue, n_persons, n_person_waves, n_clusters),
  file.path(source_dir, "fig3_formal_interaction_contrasts.tsv")
)

path_plot <- function(title, hero = FALSE) {
  d <- contrast_paths %>% filter(panel_title == title)
  s <- contrast_data %>% filter(panel_title == title)
  ymax <- 5.4
  ggplot(d, aes(x, fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = pal[["use_light"]], alpha = if (hero) 0.85 else 0.65) +
    geom_hline(yintercept = 0, linewidth = 0.36, colour = pal[["reference"]]) +
    geom_line(linewidth = if (hero) 0.9 else 0.7, colour = pal[["use_dark"]]) +
    geom_point(data = filter(d, x %in% c(-1, 1)), size = if (hero) 2.0 else 1.6, shape = 21, fill = pal[["use"]], colour = "white", stroke = 0.3) +
    annotate(
      "text", x = -0.96, y = ymax - 0.35, label = s$reference_text,
      hjust = 0, vjust = 1, size = if (hero) 2.2 else 1.75, colour = pal[["muted"]]
    ) +
    annotate(
      "text", x = 0.96, y = -ymax + 0.35,
      label = sprintf("%+.2f pp per SD\n%s", s$interaction_effect, s$p_label),
      hjust = 1, vjust = 0, size = if (hero) 2.2 else 1.75,
      colour = if (s$pvalue < 0.05) pal[["use_dark"]] else pal[["muted"]], lineheight = 0.95
    ) +
    scale_x_continuous(breaks = c(-1, 0, 1), labels = c("-1 SD", "Mean", "+1 SD"), expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = c(-ymax, ymax), breaks = c(-4, 0, 4)) +
    labs(
      title = title,
      x = if (hero) "TCM hospital bed density" else NULL,
      y = if (hero) "Additional difference in use\n(percentage points)" else NULL
    ) +
    theme_paper(if (hero) 6.8 else 5.9) +
    theme(
      axis.text.y = if (hero) element_text() else element_blank(),
      axis.ticks.y = if (hero) element_line(linewidth = 0.3) else element_blank(),
      axis.line.y = if (hero) element_line(linewidth = 0.3) else element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

p_multi <- path_plot("Multimorbidity", hero = TRUE)
p_rural <- path_plot("Rural residence")
p_education <- path_plot("Lower education")
p_income <- path_plot("Lower income")
p_adl <- path_plot("ADL limitation")

fig3_design <- "
AAAA
BCDE
"

fig3 <- p_multi + p_rural + p_education + p_income + p_adl +
  plot_layout(design = fig3_design, heights = c(1.25, 1)) +
  plot_annotation(
    title = "Interaction contrasts across prespecified groups",
    subtitle = "Lines show the additional difference relative to each reference group; shaded fans are 95% confidence intervals",
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(family = "Arial", size = 8, face = "bold", colour = pal[["ink"]]),
      plot.subtitle = element_text(family = "Arial", size = 6.2, colour = pal[["muted"]], margin = margin(b = 3))
    )
  )

save_pub(fig3, "fig3_need_and_equity_paths", height_mm = 118)

# Supplementary Figure 1: estimates plus an explicit specification matrix.
pick_main <- function(indicator, model_id) {
  main_models %>%
    filter(supply_indicator == indicator, model == model_id) %>%
    slice(1)
}

bed_name <- "TCM hospital beds per 10,000 population"
physician_name <- "TCM physicians per 10,000 population"

specs <- bind_rows(
  pick_main(bed_name, "M2_covariate_adjusted_lpm") %>% mutate(spec = "S1", label = "Primary bed-density model", exposure = "Bed density", population = "Age 60+", outcome_group = "Primary outcome", timing = "Concurrent", adjustment = "Full adjustment", uncertainty = "Province clustered"),
  pick_main(physician_name, "M2_covariate_adjusted_lpm") %>% mutate(spec = "S2", label = "Physician-density indicator", exposure = "Physician density", population = "Age 60+", outcome_group = "Primary outcome", timing = "Concurrent", adjustment = "Full adjustment", uncertainty = "Province clustered"),
  pick_main(bed_name, "M1_minimal_adjusted_lpm") %>% mutate(spec = "S3", label = "Minimal adjustment", exposure = "Bed density", population = "Age 60+", outcome_group = "Primary outcome", timing = "Concurrent", adjustment = "Minimal adjustment", uncertainty = "Province clustered"),
  robustness %>% filter(check_type == "age_threshold") %>% slice(1) %>% mutate(spec = "S4", label = "Age 45+", exposure = "Bed density", population = "Age 45+", outcome_group = "Primary outcome", timing = "Concurrent", adjustment = "Full adjustment", uncertainty = "Province clustered"),
  robustness %>% filter(check_type == "covariate_set") %>% slice(1) %>% mutate(spec = "S5", label = "Income omitted", exposure = "Bed density", population = "Age 60+", outcome_group = "Primary outcome", timing = "Concurrent", adjustment = "No income covariate", uncertainty = "Province clustered"),
  robustness %>% filter(check_type == "lag_structure", str_detect(supply_indicator, "beds")) %>% slice(1) %>% mutate(spec = "S6", label = "Lagged bed density", exposure = "Bed density", population = "Age 60+", outcome_group = "Primary outcome", timing = "Lagged", adjustment = "Full adjustment", uncertainty = "Province clustered"),
  robustness %>% filter(check_type == "alternative_outcome", str_detect(outcome, "hospital visit")) %>% slice(1) %>% mutate(spec = "S7", label = "Strict hospital visit", exposure = "Bed density", population = "Age 60+", outcome_group = "Strict outcome", timing = "Concurrent", adjustment = "Full adjustment", uncertainty = "Province clustered"),
  robustness %>% filter(check_type == "alternative_outcome", str_detect(outcome, "Broader")) %>% slice(1) %>% mutate(spec = "S8", label = "Broader use", exposure = "Bed density", population = "Age 60+", outcome_group = "Broader outcome", timing = "Concurrent", adjustment = "Full adjustment", uncertainty = "Province clustered"),
  cluster_checks %>% filter(check_type == "province_cluster_bootstrap") %>% slice(1) %>% mutate(spec = "S9", label = "Province resampling", exposure = "Bed density", population = "Age 60+", outcome_group = "Primary outcome", timing = "Concurrent", adjustment = "Full adjustment", uncertainty = "Province resampling"),
  weighted_models %>% filter(model == "W1_adjusted_cross_sectional_weight") %>% slice(1) %>% mutate(spec = "S10", label = "Adjusted respondent weight", exposure = "Bed density", population = "Age 60+", outcome_group = "Primary outcome", timing = "Concurrent", adjustment = "Full adjustment", uncertainty = "Adjusted respondent weight"),
  identification_models %>% filter(model == "I1_individual_fixed_effects") %>% slice(1) %>% mutate(spec = "S11", label = "Individual fixed effects", exposure = "Bed density", population = "Age 60+", outcome_group = "Primary outcome", timing = "Concurrent", adjustment = "Individual FE", uncertainty = "Province clustered"),
  identification_models %>% filter(model == "I2_mundlak_correlated_effects") %>% slice(1) %>% mutate(spec = "S12", label = "Mundlak correlated effects", exposure = "Bed density", population = "Age 60+", outcome_group = "Primary outcome", timing = "Concurrent", adjustment = "Mundlak", uncertainty = "Province clustered"),
  identification_models %>% filter(model == "I3_province_linear_trends") %>% slice(1) %>% mutate(spec = "S13", label = "Province-specific trends", exposure = "Bed density", population = "Age 60+", outcome_group = "Primary outcome", timing = "Concurrent", adjustment = "Province trends", uncertainty = "Province clustered")
) %>%
  transmute(spec, label, effect, ci_lower, ci_upper, exposure, population, outcome_group, timing, adjustment, uncertainty) %>%
  mutate(spec = factor(spec, levels = paste0("S", 1:13)), supported = ci_lower > 0)

feature_order <- c(
  "Bed density", "Physician density", "Age 60+", "Age 45+",
  "Primary outcome", "Strict outcome", "Broader outcome",
  "Concurrent", "Lagged", "Full adjustment", "Minimal adjustment",
  "No income covariate", "Individual FE", "Mundlak", "Province trends",
  "Province clustered", "Province resampling", "Adjusted respondent weight"
)

spec_matrix <- specs %>%
  select(spec, exposure, population, outcome_group, timing, adjustment, uncertainty) %>%
  pivot_longer(-spec, names_to = "family", values_to = "active_feature") %>%
  select(spec, active_feature) %>%
  right_join(
    expand_grid(spec = factor(paste0("S", 1:13), levels = paste0("S", 1:13)), feature = feature_order),
    by = "spec", relationship = "many-to-many"
  ) %>%
  mutate(active = feature == active_feature) %>%
  group_by(spec, feature) %>%
  summarise(active = any(active), .groups = "drop") %>%
  mutate(feature = factor(feature, levels = rev(feature_order)))

write_tsv(specs %>% mutate(spec = as.character(spec)), file.path(source_dir, "figs1_model_multiverse_estimates.tsv"))
write_tsv(spec_matrix %>% mutate(spec = as.character(spec), feature = as.character(feature)), file.path(source_dir, "figs1_model_multiverse_matrix.tsv"))

figs1a <- ggplot(specs, aes(spec, effect)) +
  geom_hline(yintercept = 0, linewidth = 0.32, colour = pal[["muted"]]) +
  geom_linerange(aes(ymin = ci_lower, ymax = ci_upper), linewidth = 0.46, colour = pal[["muted"]]) +
  geom_point(aes(fill = supported), shape = 21, size = 2.3, stroke = 0.38, colour = pal[["ink"]]) +
  geom_text(aes(label = sprintf("%.2f", effect)), nudge_y = 0.48, size = 1.9, colour = pal[["muted"]]) +
  scale_fill_manual(values = c("TRUE" = pal[["supply"]], "FALSE" = "white"), guide = "none") +
  scale_y_continuous(limits = c(-4.5, 5.5), breaks = c(-4, -2, 0, 2, 4)) +
  labs(
    title = "Sensitivity estimates",
    subtitle = "Percentage-point difference per 1-SD higher province-year supply",
    x = NULL,
    y = "Estimate (95% CI)"
  ) +
  theme_paper(6.4) +
  theme(axis.text.x = element_text(face = "bold", colour = pal[["ink"]]))

figs1b <- ggplot(spec_matrix, aes(spec, feature, fill = active)) +
  geom_tile(width = 0.82, height = 0.78, linewidth = 0.22, colour = pal[["grid"]]) +
  scale_fill_manual(values = c("TRUE" = pal[["supply_dark"]], "FALSE" = "white"), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_paper(5.8) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_blank(), axis.text.y = element_text(colour = pal[["ink"]]),
    plot.margin = margin(0, 5, 4, 5)
  )

figs1 <- figs1a / figs1b + plot_layout(heights = c(1.15, 1.35)) + plot_annotation(tag_levels = "A")
save_pub(figs1, "figs1_model_multiverse", height_mm = 150)

# Supplementary Figure 2: outcome availability, event count, and event rate.
outcome_long <- panel_raw %>%
  filter(main_model_age60) %>%
  transmute(
    year,
    `Primary treatment outcome` = as.numeric(primary_condition_tcm_any),
    `Strict TCM hospital visit` = as.numeric(strict_tcm_hospital_visit),
    `Broader TCM-related use` = as.numeric(broader_tcm_use_2011_2015)
  ) %>%
  pivot_longer(-year, names_to = "outcome", values_to = "value") %>%
  group_by(year, outcome) %>%
  summarise(
    respondents = sum(!is.na(value)),
    events = sum(value == 1, na.rm = TRUE),
    event_rate = if_else(respondents > 0, 100 * events / respondents, NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    year = factor(year, levels = c(2011, 2013, 2015, 2018)),
    outcome = factor(outcome, levels = rev(c("Primary treatment outcome", "Strict TCM hospital visit", "Broader TCM-related use"))),
    available = respondents > 0,
    rate_label = if_else(available, sprintf("%.1f%%", event_rate), "Not available")
  )

write_tsv(
  outcome_long %>% mutate(year = as.character(year), outcome = as.character(outcome)),
  file.path(source_dir, "figs2_outcome_observability.tsv")
)

available_outcomes <- outcome_long %>% filter(available)
missing_outcomes <- outcome_long %>% filter(!available)

figs2 <- ggplot() +
  geom_tile(
    data = outcome_long,
    aes(year, outcome), width = 0.82, height = 0.72,
    fill = "white", colour = pal[["grid"]], linewidth = 0.32
  ) +
  geom_point(
    data = available_outcomes,
    aes(year, outcome, size = events, fill = event_rate),
    shape = 21, colour = "white", stroke = 0.42
  ) +
  geom_text(
    data = available_outcomes,
    aes(year, outcome, label = rate_label),
    size = 2.15, colour = pal[["ink"]]
  ) +
  geom_text(
    data = missing_outcomes,
    aes(year, outcome, label = "Not available"),
    size = 2.05, colour = pal[["muted"]]
  ) +
  scale_size_area(max_size = 18, breaks = c(100, 1000, 3000), labels = comma) +
  scale_fill_gradient(low = pal[["reference_light"]], high = pal[["use"]], limits = c(0, max(available_outcomes$event_rate))) +
  guides(
    size = guide_legend(title = "Events", order = 1, override.aes = list(fill = pal[["neutral_light"]], colour = pal[["muted"]])),
    fill = guide_colorbar(title = "Event rate (%)", order = 2, barwidth = unit(28, "mm"), barheight = unit(2.4, "mm"))
  ) +
  labs(
    title = "Outcome observability across the longitudinal panel",
    subtitle = "Circle area represents event count; fill represents event rate among respondents aged 60+",
    x = "Survey wave",
    y = NULL
  ) +
  theme_paper(7.0) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.y = element_text(face = "bold", colour = pal[["ink"]]),
    legend.position = "bottom", legend.box = "horizontal",
    panel.grid = element_blank()
  )

save_pub(figs2, "figs2_outcome_observability", height_mm = 88)

cat(sprintf("Wrote publication figures to %s\n", out_dir))
cat(sprintf("Wrote aggregated figure source data to %s\n", source_dir))
