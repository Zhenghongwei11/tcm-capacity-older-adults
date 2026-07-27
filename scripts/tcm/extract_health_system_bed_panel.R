#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)

argval <- function(flag, default) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1]]
}

manual_tsv <- argval(
  "--scanned-transcription",
  "data/metadata/health_system_beds_scanned_transcription.tsv"
)
core_tsv <- argval(
  "--population-panel",
  "results/tcm/province_year_tcm_core_density.tsv"
)
out_tsv <- argval(
  "--output",
  "results/tcm/province_year_health_system_covariates.tsv"
)

dir.create(dirname(out_tsv), recursive = TRUE, showWarnings = FALSE)

province_labels <- tribble(
  ~province, ~short_label,
  "北京市", "北京",
  "天津市", "天津",
  "河北省", "河北",
  "山西省", "山西",
  "内蒙古自治区", "内蒙古",
  "辽宁省", "辽宁",
  "吉林省", "吉林",
  "黑龙江省", "黑龙江",
  "上海市", "上海",
  "江苏省", "江苏",
  "浙江省", "浙江",
  "安徽省", "安徽",
  "福建省", "福建",
  "江西省", "江西",
  "山东省", "山东",
  "河南省", "河南",
  "湖北省", "湖北",
  "湖南省", "湖南",
  "广东省", "广东",
  "广西壮族自治区", "广西",
  "海南省", "海南",
  "重庆市", "重庆",
  "四川省", "四川",
  "贵州省", "贵州",
  "云南省", "云南",
  "西藏自治区", "西藏",
  "陕西省", "陕西",
  "甘肃省", "甘肃",
  "青海省", "青海",
  "宁夏回族自治区", "宁夏",
  "新疆维吾尔自治区", "新疆"
)

sources <- tribble(
  ~year, ~edition, ~pdf_path, ~pdf_page, ~source_id, ~source_url, ~expected_total, ~expected_hospital, ~expected_comprehensive,
  2011L, "2012", "data/raw/tcm_supply/official_pdfs/china_health_statistics_yearbook_2012_nhc_full.pdf", 80L,
  "nhc_yearbook_2012_full_pdf", "https://www.nhc.gov.cn/mohwsbwstjxxzx/tjtjnj/202605/e2ea85f22ed349e78ab5a776dcced4ea.shtml",
  5159889, 3705118, 2670729,
  2013L, "2014", "data/raw/tcm_supply/official_pdfs/china_health_statistics_yearbook_2014_nhc_full.pdf", 76L,
  "nhc_yearbook_2014_full_pdf", "https://www.nhc.gov.cn/mohwsbwstjxxzx/c100228/201506/a2f73d0edf1d4e678eb5978af952a60f.shtml",
  6181891, 4578601, 3255153,
  2015L, "2016", "data/raw/tcm_supply/official_pdfs/china_health_statistics_yearbook_2016_nhc_archived.pdf", 76L,
  "nhc_yearbook_2016_pdf_archived", "https://www.nhc.gov.cn/mohwsbwstjxxzx/tjtjnj/201706/05a4a714a05e44b59a86ff271675c653.shtml",
  7015214, 5330580, 3721036,
  2018L, "2019", "data/raw/tcm_supply/official_pdfs/china_health_statistics_yearbook_2019_nhc.pdf", 78L,
  "nhc_yearbook_2019_pdf", "https://www.nhc.gov.cn/mohwsbwstjxxzx/tjtjnj/202006/456391e2294244c593e871ac254f8aaa.shtml",
  8404078, 6519749, 4378892
)

extract_pdf_text_page <- function(pdf_path, pdf_page) {
  if (!file.exists(pdf_path)) stop("Missing source PDF: ", pdf_path)
  output <- system2(
    "pdftotext",
    c("-f", pdf_page, "-l", pdf_page, "-layout", pdf_path, "-"),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0) stop("pdftotext failed for ", pdf_path)
  output
}

parse_text_layer_table <- function(year, pdf_path, pdf_page) {
  lines <- extract_pdf_text_page(pdf_path, pdf_page)
  normalized <- str_replace_all(lines, "\\s+", "")

  map_dfr(seq_len(nrow(province_labels)), function(i) {
    label <- province_labels$short_label[[i]]
    hit_index <- which(str_starts(normalized, fixed(label)))
    if (length(hit_index) != 1) {
      stop(sprintf("Expected one row for %s in %s, found %s", label, year, length(hit_index)))
    }
    hit <- lines[[hit_index]]
    values <- str_extract_all(hit, "[0-9]+")[[1]] |> as.numeric()
    if (length(values) < 3) stop("Fewer than three numeric fields for ", label, " in ", year)

    tibble(
      year = as.integer(year),
      province = province_labels$province[[i]],
      value_total_health_institution_beds = values[[1]],
      value_hospital_beds = values[[2]],
      value_comprehensive_hospital_beds = values[[3]],
      source_pdf_page = as.integer(pdf_page),
      extraction_method = "official_pdf_text_layer_parsed"
    )
  })
}

manual <- read_tsv(manual_tsv, show_col_types = FALSE) %>%
  mutate(year = as.integer(year))

if (nrow(manual) != 62 || !setequal(unique(manual$year), c(2011L, 2013L))) {
  stop("Scanned transcription must contain 31 rows for both 2011 and 2013")
}

text_layer <- map_dfr(c(2015L, 2018L), function(y) {
  src <- filter(sources, year == y)
  parse_text_layer_table(y, src$pdf_path[[1]], src$pdf_page[[1]])
})

beds <- bind_rows(manual, text_layer) %>%
  arrange(year, match(province, province_labels$province))

if (nrow(beds) != 124 || anyDuplicated(beds[c("province", "year")])) {
  stop("Expected a unique 31-province panel for four years")
}

for (i in seq_len(nrow(sources))) {
  src <- sources[i, ]
  part <- filter(beds, year == src$year)
  checks <- c(
    sum(part$value_total_health_institution_beds) == src$expected_total,
    sum(part$value_hospital_beds) == src$expected_hospital,
    sum(part$value_comprehensive_hospital_beds) == src$expected_comprehensive
  )
  if (!all(checks)) stop("Published national totals do not match province sums for ", src$year)
}

population <- read_tsv(core_tsv, show_col_types = FALSE) %>%
  select(province, year, denominator_population, denominator_source_id) %>%
  mutate(year = as.integer(year))

panel <- beds %>%
  left_join(sources, by = "year") %>%
  left_join(population, by = c("province", "year")) %>%
  mutate(
    source_name = paste0("China Health Statistical Yearbook ", edition),
    table_id = "3-1-3",
    table_title = paste0(year, "年各地区医疗卫生机构床位数"),
    resource_definition = paste(
      "Comprehensive hospitals are a mutually exclusive hospital category in yearbook table 3-1-3;",
      "they are used as a general health-system capacity co-exposure, not as a pure Western-medicine or placebo exposure."
    ),
    denominator_unit = "persons",
    denominator_notes = "Uses the same official NATCM province-year population denominator as the main TCM core-density panel.",
    value_per_10000_population_total_health_institution_beds =
      value_total_health_institution_beds / denominator_population * 10000,
    value_per_10000_population_hospital_beds =
      value_hospital_beds / denominator_population * 10000,
    value_per_10000_population_comprehensive_hospital_beds =
      value_comprehensive_hospital_beds / denominator_population * 10000,
    extraction_status = extraction_method,
    quality_flag = if_else(
      extraction_method == "official_pdf_text_layer_parsed",
      "official_pdf_text_parsed_national_total_verified",
      "scan_transcribed_twice_national_total_verified"
    ),
    unit = "beds",
    notes = paste0(
      "Official yearbook PDF page ", source_pdf_page,
      "; 31 province rows and all three published national totals verified."
    )
  ) %>%
  select(
    province, year, source_id, source_name, source_url, source_edition = edition,
    table_id, table_title, source_pdf_page, extraction_status, quality_flag,
    resource_definition, denominator_population, denominator_unit,
    denominator_source_id, denominator_notes,
    value_total_health_institution_beds, value_hospital_beds,
    value_comprehensive_hospital_beds,
    value_per_10000_population_total_health_institution_beds,
    value_per_10000_population_hospital_beds,
    value_per_10000_population_comprehensive_hospital_beds,
    unit, notes
  ) %>%
  arrange(year, match(province, province_labels$province))

if (anyNA(panel$denominator_population)) stop("Missing population denominator after linkage")

write_tsv(panel, out_tsv)
cat("Wrote ", out_tsv, "\n", sep = "")
print(panel %>% group_by(year) %>% summarise(
  provinces = n(),
  total_beds = sum(value_total_health_institution_beds),
  hospital_beds = sum(value_hospital_beds),
  comprehensive_hospital_beds = sum(value_comprehensive_hospital_beds),
  .groups = "drop"
))
