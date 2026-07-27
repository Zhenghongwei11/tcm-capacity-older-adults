#!/usr/bin/env Rscript

sources <- data.frame(
  source_id = c(
    "natcm_2011_a86_xlsx",
    "natcm_2013_a95_xls",
    "natcm_2015_a95_xls",
    "natcm_2018_a95_xlsx"
  ),
  url = c(
    "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2011/A86.xlsx",
    "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2013/A95.xls",
    "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2015/A95.xls",
    "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2018/A95.xlsx"
  ),
  local_path = c(
    "data/raw/tcm_supply/official_excels/natcm_2011_extract/A86.xlsx",
    "data/raw/tcm_supply/official_excels/natcm_core_density/A95_2013.xls",
    "data/raw/tcm_supply/official_excels/natcm_core_density/A95_2015.xls",
    "data/raw/tcm_supply/official_excels/natcm_core_density/A95_2018.xlsx"
  ),
  sha256 = c(
    "e0972df9e40e7c23fb44627f7662c0221aff5d003abf82cc3320f7891d6031e3",
    "bfff703be5a8fdc79244626dffbe6597617f5e0b3b6343d648ed29711a7bb995",
    "17e724a3a6943e07d1737bdec70eb132afbb4ecfc80d16c9f90ca2e40ec545b1",
    "448dd8218c12ae24cce60292e325a5e4135681942cd88dbc88424f51051cb28a"
  ),
  stringsAsFactors = FALSE
)

download_one <- function(row) {
  destination <- row[["local_path"]]
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(destination)) {
    download.file(row[["url"]], destination, mode = "wb", quiet = TRUE)
  }

  observed <- unname(tools::sha256sum(destination))
  if (!identical(observed, row[["sha256"]])) {
    stop(sprintf(
      "Checksum mismatch for %s\nexpected: %s\nobserved: %s",
      row[["source_id"]], row[["sha256"]], observed
    ))
  }

  data.frame(
    source_id = row[["source_id"]],
    local_path = destination,
    sha256 = observed,
    status = "ready",
    stringsAsFactors = FALSE
  )
}

results <- do.call(
  rbind,
  lapply(seq_len(nrow(sources)), function(i) download_one(sources[i, ]))
)

print(results, row.names = FALSE)
