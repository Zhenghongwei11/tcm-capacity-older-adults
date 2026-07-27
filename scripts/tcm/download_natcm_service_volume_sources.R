#!/usr/bin/env Rscript

sources <- data.frame(
  source_id = c(
    "natcm_2011_b35_xls",
    "natcm_2013_b38_xls",
    "natcm_2015_b38_xls",
    "natcm_2018_b38_xls"
  ),
  url = c(
    "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2011/B35.xls",
    "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2013/B38.xls",
    "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2015/B38.xls",
    "http://www.natcm.gov.cn/2020tjzb/%E5%85%A8%E5%9B%BD%E4%B8%AD%E5%8C%BB%E8%8D%AF%E7%BB%9F%E8%AE%A1%E6%91%98%E7%BC%96/atog/2018/B38.xls"
  ),
  local_path = c(
    "data/raw/tcm_supply/official_excels/natcm_service_volume/B35_2011.xls",
    "data/raw/tcm_supply/official_excels/natcm_service_volume/B38_2013.xls",
    "data/raw/tcm_supply/official_excels/natcm_service_volume/B38_2015.xls",
    "data/raw/tcm_supply/official_excels/natcm_service_volume/B38_2018.xls"
  ),
  sha256 = c(
    "9210c94db108b5818e1b7dae53e6504aee75663eca10e5aba30d605337443dde",
    "ef58011bd565b6eb8e3e436edadbc10b5dedb62770e64e2033a1d83426288f2b",
    "be944f431181068712fee5ebd075909e57390fc82415ccaa78169f356391eb7e",
    "1bd5595042d70480a0fea18b75ccba24da39e14d8d2b6d5597d814e8f617c49c"
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
