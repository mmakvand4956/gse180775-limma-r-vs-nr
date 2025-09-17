# 01_download.R — download & unpack GSE180775 raw files
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
options(timeout = max(3600L, getOption("timeout")))
dir.create("data/raw/GSE180775", recursive = TRUE, showWarnings = FALSE)

ok <- FALSE
try({
  if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery", ask = FALSE, update = FALSE)
  library(GEOquery)
  getGEOSuppFiles("GSE180775", baseDir = "data/raw", makeDirectory = TRUE, fetch_files = TRUE)
  ok <- file.exists("data/raw/GSE180775/GSE180775_RAW.tar")
}, silent = TRUE)

if (!ok) {
  if (!requireNamespace("curl", quietly = TRUE)) install.packages("curl")
  curl::curl_download(
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE180775&format=file",
    destfile = "data/raw/GSE180775/GSE180775_RAW.tar", mode = "wb"
  )
  ok <- file.exists("data/raw/GSE180775/GSE180775_RAW.tar")
}
stopifnot(ok)
utils::untar("data/raw/GSE180775/GSE180775_RAW.tar", exdir = "data/raw/GSE180775")
