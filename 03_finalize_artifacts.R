# 03_finalize_artifacts.R — build the submission ZIP
if (!requireNamespace("zip", quietly = TRUE)) install.packages("zip")
dir.create("artifact", showWarnings = FALSE)

files_to_zip <- c(
  list.files(pattern="\\.(R|Rmd|qmd)$"),
  "metadata/targets.tsv",
  list.files("results", recursive = TRUE, full.names = TRUE),
  "renv.lock", "README.md", "LICENSE", "CITATION.cff"
)
files_to_zip <- files_to_zip[file.exists(files_to_zip)]
zip::zip(zipfile = "artifact/code_and_results_v1.0.zip", files = files_to_zip)
