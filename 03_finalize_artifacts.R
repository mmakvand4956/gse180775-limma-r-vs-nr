 # 03_finalize_artifacts.R — make PCA_samples.png and the ZIP artifact
> options(stringsAsFactors = FALSE)
> 
> dir.create("results",  showWarnings = FALSE)
> dir.create("artifact", showWarnings = FALSE)
> 
> # --- Load normalized matrix ----------------------------------------------------
> E <- read.delim("results/normalized_log2_expression.tsv",
+                 check.names = FALSE, row.names = 1)
> E <- as.matrix(E)
> E <- E[is.finite(rowSums(E)), , drop = FALSE]
> 
> # --- Load targets and map sample -> group -------------------------------------
> targets <- read.delim("metadata/targets.tsv", check.names = FALSE)
> 
> # Normalize column names from older sheets if needed
> if (!"Group" %in% names(targets) && "Type" %in% names(targets)) targets$Group <- targets$Type
> if (!"File"  %in% names(targets)) {
+   if ("FileName" %in% names(targets)) targets$File <- targets$FileName else targets$File <- NA_character_
+ }
> 
> # 1) primary mapping by basename of raw file
> map1 <- setNames(targets$Group, basename(targets$File))
> grp  <- unname(map1[colnames(E)])
> 
> # 2) fallback mapping by GSM if needed
> if (any(is.na(grp)) && "GSM" %in% names(targets)) {
+   gsm_in_col <- regmatches(colnames(E), regexpr("GSM\\d+", colnames(E)))
+   map2 <- setNames(targets$Group, targets$GSM)
+   grp2 <- unname(map2[gsm_in_col])
+   idx <- is.na(grp) & !is.na(grp2)
+   grp[idx] <- grp2[idx]
+ }
> 
> if (any(is.na(grp))) {
+   stop("Could not map these columns to a group: ",
+        paste(colnames(E)[is.na(grp)], collapse = ", "))
+ }
> group <- factor(grp)
> 
> # --- PCA figure ----------------------------------------------------------------
> pc <- prcomp(t(E), scale. = TRUE)
> png("results/PCA_samples.png", width = 1600, height = 1200, res = 200)
> plot(pc$x[,1], pc$x[,2], pch = 19, xlab = "PC1", ylab = "PC2",
+      col = as.integer(group), main = "PCA (log2 normalized)")
> legend("topright", legend = levels(group), col = seq_along(levels(group)), pch = 19)
> dev.off()
null device 
          1 
> 
> # --- ZIP artifact --------------------------------------------------------------
> if (!requireNamespace("zip", quietly = TRUE)) install.packages("zip")
> files_to_zip <- c(
+   list.files(pattern = "\\.(R|Rmd|qmd)$"),
+   "metadata/targets.tsv",
+   list.files("results", recursive = TRUE, full.names = TRUE),
+   "renv.lock", "README.md", "LICENSE", "CITATION.cff"
+ )
> files_to_zip <- files_to_zip[file.exists(files_to_zip)]
> zip::zip(zipfile = "artifact/code_and_results_v1.0.zip", files = files_to_zip)
> 
> # --- Final verification --------------------------------------------------------
> expected <- c(
+   "metadata/targets.tsv",
+   "results/deg_RvsNR.tsv",
+   "results/volcano_RvsNR.png",
+   "results/PCA_samples.png",
+   "results/heatmap_top100_RvsNR.png",
+   "results/normalized_log2_expression.tsv",
+   "results/sessionInfo.txt",
+   "artifact/code_and_results_v1.0.zip",
+   "renv.lock"
+ )
> print(all(file.exists(expected)))
[1] TRUE
> print(data.frame(file = expected, exists = file.exists(expected)))
                                    file exists
1                   metadata/targets.tsv   TRUE
2                  results/deg_RvsNR.tsv   TRUE
3              results/volcano_RvsNR.png   TRUE
4                results/PCA_samples.png   TRUE
5       results/heatmap_top100_RvsNR.png   TRUE
6 results/normalized_log2_expression.tsv   TRUE
7                results/sessionInfo.txt   TRUE
8     artifact/code_and_results_v1.0.zip   TRUE
9                              renv.lock   TRUE
> 
