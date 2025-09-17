# 02_analysis.R — Agilent one-color limma workflow (R vs NR) for 8 samples

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma","pheatmap","EnhancedVolcano"), ask = FALSE, update = FALSE)
library(limma); library(pheatmap); library(EnhancedVolcano)

dir.create("metadata", showWarnings = FALSE)
dir.create("results",  showWarnings = FALSE)

# --- map GSM IDs -> actual raw files ------------------------------------------
gsm_ids <- c("GSM5470636","GSM5470637","GSM5470642","GSM5470643","GSM5470644","GSM5470648","GSM5470654","GSM5470667")
group_vec <- c("R","NR","NR","NR","NR","NR","R","NR")

raw_dir <- "data/raw/GSE180775"
files_all <- list.files(raw_dir, recursive = TRUE, pattern = "\\.(txt|TXT)(\\.gz)?$", full.names = TRUE)
pick_one_file <- function(gsm){
  cand <- files_all[grepl(paste0("(^|/)", gsm, "(_|\\.)"), basename(files_all), ignore.case = TRUE)]
  if (!length(cand)) cand <- files_all[grepl(gsm, files_all, ignore.case = TRUE)]
  if (!length(cand)) stop("No file for ", gsm)
  cand[which.max(file.info(cand)$size)]
}
picked <- vapply(gsm_ids, pick_one_file, character(1))
targets <- data.frame(GSM=gsm_ids, Group=group_vec, File=picked, check.names=FALSE)
write.table(targets, "metadata/targets.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# --- read Agilent Feature Extraction files ------------------------------------
x <- read.maimages(files = targets$File, source = "agilent", green.only = TRUE)

# optional: remove control probes
if ("ControlType" %in% colnames(x$genes)) {
  keep_ctl <- is.na(x$genes$ControlType) | x$genes$ControlType == 0
  x <- x[keep_ctl, ]
}

# --- background correction & normalization ------------------------------------
x <- backgroundCorrect(x, method="normexp", offset=50)
x <- normalizeBetweenArrays(x, method="quantile")

# --- filter & collapse duplicates ---------------------------------------------
E <- x$E
keep <- rowSums(E > log2(32)) >= 2
x <- x[keep, ]
if (!"ProbeName" %in% colnames(x$genes)) x$genes$ProbeName <- if (!is.null(rownames(x$E))) rownames(x$E) else seq_len(nrow(x$E))
x <- avereps(x, ID = x$genes$ProbeName)

# --- design, contrast, fit ----------------------------------------------------
group  <- factor(targets$Group, levels = sort(unique(targets$Group)))
design <- model.matrix(~0 + group); colnames(design) <- levels(group)
contrast <- makeContrasts(RvsNR = R - NR, levels = design)

fit  <- lmFit(x, design)
fit2 <- contrasts.fit(fit, contrast)
fit3 <- eBayes(fit2, trend = TRUE, robust = TRUE)

# --- results & figures --------------------------------------------------------
res <- topTable(fit3, coef="RvsNR", number=Inf, adjust.method="BH")
sym_cols <- c("GeneName","GENE_SYMBOL","GeneSymbol","Symbol","Gene.symbol","Gene_Symbol")
sym_col  <- sym_cols[sym_cols %in% colnames(x$genes)][1]
if (!is.na(sym_col)) { sym <- x$genes[,sym_col]; names(sym) <- rownames(x$E); res$Symbol <- sym[rownames(res)] }
write.table(res, "results/deg_RvsNR.tsv", sep="\t", quote=FALSE)

png("results/MAplot_RvsNR.png", 1600, 1200, res=200); plotMA(fit2, coef=1, main="MA: R vs NR"); abline(h=0,lty=2); dev.off()
png("results/SA_trend.png", 1600, 1200, res=200); plotSA(fit3, main="Mean–variance trend"); dev.off()

lab_vec <- if ("Symbol" %in% colnames(res)) res$Symbol else rownames(res)
png("results/volcano_RvsNR.png", 1800, 1400, res=200)
print(EnhancedVolcano(res, lab=lab_vec, x="logFC", y="adj.P.Val", pCutoff=0.05, FCcutoff=1, title="R vs NR"))
dev.off()

top_idx <- head(order(res$adj.P.Val), 100)
ann <- data.frame(Group=group); rownames(ann) <- colnames(x$E)
png("results/heatmap_top100_RvsNR.png", 1600, 2000, res=200)
pheatmap(x$E[top_idx,], cluster_rows=TRUE, cluster_cols=TRUE, scale="row", fontsize=8,
         annotation_col=ann, main="Top 100 DE: R vs NR")
dev.off()

pc <- prcomp(t(x$E), scale.=TRUE)
png("results/PCA_samples.png", 1600, 1200, res=200)
plot(pc$x[,1], pc$x[,2], pch=19, xlab="PC1", ylab="PC2", col=as.integer(group), main="PCA (log2 normalized)")
legend("topright", legend=levels(group), col=seq_along(levels(group)), pch=19)
dev.off()

write.table(x$E, "results/normalized_log2_expression.tsv", sep="\t", quote=FALSE)
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
