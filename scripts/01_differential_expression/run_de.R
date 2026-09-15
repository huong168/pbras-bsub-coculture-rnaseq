## =============================================================================
## Differential expression, limma-voom + TMM
##
## Reproduces the pipeline used throughout the project, which follows the
## workflow the laboratory applies to its Penicillium co-culture RNA-seq data:
##
##   readDGE -> filter rowSums(cpm > 1) >= 3 -> calcNormFactors(method = "TMM")
##        -> voom -> lmFit -> contrasts.fit -> treat(lfc = 1) -> topTreat
##
## The contrast is treatment - control, i.e. (co-culture with DK1042 wild type)
## minus (co-culture with the dhbF knock-out). A negative logFC therefore means
## the gene is expressed at a LOWER level when the bacillibactin-producing wild
## type is present.
##
## Usage
##   Rscript scripts/01_differential_expression/run_de.R
##
## Outputs, written to results/de_tables/
##   DE_<short>_treat.tsv        raw topTreat output
##   DE_<short>_annotated.csv    same table plus fold change, BGC class,
##                               ftm enzyme name and product description
##   DE_<short>_eBayes.tsv       standard eBayes test, no fold-change threshold
##   cpm_<short>.tsv, rpkm_<short>.tsv
## =============================================================================

suppressPackageStartupMessages({library(limma); library(edgeR)})
## Find config.R relative to this script, so the script runs from any directory.
local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  d <- if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else getwd()
  repeat {
    if (file.exists(file.path(d, ".repo-root"))) break
    up <- dirname(d); if (identical(up, d)) { d <- getwd(); break }; d <- up
  }
  cfg <- file.path(Sys.getenv("PBRAS_ROOT", unset = d), "scripts", "config.R")
  if (!file.exists(cfg)) stop("Cannot find scripts/config.R; set PBRAS_ROOT.")
  source(cfg, local = FALSE)
})

## ---- annotation shared by every design -------------------------------------
gene_len <- read.delim(file.path(ANNOT, "gene_length_Penbrz1.txt"),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
bgc <- read.csv(file.path(ANNOT, "bgc_genes.csv"), stringsAsFactors = FALSE)
prod <- read.delim(file.path(ANNOT, "gene_annot.tsv"), sep = "\t",
                   stringsAsFactors = FALSE)

run_one_design <- function(d) {
  message("\n=== ", d$label, " (", paste(d$samples, collapse = ", "), ") ===")
  s <- SAMPLES[match(d$samples, SAMPLES$sample), ]
  ## readDGE expects control samples first so that treatment - control is
  ## positive for genes higher in the wild-type co-culture.
  s <- s[order(factor(s$group, levels = c("control", "treatment"))), ]

  Fg <- readDGE(file.path(COUNTS, s$file), columns = c(1, 3))
  colnames(Fg) <- s$sample
  group <- factor(s$group, levels = c("control", "treatment"))
  Fg$samples$group <- group
  Fg$genes <- gene_len[match(rownames(Fg), gene_len$gene_ID), ]

  ## Expression filter. Requiring three of the four libraries above 1 CPM is
  ## stricter than edgeR's filterByExpr(), which would size its threshold to
  ## the two-sample groups; scripts/03_qc/filter_sensitivity.R runs both and
  ## shows that the verruculogen result is the same either way.
  keep <- rowSums(cpm(Fg) > 1) >= 3
  Fg <- Fg[keep, , keep.lib.sizes = FALSE]
  Fg <- calcNormFactors(Fg, method = "TMM")
  message("Genes passing the expression filter: ", nrow(Fg))

  ## Normalised expression, computed AFTER filtering and TMM so that the MDS
  ## plot and the exported CPM table describe the same data the model sees.
  cpm_all  <- cpm(Fg)
  lcpm_all <- cpm(Fg, log = TRUE)

  design <- model.matrix(~ 0 + group)
  colnames(design) <- gsub("group", "", colnames(design))
  cm <- makeContrasts(controlvstreatment = treatment - control,
                      levels = colnames(design))

  v    <- voom(Fg, design)
  vfit <- contrasts.fit(lmFit(v, design), contrasts = cm)
  efit <- eBayes(vfit)                       # standard moderated t test
  tfit <- treat(vfit, lfc = LFC_TREAT)       # compound hypothesis |logFC| > 1

  res_t <- topTreat(tfit, coef = 1, n = Inf)
  res_e <- topTable(efit, coef = 1, n = Inf)

  message("decideTests(treat):");  print(summary(decideTests(tfit)))
  message("decideTests(eBayes):"); print(summary(decideTests(efit)))

  ## ---- annotated table ------------------------------------------------------
  ann <- data.frame(
    gene_ID     = res_t$gene_ID,
    gene_length = res_t$gene_length,
    logFC       = res_t$logFC,
    fold_change = 2^res_t$logFC,            # ratio treatment / control
    direction   = ifelse(res_t$logFC >= 0, "up in DK1042 WT", "down in DK1042 WT"),
    AveExpr     = res_t$AveExpr,
    t           = res_t$t,
    P.Value     = res_t$P.Value,
    adj.P.Val   = res_t$adj.P.Val,
    BGC_class   = bgc$BGC_class[match(res_t$gene_ID, bgc$locus_tag)],
    enzyme      = unname(FTM[res_t$gene_ID]),
    product     = prod$product[match(res_t$gene_ID, prod$locus_tag)],
    stringsAsFactors = FALSE)
  ann$BGC_class[is.na(ann$BGC_class)] <- ""
  ann$enzyme[is.na(ann$enzyme)] <- ""
  ann$product[is.na(ann$product)] <- ""

  write.table(res_t, file.path(DE_DIR, sprintf("DE_%s_treat.tsv", d$short)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(res_e, file.path(DE_DIR, sprintf("DE_%s_eBayes.tsv", d$short)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(ann, file.path(DE_DIR, sprintf("DE_%s_annotated.csv", d$short)),
            row.names = FALSE)
  write.table(cpm_all, file.path(DE_DIR, sprintf("cpm_%s.tsv", d$short)),
              sep = "\t", quote = FALSE)
  write.table(rpkm(Fg, Fg$genes$gene_length),
              file.path(DE_DIR, sprintf("rpkm_%s.tsv", d$short)),
              sep = "\t", quote = FALSE)

  ## ---- QC figures, as in the laboratory workflow ----------------------------
  png(file.path(FIG_DIR, sprintf("qc_voom_%s.png", d$short)),
      width = 1600, height = 1400, res = 200)
  voom(Fg, design, plot = TRUE); dev.off()

  png(file.path(FIG_DIR, sprintf("qc_MDS_%s.png", d$short)),
      width = 1800, height = 1500, res = 220)
  plotMDS(lcpm_all[, s$sample], col = ifelse(s$group == "control", "grey45", "black"),
          pch = 19, cex = 1.5,
          xlab = "Leading log2 fold change, dimension 1",
          ylab = "Leading log2 fold change, dimension 2")
  title(main = d$label); dev.off()

  png(file.path(FIG_DIR, sprintf("qc_MD_%s.png", d$short)),
      width = 1800, height = 1500, res = 220)
  plotMD(tfit, column = 1, status = decideTests(tfit)[, 1], main = d$label,
         xlim = c(-3, 13), ylim = c(-12, 12), legend = FALSE, las = 1, cex = .7)
  dev.off()

  ## ---- the pre-specified cluster --------------------------------------------
  ftm <- ann[ann$enzyme != "", c("enzyme", "gene_ID", "logFC", "P.Value", "adj.P.Val")]
  ftm <- ftm[order(match(ftm$enzyme, FTM_ORDER)), ]
  ## Benjamini-Hochberg within the family of the eight pre-specified genes
  ftm$adj.P.within_cluster <- p.adjust(ftm$P.Value, method = "BH")
  write.csv(ftm, file.path(DE_DIR, sprintf("ftm_cluster_%s.csv", d$short)),
            row.names = FALSE)
  message("ftm cluster:")
  print(format(ftm, digits = 3), row.names = FALSE)

  invisible(ftm)
}

args <- commandArgs(trailingOnly = TRUE)
wanted <- if (length(args)) args else names(DESIGNS)
unknown <- setdiff(wanted, names(DESIGNS))
if (length(unknown)) stop("Unknown design(s): ", paste(unknown, collapse = ", "),
                          ". Available: ", paste(names(DESIGNS), collapse = ", "))
invisible(lapply(DESIGNS[wanted], run_one_design))
message("\nDone. Tables written to ", DE_DIR)
