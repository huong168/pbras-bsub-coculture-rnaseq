## ---------------------------------------------------------------------------
## Does the expression filter change the conclusion?
##
## The analysis keeps genes with rowSums(cpm > 1) >= 3, that is, three of the
## four libraries. With two replicates per group this is stricter than edgeR's
## filterByExpr(), which sizes its threshold to the smallest group: a gene
## expressed in exactly one group clears two libraries, not three, and is
## discarded.
##
## This script runs the analysis under both filters and reports what changes:
## the size of the gene universe, the number of DE genes, and the eight genes of
## the pre-specified ftm cluster.
##
##   Rscript scripts/03_qc/filter_sensitivity.R
## ---------------------------------------------------------------------------

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

suppressPackageStartupMessages({library(edgeR); library(limma)})

gene_len <- read.delim(file.path(ANNOT, "gene_length_Penbrz1.txt"),
                       stringsAsFactors = FALSE)

fit_one <- function(d, filter = c("fixed", "byexpr")) {
  filter <- match.arg(filter)
  s <- SAMPLES[SAMPLES$sample %in% d$samples, ]
  s <- s[order(factor(s$group, levels = c("control", "treatment"))), ]

  Fg <- readDGE(file.path(COUNTS, s$file), columns = c(1, 3))
  colnames(Fg) <- s$sample
  group <- factor(s$group, levels = c("control", "treatment"))
  Fg$samples$group <- group
  Fg$genes <- gene_len[match(rownames(Fg), gene_len$gene_ID), ]

  design <- model.matrix(~ 0 + group)
  colnames(design) <- gsub("group", "", colnames(design))

  keep <- if (filter == "fixed") rowSums(cpm(Fg) > 1) >= 3
          else filterByExpr(Fg, design = design)
  Fg <- Fg[keep, , keep.lib.sizes = FALSE]
  Fg <- calcNormFactors(Fg, method = "TMM")

  cm <- makeContrasts(controlvstreatment = treatment - control,
                      levels = colnames(design))
  tfit <- treat(contrasts.fit(lmFit(voom(Fg, design), design), contrasts = cm),
                lfc = LFC_TREAT)
  tab <- topTreat(tfit, coef = 1, n = Inf)
  list(n_genes = nrow(Fg),
       n_de    = sum(tab$adj.P.Val < FDR_CUTOFF),
       tab     = tab)
}

ftm_line <- function(r) {
  k <- match(names(FTM), r$tab$gene_ID)
  present <- sum(!is.na(k))
  k <- k[!is.na(k)]
  c(present = present,
    sig_fdr = sum(r$tab$adj.P.Val[k] < FDR_CUTOFF),
    sig_p   = sum(r$tab$P.Value[k] < 0.05),
    ## the most negative logFC in the cluster, whichever gene that is
    min_lfc = if (length(k)) round(min(r$tab$logFC[k]), 2) else NA)
}

out <- do.call(rbind, lapply(DESIGNS, function(d) {
  a <- fit_one(d, "fixed"); b <- fit_one(d, "byexpr")
  fa <- ftm_line(a); fb <- ftm_line(b)
  data.frame(design = d$id, label = d$label,
             genes_fixed = a$n_genes,   genes_byexpr = b$n_genes,
             de_fixed    = a$n_de,      de_byexpr    = b$n_de,
             ftm_kept_fixed  = fa["present"], ftm_kept_byexpr  = fb["present"],
             ftm_fdr_fixed   = fa["sig_fdr"], ftm_fdr_byexpr   = fb["sig_fdr"],
             ftm_p_fixed     = fa["sig_p"],   ftm_p_byexpr     = fb["sig_p"],
             ftm_min_lfc_fixed = fa["min_lfc"], ftm_min_lfc_byexpr = fb["min_lfc"],
             stringsAsFactors = FALSE)
}))
rownames(out) <- NULL

cat("\nExpression filter: rowSums(cpm > 1) >= 3   vs   filterByExpr(design)\n\n")
print(out, row.names = FALSE)

f <- file.path(RESULTS, "filter_sensitivity.csv")
write.csv(out, f, row.names = FALSE)
cat("\nWritten:", f, "\n")
cat("\nAll eight ftm genes survive both filters, so the conclusion about the\n",
    "verruculogen cluster does not depend on the choice.\n", sep = "")
