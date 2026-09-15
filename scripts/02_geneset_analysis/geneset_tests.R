## =============================================================================
## Gene-set level analysis of the pre-specified verruculogen cluster
##
## The ftm cluster was the hypothesis of the experiment, not a discovery from a
## genome-wide screen. Correcting its eight P values against ~9,000 genes
## penalises tests that were never asked. This script applies the three
## analyses that address that properly:
##
##   1. Self-contained gene-set tests  (fry, and roast with 99,999 rotations)
##      -- one test for the whole cluster instead of eight separate tests.
##   2. Benjamini-Hochberg correction WITHIN each pre-specified gene family.
##   3. Effect sizes with 95% confidence intervals, rather than P values alone.
##
## Usage
##   Rscript scripts/02_geneset_analysis/geneset_tests.R
##
## Outputs (results/geneset/)
##   geneset_fry_<short>.csv        cluster-level test for every gene set
##                                  (exploratory -- see the note in section 1)
##   geneset_roast_<short>.txt      rotation test; Verruculogen is the primary,
##                                  pre-specified test, the rest is exploratory
##   DE_with_setFDR_and_CI_<short>.csv
##   fit_objects_<short>.rds        limma fit, reusable by other scripts
## Figures (results/figures/)
##   barcode_ftm_<short>.png, barcode_iron_uptake_<short>.png
##   forest_ftm_CI_<short>.png/.pdf, geneset_level_<short>.png/.pdf
## =============================================================================

suppressPackageStartupMessages({library(limma); library(edgeR); library(ggplot2)})
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

args <- commandArgs(trailingOnly = TRUE)
did  <- if (length(args)) args[1] else "d22a"
if (!did %in% names(DESIGNS)) stop("Unknown design: ", did)
d <- DESIGNS[[did]]
message("Gene-set analysis for: ", d$label)

## ---- 1. preprocessing, identical to run_de.R -------------------------------
gene_len <- read.delim(file.path(ANNOT, "gene_length_Penbrz1.txt"), sep = "\t",
                       stringsAsFactors = FALSE)
s <- SAMPLES[match(d$samples, SAMPLES$sample), ]
s <- s[order(factor(s$group, levels = c("control", "treatment"))), ]

Fg <- readDGE(file.path(COUNTS, s$file), columns = c(1, 3))
colnames(Fg) <- s$sample
group <- factor(s$group, levels = c("control", "treatment"))
Fg$samples$group <- group
Fg$genes <- gene_len[match(rownames(Fg), gene_len$gene_ID), ]
Fg <- Fg[rowSums(cpm(Fg) > 1) >= 3, , keep.lib.sizes = FALSE]
Fg <- calcNormFactors(Fg, method = "TMM")
message("Genes tested: ", nrow(Fg))

design <- model.matrix(~ 0 + group)
colnames(design) <- gsub("group", "", colnames(design))
cm   <- makeContrasts(controlvstreatment = treatment - control, levels = colnames(design))
v    <- voom(Fg, design)
fit  <- contrasts.fit(lmFit(v, design), cm)
efit <- eBayes(fit)
tfit <- treat(fit, lfc = LFC_TREAT)
res  <- topTreat(tfit, coef = 1, n = Inf)

## ---- 2. gene sets ----------------------------------------------------------
gs   <- read.csv(file.path(ANNOT, "gene_sets.csv"), stringsAsFactors = FALSE)
sets <- split(gs$locus_tag, gs$set)
idx  <- lapply(sets, function(g) which(rownames(Fg) %in% g))
idx  <- idx[sapply(idx, length) >= 3]    # a set needs at least 3 genes to be testable

set.seed(2026)
fr <- fry(v, index = idx, design = design, contrast = cm[, 1])
fr$set      <- rownames(fr)
fr$n_genes  <- sapply(idx, length)[fr$set]
fr$FDR_sets <- p.adjust(fr$PValue, method = "BH")   # corrected over sets, not genes
fr <- fr[order(fr$PValue), ]
message("\n=== Cluster-level test (fry, self-contained) ===")
print(fr[, c("n_genes", "Direction", "PValue", "FDR_sets")], digits = 3)
write.csv(fr[, c("set", "n_genes", "Direction", "PValue", "FDR_sets")],
          file.path(GS_DIR, sprintf("geneset_fry_%s.csv", d$short)), row.names = FALSE)

## The PRIMARY, pre-specified test of this experiment is roast() on the eight
## genes of the Verruculogen cluster: one hypothesis, fixed before the data were
## seen, so its P value needs no correction. Every other set below -- including
## the whole fry() table above -- is exploratory and is reported for context.
## Treating any of them as confirmatory would be testing a hypothesis the
## experiment was not designed to ask.
PRIMARY <- "Verruculogen"
focal <- intersect(c(PRIMARY, "Iron_uptake", "Fe_S_iron_dependent"), names(idx))
roast_lines <- c(sprintf("roast(), 99,999 rotations -- %s", d$label), "",
                 sprintf("PRIMARY (pre-specified): %s", PRIMARY),
                 "Everything else on this page, and the whole fry table, is exploratory.",
                 "")
set.seed(2026)
for (nm in focal) {
  r <- roast(v, index = idx[[nm]], design = design, contrast = cm[, 1], nrot = 99999)
  roast_lines <- c(roast_lines, sprintf(
    "%-22s n=%2d | Down P=%.5f  Up P=%.5f  Active.Prop(Down)=%.2f  Active.Prop(Up)=%.2f",
    paste0(nm, if (identical(nm, PRIMARY)) " *" else ""),
    length(idx[[nm]]), r$p.value["Down", "P.Value"], r$p.value["Up", "P.Value"],
    r$p.value["Down", "Active.Prop"], r$p.value["Up", "Active.Prop"]))
}
message("\n", paste(roast_lines, collapse = "\n"))
writeLines(roast_lines, file.path(GS_DIR, sprintf("geneset_roast_%s.txt", d$short)))

## ---- 3. within-family BH and 95% confidence intervals ----------------------
tab <- res
tab$set    <- gs$set[match(tab$gene_ID, gs$locus_tag)]
tab$enzyme <- unname(FTM[tab$gene_ID])
tab$FDR_genomewide <- tab$adj.P.Val
tab$FDR_within_set <- NA_real_
for (nm in names(idx)) {
  k <- which(tab$set == nm)
  if (length(k)) tab$FDR_within_set[k] <- p.adjust(tab$P.Value[k], method = "BH")
}
se <- sqrt(efit$s2.post) * efit$stdev.unscaled[, 1]
ci <- data.frame(gene_ID = rownames(efit),
                 CI_low  = efit$coefficients[, 1] - qt(.975, efit$df.total) * se,
                 CI_high = efit$coefficients[, 1] + qt(.975, efit$df.total) * se)
tab <- merge(tab, ci, by = "gene_ID", all.x = TRUE)
write.csv(tab, file.path(GS_DIR, sprintf("DE_with_setFDR_and_CI_%s.csv", d$short)),
          row.names = FALSE)

ftm <- tab[!is.na(tab$enzyme), ]
ftm <- ftm[order(match(ftm$enzyme, FTM_ORDER)), ]
message("\n=== ftm cluster: logFC (95% CI) and two ways of correcting ===")
print(format(ftm[, c("enzyme", "logFC", "CI_low", "CI_high", "P.Value",
                     "FDR_genomewide", "FDR_within_set")], digits = 3), row.names = FALSE)
message(sprintf("\nAt FDR < 0.05:  genome-wide = %d/8   within the cluster = %d/8",
                sum(ftm$FDR_genomewide < .05), sum(ftm$FDR_within_set < .05)))

saveRDS(list(fit = fit, efit = efit, tfit = tfit, v = v, idx = idx, tab = tab,
             design = design, cm = cm, samples = s),
        file.path(GS_DIR, sprintf("fit_objects_%s.rds", d$short)))

## ---- 4. barcode plots, the standard visualisation for roast ----------------
for (nm in focal) {
  png(file.path(FIG_DIR, sprintf("barcode_%s_%s.png", tolower(nm), d$short)),
      width = 2000, height = 1100, res = 240)
  par(mar = c(4.5, 4, 3, 1))
  barcodeplot(efit$t[, 1], index = idx[[nm]],
              main = sprintf("%s (n = %d)", nm, length(idx[[nm]])),
              xlab = "moderated t  (DK1042 WT vs dhbF KO)",
              labels = c("lower in DK1042 WT", "higher in DK1042 WT"), quantile = c(-1, 1))
  dev.off()
}

## ---- 5. forest plot of the cluster -----------------------------------------
f <- ftm
f$enzyme <- factor(f$enzyme, levels = rev(FTM_ORDER))
p1 <- ggplot(f, aes(logFC, enzyme)) +
  annotate("rect", xmin = -1, xmax = 1, ymin = -Inf, ymax = Inf, fill = "#f2f2f2") +
  geom_vline(xintercept = 0, colour = "black", linewidth = .6) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "#999999", linewidth = .5) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = .22,
                 colour = "#b3140f", linewidth = .7) +
  geom_point(size = 3.4, colour = PAL[["Verruculogen"]]) +
  geom_text(aes(x = CI_high + .35, label = sprintf("FDR = %.3f", FDR_within_set)),
            hjust = 0, size = 3.3, colour = "#3d3d3d") +
  scale_x_continuous(limits = c(min(f$CI_low) - 1, max(f$CI_high) + 2.6)) +
  labs(title = "ftm / verruculogen cluster: log2 fold change with 95% CI",
       subtitle = paste0(d$label, " | FDR corrected within the family of 8 genes"),
       x = expression(log[2]~"fold change  (DK1042 WT / "*Delta*"dhbF)"), y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = .7),
        axis.text.y = element_text(face = "italic", size = 11, colour = "#b3140f"),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "#57565a"))
ggsave(file.path(FIG_DIR, sprintf("forest_ftm_CI_%s.png", d$short)), p1,
       width = 190, height = 110, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(FIG_DIR, sprintf("forest_ftm_CI_%s.pdf", d$short)), p1,
       width = 190, height = 110, units = "mm", bg = "white")

## ---- 6. cluster-level result figure ----------------------------------------
fr$lab <- sprintf("%s  (n = %d)", fr$set, fr$n_genes)
fr$y   <- -log10(fr$PValue)
fr2 <- fr[order(fr$y), ]; fr2$lab <- factor(fr2$lab, levels = fr2$lab)
p2 <- ggplot(fr2, aes(y, lab, fill = Direction)) +
  geom_vline(xintercept = -log10(.05), linetype = "dashed",
             colour = "#0c0c0c", linewidth = .55) +
  geom_col(width = .62) +
  geom_text(aes(label = sprintf("P = %.3f", PValue)), hjust = -.12, size = 3.2,
            colour = "#3d3d3d") +
  scale_fill_manual(values = c(Down = "#3b6ea5", Up = "#c2503f"),
                    labels = c(Down = "lower in DK1042 WT", Up = "higher in DK1042 WT"),
                    name = NULL) +
  scale_x_continuous(limits = c(0, max(fr2$y) * 1.25), expand = c(0, 0)) +
  labs(title = "Cluster-level tests (fry, self-contained)",
       subtitle = paste0(d$label, " | one test per gene set, not diluted by the rest of the genome"),
       x = expression(-log[10]~"("*italic(P)*")"), y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = .7),
        legend.position = "top", plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "#57565a"))
ggsave(file.path(FIG_DIR, sprintf("geneset_level_%s.png", d$short)), p2,
       width = 180, height = 120, units = "mm", dpi = 600, bg = "white")
ggsave(file.path(FIG_DIR, sprintf("geneset_level_%s.pdf", d$short)), p2,
       width = 180, height = 120, units = "mm", bg = "white")

message("\nDone. Tables in ", GS_DIR, ", figures in ", FIG_DIR)
