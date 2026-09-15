## =============================================================================
## Quality control: is the signal confounded by the fungal / bacterial ratio?
##
## The libraries are mixed: each sample contains transcripts from both the
## fungus and the bacterium, and the proportion of reads that map to the fungal
## genome differs between the two groups. TMM normalises the composition of the
## fungal library; it cannot correct a difference in the fungal-to-bacterial
## biomass ratio. This script quantifies the problem and shows why it cannot be
## fixed by adding a covariate.
##
## Usage
##   Rscript scripts/03_qc/qc_bacterial_load.R
##
## Outputs (results/figures/)
##   qc_fungal_fraction_<short>.png   the fraction per sample, by group
##   qc_confound_ftm_<short>.png      how tightly the ftm genes track that fraction
## Console / results/geneset/qc_bacterial_load_<short>.txt
##   the covariate-adjusted model and why it is over-adjusted
## =============================================================================

suppressPackageStartupMessages({
  library(limma); library(edgeR); library(ggplot2); library(ggrepel)})
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
d    <- DESIGNS[[did]]
s    <- SAMPLES[match(d$samples, SAMPLES$sample), ]
s    <- s[order(factor(s$group, levels = c("control", "treatment"))), ]
frac <- setNames(s$fungal_pct, s$sample)
grp  <- factor(s$group, levels = c("control", "treatment"))
log  <- c()
say  <- function(...) { m <- sprintf(...); message(m); log <<- c(log, m) }

say("=== Fraction of reads mapping to the fungal genome -- %s ===", d$label)
print(data.frame(sample = s$label, group = ifelse(grp == "control", "KO", "WT"),
                 fungal_pct = frac), row.names = FALSE)
say("Mean: WT = %.1f%%   KO = %.1f%%", mean(frac[grp == "treatment"]),
    mean(frac[grp == "control"]))
tt <- t.test(frac ~ grp)
say("t-test between groups: t = %.2f, P = %.3f", tt$statistic, tt$p.value)
say("Correlation between fungal fraction and group (KO=0/WT=1): r = %.2f",
    cor(frac, as.numeric(grp) - 1))

## ---- Figure A: the fraction itself -----------------------------------------
dd <- data.frame(sample = s$label, pct = frac,
                 group = factor(ifelse(grp == "control", "dhbF KO", "DK1042 WT"),
                                levels = c("DK1042 WT", "dhbF KO")))
pA <- ggplot(dd, aes(group, pct, fill = group)) +
  geom_boxplot(width = .45, alpha = .35, outlier.shape = NA, colour = "#3d3d3d") +
  geom_point(size = 3.4, shape = 21, colour = "#1a1a1a", stroke = .5) +
  geom_text_repel(aes(label = sample), size = 3.2, colour = "#3d3d3d", nudge_x = .28,
                  direction = "y", min.segment.length = 0,
                  segment.colour = "#b0b0b0", seed = 3) +
  scale_fill_manual(values = c("DK1042 WT" = "#c2503f", "dhbF KO" = "#3b6ea5"),
                    guide = "none") +
  scale_y_continuous(limits = c(0, 90)) +
  labs(title = "A. Reads mapping to the fungal genome",
       subtitle = sprintf("WT = %.1f%% vs KO = %.1f%%  (t-test P = %.3f)",
                          mean(frac[grp == "treatment"]), mean(frac[grp == "control"]),
                          tt$p.value),
       x = NULL, y = "% of reads pseudo-aligned to P. brasilianum") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = .7),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "#57565a"))
ggsave(file.path(FIG_DIR, sprintf("qc_fungal_fraction_%s.png", d$short)), pA,
       width = 150, height = 110, units = "mm", dpi = 600, bg = "white")

## ---- the fitted model, reused from the gene-set step if available ----------
rds <- file.path(GS_DIR, sprintf("fit_objects_%s.rds", d$short))
if (!file.exists(rds))
  stop("Missing ", rds, " -- run scripts/02_geneset_analysis/geneset_tests.R ", did, " first.")
o <- readRDS(rds); v <- o$v; idx <- o$idx; tab <- o$tab
E <- v$E[, s$sample, drop = FALSE]

## ---- Figure B: do the ftm genes track the fraction more than other genes? --
r   <- apply(E, 1, function(x) suppressWarnings(cor(x, frac[colnames(E)])))
lfc <- tab$logFC[match(rownames(E), tab$gene_ID)]
k   <- idx[["Verruculogen"]]
lo  <- min(abs(lfc[k]), na.rm = TRUE) * .95
hi  <- max(abs(lfc[k]), na.rm = TRUE) * 1.05
matched <- which(abs(lfc) >= lo & abs(lfc) <= hi)
matched <- setdiff(matched, k)

say("\nftm cluster: mean |logFC| = %.2f, median |r| with the fungal fraction = %.2f",
    mean(abs(lfc[k]), na.rm = TRUE), median(abs(r[k])))
say("Other genes of comparable |logFC| (%.1f-%.1f, n = %d): median |r| = %.2f",
    lo, hi, length(matched), median(abs(r[matched])))
w <- wilcox.test(abs(r[k]), abs(r[matched]))
say("Wilcoxon P = %.2e  (EXPLORATORY -- see the caveats below)", w$p.value)
say("  -> the cluster tracks the fungal read fraction %s than the background",
    ifelse(median(abs(r[k])) > median(abs(r[matched])), "MORE closely", "less closely"))
say("")
say("Read this comparison as descriptive, not as a hypothesis test. Three")
say("things make the P value far smaller than the evidence warrants:")
say("  1. The fungal read fraction is not a measurement of biomass. It is the")
say("     share of a mixed library that maps to the fungal genome, which also")
say("     moves with RNA yield, rRNA depletion and bacterial transcriptional")
say("     activity.")
say("  2. Wilcoxon treats every gene as an independent observation. The eight")
say("     ftm genes sit in one co-regulated cluster and are anything but")
say("     independent, so the effective sample size is far below n = 8.")
say("  3. The background set is chosen to match the OBSERVED |logFC| of the")
say("     ftm genes, and |logFC| is itself a function of the group difference")
say("     that the fungal fraction tracks. The matching is therefore partly")
say("     circular.")
say("")
say("With %d libraries, the honest conclusion is that a contribution of", ncol(E))
say("differing fungal biomass CANNOT BE EXCLUDED from RNA-seq alone. Settling")
say("it needs qRT-PCR against a fungal housekeeping gene, or qPCR of fungal")
say("ITS against bacterial 16S, or dry weight.")

dd2 <- rbind(data.frame(g = "Other genes, comparable |logFC|", r = abs(r[matched])),
             data.frame(g = sprintf("ftm cluster (n = %d)", length(k)), r = abs(r[k])))
dd2$g <- factor(dd2$g, levels = unique(dd2$g))
pB <- ggplot(dd2, aes(g, r, fill = g)) +
  geom_violin(alpha = .3, colour = NA) +
  geom_boxplot(width = .22, alpha = .6, outlier.shape = NA, colour = "#3d3d3d") +
  scale_fill_manual(values = c("#9aa0a6", PAL[["Verruculogen"]]), guide = "none") +
  labs(title = "B. Association with the fungal read fraction",
       subtitle = sprintf("Exploratory; Wilcoxon P = %.1e, genes are not independent",
                          w$p.value),
       x = NULL, y = "|r| between expression and fungal read fraction") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = .7),
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, colour = "#57565a"))
ggsave(file.path(FIG_DIR, sprintf("qc_confound_ftm_%s.png", d$short)), pB,
       width = 150, height = 110, units = "mm", dpi = 600, bg = "white")

## ---- what happens if the fraction is added as a covariate ------------------
gene_len <- read.delim(file.path(ANNOT, "gene_length_Penbrz1.txt"), sep = "\t",
                       stringsAsFactors = FALSE)
Fg <- readDGE(file.path(COUNTS, s$file), columns = c(1, 3)); colnames(Fg) <- s$sample
Fg$samples$group <- grp
Fg$genes <- gene_len[match(rownames(Fg), gene_len$gene_ID), ]
Fg <- Fg[rowSums(cpm(Fg) > 1) >= 3, , keep.lib.sizes = FALSE]
Fg <- calcNormFactors(Fg, method = "TMM")
d2 <- model.matrix(~ 0 + grp + scale(frac))
colnames(d2) <- c("control", "treatment", "fungal_frac")
v2  <- voom(Fg, d2)
cm2 <- makeContrasts(treatment - control, levels = colnames(d2))
r2  <- topTable(eBayes(contrasts.fit(lmFit(v2, d2), cm2)), coef = 1, n = Inf)
f2  <- r2[r2$gene_ID %in% names(FTM), ]
say("\n=== ftm cluster AFTER adjusting for the fungal fraction ===")
say("mean logFC = %.2f   (unadjusted: %.2f)", mean(f2$logFC), mean(lfc[k], na.rm = TRUE))
set.seed(1)
fr2 <- fry(v2, index = list(Verruculogen = which(rownames(Fg) %in% names(FTM))),
           design = d2, contrast = cm2[, 1])
say("Cluster-level test after adjustment: P = %.3f", fr2$PValue[1])
say("")
say("INTERPRETATION. The adjusted model must not be read as evidence that the")
say("effect disappears. The covariate correlates with the group at r = %.2f, so the",
    cor(frac, as.numeric(grp) - 1))
say("two are nearly collinear and the group effect is barely estimable. More")
say("importantly, if bacillibactin restricts fungal growth then the lower fungal")
say("fraction in the wild-type co-culture is part of the biological effect, i.e. a")
say("mediator rather than a confounder, and adjusting for a mediator removes the")
say("very quantity the experiment set out to measure. This limitation is resolved")
say("experimentally -- by qRT-PCR normalised to a fungal housekeeping gene, by")
say("direct biomass measurement, and by LC-MS quantification of verruculogen --")
say("not statistically.")

writeLines(log, file.path(GS_DIR, sprintf("qc_bacterial_load_%s.txt", d$short)))
message("\nWritten: ", file.path(GS_DIR, sprintf("qc_bacterial_load_%s.txt", d$short)))
