## =============================================================================
## Publication volcano plot for a single design, in the style used by the
## laboratory (200 x 160 mm, 600 dpi, Annotation legend, red verruculogen genes).
##
## Usage
##   Rscript scripts/04_figures/volcano_labstyle.R [design_id] [stat] [y_max]
##     design_id  d22a                                (the only design)
##     stat       fdr | p                             (default fdr)
##     y_max      upper limit of the y axis           (default: auto)
##
##   # the figure used in the manuscript
##   Rscript scripts/04_figures/volcano_labstyle.R d22a fdr 2.5
##
## Options set in the block below control which extra genes are labelled:
##   LABEL_NRPS_BY  "none"        label only FtmA-H            (default)
##                  "significant" also label NRPS genes passing both thresholds
##                  "largeFC"     also label NRPS genes with |log2FC| >= 2
##
## Output  results/figures/volcano_labstyle_<short>_<stat>.png / .pdf
## =============================================================================

suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(ggrepel)})
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
did   <- if (length(args) >= 1) args[1] else "d22a"
STAT  <- if (length(args) >= 2 && tolower(args[2]) %in% c("p","pvalue")) "p" else "fdr"
Y_ARG <- if (length(args) >= 3) as.numeric(args[3]) else NA_real_
LABEL_NRPS_BY <- "none"

d <- DESIGNS[[did]]; if (is.null(d)) stop("Unknown design: ", did)
STAT_COL   <- if (STAT == "p") "P.Value" else "adj.P.Val"
STAT_LABEL <- if (STAT == "p") "P" else "FDR"

df <- read.csv(file.path(DE_DIR, sprintf("DE_%s_annotated.csv", d$short)),
               stringsAsFactors = FALSE) %>%
  filter(!is.na(logFC), !is.na(.data[[STAT_COL]])) %>%
  mutate(stat = .data[[STAT_COL]], y = -log10(stat),
         enzyme = trimws(ifelse(is.na(enzyme), "", enzyme)),
         bgc    = trimws(ifelse(is.na(BGC_class), "", BGC_class)),
         is_up   = stat < FDR_CUTOFF & logFC >=  FC_CUTOFF,
         is_down = stat < FDR_CUTOFF & logFC <= -FC_CUTOFF)

nrps <- switch(LABEL_NRPS_BY,
  none        = df[0, ],
  significant = df %>% filter(bgc == "NRPS", is_up | is_down),
  largeFC     = df %>% filter(bgc == "NRPS", abs(logFC) >= FC_CUTOFF),
  stop("LABEL_NRPS_BY must be none, significant or largeFC"))

df <- df %>% mutate(plot_group = case_when(
  nzchar(enzyme)                      ~ "Verruculogen",
  gene_ID %in% nrps$gene_ID           ~ "NRPS",
  is_up                               ~ "Upregulated",
  is_down                             ~ "Downregulated",
  TRUE ~ NA_character_))

Y_MAX <- if (!is.na(Y_ARG)) Y_ARG else ceiling(max(df$y) * 2) / 2
n_out <- sum(df$y > Y_MAX)
if (n_out > 0) {
  message(sprintf("Note: %d point(s) lie above y = %.1f and are not drawn.", n_out, Y_MAX))
  df <- df %>% filter(y <= Y_MAX)
}
X_LIM <- max(abs(df$logFC)) * 1.18

eg <- df %>% filter(nzchar(enzyme))
bg <- df %>% filter(is.na(plot_group))
hl <- df %>% filter(!is.na(plot_group))
lv <- names(PAL)[names(PAL) %in% hl$plot_group]
hl$plot_group <- factor(hl$plot_group, levels = lv)

p <- ggplot(df, aes(logFC, y)) +
  geom_point(data = bg, colour = "#c3c2b7", alpha = .45, size = 2.1, stroke = 0) +
  geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF), linetype = "dashed",
             colour = "#0c0c0c", linewidth = .8) +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed",
             colour = "#0c0c0c", linewidth = .8) +
  geom_point(data = subset(hl, plot_group != "Verruculogen"),
             aes(colour = plot_group), size = 2.8) +
  geom_point(data = subset(hl, plot_group == "Verruculogen"),
             aes(colour = plot_group), size = 2.8) +
  geom_text_repel(data = eg, aes(label = enzyme), colour = "#f10808", size = 6.0,
                  family = FIG_FONT, box.padding = .75, point.padding = .55,
                  min.segment.length = 0, segment.colour = "#555555",
                  segment.size = .65, max.overlaps = Inf, seed = 123,
                  xlim = c(-X_LIM * .97, X_LIM * .97),
                  ylim = c(Y_MAX * .37, Y_MAX * .99), show.legend = FALSE) +
  scale_colour_manual(values = PAL, breaks = lv, name = "Annotation", drop = TRUE) +
  scale_x_continuous(breaks = seq(-10, 10, 5), labels = function(x) sprintf("%.1f", x)) +
  scale_y_continuous(breaks = seq(0, Y_MAX, if (Y_MAX > 4) 1 else 0.5),
                     labels = function(y) sprintf(if (Y_MAX > 4) "%.0f" else "%.1f", y)) +
  coord_cartesian(xlim = c(-X_LIM, X_LIM), ylim = c(0, Y_MAX)) +
  labs(title = "Volcano Plot: Upregulated and Downregulated Genes",
       subtitle = sprintf("%s < %s and absolute log2 fold change >= %s; Verruculogen is red",
                          STAT_LABEL, FDR_CUTOFF, FC_CUTOFF),
       x = expression(log[2]~"(Fold Change)"),
       y = if (STAT == "p") expression(-log[10]~"("*italic(P)*")")
           else expression(-log[10]~"(FDR)")) +
  theme_minimal(base_size = 12, base_family = FIG_FONT) +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = .8),
        plot.title = element_text(face = "bold", hjust = .5, size = 15),
        plot.subtitle = element_text(hjust = .5, colour = "#52514e", size = 12),
        axis.title = element_text(size = 15, face = "bold", colour = "black"),
        axis.text  = element_text(size = 14, face = "bold", colour = "black"),
        axis.line  = element_line(colour = "black", linewidth = .8),
        axis.ticks = element_line(colour = "black", linewidth = .8),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11), legend.position = "right",
        legend.key.size = grid::unit(7, "mm"),
        plot.margin = margin(4, 4, 4, 4, unit = "mm")) +
  guides(colour = guide_legend(override.aes = list(size = 3.6)))

## gene-ID labels for the NRPS genes, placed in the empty lower-left corner so
## that two independent repel layers cannot collide
if (nrow(nrps)) {
  idl <- nrps %>% filter(y <= Y_MAX) %>% arrange(desc(y)) %>%
    mutate(lab_x = -X_LIM * 0.86,
           lab_y = seq(from = Y_MAX * 0.34, by = -Y_MAX * 0.125, length.out = n()))
  p <- p +
    geom_segment(data = idl, aes(x = lab_x + 2.75, y = lab_y, xend = logFC, yend = y),
                 colour = "#444444", linewidth = .32) +
    geom_text(data = idl, aes(lab_x, lab_y, label = gene_ID), colour = "#111111",
              size = 5.3, family = FIG_FONT, hjust = 0)
}

out <- file.path(FIG_DIR, sprintf("volcano_labstyle_%s_%s", d$short, STAT))
ggsave(paste0(out, ".png"), p, width = 200, height = 160, units = "mm",
       dpi = 600, bg = "white")
ggsave(paste0(out, ".pdf"), p, width = 200, height = 160, units = "mm", bg = "white")
message("Saved ", out, ".png / .pdf")
message(sprintf("Up: %d   Down: %d   ftm above the line: %d/%d   y axis 0-%.1f",
                sum(df$is_up), sum(df$is_down), sum(eg$stat < FDR_CUTOFF), nrow(eg), Y_MAX))
