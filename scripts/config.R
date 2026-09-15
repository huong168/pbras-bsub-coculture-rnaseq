## ---------------------------------------------------------------------------
## Shared configuration. Every R script in this repository starts with
##   source("<repo>/scripts/config.R")
## The repository root is found from the .repo-root marker, so the scripts run
## from any working directory. No absolute paths appear anywhere in the code.
## ---------------------------------------------------------------------------

## Locate the repository root by walking up from a starting directory until the
## .repo-root marker is found. The starting point is the script's own location
## when one is available (Rscript passes it as --file=), so the scripts work
## from any working directory; interactive sessions fall back to getwd(), which
## must then be somewhere inside the repository. PBRAS_ROOT overrides both.
script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f[1], mustWork = FALSE)))
  if (!is.null(sys.frames()[[1]]$ofile))
    return(dirname(normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE)))
  NULL
}

walk_up <- function(start) {
  p <- tryCatch(normalizePath(start, mustWork = TRUE), error = function(e) NULL)
  if (is.null(p)) return(NULL)
  repeat {
    if (file.exists(file.path(p, ".repo-root"))) return(p)
    up <- dirname(p)
    if (identical(up, p)) return(NULL)
    p <- up
  }
}

find_repo_root <- function() {
  for (start in c(script_dir(), getwd())) {
    r <- walk_up(start)
    if (!is.null(r)) return(r)
  }
  stop("Could not locate the repository root (.repo-root marker). ",
       "Set PBRAS_ROOT to the repository directory.")
}

ROOT    <- Sys.getenv("PBRAS_ROOT", unset = find_repo_root())
DATA    <- file.path(ROOT, "data")
COUNTS  <- file.path(DATA, "counts")
ANNOT   <- file.path(DATA, "annotation")
RESULTS <- file.path(ROOT, "results")
DE_DIR  <- file.path(RESULTS, "de_tables")
GS_DIR  <- file.path(RESULTS, "geneset")
FIG_DIR <- file.path(RESULTS, "figures")
for (d in c(RESULTS, DE_DIR, GS_DIR, FIG_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## --- Sample metadata --------------------------------------------------------
## control   = co-culture with the bacillibactin-deficient dhbF knock-out
## treatment = co-culture with B. subtilis DK1042 wild type
## Contrast is treatment - control, so logFC < 0 means "reduced in the presence
## of the wild type".
SAMPLES <- data.frame(
  sample = c("KO1", "KO3", "WT1", "WT2"),
  file   = c("htseq_KO1.tab", "htseq_KO3.tab", "htseq_WT1.tab", "htseq_WT2.tab"),
  group  = c("control", "control", "treatment", "treatment"),
  label  = c("KO-1", "KO-3", "3610-1", "3610-2"),
  ## Percentage of reads pseudo-aligning to the P. brasilianum genome, as
  ## measured for the quantification distributed in data/counts/. These are
  ## overwritten below if data/annotation/fungal_read_fraction.csv exists,
  ## which kallisto_quant.py writes from kallisto's own run_info.json -- so a
  ## re-quantification at full depth, or with a strandedness setting,
  ## updates the QC automatically instead of silently reusing these numbers.
  fungal_pct = c(19.28, 18.58, 57.80, 77.75),
  stringsAsFactors = FALSE
)

## --- The design ------------------------------------------------------------
DESIGNS <- list(
  d22a = list(id = "d22a", samples = c("KO1", "KO3", "WT1", "WT2"),
              label = "2 WT vs 2 KO", short = "2WTx2KO")
)
DESIGN <- DESIGNS$d22a          # the single design analysed in this repository

## Replace the recorded fungal read fractions with measured ones when a
## re-quantification has produced them.
local({
  f <- file.path(ANNOT, "fungal_read_fraction.csv")
  if (!file.exists(f)) return(invisible(NULL))
  m <- read.csv(f, stringsAsFactors = FALSE)
  if (!all(c("sample", "fungal_pct") %in% names(m)))
    stop(f, " must have columns 'sample' and 'fungal_pct'.")
  hit <- match(SAMPLES$sample, m$sample)
  if (anyNA(hit)) {
    warning("fungal_read_fraction.csv is missing: ",
            paste(SAMPLES$sample[is.na(hit)], collapse = ", "),
            " -- keeping the recorded values for those.")
    hit <- hit[!is.na(hit)]
  }
  SAMPLES$fungal_pct[!is.na(match(SAMPLES$sample, m$sample))] <<- m$fungal_pct[hit]
  message("fungal_pct read from ", f)
})

## --- The pre-specified hypothesis: the verruculogen (ftm) cluster -----------
FTM <- c(PMG11_03146 = "FtmA", PMG11_03147 = "FtmC", PMG11_03148 = "FtmD",
         PMG11_03149 = "FtmB", PMG11_03150 = "FtmE", PMG11_03151 = "FtmF",
         PMG11_03152 = "FtmG", PMG11_03153 = "FtmH")
FTM_ORDER <- c("FtmA","FtmB","FtmC","FtmD","FtmE","FtmF","FtmG","FtmH")

## --- Thresholds used throughout ---------------------------------------------
LFC_TREAT  <- 1     # null hypothesis of treat(): |log2FC| <= 1
FC_CUTOFF  <- 2     # fold-change cutoff used for colouring the volcano plots
FDR_CUTOFF <- 0.05

## --- Shared figure style ----------------------------------------------------
FIG_FONT <- if ("Arial" %in% names(grDevices::postscriptFonts())) "Arial" else "sans"
PAL <- c(Verruculogen = "#d92b2b", Upregulated = "#f5cfc8",
         Downregulated = "#bed3e4", NRPS = "#7207ec", PKS = "#8e44ad",
         Terpene = "#e67e22", RiPP = "#00a6a6", Betalactone = "#8c6d31",
         Siderophore = "#d4a017")
