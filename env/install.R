## Install the R packages the analysis needs, pinned to the versions it was run
## with.
##
##   Rscript env/install.R              # pinned versions (reproducible)
##   Rscript env/install.R --latest     # whatever is current (convenient)
##
## The pins matter. limma's treat() and edgeR's filterByExpr() have both changed
## behaviour between releases, so an unpinned install can reproduce the code
## without reproducing the numbers. env/sessionInfo.txt records the full session
## the published results came from.

PINS <- c(limma = "3.58.1", edgeR = "4.0.16", statmod = "1.5.0",
          ggplot2 = "3.4.4", ggrepel = "0.9.5", dplyr = "1.1.4")
BIOC_VERSION <- "3.18"          # the Bioconductor release limma 3.58.1 belongs to

latest <- "--latest" %in% commandArgs(trailingOnly = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

if (latest) {
  BiocManager::install(c("limma", "edgeR"), ask = FALSE, update = FALSE)
  cran <- c("statmod", "ggplot2", "ggrepel", "dplyr")
  missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) install.packages(missing)
} else {
  BiocManager::install(version = BIOC_VERSION, ask = FALSE, update = FALSE)
  BiocManager::install(c("limma", "edgeR"), ask = FALSE, update = FALSE)
  for (p in c("statmod", "ggplot2", "ggrepel", "dplyr")) {
    have <- requireNamespace(p, quietly = TRUE) &&
            as.character(packageVersion(p)) == PINS[[p]]
    if (!have) remotes::install_version(p, version = PINS[[p]], upgrade = "never")
  }
}

cat("\nInstalled. Versions in use (pinned version in brackets):\n")
ok <- TRUE
for (p in names(PINS)) {
  v <- tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING")
  flag <- if (v == PINS[[p]]) " " else "*"
  if (v != PINS[[p]]) ok <- FALSE
  cat(sprintf(" %s %-10s %-10s [%s]\n", flag, p, v, PINS[[p]]))
}
if (!ok)
  cat("\n* marks a version that differs from the one the published results were\n",
      "  produced with. The analysis will still run; small numerical differences\n",
      "  are possible. See env/sessionInfo.txt.\n", sep = "")
