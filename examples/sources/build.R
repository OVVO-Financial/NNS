#!/usr/bin/env Rscript
# Rebuild the example write-ups and pages from their R Markdown sources.
#
#   Rscript examples/sources/build.R                  # every source
#   Rscript examples/sources/build.R Normalization    # one or more sources by name
#   Rscript examples/sources/build.R --html           # HTML previews instead of PDFs
#   Rscript examples/sources/build.R --out DIR ...    # write somewhere other than examples/
#
# Each source renders to its published file in examples/ (the name in `outputs`
# below). PDF write-ups need a LaTeX installation (for example
# tinytex::install_tinytex()); without one, or with --html, they render to HTML
# next to the PDF instead. Run from the repository root with NNS installed.

args <- commandArgs(trailingOnly = TRUE)
html_only <- "--html" %in% args
args <- setdiff(args, "--html")
out_dir <- NULL
if ("--out" %in% args) {
  i <- match("--out", args)
  out_dir <- args[i + 1L]
  args <- args[-c(i, i + 1L)]
}

here <- if (file.exists("examples/sources/build.R")) "examples/sources" else "."
if (is.null(out_dir)) out_dir <- normalizePath(file.path(here, ".."))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Source name -> published file.
outputs <- c(
  "5_way_rank_probabilities_using_R"                         = "5_way_rank_probabilities_using_R.pdf",
  "Bayes_Theorem_From_Partial_Moments"                       = "Bayes' Theorem From Partial Moments.pdf",
  "Causal_Inference_with_NNS_stack"                          = "Causal_Inference_with_NNS_stack.pdf",
  "Continuous_CDFs_and_ANOVA_with_NNS"                       = "Continuous_CDFs_and_ANOVA_with_NNS.pdf",
  "Machine_Learning"                                         = "Machine_Learning.pdf",
  "NNS_ARMA"                                                 = "NNS ARMA.pdf",
  "NNS_vs_KNN_MNIST"                                         = "NNS vs KNN MNIST dataset.pdf",
  "Nonlinear_Correlation_and_Dependence_Using_NNS"           = "Nonlinear_Correlation_and_Dependence_Using_NNS.pdf",
  "Nonlinear_FWL"                                            = "Nonlinear_FWL.pdf",
  "Normalization"                                            = "Normalization.pdf",
  "Partitional_Estimation_Using_Partial_Moments"             = "Partitional_Estimation_Using_Partial_Moments.pdf",
  "Spiral_Matching_Example"                                  = "Spiral Matching Example.pdf",
  "tides"                                                    = "tides.html",
  "Sunspots_example"                                         = "Sunspots_example.html",
  "VAR_example"                                              = "VAR_example.html",
  "Causal_Inference_Amongst_Macroeconomic_Variables_Using_NNS" = "Causal_Inference_Amongst_Macroeconomic_Variables_Using_NNS.html",
  "prophet_NNS_comparison"                                   = "prophet_NNS_comparison.html"
)

todo <- if (length(args)) args else names(outputs)
unknown <- setdiff(todo, names(outputs))
if (length(unknown)) stop("Unknown source(s): ", paste(unknown, collapse = ", "))

has_latex <- nzchar(Sys.which("pdflatex")) ||
  (requireNamespace("tinytex", quietly = TRUE) && tinytex::is_tinytex())

for (name in todo) {
  input <- file.path(here, paste0(name, ".Rmd"))
  target <- outputs[[name]]
  pdf <- grepl("\\.pdf$", target)
  if (pdf && (html_only || !has_latex)) {
    if (!html_only) message("No LaTeX found: rendering ", name, " to HTML.")
    format <- "html_document"
    target <- sub("\\.pdf$", ".html", target)
  } else {
    format <- if (pdf) "pdf_document" else "html_document"
  }
  message("Rendering ", input, " -> ", file.path(out_dir, target))
  # Each source renders in its own R session, so packages one source attaches
  # cannot affect another or appear in its sessionInfo().
  script <- tempfile(fileext = ".R")
  writeLines(sprintf(paste("rmarkdown::render(%s, output_format = %s, output_file = %s,",
                           "output_dir = %s, intermediates_dir = tempdir(), quiet = TRUE)"),
                     deparse(normalizePath(input)), deparse(format), deparse(target),
                     deparse(normalizePath(out_dir))), script)
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script))
  unlink(script)
  if (!identical(status, 0L)) stop("Rendering ", name, " failed.")
}
