# Run from the repository root. Generated HTML is written beside its source.
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
    !requireNamespace("knitr", quietly = TRUE)) {
  stop("Install the suggested knitr and rmarkdown packages to render the companion.")
}
if (!rmarkdown::pandoc_available()) stop("Pandoc is required for the HTML build.")
local({
  old <- setwd("examples")
  on.exit(setwd(old), add = TRUE)
  knitr::opts_chunk$set(fig.path = "revised-companions_files/figure-html/")
  knitr::knit("revised-companions.Rmd", output = "revised-companions.md",
              quiet = TRUE, envir = new.env(parent = globalenv()))
  rmarkdown::pandoc_convert("revised-companions.md", to = "html",
                            output = "revised-companions.html",
                            options = c("--standalone", "--toc", "--mathjax"))
})
