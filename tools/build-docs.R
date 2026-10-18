# Build website sources separately from the short repository README.
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
stage <- tempfile("frpp-docs-")
dir.create(stage)
files <- c("DESCRIPTION", "NAMESPACE", "README.md", "NEWS.md", "_pkgdown.yml",
           "R", "man", "vignettes", "pkgdown")
if (!all(file.copy(file.path(root, files), stage, recursive = TRUE)))
  stop("Could not stage the documentation sources.")
articles <- file.path(stage, "vignettes", "articles")
dir.create(articles)
article_sources <- list.files(file.path(root, "vignettes_ext"),
                              pattern = "[.]Rmd$", full.names = TRUE)
if (!length(article_sources) || !all(file.copy(article_sources, articles)))
  stop("Could not stage the website guides.")
destination <- Sys.getenv("FRPP_DOCS_DEST", file.path(root, "docs"))
pkgdown::build_site(pkg = stage, override = list(destination = destination),
  preview = FALSE, new_process = FALSE, install = FALSE)
for (article in basename(article_sources)) {
  page <- file.path(destination, "articles", sub("[.]Rmd$", ".html", article))
  html <- readLines(page, warn = FALSE, encoding = "UTF-8")
  html <- gsub("/vignettes/articles/", "/vignettes_ext/", html, fixed = TRUE)
  writeLines(html, page, useBytes = TRUE)
}
