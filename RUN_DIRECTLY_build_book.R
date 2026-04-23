
# Make sure the R environment w/ requisite
# packages and versions is loaded
if (is.null(renv::project())) {
  renv::load()
}

# Render the book within this environment
bookdown::render_book("index.Rmd", "bookdown::gitbook")

# Remove vulnerable jquery packaged with gitbook
# (need to eventually migrate to bs4).
# First, let's get path to older (vulnerable) versions:
old_paths <- Sys.glob("docs/libs/*/jquery*.js")
replace <- c()
for (p in old_paths) {

  # Check whether version number in first line is 1 or 2.
  # (OK if version starts with 3).
  check <- grepl(" jQuery v2\\.| jQuery v1\\.", readLines(con = p, n = 1))

  # Add any with an old version to the replace list
  if (check) {
    replace <- c(replace, p)
  }

}

# Replace with a modern version (no vulnerability).
# Bandaid over the real problem (need to migrate from gitbook).
tmp <- tempfile(fileext = ".js")
download.file(
  "https://code.jquery.com/jquery-3.7.1.min.js",
  destfile = tmp,
  quiet = TRUE
  )
for (r in replace) {
  file.copy(tmp, old_paths, overwrite = TRUE)
}
