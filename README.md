This is the repository for statistical standard operating procedures for the
OES. Currently in R Markdown and **bookdown** (https://github.com/rstudio/bookdown).

You can find the preview of this example at <http://gsa-oes.github.com/sop>.

# To make changes and build

Right now this project requires the bookdown package for R. We have to use version .7 for now.

```
library(devtools)
install_version("bookdown",version=0.7)
```

We also need to use version 1.33 of knitr.

```
install_version("knitr", version = 1.33)
```

# To update the web preview:

First, commit and push any changes to the _book subdirectory within the Book
directory. (That is where `render_bookdown()` creates the html files from the
Rmd files).

Second, move to the root directory (i.e. `cd ..` on a Unix machine at the
command line).

Third, the following steps work on Unix based machines:

```
rm -rf book-output
git clone -b gh-pages git@github.com:gsa-oes/sop.git book-output
cd book-output
cp -r ../Book/_book/* ./
git add --all *
git commit -m"Update the book" || true
git push -q origin gh-pages

```
