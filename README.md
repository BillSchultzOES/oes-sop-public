This is the repository for statistical standard operating procedures for the
OES. Currently in R Markdown and **bookdown** (https://github.com/rstudio/bookdown).

You can find the preview of this example at <http://gsa-oes.github.com/sop>.

# Package management

Maintaining the SOP has been made difficult in the past by inconsistent package versions across the computers of those building it, or by occasional package updates breaking things. To fix the packages used in time and make the document more stable, we'll try using [renv](https://rstudio.github.io/renv/articles/renv.html). See that vignette for more details on what it does and how it works.

In brief, `renv` helps us maintain a curated set of installed packages that is local to this R project. The file "renv.lock" (the "lock file") contains a snapshot of the packages that currently should be used to compile the SOP, including the version of each that should be installed. You'll have to work with `renv` to install the right packages and build the SOP on your machine, or if you want to add any new packages. But it doesn't actually make updating the SOP much more complicated.

After cloning the most recent version of this repository onto your machine, opening the project should prompt `renv` to automatically install itself on your computer. After this, running `renv::restore` while the R project is open will automatically prepare the right set of packages/versions in the project directory on your machine, based on the instructions from the lock file. If you run into a message saying `renv` isn't activated yet when trying to restore from the lock file, run `renv::activate` or `renv::load`. This should happen automatically according to the vignette, but I had to run `renv::load` to get it to activate properly (not sure if this is related to our security setup).

If you ever install a new package to include in a chapter of the SOP, run `renv::snapshot()` to update the lock file, and make sure you commit the new lock file to Github alonside your other changes (also be sure to commit ".Rprofile," "renv/settings.json," and "renv/activate.R"). Then, others can use `renv::restore` to make sure this package/version is used when they build the updated SOP on their own machines.

# To make changes and build

This project requires the `bookdown` package for R. We have to use version .7 for now.

```
library(devtools)
install_version("bookdown",version=0.7)
```

We also need to use version 1.33 of `knitr`.

```
install_version("knitr", version = 1.33)
```

Finally, it may be necessary to install a development version of the `bfe` package, depending on the version of R you have loaded.

```
devtools::install_github('gibbonscharlie/bfe')
```

The current reliance `renv` should take care of all of those version control issues for us, and help prevent new ones if any packages change in the future in a way that breaks our code. But this section of the README has been retained just in case we decide to drop `renv` in the future.

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
