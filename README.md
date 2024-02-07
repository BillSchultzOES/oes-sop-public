This repository contains the Standard Operating Procedures for statistical analyses done by the Office of Evaluation Sciences in the US General Services Administration. It is currently built using R Markdown and [bookdown](https://github.com/rstudio/bookdown). Specifically, this version of the OES SOP was built using `R` 4.2.1 (2022-06-23 ucrt), `RStudio` 1.2.1335, and `pandoc` 2.6. Math is rendered using the `katex` package 1.4.1.

# Package management

To control package versions used in the SOP, we're relying on [renv](https://rstudio.github.io/renv/articles/renv.html). See that vignette for more details on what `renv` does and how it works.

Briefly, `renv` maintains a curated set packages local to this R project (in the "renv" folder of this repository). The file "renv.lock" (the "lock file") in this repository contains a snapshot of the packages that should be used to compile the SOP, along with the version we want of each. You'll need to use `renv` to install the right packages and build the SOP on your machine, or add any new packages for use in one of the chapters. But this doesn't make updating the SOP much more complicated.

After cloning the most recent version of this repository onto your machine, opening the R project should prompt `renv` to automatically install itself. Once installed, if it doesn't load automatically (you'll see a message in the console), run `renv::load()`. After this, run `renv::restore` while the R project is open to automatically prepare the right set of packages/versions in the project directory on your machine, based on the instructions from the lock file. If you run into a message saying `renv` isn't activated yet when trying to restore from the lock file, make sure youv'e run `renv::activate` and `renv::load` first.

Once you've opened the R project and loaded `renv`, you can install packages as normal, and they'll be installed in the project directory rather than your machine's main R package directory. After you complete the install, add this package to the list in "index.Rmd" and run `renv::snapshot()` to update the lock file. Make sure you commit the new lock file to Github alonside your other changes (also be sure to commit ".Rprofile," "renv/settings.json," and "renv/activate.R"). Committing the entire repository when you want to update Github with your changes will accomplish this.

# Important package versions

This project requires the `bookdown` package for R. We are using version .7 for now.

```
library(devtools)
install_version("bookdown",version=0.7)
```

We also need to use version 1.33 of `knitr`.

```
install_version("knitr", version = 1.33)
```

Our current use of `renv` should take care of all of these version control issues for us, and help prevent new ones if any packages change in the future in a way that would otherwise break our code. But this section of the README has been retained just in case we decide to drop `renv` in the future.

# Adding or removing chapters

To add or remove chapters, once the `.rmd` file for that chapter is finished, open "_bookdown.yml"" and add its name to the list (or take it's name out of the list). Be sure to update the numbers at the beginning of chapter names accordingly. To help keep things clear, please only assign numbers to the names of *live* chapters. Also, be sure to keep the glossary, appendix, and references at the end of this list.

# How can I update content?

* If you want to update these instructions, open "README.md".

* If you want to update the introductory page, open "index.Rmd".

* If you want to update the content of a particular chapter, open one of the numbered .Rmd files (e.g.: "01-causalinference.Rmd", "02-basics.Rmd", etc.). You can see a list of live chapters in "_bookdown.yml".

  + To comment-out text OUTSIDE of a code chunk, highlight it press `cntrl+shift+c`
  
  + To hide a code chunk from the SOP, set the following options: `echo = F, eval = F`
  
  + There are a number of work-in-prgress sections or code chunks from earlier drafts that are currently hidden in this way in case we want to return to them later 
  
* If you want to update how the tabbed code chunks look, open "style.css". This requires working with CSS, though.

# More detail on code chunks

To include the "copy-to-clipboard" button in a new chapter (using the `klippy` package), include the following in an R code chunk near the beginning. This is already in all existing chapters.

```
klippy::klippy(all_precode = T, position = c("top", "right"))
```

## R code

To add a normal R code chunk (without the Stata option), just add a chunk as in any other R Markdown file. There is CSS/Javascript/HTML on the back-end that automatically creates a "Code/Hide" button for all R code chunks. You don't need to do anything to generate this.

In case this is useful context for others: that back-end code finds the R code chunks based on HTML tags that appear to be linked to the `highlight = T` option in R code chunks (defaults to `T`). Setting this to `F` instead for a particular code chunk, therefore, turns off the "Code/Hide" button for that chunk.

## Tabs with R and Stata options

Adding a code chunk with tabs for R/Stata options is a bit tricker (take a look at existing chapters first to get a sense of what this looks like). Note that these tabbed chunks are *just for illustrative purposes* and aren't actually being run (they all have `eval = F`). After every tabbed chunk, you'll see another hidden (`echo = F`) R code chunk that is actually being run and generating results/figures shown in the SOP (`eval = T`). These hidden chunks just repeat the R code from the tabbed chunks.

The primary challenge to adding tabbed chunks is making sure each tab has a unique HTML ID (e.g., "ch3R1," "ch3Stata1," and "ch3Hide1" in the first tabbed chunk of chapter 3). If they don't, only the first instance of a duplicated ID will work properly; the second with this ID just won't open and show the code. Luckily, we've developed a keyboard macro that automates HTML ID generation without you needing to keep track manually. This macro also a inserts template, so there's no need to do any copy-pasting (except copying the right R and Stata code into that template)!

There are two steps to set this up. First, around when we first start using R code in any existing chapter, you'll notice the following code chunk:

```
# cnum is modified automatically to iteratively count chunks
# when using the oes_code_tab markdown snippet. each use of
# the snippet adds a value of 1.
ch <- 3
cnum <- 0
```

This creates two R objects in the environment: a chapter number and a running chunk number, the latter initialized at 0 at the beginning of the chapter. If you're creating a new chapter, add this chunk near the top. If you're modifying an existing chapter, run these lines of code with `cnum` set higher than any existing chunk numbers, instead of setting it to 0 (it just needs to be unique, not actually sequential). Scroll to the last tabbed chunk in that chapter to see what the current count is. E.g., in chapter 3, the tabs in the last tabbed chunk are numbered 12 ("ch3R12").

Second, in RStudio, go to "Tools/Global Options/Code/Edit Snippets/Markdown". At the bottom of the list of Markdown snippets, copy the following:

```
snippet {oes_code_tab}
	`r cnum <- cnum + 1`
	::: {.tab} 
	<button class="tablinks" onclick="unrolltab(event, '`r paste0("ch", ch, "R", cnum)`')">R code</button>
	<button class="tablinks" onclick="unrolltab(event, '`r paste0("ch", ch, "Stata", cnum)`')">Stata code</button>
	<button class="tablinks" onclick="unrolltab(event, '`r paste0("ch", ch, "Hide", cnum)`')">Hide</button>
	::: {#`r paste0("ch", ch, "R", cnum)` .tabcontent} 
	<br />
	```{r, highlight = F, eval = F}
	# R code here
	```
	:::
	::: {#`r paste0("ch", ch, "Stata", cnum)` .tabcontent} 
	<br />
	```{stata, highlight = F, eval = F}
	* Stata code here
	```
	::: 
	::: {#`r paste0("ch", ch, "Hide", cnum)` .tabcontent}
	::: 
	:::
```
Code snippets are keyboard macros you can use to more easily insert a common chunk of code without having to re-type it every time. Once `ch` and `cnum` exist in R, this macro will automatically generate a template tabbed chunk with a new, sequential HTML ID (ch = current ch, cnum = current cnum + 1). All you have to do now is type:

```
{oes_code_tab}
```

And then hit `shift+tab` when your cursor is to the right of `}`. Then copy in R and Stata code where the comments indicate (# R code here, * Stata code here).

# Building the book

In theory, there is more than one way to compile this book to .html, .pdf, and .epub. But right now, the first option works, while the second option runs into errors that prevent compling. This may be an artifact of our need to rely on older versions of certain packages. This may be worth addressing once content is smoothed out and we have the new first online.

- (**Working**) Open the sop.Rproj file in RStudio. There should be a "Build" tab in the same pane as the "Environment" tab. Under "Build Book" in this tab, click on `bookdown::gitbook`, which is the format we currently have displayed online as well.

- (**Currently not working**) You can use `render_book()` as described [here](https://bookdown.org/yihui/bookdown/build-the-book.html)

- RStudio also has some [add-ins](https://bookdown.org/yihui/bookdown/rstudio-ide.html) that make bookdown easier to use

# Stata code examples

The "Stata code" folder in this repository contains .do files you can use to actually run the Stata code in different chapters (with all the necessary input data saved as .csv files). It may be useful to provide similar R scripts for each chapter in the future. Available for Fellows who are interested in running the Stata code themselves.
