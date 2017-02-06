## Installation

On Windows systems, you must first have 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)
installed and correctly configured.

The package and its dependencies are installed with the following
commands:

```r
install.packages(c('Rcpp', 'numDeriv', 'devtools', 'knitr'), repos = 'https://cloud.r-project.org/')
devtools::install_github('usnistgov/potMax', build_vignettes = TRUE)
```

Building the vignette takes a while so be patient. It is worth the
time because the vignette contains information about the underlying
statistical methods as well as instructions for using the functions
provided by the package.

## Use

After installation, see the vignette called *potMax_details*.

```r
library(potMax)
vignette('potMax_details')
```