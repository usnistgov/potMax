## Installation

```r
install.packages(c('Rcpp', 'numDeriv', 'devtools'))
devtools::install_github('usnistgov/potMax', build_vignettes = TRUE)
```

Building the vignette takes a while so be patient. It is worth the
time because the vignette contains information about the underlying
statistical methods as well as instructions for using the functions
provided by the package.

### Windows Note

On Windows systems, you must first have 
[Rtools](https://cran.r-project.org/bin/windows/Rtools/)
installed and correctly configured.

## Use

After installation, see the vignette called *potMax_details*.

```r
library(potMax)
vignette('potMax_details')
```