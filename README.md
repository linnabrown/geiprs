# GEiPRS - A Fast and Powerful Machine Learning Method for Polygenic Risk Score Prediction by Leveraging Genotype-Environment Interactions

License: GPL-2

### References: 
  - Le Huang#, Wujuan Zhong#, Song Zhai, and Judong Shen. GEiPRS: A Fast and Powerful Machine Learning Method for Polygenic Risk Score Prediction by Leveraging Genotype-Environment Interactions. 

### Installation:
Most of the requirements of geiprs are available from CRAN. It also depends on the `pgenlibr`, `glmnet/glmnetPlus` and `cindex` (for survival analysis) packages. One can install them by running the following commands in R. Notice that the installation of `pgenlibr` requires [zstd(>=1.4.4)](https://github.com/facebook/zstd). It can be built from source or simply available from [conda](https://anaconda.org/conda-forge/zstd), [pip](https://pypi.org/project/zstd/) or [brew](https://formulae.brew.sh/formula/zstd).

```r
library(devtools)

devtools::install_github("junyangq/glmnetPlus")
devtools::install_github("chrchang/plink-ng", subdir="/2.0/cindex")
devtools::install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
devtools::install_github("dajmcdon/sparsegl")

# install geiprs
# devtools::install_github("linnabrown/geiprs")
devtools::install("/nas/longleaf/home/lehuang/geiprs")


# longleaf environment in UNC

module load r/4.1.3


```
We assume the users already have PLINK 2.0. Otherwise it can be installed from https://www.cog-genomics.org/plink/2.0/.


# Example codes
library(stringr)
library(geiprs)
library(dplyr)
library(data.table)
configs <- list(
  verbose = TRUE,
  num.groups.batch = 100,
  nCores = 8,
  memory = 16000,
  save = TRUE,
  results.dir = "results",
  save.dir = "intermediate",
  plink2.path = "/nas/longleaf/home/lehuang/geiprs/other_files/plink2"
)

output <- geiprs(
  genotype.pfile = "/nas/longleaf/home/lehuang/geiprs/inst/extdata/sample",
  phenotype.file = "/nas/longleaf/home/lehuang/geiprs/inst/extdata/sample.phe",
  phenotype = "QPHE",
  env = "sex",
  tau = 0.9,
  family = "gaussian",
  configs = configs
)

## Performance

```
[2025-01-31 11:57:55 snpnet-ge] End snpnet. Time elapsed: 3.1815 mins
           used  (Mb) gc trigger  (Mb) max used  (Mb)
Ncells  2012616 107.5    3049137 162.9  3049137 162.9
Vcells 12238547  93.4   31174882 237.9 25911727 197.7
```