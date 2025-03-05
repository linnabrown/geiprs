# GEiPRS - A Fast and Powerful Machine Learning Method for Polygenic Risk Score Prediction by Leveraging Genotype-Environment Interactions

License: GPL-2

### References: 
  - Le Huang#, Wujuan Zhong#, Song Zhai, and Judong Shen. GEiPRS: A Fast and Powerful Machine Learning Method for Polygenic Risk Score Prediction by Leveraging Genotype-Environment Interactions. 

### Installation:
GEiPRS can be installed by running the following scripts in R. Notice that the installation of `pgenlibr` requires [zstd(>=1.4.4)](https://github.com/facebook/zstd). It can be built from source or simply available from [conda](https://anaconda.org/conda-forge/zstd), [pip](https://pypi.org/project/zstd/) or [brew](https://formulae.brew.sh/formula/zstd).

```r
library(devtools)

devtools::install_github("junyangq/glmnetPlus")
devtools::install_github("chrchang/plink-ng", subdir="/2.0/cindex")
devtools::install_github("chrchang/plink-ng", subdir="/2.0/pgenlibr")
devtools::install_github("dajmcdon/sparsegl")

# install GEiPRS
devtools::install_github("linnabrown/geiprs")
```
We assume the users already have PLINK 2.0. Otherwise it can be installed from https://www.cog-genomics.org/plink/2.0/.


Tutorial of using R package GEiPRS is available at `vignettes/vignette.pdf` and `vignettes/vignette.html`.


## Acknowledgments
This package builds upon the [`snpnet`](https://github.com/junyangq/snpnet.git) package developed by Junyang Qian, Trevor Hastie, Yosuke Tanigawa, Ruilin Li, Manuel A Rivas, and Christopher Chang. We thank them for their contributions. 
