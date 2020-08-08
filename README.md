# Phenotype

Large-scale phenotypic data processing is essential in research. Researchers need to eliminate outliers from the data in order to obtain true and reliable results. Best linear unbiased prediction (BLUP) is a standard method for estimating random effects of a mixed model. This method can be used to process phenotypic data under different conditions and is widely used in animal and plant breeding. The 'Phenotype' can remove outliers from phenotypic data and performs the best linear unbiased prediction (BLUP), help researchers quickly complete phenotypic data analysis.

## Installation

``` r
install.packages("Phenotype")
```

## Usage

### Remove outliers from phenotypic data

``` r
library("Phenotype")
data("wheatds")
inlier <- outlier(wheatds, sample = "Line", loc = "Env", rep = "Rep", phe = "DS", mode = "blup")
```

### Calculate statistical indicators of phenotypic data

``` r
data("wheatds")
inlier <- outlier(wheatds, sample = "Line", loc = "Env", rep = "Rep", phe = "DS", mode = "blup")
stat_out <- stat(x = inlier, sample = "Sample", phe = "inlier")
```

### Histogram drawing

``` r
data("wheatds")
inlier <- outlier(wheatds, sample = "Line", loc = "Env", rep = "Rep", phe = "DS", mode = "blup")
stat_out <- stat(x = inlier, sample = "Sample", phe = "inlier")
histplot(x = stat_out$mean)
```

### Performs the Best Linear Unbiased Prediction (BLUP)

``` r
data("wheatds")
blup_out <- blup(wheatds, sample = "Line", loc = "Env", rep = "Rep", phe = "DS")
```

## Contact
For any bugs/issues/suggestions, please send emails to: Peng Zhao pengzhao@nwafu.edu.cn.
