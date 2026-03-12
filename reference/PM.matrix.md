# Partial Moment Matrix

Builds a list containing CUPM, DUPM, DLPM, CLPM and the overall
covariance matrix.

## Usage

``` r
PM.matrix(LPM_degree, UPM_degree, target, variable, pop_adj, norm = FALSE)
```

## Arguments

- LPM_degree:

  numeric; lower partial moment degree (0 = freq, 1 = area).

- UPM_degree:

  numeric; upper partial moment degree (0 = freq, 1 = area).

- target:

  numeric vector; thresholds for each column (defaults to colMeans).

- variable:

  numeric matrix or data.frame.

- pop_adj:

  logical; TRUE adjusts population vs. sample moments.

- norm:

  logical; default FALSE. If TRUE, each quadrant matrix is cell-wise
  normalized so their sum is 1 at each (i,j).

## Value

A list: \$cupm, \$dupm, \$dlpm, \$clpm, \$cov.matrix.

## Details

Partial Moment Matrix

## Note

When `norm = TRUE`, each cell (i,j) of the four quadrant matrices is
normalized so that their sum equals 1. In this case, `$cov.matrix` is
computed as `$cupm + $clpm - $dupm - $dlpm`, yielding a dimensionless,
signed dependence measure bounded between -1 and 1. This representation
discards magnitude information and is therefore a lossy nonlinear
correlation matrix. A higher fidelity nonlinear correlation matrix is
available via the `NNS.dep` function.

## Examples

``` r
set.seed(123)
A <- cbind(rnorm(100), rnorm(100), rnorm(100))

# Uses norm = FALSE by default
PM.matrix(1, 1, target = NULL, variable = A, pop_adj = TRUE)
#> $cupm
#>           [,1]      [,2]      [,3]
#> [1,] 0.4299250 0.1033601 0.1338273
#> [2,] 0.1033601 0.5411626 0.1339439
#> [3,] 0.1338273 0.1339439 0.5065920
#> 
#> $dupm
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.1469182 0.2031944
#> [2,] 0.1560924 0.0000000 0.1345625
#> [3,] 0.1504256 0.1492284 0.0000000
#> 
#> $dlpm
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.1560924 0.1504256
#> [2,] 0.1469182 0.0000000 0.1492284
#> [3,] 0.2031944 0.1345625 0.0000000
#> 
#> $clpm
#>           [,1]      [,2]      [,3]
#> [1,] 0.4033078 0.1559295 0.1077888
#> [2,] 0.1559295 0.3939005 0.1779344
#> [3,] 0.1077888 0.1779344 0.3956781
#> 
#> $cov.matrix
#>             [,1]        [,2]        [,3]
#> [1,]  0.83323283 -0.04372107 -0.11200394
#> [2,] -0.04372107  0.93506310  0.02808746
#> [3,] -0.11200394  0.02808746  0.90227005
#> 

# Enable normalization
PM.matrix(1, 1, target = NULL, variable = A, pop_adj = TRUE, norm = TRUE)
#> $cupm
#>           [,1]      [,2]      [,3]
#> [1,] 0.5159723 0.1838166 0.2248306
#> [2,] 0.1838166 0.5787444 0.2248629
#> [3,] 0.2248306 0.2248629 0.5614638
#> 
#> $dupm
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.2612807 0.3413677
#> [2,] 0.2775962 0.0000000 0.2259014
#> [3,] 0.2527159 0.2505223 0.0000000
#> 
#> $dlpm
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.2775962 0.2527159
#> [2,] 0.2612807 0.0000000 0.2505223
#> [3,] 0.3413677 0.2259014 0.0000000
#> 
#> $clpm
#>           [,1]      [,2]      [,3]
#> [1,] 0.4840277 0.2773064 0.1810858
#> [2,] 0.2773064 0.4212556 0.2987135
#> [3,] 0.1810858 0.2987135 0.4385362
#> 
#> $cov.matrix
#>             [,1]        [,2]        [,3]
#> [1,]  1.00000000 -0.07775396 -0.18816729
#> [2,] -0.07775396  1.00000000  0.04715278
#> [3,] -0.18816729  0.04715278  1.00000000
#> 

# Use 0's for targets
PM.matrix(1, 1, target = rep(0, ncol(A)), variable = A, pop_adj = TRUE)        
#> $cupm
#>            [,1]       [,2]      [,3]
#> [1,] 0.50064636 0.09886398 0.1716709
#> [2,] 0.09886398 0.46466892 0.1341120
#> [3,] 0.17167087 0.13411205 0.6060483
#> 
#> $dupm
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.1142550 0.2066309
#> [2,] 0.1966875 0.0000000 0.1792208
#> [3,] 0.1423007 0.1102007 0.0000000
#> 
#> $dlpm
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.1966875 0.1423007
#> [2,] 0.1142550 0.0000000 0.1102007
#> [3,] 0.2066309 0.1792208 0.0000000
#> 
#> $clpm
#>            [,1]      [,2]       [,3]
#> [1,] 0.34084226 0.1585364 0.07625754
#> [2,] 0.15853639 0.4820773 0.17031044
#> [3,] 0.07625754 0.1703104 0.31088014
#> 
#> $cov.matrix
#>             [,1]        [,2]        [,3]
#> [1,]  0.84148862 -0.05354215 -0.10100318
#> [2,] -0.05354215  0.94674624  0.01500096
#> [3,] -0.10100318  0.01500096  0.91692847
#> 

# Use variable medians as targets
PM.matrix(1, 1, target = apply(A, 2, "median"), variable = A, pop_adj = TRUE)  
#> $cupm
#>           [,1]      [,2]      [,3]
#> [1,] 0.4514367 0.1309246 0.1529324
#> [2,] 0.1309246 0.6379715 0.1720977
#> [3,] 0.1529324 0.1720977 0.5749040
#> 
#> $dupm
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.1581679 0.2136017
#> [2,] 0.1384553 0.0000000 0.1201787
#> [3,] 0.1381406 0.1492058 0.0000000
#> 
#> $dlpm
#>           [,1]      [,2]      [,3]
#> [1,] 0.0000000 0.1384553 0.1381406
#> [2,] 0.1581679 0.0000000 0.1492058
#> [3,] 0.2136017 0.1201787 0.0000000
#> 
#> $clpm
#>            [,1]      [,2]       [,3]
#> [1,] 0.38262527 0.1254005 0.08925276
#> [2,] 0.12540048 0.3112238 0.13547616
#> [3,] 0.08925276 0.1354762 0.33458706
#> 
#> $cov.matrix
#>             [,1]        [,2]        [,3]
#> [1,]  0.83406192 -0.04029807 -0.10955714
#> [2,] -0.04029807  0.94919535  0.03818937
#> [3,] -0.10955714  0.03818937  0.90949103
#> 
```
