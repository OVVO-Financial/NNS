# NNS Term Matrix

Generates a term matrix for text classification use in
[NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md).

## Usage

``` r
NNS.term.matrix(x, oos = NULL)
```

## Arguments

- x:

  mixed data.frame; character/numeric; A two column dataset should be
  used. Concatenate text from original sources to comply with format.
  Also note the possibility of factors in `"DV"`, so
  `"as.numeric(as.character(...))"` is used to avoid issues.

- oos:

  mixed data.frame; character/numeric; Out-of-sample text dataset to be
  classified.

## Value

Returns the text as independent variables `"IV"` and the classification
as the dependent variable `"DV"`. Out-of-sample independent variables
are returned with `"OOS"`.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- data.frame(cbind(c("sunny", "rainy"), c(1, -1)))
NNS.term.matrix(x)

### Concatenate Text with space separator, cbind with "DV"
x <- data.frame(cbind(c("sunny", "rainy"), c("windy", "cloudy"), c(1, -1)))
x <- data.frame(cbind(paste(x[ , 1], x[ , 2], sep = " "), as.numeric(as.character(x[ , 3]))))
NNS.term.matrix(x)

### NYT Example
require(RTextTools)
data(NYTimes)

### Concatenate Columns 3 and 4 containing text, with column 5 as DV
NYT <- data.frame(cbind(paste(NYTimes[ , 3], NYTimes[ , 4], sep = " "),
                     as.numeric(as.character(NYTimes[ , 5]))))
NNS.term.matrix(NYT)
} # }
```
