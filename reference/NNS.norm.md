# NNS Normalization

Normalizes a matrix of variables based on nonlinear scaling
normalization method.

## Usage

``` r
NNS.norm(X, linear = FALSE, chart.type = NULL, location = "topleft")
```

## Arguments

- X:

  a numeric matrix or data frame, or a list.

- linear:

  logical; `FALSE` (default) Performs a linear scaling normalization,
  resulting in equal means for all variables.

- chart.type:

  options: ("l", "b"); `NULL` (default). Set `(chart.type = "l")` for
  line, `(chart.type = "b")` for boxplot.

- location:

  Sets the legend location within the plot, per the `x` and `y`
  co-ordinates used in base graphics
  [legend](https://rdrr.io/r/graphics/legend.html).

## Value

Returns a [data.frame](https://rdrr.io/r/base/data.frame.html) of
normalized values.

## Note

Unequal vectors provided in a list will only generate `linear=TRUE`
normalized values.

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
x <- rnorm(100) ; y <- rnorm(100)
A <- cbind(x, y)
NNS.norm(A)

### Normalize list of unequal vector lengths

vec1 <- c(1, 2, 3, 4, 5, 6, 7)
vec2 <- c(10, 20, 30, 40, 50, 60)
vec3 <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3)

vec_list <- list(vec1, vec2, vec3)
NNS.norm(vec_list)
} # }
```
