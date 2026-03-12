# NNS Nowcast

Wrapper function for NNS nowcasting method using the nonparametric
vector autoregression
[NNS.VAR](https://OVVO-Financial.github.io/NNS/reference/NNS.VAR.md),
and Federal Reserve Nowcasting variables.

## Usage

``` r
NNS.nowcast(
  h = 1,
  additional.regressors = NULL,
  additional.sources = NULL,
  naive.weights = FALSE,
  specific.regressors = NULL,
  start.date = "2000-01-03",
  keep.data = FALSE,
  status = TRUE,
  ncores = NULL
)
```

## Arguments

- h:

  integer; `(h = 1)` (default) Number of periods to forecast. `(h = 0)`
  will return just the interpolated and extrapolated values up to the
  current month.

- additional.regressors:

  character; `NULL` (default) add more regressors to the base model. The
  format must utilize the
  [`getSymbols`](https://rdrr.io/pkg/quantmod/man/getSymbols.html)
  format for FRED data, else specify the source.

- additional.sources:

  character; `NULL` (default) specify the `source` argument per
  [`getSymbols`](https://rdrr.io/pkg/quantmod/man/getSymbols.html) for
  each `additional.regressors` specified.

- naive.weights:

  logical; `TRUE` Equal weights applied to univariate and multivariate
  outputs in ensemble. `FALSE` (default) will apply weights based on the
  number of relevant variables detected.

- specific.regressors:

  integer; `NULL` (default) Select individual regressors from the base
  model per Viole (2020) listed in the `Note` below.

- start.date:

  character; `"2000-01-03"` (default) Starting date for all data series
  download.

- keep.data:

  logical; `FALSE` (default) Keeps downloaded variables in a new
  environment `NNSdata`.

- status:

  logical; `TRUE` (default) Prints status update message in console.

- ncores:

  integer; value specifying the number of cores to be used in the
  parallelized subroutine
  [NNS.ARMA.optim](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.optim.md).
  If NULL (default), the number of cores to be used is equal to the
  number of cores of the machine - 1.

## Value

Returns the following matrices of forecasted variables:

- `"interpolated_and_extrapolated"` Returns a `data.frame` of the linear
  interpolated and
  [NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md)
  extrapolated values to replace `NA` values in the original `variables`
  argument. This is required for working with variables containing
  different frequencies, e.g. where `NA` would be reported for
  intra-quarterly data when indexed with monthly periods.

- `"relevant_variables"` Returns the relevant variables from the
  dimension reduction step.

- `"univariate"` Returns the univariate
  [NNS.ARMA](https://OVVO-Financial.github.io/NNS/reference/NNS.ARMA.md)
  forecasts.

- `"multivariate"` Returns the multi-variate
  [NNS.reg](https://OVVO-Financial.github.io/NNS/reference/NNS.reg.md)
  forecasts.

- `"ensemble"` Returns the ensemble of both `"univariate"` and
  `"multivariate"` forecasts.

## Note

Specific regressors include:

1.  `PAYEMS` – Payroll Employment

2.  `JTSJOL` – Job Openings

3.  `CPIAUCSL` – Consumer Price Index

4.  `DGORDER` – Durable Goods Orders

5.  `RSAFS` – Retail Sales

6.  `UNRATE` – Unemployment Rate

7.  `HOUST` – Housing Starts

8.  `INDPRO` – Industrial Production

9.  `DSPIC96` – Personal Income

10. `BOPTEXP` – Exports

11. `BOPTIMP` – Imports

12. `TTLCONS` – Construction Spending

13. `IR` – Import Price Index

14. `CPILFESL` – Core Consumer Price Index

15. `PCEPILFE` – Core PCE Price Index

16. `PCEPI` – PCE Price Index

17. `PERMIT` – Building Permits

18. `TCU` – Capacity Utilization Rate

19. `BUSINV` – Business Inventories

20. `ULCNFB` – Unit Labor Cost

21. `IQ` – Export Price Index

22. `GACDISA066MSFRBNY` – Empire State Mfg Index

23. `GACDFSA066MSFRBPHI` – Philadelphia Fed Mfg Index

24. `PCEC96` – Real Consumption Spending

25. `GDPC1` – Real Gross Domestic Product

26. `ICSA` – Weekly Unemployment Claims

27. `DGS10` – 10-year Treasury rates

28. `T10Y2Y` – 2-10 year Treasury rate spread

29. `WALCL` – Total Assets

30. `PALLFNFINDEXM` – Global Price Index of All Commodities

31. `FEDFUNDS` – Federal Funds Effective Rate

32. `PPIACO` – Producer Price Index All Commodities

33. `CIVPART` – Labor Force Participation Rate

34. `M2NS` – M2 Money Supply

35. `ADPMNUSNERNSA` – ADP Payrolls

## References

Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics:
Using Partial Moments" (ISBN: 1490523995)

Viole, F. (2019) "Multi-variate Time-Series Forecasting: Nonparametric
Vector Autoregression Using NNS"
[doi:10.2139/ssrn.3489550](https://doi.org/10.2139/ssrn.3489550)

Viole, F. (2020) "NOWCASTING with NNS"
[doi:10.2139/ssrn.3589816](https://doi.org/10.2139/ssrn.3589816)

## Author

Fred Viole, OVVO Financial Systems

## Examples

``` r
 if (FALSE) { # \dontrun{
 ## Interpolates / Extrapolates all variables to current month
 NNS.nowcast(h = 0)
 
 ## Additional regressors and sources specified
 NNS.nowcast(h = 0, additional.regressors = c("SPY", "USO"), 
             additional.sources = c("yahoo", "yahoo"))
             
              
 ### PREDICTION INTERVALS 
 ## Store NNS.nowcast output
 nns_estimates <- NNS.nowcast(h = 12)           
 
 # Create bootstrap replicates using NNS.meboot (GDP Variable)
 gdp_replicates <- NNS.meboot(nns_estimates$ensemble$GDPC1, 
                              rho = seq(0,1,.25), 
                              reps = 100)["replicates",]
                              
 replicates <- do.call(cbind, gdp_replicates)
 
 # Apply UPM.VaR and LPM.VaR for desired prediction interval...95 percent illustrated
 # Tail percentage used in first argument per {LPM.VaR} and {UPM.VaR} functions
 lower_GDP_CIs <- apply(replicates, 1, function(z) LPM.VaR(0.025, 0, z))
 upper_GDP_CIs <- apply(replicates, 1, function(z) UPM.VaR(0.025, 0, z))
 
 # View results
 cbind(nns_estimates$ensemble$GDPC1, lower_GDP_CIs, upper_GDP_CIs)
 } # }
```
