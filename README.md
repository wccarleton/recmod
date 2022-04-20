
<!-- README.md is generated from README.Rmd. Please edit that file -->
# recm

<!-- badges: start -->
<!-- badges: end -->
The recm package provides tools for estimating parameters for Radiocarbon-dated Event Count (REC) Models and some relevant plotting functions. REC models are count-based regression models that account for radiocarbon-dating uncertainty in event times. They can be used in cases where abundances of radiocarbon dates are used as a proxy for some process that occurred in the past. Like any regression model, covariates are used to test hypotheses about potential causal relationships between those covariates and through-time variation in radiocarbon-dated sample abundances. Unlike previously developed REC model packages, this package employs an 'end-to-end' Bayesian framework, which means that the likelihoods of individual event times (given radiocarbon measurement and calibration uncertainties) are included directly in the overall likelihood for the model.

For now, the recm package only includes a simple Poisson REC model. It includes no autocorrelation terms and accounts only for radiocarbon-dating uncertainty in the dependent variable. It does not account for chronological or other uncertainty in the included independent variables. In time, the package will be extended to include covariate uncertainties (both chronological and measurement) and more complicated count-based time-series regression models.

## Installation

You can install the released version of recm from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("recm")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("wccarleton/recmod")
```

## Example

Imagine a set of true event times (in calendar years BP) generated from a simple exponentially increasing process. This could be like a simple exponential population growth model with a positive growth rate of 0.002 over an arbitrary finite timespan (say, 4000--3001 BP) sampled at an annual resolution:

``` r
library(recm)

process <- exp(c(1:1000)*0.002)

nevents <- 100

event_times <- sample(x = 4000:3001,
                    size = nevents,
                    replace = T,
                    prob = process)
```

This sample of true calendar dates (years) can then be used to produce a set of plausible radiocarbon determinations using the IntCal package:

``` r
library(IntCal)

dates <- do.call(rbind,lapply(event_times,IntCal::calBP.14C))
```

We can then run a REC model analysis in an attempt to reconstruct the true growth rate for the exponential process (0.002, from above). To do so, we would need to include the passage of time as the sole covariate, since that's what defined the inputs for the target process. To include the covariate in the REC model, we create an X matrix with a column of ones representing the model intercept and a column containing a sequence that represents the passage of time scaled to whatever temporal resolution we want to use for the analysis. The temporal resolution won't affect the results as long as we have enough observations to estimate the model's parameters. The resolution just scales the regression coefficient in the model. In this case, we can speed up our estimation process by using a decadal resolution:

``` r
X <- cbind(rep(1, 1000), 1:1000 * 10)
```

Next, we need to establish a temporal grid that will be used within the recm function to bin the event times, creating probable count sequences to be used in the core Poisson regression model:

``` r
t_edges <- seq(4000, 3000,-10)
t_mids <- (t_edges - 5)[-length(t_edges)]
```

The recm function runs a Metropolis-Hastings MCMC to sample the Poisson regression model's posterior distribution. The simulation takes several arguments. The key ones to note here are

-   *niter*, which determines the number of iterations the MCMC is run for;
-   *adapt*, which is a logical (T or F) that determines whether the engage an adaptive algorithm for finding optimal proposal distribution scales for the main regression parameters (i.e., excluding the scales for proposal distributions pertaining to the event times);
-   *adapt\_amount*, *adapt\_interval*, and *adapt\_window*, which are all related to the adaptive algorithm;
-   *scales*, which is a vector of starting values for the proposal distribution scales for the main regression parameter.

The recm function should in most cases be called the first time with `adapt = T` in order to find good values for proposal distribution scales. It will return a list including `$scales`, which will be a matrix of scale values that have been adapted throughout the simulation. Averages of these values (one for each of the adapted parameter proposal distributions) can be passed to the recm function again as the *scales* parameter in order to run the MCMC for longer and get good samples for the posteriors of the regression parameters.

Unlike some other regression functions, the dependent variable is passed into recm as a matrix of radiocarbon dates (simulated from above) and the covariates are passed in as the `X` matrix. The third parameter is a vector of temporal bin edges that define a temporal grid for the analysis.

``` r
mysamples <- recm(dates,
                    X,
                    t_edges,
                    niter = 2000,
                    adapt = T,
                    adapt_amount = 0.1,
                    adapt_interval = 100,
                    adapt_window = c(0.21, 0.25),
                    scales = c(0.005,0.0005))
#> No starting values provided ('startvals == NULL'). Using defaults.
#> No calibration curve provided. Using IntCal.
#> Determining calibrated date ranges...
#> Running MCMC...

newscales = colMeans(mysamples$scales[-c(1:10),])

mysamples_adapted <- recm(dates,
                    X,
                    t_edges,
                    niter = 20000,
                    adapt = F,
                    scales = newscales)
#> No starting values provided ('startvals == NULL'). Using defaults.
#> No calibration curve provided. Using IntCal.
#> Determining calibrated date ranges...
#> Running MCMC...
```

The second run of the recm function above returned only the MCMC chains as a matrix with the first columns containing samples for the posteriors of the regression coefficients (two in this case corresponding to the two columns of the `X` matrix). The remaining columns contain samples of the event times, which will be samples from the calibrated date distributions---these are ordered the same as the radiocarbon dates in the matrix passed as the `dates` argument. These MCMC samples can then be plotted in order to check for mixing/convergence and then used to estimate posterior distributions and perform inference:

<img src="man/figures/README-plots-1.png" width="100%" /><img src="man/figures/README-plots-2.png" width="100%" />