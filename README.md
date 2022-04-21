
<!-- README.md is generated from README.Rmd. Please edit that file -->
# recm

<!-- badges: start -->
<!-- badges: end -->
The R::recm package provides tools for estimating parameters for Radiocarbon-dated Event Count (REC) Models and some relevant plotting functions. REC models are count-based regression models that account for radiocarbon-dating uncertainty in event times. They can be used in cases where abundances of radiocarbon dates are used as a proxy for some process that occurred in the past. Like any regression model, covariates are used to test hypotheses about potential causal relationships between those covariates and through-time variation in radiocarbon-dated sample abundances. Unlike previously developed REC model packages, this package employs an 'end-to-end' Bayesian framework, which means that the likelihoods of individual event times (given radiocarbon measurement and calibration uncertainties) are included directly in the overall likelihood for the model.

For now, the recm package only includes a simple Poisson REC model. It includes no autocorrelation terms and accounts only for radiocarbon-dating uncertainty in the dependent variable. It does not account for chronological or other uncertainty in the included independent variables. In time, the package will be extended to include covariate uncertainties (both chronological and measurement) and more complicated count-based time-series regression models.

## Installation

<!-- not on CRAN yet
You can install the released version of recm from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("recm")
```
-->
You can install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("wccarleton/recmod")
```

## Example

Imagine a set of true event times (in calendar years BP) generated from a simple exponentially increasing process. This could be like a simple exponential population growth model with a positive growth rate of 0.002 over an arbitrary finite timespan (say, 4000--3001 BP) sampled at an annual resolution:

``` r
library(recm)

process <- exp(c(1:1000) * 0.002)

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

We can then run a REC model analysis in an attempt to reconstruct the true growth rate for the exponential process (0.002, from above). To do so, we would need to include the passage of time as the sole covariate, since that's what defined the inputs for the target process. To include the covariate in the REC model, we create an `X` matrix with a column of ones representing the model intercept and a column containing a sequence that represents the passage of time scaled to whatever temporal resolution we want to use for the analysis. The temporal resolution won't affect the results as long as we have enough observations to estimate the model's parameters. The resolution just scales the regression coefficient in the model. In this case, we can speed up our estimation process by using a decadal resolution, meaning that we will have 100 decade-long time-bins spanning 1000 years:

``` r
X <- cbind(rep(1, 100), 1:100 * 10)
```

Next, we need to establish a temporal grid defining the bins that will be used within the recm function to bin the event times. Binning the event times creates probable count sequences to be used in the core Poisson regression model:

``` r
t_edges <- seq(4000, 3000,-10)
t_mids <- (t_edges - 5)[-length(t_edges)]
```

The `recm::recm()` function runs a Metropolis-Hastings MCMC to sample the Poisson regression model's posterior distribution. The simulation takes several arguments. The key ones to note here are

-   *niter*, which determines the number of iterations the MCMC is run for;
-   *adapt*, which is a logical (T or F) that determines whether the engage an adaptive algorithm for finding optimal proposal distribution scales for the main regression parameters (i.e., excluding the scales for proposal distributions pertaining to the event times);
-   *adapt\_amount*, *adapt\_interval*, and *adapt\_window*, which are all related to the adaptive algorithm;
-   *scales*, which is a vector of starting values for the proposal distribution scales for the main regression parameter.

The `recm()` function should in most cases be called the first time with `adapt = T` in order to find good values for proposal distribution scales. It will return a list including `$scales`, which will be a matrix of scale values that have been adapted throughout the simulation. Averages of these values (one for each of the adapted parameter proposal distributions) can be passed to the recm function again as the `scales` parameter in order to run the MCMC for longer and get good samples of the model's posteriors.

Unlike some other regression functions, the dependent variable is passed into `recm()` as a matrix of radiocarbon dates (simulated from above) and the covariates are passed in as the `X` matrix. The third parameter is a vector of temporal bin edges that define a temporal grid for the analysis.

``` r
mysamples_adapt <- recm(dates,
                    X,
                    t_edges,
                    niter = 10000,
                    adapt = T,
                    adapt_amount = 0.1,
                    adapt_interval = 50,
                    adapt_window = c(0.21, 0.25),
                    scales = c(0.1, 0.01))
#> No starting values provided ('startvals == NULL'). Using defaults.
#> No calibration curve provided. Using IntCal.
#> Determining calibrated date ranges...
#> Running MCMC...

burnin = 1:floor(dim(mysamples_adapt$scales)[1] * 0.25)

newscales <- colMeans(mysamples_adapt$scales[-burnin,])

mysamples <- recm(dates,
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

The second run of the recm function above returned only the MCMC chains as a matrix with the first columns containing samples for the posteriors of the regression coefficients (two in this case corresponding to the two columns of the `X` matrix). The remaining columns contain samples of the event times, which will be samples from the calibrated date distributions---these are ordered the same as the radiocarbon dates in the matrix passed as the `dates` argument. These MCMC samples can then be plotted in order to check for mixing/convergence and then used to estimate posterior distributions and perform inference.

### MCMC Chain for target regression coefficient

<img src="man/figures/README-plots-1.png" width="100%" />

### Density estimate for posterior of target regression coefficient

<img src="man/figures/README-plots-2-1.png" width="100%" />

## Note

In a practical, real-world analysis you probably won't be trying to find parameters for a model defined by a simple exponential function of time. Rather, you are likely to be interested in testing one or more hypotheses about the nature of through-time changes in the abundance of radiocarbon-dated samples (events, like settlement occupations, or burials). As a result, the `X` parameter passed to `recm()` should be a matrix containing relevant covariates---e.g., a temperature reconstruction or proxy for past precipitation levels. The first column of `X` should be a column of 1's in most cases (representing the inclusion of an intercept in the regression) followed by one column for each relevant covariate/predictor variable.

Those predictor variables will need to be put onto the same temporal grid used to define the whole analysis (defined in the `recm()` function by the `t_edges` parameter). Be sure to prepare the `X` matrix accordingly. If, for instance, you wanted to test whether past temperature levels affected population sizes (using radiocarbon-dated events as a proxy for the latter), then you will have to align the temperature reconstruction onto the temporal grid defined by `t_edges`. This could be as simple as taking a running average of the temperature reconstruction and sampling it at the midpoints of the grid bins to produce a sequence of temperature averages at the appropriate temporal resolution (e.g., decadal averages, or whatever is appropriate).

Of course, almost any potential covariate will have chronological and measurement uncertainties and the `recm()` function currently does not include those. Future updates will include options for representing covariate uncertainties ('errors-in-variables').
