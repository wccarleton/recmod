#' Bayesian Radiocarbon-dated Event Count Model.
#'
#' @description This function runs a Metropolis-Hastings MCMC simulation to
#'  estimate the parameters of a count-based (Poisson) regression model
#'  accounting for radiocarbon-dating chronological uncertainty in event times.
#'
#' @param dates Matrix containing radiocarbon date means and errors in the
#'  first and second columns, respectively.
#' @param X A matrix containing the covariates for the model (independent
#'  variables) including an intercept if desired. Each column should contain
#'  one variable with the leading column being the intercept (a column of 1's).
#' @param t_edges A vector of temporal bin edges that determine the temporal
#'  scale and span of the event count series being modelled. These edges create
#'  intervals that are assumed to be left open and right closed when events are
#'  counted in them.
#' @param model A string specifying the type of model to use. Currently only
#'  Poisson ('pois') models are supported.
#' @param startvals A numeric vector of starting values for the MCMC. The
#'  leading n elements will be the starting values for the regression
#'  parameters (coefficients associated with the columns of X) while the
#'  remaining columns should contain initial values for the dates (in years BP
#'  1950).
#' @param niter Integer indicating the number of MCMC iterations.
#' @param adapt Logical (default = True), indicating whether to use an adaptive
#'  MCMC algorithm for finding optimal variances for proposal distributions.
#' @param adapt_amount A scalar containing the ratio by which the proposal
#'  variance is changed during the MCMC adapt step when adapt = True---e.g.,
#'  adapt_amount = 0.1 (default) grows/shrinks the variance by 10%.
#' @param adapt_interval (default = 100) A numeric scalar indicating the the n-
#'  th mcmc iteration during which the simulation will attempt to adapt the
#'  proposal variances.
#' @param adapt_window A numeric vector with two elements that define the lower
#'  and upper bounds of the target acceptance value for the MCMC. This is only
#'  used when adapt = T. Default is [0.21, 0.25].
#' @param scales A numeric vector of scales for the proposal distributions. The
#'  order is important. Each element refers to the scale of a proposal function
#'  for the model parameters in the following order: (b_1, b_2, ..., b_nX),
#'  where `b` refers to a regression coefficient and `nX' is the number of
#'  columns in X. If Null, naive estimates are used.
#' @param priors A numeric vector of parameter values for the model priors---the
#'  order is important and corresponds to order of appearance of the prior
#'  density functions in the source code for the `prior` function. If Null,
#'  very wide priors are used--e.g., N(0, 1000).
#' @param BP Logical (default = True), indicating whether to assume dates
#'  provided are in BP. If so, dates (ages) increase into the past. It is up to
#'  the user to ensure that the X columns are oriented appropriately.
#' @param calcurve Matrix containing a calibration curve. First column should
#'  be calibrated ages/dates; second column should contain radiocarbon years or
#'  fraction modern values; third column should contain errors (standard
#'  deviation). If Null, IntCal::ccurve() is called, meaning that the latest
#'  IntCal curve available from that package is used.
#'
#' @return If adapt = True, returns a named list including: the mcmc
#'  $samples; $acceptance rates for all model parameters; and $scales found
#'  with the adaptive algorithm. Otherwise, returns only the mcmc samples (mcmc
#'  chains) as a matrix where each column contains samples for one of the
#'  model's parameters. The mcmc samples matrix will be sorted such that the
#'  first 'n' columns contain regression coefficients (in the same order as the
#'  columns in the X argument) where 'n' is the number of covariates included;
#'  the next 'm' columns will contain any other relevant model parameters as
#'  specified by the model itself (for now only a simple Poisson regression is
#'  supported and so 'm' = 0); and the remaining 'l' columns will contain
#'  samples of 'l' event ages, where 'l' is the number of dated events included
#'  in the 'dates' argument.
#'
#' @import progress stats utils IntCal
#' @export

recm <- function(dates,
                X,
                t_edges,
                model = "pois",
                startvals = NULL,
                niter = 10000,
                adapt = T,
                adapt_amount = 0.1,
                adapt_interval = 100,
                adapt_window = c(0.21, 0.25),
                scales = NULL,
                priors = NULL,
                BP = T,
                calcurve = NULL){

    # set up variables and initialize defaults if needed

    N <- dim(dates)[1] # number of events (dates)
    nX <- dim(X)[2] # number of regression coefficients
    nparams <- nX # number of model parameters excluding event times
    nedges <- length(t_edges) # number of temporal grid edges

    if (is.null(startvals)){
        alert <- paste("No starting values provided ('startvals == NULL').",
                        "Using defaults.",
                        sep = " ")
        message(alert)
        startvals <- c(rep(0, nX), dates[, 1])
    }

    if (is.null(scales)){
        alert <- paste("No scales provided ('scales == NULL').",
                        "Using defaults.",
                        sep = " ")
        message(alert)
        scales <- rep(0.001, nX)
    }

    if (is.null(calcurve)){
        alert <- paste("No calibration curve provided. Using IntCal.",
                    sep = "")
        message(alert)
        calcurve = IntCal::ccurve()
    }

    # get calibrated date approximations for use with proposals
    message("Determining calibrated date ranges...")

    if (BP){
        curve_limits <- rev(range(calcurve[, 1]))
    }else{
        curve_limits <- range(calcurve[, 1])
    }
    cal_dates <- curve_limits[1]:curve_limits[2]

    cal_matrix <- array(dim = c(length(cal_dates), N))

    pb <- progress::progress_bar$new(total = N)
    for(j in 1:N){
        pb$tick()
        cal_matrix[,j] <- exp(cal_likelihood(dates[j, 1],
                                        dates[j, 2],
                                        cal_dates,
                                        calcurve))
    }
    t_range <- t(apply(cal_matrix,2,function(x){
                                        range(cal_dates[-which(x == 0)])
                                        }))

    # setup the mcmc chain container

    chain <- array(dim = c(niter + 1, nparams + N))

    chain[1, ] <- startvals

    # initialize containers for adapt procedure if needed

    if (adapt){
        n_adapts <- floor(niter / adapt_interval)
        acceptance <- array(dim = c(n_adapts, nparams))
        scales_matrix <- array(dim = c(n_adapts, nparams))
    }

    # initial posterior density

    proposal <- chain[1, ]

    t_sample <- which(proposal[-c(1:nparams)] <= t_edges[1] &
                    proposal[-c(1:nparams)] > t_edges[nedges])

    Y <- count_events(x = proposal[-c(1:nparams)][t_sample],
                    breaks = t_edges,
                    BP = BP)

    pd_previous <- posterior(dates,
                        Y,
                        X,
                        proposal,
                        nX,
                        priors,
                        calcurve)

    message("Running MCMC...")

    pb <- progress::progress_bar$new(total = niter)

    # main loop

    for (j in 1:niter){
        pb$tick()

        # adapt step

        if(adapt & (j %% adapt_interval == 0)){
            interval_index <- c(j - adapt_interval + 1):j
            diffs <- apply(as.matrix(chain[interval_index, 1:nparams]), 2, diff)
            acceptance_rate <- 1 - colMeans(diffs == 0)
            acceptance[j / adapt_interval, ] <- acceptance_rate
            scales <- adaptScale(acceptance_rate,
                                scales,
                                adapt_amount,
                                adapt_window)
            scales_matrix[j / adapt_interval, ] <- scales
        }

        chain[j + 1, ] <- chain[j, ]

        # each regression param separately
        for (l in 1:nX){

            # propose step

            proposal <- chain[j + 1, ]
            proposal[l] <- propose_reg(chain[j + 1, l], scales[l])

            pd_proposal <- posterior(dates,
                                Y,
                                X,
                                proposal,
                                nX,
                                priors,
                                calcurve)

            # accept step

            accept <- exp(pd_proposal - pd_previous)

            if (runif(1) < accept){
                chain[j + 1, ] <- proposal
                pd_previous <- pd_proposal
            }
        }

        # each date separately

        for (l in 1:N){

            # propose step

            proposal_d <- chain[j + 1, ]
            proposal_d[nparams + l] <- propose_t(t_range[l, ])

            t_sample <- which(proposal_d[-c(1:nparams)] <= t_edges[1] &
                            proposal_d[-c(1:nparams)] > t_edges[nedges])

            Y <- count_events(x = proposal_d[-c(1:nparams)][t_sample],
                            breaks = t_edges,
                            BP = BP)

            pd_proposal <- posterior(dates,
                                    Y,
                                    X,
                                    proposal_d,
                                    nX,
                                    priors,
                                    calcurve)

            # accept step

            accept <- exp(pd_proposal - pd_previous)
            if(runif(1) < accept){
                chain[j + 1, ] <- proposal_d
                pd_previous <- pd_proposal
            }
        }
    }
    if (adapt){
        return(list(
                samples = chain,
                acceptance = acceptance,
                scales = scales_matrix))
    }else{
        return(chain)
    }
}

#' Calibration likelihood
#'
#' @param c14_mean Vector containing one or more uncalibrated radiocarbon-date
#'  means or Fraction Modern values (should correspond to the calibration curve
#'  provided to 'calcurve').
#' @param c14_err Vector containing one or more uncalibrated radiocarbon-date
#'  errors---must be the same length as c14_mean, but this is not checked
#'  internally
#' @param caldate Vector of dates on the same calendar as the provided
#'  calibration curve---typically this will be a calendar year BP (where
#'  present is 1950)
#' @param calcurve Matrix containing a calibration curve. First column should
#'  be calibrated ages/dates; second column should contain radiocarbon years or
#'  fraction modern values; third column should contain errors (standard
#'  deviation).
#'
#' @return Vector of log-likelihoods for the provided caldate(s):
#'  p(caldate|c14_mean, c14_err, calcurve)

cal_likelihood <- function(c14_mean,
                    c14_err,
                    caldate,
                    calcurve){
    curve_c14_mean <- approx(x = calcurve[, 1],
                            y = calcurve[, 2],
                            xout = caldate)$y
    curve_err <- approx(x = calcurve[, 1],
                        y = calcurve[, 3],
                        xout = caldate)$y
    ll <- dnorm(x = curve_c14_mean,
                mean = c14_mean,
                sd = sqrt(c14_err^2 + curve_err^2),
                log = T)
    return(ll)
}

#' Count sequence likelihood
#'
#' @param t_sample Vector containing a sample of potential ages, one for each
#'  event in the dataset.
#' @param dates Matrix containing radiocarbon date means and errors in the
#'  first and second columns, respectively.
#' @param calcurve Matrix containing a calibration curve. First column should
#'  be calibrated ages/dates; second column should contain radiocarbon years or
#'  fraction modern values; third column should contain errors (standard
#'  deviation).
#'
#' @return Scalar summed log-likelihood of t_sample given the uncertainties
#'  defined by the dates and calibration curve. This log-likelihood indicates
#'  the likelihood for the sample of ages as a group (i.e., the p(t_sample|
#'  c14_mean, c14_err, calcurve)) and, so, corresponds to the likelihood for
#'  the corresponding count sequence.

count_likelihood <- function(t_sample,
                            dates,
                            calcurve){
    ll_vector <- cal_likelihood(dates[, 1],
                            dates[, 2],
                            t_sample,
                            calcurve)
    sll <- sum(ll_vector)
    return(sll)
}

#' Regression model likelihood
#'
#' @param Y Vector containing a probable count sequence.
#' @param X Vector containing covariates, including an intercept as the first
#'  element if desired.
#' @param params Vector containing regression coefficients.
#'
#' @return Scalar log-likelihood of the proposed regression parameters, given Y
#'  and X.

reg_likelihood <- function(Y,
                        X,
                        params){
    pred <- exp(as.matrix(X) %*% params)
    loglike <- dpois(x = Y,
                    lambda = pred,
                    log = T)
    sll <- sum(loglike)
    return(sll)
}

#' Prior probability for regression parameters
#'
#' @param params Vector of proposed regression parameter values.
#' @param priors Vector (Matrix) containing two elements (columns), one for the
#'  prior mean(s) and the other for the standard deviation(s) of Gaussian
#'  distribution(s)
#'
#' @return Scalar summed log-likelihood for the proposed regression parameter
#'  values given the priors.

prior <- function(params,
                priors = NULL){
    if(is.null(priors)){
        priors <- c(0, 1000)
    }
    B_priors <- dnorm(x = params,
                    mean = priors[1],
                    sd = priors[2],
                    log = T)
    sll <- sum(B_priors)
    return(sll)
}

#' Posterior for REC model regression parameters
#'
#' @param c14_mean Vector containing one or more uncalibrated radiocarbon-date
#'  means or Fraction Modern values (should correspond to the calibration curve
#'  provided to 'calcurve').
#' @param c14_err Vector containing one or more uncalibrated radiocarbon-date
#'  errors---must be the same length as c14_mean, but this is not checked
#'  internally
#' @param caldate Vector of dates on the same calendar as the provided
#'  calibration curve---typically this will be a calendar year BP (where
#'  present is 1950)
#' @param calcurve Matrix containing a calibration curve. First column should
#'  be calibrated ages/dates; second column should contain radiocarbon years or
#'  fraction modern values; third column should contain errors (standard
#'  deviation).
#'
#' @return Vector of log-likelihoods for the provided caldate(s):
#'  p(caldate|c14_mean, c14_err, calcurve)

posterior <- function(dates,
                    Y,
                    X,
                    params,
                    nX,
                    priors,
                    calcurve){
    post <- count_likelihood(params[-c(1:nX)], dates, calcurve) +
            reg_likelihood(Y, X, params[1:nX]) +
            prior(params[1:nX], priors)
    return(post)
}

#' Proposal function for regression parameters
#'
#' @param Bs Vector containing the current regression coefficient value(s) used
#'  in the MCMC.
#' @param v scale for Gaussian proposal distribution.
#'
#' @return Vector proposed regression coefficient parameter values.

propose_reg <- function(Bs,
                        v){
    nX <- length(Bs)
    Bs <- rnorm(nX, mean = Bs, sd = v)
    return(Bs)
}

#' Proposal function for event times
#'
#' @param t_range Vector containing two elements referring to the start
#'  (t_range[1]) and end (t_range[2]) of a span of time from which to sample a
#'  probable event age.
#'
#' @return Scalar containing a possible event time (age) from a uniform
#'  distribution with limits defined by the t_range parameter.

propose_t <- function(t_range){
    t_sample <- runif(1, t_range[1], t_range[2])
    return(t_sample)
}

#' Adapt proposal distribution scales
#'
#' @param acceptance_rate Vector of current acceptance rate for MCMC parameters.
#' @param s Vector of current proposal distribution scales.
#' @param adapt_amount Scalar indicating the factor to grow/shrink the current
#'  scale(s) by---e.g., 0.1 (10% change) would cause a given scale to grow to
#'  1.1x or shrink to 0.9x its previous value.
#' @param adapt_window Vector containing two elements that together define a
#'  target acceptance rate window for the scale adaptation algorithm. If the
#'  current acceptance rate for a given parameter falls outside the window
#'  (i.e., is greater/less than the first/second element of the adapt_window
#'  vector), the relevant proposal scale is decreased/increased by
#'  adapt_amount.
#'
#' @return Vector of adapted scales.

adaptScale <- function(acceptance_rate,
                    s,
                    adapt_amount,
                    adapt_window){
    s <- s - (s * adapt_amount * (acceptance_rate < adapt_window[1]))
    s <- s + (s * adapt_amount * (acceptance_rate > adapt_window[2]))
    return(s)
}
