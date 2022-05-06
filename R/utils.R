#' Count events
#'
#' @description Bins event times given in 'x' into time grid defined by
#'  'breaks' using graphics::hist.
#'
#' @param x Vector of event times.
#' @param breaks Vector of temporal bin edges.
#' @param BP Logical (default = True), indicating whether to assume dates
#'  provided are in BP. If so, dates (ages) increase into the past. The
#'  returned count sequence will be re-sorted appropriately so that age
#'  increases down the length of the vector.
#'
#' @return Vector containing (temporally binned) count sequence

count_events <- function(x,
                        breaks,
                        BP = TRUE){
    if(BP){
        counts <- graphics::hist(
                        x,
                        breaks = breaks,
                        include.lowest = FALSE,
                        plot = FALSE)$counts
        return(rev(counts))
    }else{
        counts <- graphics::hist(
                        x,
                        breaks = breaks,
                        right = FALSE,
                        include.lowest = FALSE,
                        plot = F)$counts
        return(counts)
    }
}

#' Get number of radiocarbon dated events included as one or more covariates
#'
#' @description Counts the number of radiocarbon-dated events included as
#'  covariates in the REC model. This is an internal utility function.
#'
#' @param x List where each element is a matrix or dataframe containing
#'  covariate data.
#' @param xdistr Character Vector where each element is the name of a data type
#'  supported by the recm package. It is used in this function to isolate the
#'  elements of x containing radiocarbon data.
#'
#' @return Vector each element is the number of radiocarbon dated events
#'  included as each covariate in the REC model.

get_n_xevents <- function(x, xdistr){

}

name_events <- function(varname, n){
    return(paste(varname, 1:n, sep = ""))
}

## get t_ranges

get_cal_range <- function(dates, calcurve, BP){
    ndates <- dim(dates)[1]
    if (BP){
        curve_limits <- rev(range(calcurve[, 1]))
    }else{
        curve_limits <- range(calcurve[, 1])
    }
    cal_dates <- curve_limits[1]:curve_limits[2]
    cal_matrix <- array(dim = c(length(cal_dates), ndates))

    pb <- progress::progress_bar$new(total = ndates)
    for(j in 1:nyevents){
        pb$tick()
        cal_matrix[,j] <- exp(cal_likelihood(dates[j, 1],
                                        dates[j, 2],
                                        cal_dates,
                                        calcurve))
    }
    t_range <- t(apply(cal_matrix,2,function(x){
                                        range(cal_dates[-which(x == 0)])
                                        }))
    return(t_range)
}
