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
