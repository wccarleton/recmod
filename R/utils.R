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
