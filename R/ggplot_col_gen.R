#' Generate ggplot colors
#' 
#' Accepts a number of colors to generate and generates a ggplot color spectrum.
#' 
#' @param n Number of colors to generate
#' @return A vector of colors according to ggplot color generation.
#' 
ggplot_col_gen = function(n){
    hues = seq(15, 375, length = n + 1)
    return(hcl(h = hues, l = 65, c = 100)[1:n])
}