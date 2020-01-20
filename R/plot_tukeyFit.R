#' Plot results of Tukey's sensitivity analysis for either the ATE or QTEs.
#'
#' @param x 
#' @param type either "ate" (for a heatmap of average treatment effects) or "qte" for ribbon plots showing quantile treatment effects 
#' @param ... additional plotting parameters for either \code{ate_plot} or \code{qte_plot}
#'
#'
#' @examples
#' 
#' @export
#' 
plot.tukeyFit <- function(tukeyFit, type = "ate", ...) {
  
  stopifnot(class(tukeyFit)=="tukeyFit")
  
  if(type == "ate") {
       
     ate_plot(tukeyFit, ...)
     
  } else if(type == "qte") {
     qte_plot(tukeyFit, ...)
  } else {
    stop("Type must be 'ate' or 'qte'.")    
  }

}
