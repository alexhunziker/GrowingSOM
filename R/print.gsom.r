print.gsom <- function(x, ...){
  
  if(!is.null(x$mapped)){
    observations = length(x$mapped$bmn)
    return(cat("GrowingSOM object.", observations, "observations are mapped onto an existing GSOM map.\n"))
  }
  
  if(!is.null(x$prediction)){
    observations = length(x$prediction$bmn)
    return(cat("GrowingSOM object. For", observations, "observations Y values were predicted according to an existing gsom object.\n"))
  }
  
  if(x$nhood == "rect") topology = "rectangular"
  else topology = "hexagonal"
  
  cat("GrowingSOM (gsom) map with", nrow(x$nodes$position), "nodes and a", topology, "topology.\n")
  
  if(!is.null(x$nodes$predict)) cat("Model contains predictions for dependent variables.\n")
  
  if(!exists("x$data")) included = "is not included.\n"
  else included = "is included.\n"
  
  cat("Training data", included)
}