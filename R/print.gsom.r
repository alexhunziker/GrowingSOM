print.gsom <- function(gsom_object){
  
  if(!is.null(gsom_object$mapped)){
    observations = length(gsom_object$mapped$bmn)
    return(cat("GrowingSOM object.", observations, "observations are mapped onto an existing GSOM map.\n"))
  }
  
  if(!is.null(gsom_object$prediction)){
    observations = length(gsom_object$prediction$bmn)
    return(cat("GrowingSOM object. For", observations, "observations Y values were predicted according to an existing gsom object.\n"))
  }
  
  if(gsom_object$nhood == "rect") topology = "rectangular"
  else topology = "hexagonal"
  
  cat("GrowingSOM (gsom) map with", nrow(gsom_object$nodes$position), "nodes and a", topology, "topology.\n")
  
  if(!is.null(gsom_object$nodes$predict)) cat("Model contains predictions for dependent variables.\n")
  
  if(!exists("gsom_object$data")) included = "is not included.\n"
  else included = "is included.\n"
  
  cat("Training data", included)
}