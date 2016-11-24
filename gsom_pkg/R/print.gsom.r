print.gsom <- function(gsom_object){
  
  print(gsom_object)
  if(gsom_object$nhood == "rect") topology = "rectangular"
  else topology = "hexagonal"
  
  cat("GrowingSOM (gsom) object with", nrow(gsom_object$nodes$position), "nodes and a", topology, "topology.\n")
  
  if(exists(gsom_object$nodes$prediction)) cat("Model contains predictions for dependent variables.\n")
  
  if(!exists("gsom_object$data")) included = "is not included.\n"
  else included = "is included.\n"
  
  if(exists("gsom_object$mapped")) contains = "This object contains data mapped to an existing gsom object."
  else contains = "This object contains a trained gsom map."
  
  cat(contains, "Training data", included)
}