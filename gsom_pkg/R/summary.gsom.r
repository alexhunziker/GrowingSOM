#######################################
#GSOM - Growing Self Organizing Maps
#summary.r
#11/10/16 - Alex Hunziker
#######################################

# This function gives a (text) summary of a gsom_model
# Requires: trained gsom_model
# Returns: nothing

summary.gsom <- function(gsom_model){

  # Information about of Training Dataset
  observations <- sum(gsom_model$nodes$freq)
  dimenstions <- ncol(gsom_model$nodes$codes)
  
  if(exists("gsom_object$predict")) predict = paste("and", ncol(gsom_object$predict), "predicted variables")
  else predict = ""
  
  # Total number of nodes
  nodes <- nrow(gsom_model$nodes$position)
  
  # Average Distance, No of Iterations
  iterations <- nrow(gsom_model$training)
  distance <- gsom_model$training$meandist[iterations]
  
  if(is.null(gsom_model[["data"]])) sorstat = "is not"
  else sorstat = "is"
  
  # Print
  print(paste("Growing SOM map with", nodes, "nodes."))
  print(paste("Training data used:", observations, "Observations with", dimenstions, "Dimensions", predict))
  print(paste("Training data", sorstat, "stored in the model."))
  print(paste("Mean Distance to the closest unit in the map is:", distance, "(after", iterations, "iterations)"))
}