#######################################
#GSOM - Growing Self Organizing Maps
#summary.r
#11/10/16 - Alex Hunziker
#######################################

# This function gives a (text) summary of a gsom_model
# Requires: trained gsom_model
# Returns: nothing

gsom.summary <- function(gsom_model){

  # Information about of Training Dataset
  observations <- sum(gsom_model$nodes$freq)
  dimenstions <- ncol(gsom_model$nodes$weights)
  
  # Total number of nodes
  nodes <- ncol(gsom_model$nodes$position)
  
  # Average Distance, No of Iterations
  iterations <- nrow(gsom_model$training$error) #???
  distance <- gsom_model$training$error[iterations]   #???
  
  # Print
  print(paste("Growing SOM map with", nodes, "nodes."))
  print(paste("Training data used:", observations, "Observations with", dimenstions, "Dimensions"))
  print("Training is not stored in the model.")
  print(paste("Mean Distance to the closest unit in the map is:", distance))
}