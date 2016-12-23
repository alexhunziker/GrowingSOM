#######################################
#GSOM - Growing Self Organizing Maps
#summary.r
#11/10/16 - Alex Hunziker
#######################################

# This function gives a (text) summary of a gsom_object
# Requires: trained gsom_object
# Returns: nothing

summary.gsom <- function(object, ...){
  
  if(!is.null(object$training)){
    
    # Information about of Training Dataset
    observations <- sum(object$nodes$freq)
    dimenstions <- ncol(object$nodes$codes)
    
    if(!is.null(object$nodes$predict)) predict = paste("and", ncol(object$nodes$predict), "predicted variable(s)")
    else predict = ""
    
    # Total number of nodes
    nodes <- nrow(object$nodes$position)
    
    # Average Distance, No of Iterations
    iterations <- nrow(object$training)
    distance <- object$training$meandist[iterations]
    
    if(is.null(object[["data"]])) sorstat = "is not"
    else sorstat = "is"
    
    # Print
    cat("Growing SOM map with", nodes, "nodes.\n")
    cat("Training data used:", observations, "Observations with", dimenstions, "Dimensions", predict, "\n")
    cat("Training data", sorstat, "stored in the model.\n")
    cat("Mean Distance to the closest unit in the map is:", distance, "(after", iterations, "iterations)\n")
  
  } else {
    
    if(!is.null(object$nodes$predict)){
      distance <- mean(object$prediction$dist)
      observations <- length(object$prediction$dist)
    } else {
      distance <- mean(object$mapped$dist)
      observations <- length(object$mapped$dist)
    }
    
    nodes <- nrow(object$nodes$position)
    cat("Growing SOM map with", nodes, "nodes.\n")
    
    cat(observations, "have been mapped onto a trained gsom map.\n")
    
    if(!is.null(object$nodes$predict)){
      depvarno <- ncol(object$nodes$predict)
      cat("Predictions for", depvarno, "variables are stored.\n")
    }
   
    if(is.null(object[["data"]])) sorstat = "is not"
    else sorstat = "is"
    cat("Training data", sorstat, "stored in the model.\n")
    
    cat("Mean Distance to the closest unit in the map is:", distance, "(for the mapped observations)\n")
    
  }

  
}
