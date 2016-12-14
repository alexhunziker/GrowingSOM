#######################################
#GSOM - Growing Self Organizing Maps
#summary.r
#11/10/16 - Alex Hunziker
#######################################

# This function gives a (text) summary of a gsom_object
# Requires: trained gsom_object
# Returns: nothing

summary.gsom <- function(gsom_object){
  
  if(!is.null(gsom_object$training)){
    
    # Information about of Training Dataset
    observations <- sum(gsom_object$nodes$freq)
    dimenstions <- ncol(gsom_object$nodes$codes)
    
    if(!is.null(gsom_object$nodes$predict)) predict = paste("and", ncol(gsom_object$nodes$predict), "predicted variable(s)")
    else predict = ""
    
    # Total number of nodes
    nodes <- nrow(gsom_object$nodes$position)
    
    # Average Distance, No of Iterations
    iterations <- nrow(gsom_object$training)
    distance <- gsom_object$training$meandist[iterations]
    
    if(is.null(gsom_object[["data"]])) sorstat = "is not"
    else sorstat = "is"
    
    # Print
    cat("Growing SOM map with", nodes, "nodes.\n")
    cat("Training data used:", observations, "Observations with", dimenstions, "Dimensions", predict, "\n")
    cat("Training data", sorstat, "stored in the model.\n")
    cat("Mean Distance to the closest unit in the map is:", distance, "(after", iterations, "iterations)\n")
  
  } else {
    
    if(!is.null(gsom_object$nodes$predict)){
      distance <- mean(gsom_object$prediction$dist)
      observations <- length(gsom_object$prediction$dist)
    } else {
      distance <- mean(gsom_object$mapped$dist)
      observations <- length(gsom_object$mapped$dist)
    }
    
    nodes <- nrow(gsom_object$nodes$position)
    cat("Growing SOM map with", nodes, "nodes.\n")
    
    cat(observations, "have been mapped onto a trained gsom map.\n")
    
    if(!is.null(gsom_object$nodes$predict)){
      depvarno <- ncol(gsom_object$nodes$predict)
      cat("Predictions for", depvarno, "variables are stored.\n")
    }
   
    if(is.null(gsom_object[["data"]])) sorstat = "is not"
    else sorstat = "is"
    cat("Training data", sorstat, "stored in the model.\n")
    
    cat("Mean Distance to the closest unit in the map is:", distance, "(for the mapped observations)\n")
    
  }

  
}
