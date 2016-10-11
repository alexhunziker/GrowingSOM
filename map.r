#######################################
#GSOM - Growing Self Organizing Maps
#map.r
#11/10/16 - Alex Hunziker
#######################################

# This Function maps new data onto a trained gsom_model without adjusting 
# the gsom_model itself.
# Requires: trained gsom_model and testdata (DataFrame)
# Returns: mapped_data, which includes the nodes with position of nodes, frequency and average errors
#   as well as the error and winning node for each node of the testdata

gsom.map <- function(gsom_model, newdf) {
  
  nodes <- list(
    position = gsom_model$nodes$position,
    weight = gsom_model$nodes$weight,
    error = 0,
    freq = 0
  )
  
  observation <- list(
    winner = numeric(),
    error = numeric()
  )
  
  mapped_data <- list(nodes = nodes, observation = observation)
  
  for(j in 1:nrow(newdf)){
    errors <- sqrt(rowSums(sweep(gsom_model$nodes$weight, MARGIN = 2, newdf[j,], FUN="-")^2, dims=1))
    
    minError=min(errors)
    winner <- which(grepl(minError, errors))
    if(length(winner)>1) winner <- winner[1]
    
    mapped_data$nodes$freq[winner] <- gsom_model$nodes$freq[winner] + 1
    mapped_data$nodes$avgerror[winner] <- mapped_data$nodes$avgerror[winner] + minError
    
    mapped_data$observation[winner,] <- c(winner = winner, error = minError)
  }
  
  mapped_data$nodes$avgerror <- mapped_data$nodes$avgerror / mapped_data$nodes$freq
  
  return(mapped_data)
}