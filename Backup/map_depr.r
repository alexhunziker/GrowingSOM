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

gsom.map <- function(newdf, gsom_model) {
  
  nodes <- list(
    position = gsom_model$nodes$position,
    weight = gsom_model$nodes$weight,
    avgerror = rep(0, nrow(gsom_model$nodes$position)),
    freq = rep(0, nrow(gsom_model$nodes$position))
  )
  
  observation <- data.frame(
    winner = numeric(),
    error = numeric()
  )
  
  mapped_data <- list(nodes = nodes, observation = observation)
  
  newdf <- t(apply(newdf, 1, function(x){(x-gsom_model$norm_param$min)/(gsom_model$norm_param$max-gsom_model$norm_param$min)}))
  
  for(j in 1:nrow(newdf)){
    errors <- sqrt(rowSums(sweep(gsom_model$nodes$weight, MARGIN = 2, newdf[j,], FUN="-")^2, dims=1))
    
    minError=min(errors)
    winner <- which(grepl(minError, errors))
    #print(paste("Map:", winner))
    if(length(winner)>1) winner <- winner[1]
    
    mapped_data$nodes$freq[winner] <- mapped_data$nodes$freq[winner] + 1
    mapped_data$nodes$avgerror[winner] <- mapped_data$nodes$avgerror[winner] + minError
    
    mapped_data$observation[j,] <- c(winner = winner, error = minError)
  }
  
  mapped_data$nodes$avgerror <- mapped_data$nodes$avgerror / mapped_data$nodes$freq
  
  return(mapped_data)
}