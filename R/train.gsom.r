#######################################
# train.gsom - GrowingSOM
# Alex Hunziker - 2017
#######################################

# The Functions in this File are required in order to train the gsom model.

train.gsom <- function(data, spreadFactor=0.8, keepdata=FALSE, iterations=50, alpha=0.9, beta=0.5, gridsize = FALSE, nhood= "rect", initrad = NULL, ...){
  
  # Normalizing the training or testdata (mean/sd) in order to balance the impact
  # of the different properties of the dataframe
  mean <- apply(data, 2, function(x){mean(x)})
  sd <- apply(data, 2, function(x){sd(x)})
  df <- t(apply(data, 1, function(x){(x-mean)/ifelse(sd==0,1,sd)}))
  
  if(gridsize==FALSE) grow=1
  else grow=2
  
  if(gridsize == FALSE) gridsize=2
  else if(!is.numeric(gridsize)){
    stop("Grid size must be nummeric (for classical kohonen map) or FALSE (for Growing SOM).")
  }

  # Call Function grow.gsom
  gsom_object <- grow.gsom(gsom_object, df, iterations, spreadFactor, alpha, beta, gridsize = gridsize, nhood=nhood, grow=grow, initrad = initrad)

  gsom_object$nodes$codes <- t(apply(gsom_object$nodes$codes, 1, function(x){(x*sd+mean)}))
  
  norm_param <- data.frame(mean = mean, sd = sd)
  gsom_object[["norm_param"]] <- norm_param
  
  if(keepdata==TRUE){
    gsom_object[["data"]] = data
  }
  
  class(gsom_object) = "gsom"
  
  return(gsom_object)
  
}
