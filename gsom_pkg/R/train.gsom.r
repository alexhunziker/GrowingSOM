#######################################
#GSOM - Growing Self Organizing Maps
#train.r
#26/10/16 - Alex Hunziker
#######################################

# The Functions in this File are required in order to train the gsom model.
# gsom.train() is the main function, which should be called by the user.
# The performance intensive loop has been outsourced to C for performance reasons.

train.gsom <- function(data, spreadFactor=0.8, keepdata=FALSE, iterations=50, alpha=0.5, gridsize = FALSE, nhood= "rect", ...){
  
  # Normalizing the training or testdata (min/max) in order to balance the impact
  # of the different properties of the dataframe
  min <- apply(data, 2, function(x){min(x)})
  max <- apply(data, 2, function(x){max(x)})
  df <- t(apply(data, 1, function(x){(x-min)/ifelse(max==min,1,(max-min))}))
  
  if(gridsize==FALSE) grow=1
  else grow=2
  
  if(gridsize == FALSE) gridsize=2
  else if(!is.numeric(gridsize)){
    error("Grid size must be nummeric (for classical kohonen map) or FALSE (for Growing SOM).")
  }

  t1 <- Sys.time()
  gsom_model <- grow.gsom(gsom_model, df, iterations, spreadFactor, alpha, gridsize = gridsize, nhood=nhood, grow=grow)
  t2 <- Sys.time()
  print(t2-t1)
  
  norm_param <- data.frame(min = min, max = max)
  gsom_model[["norm_param"]] <- norm_param
  
  if(keepdata==TRUE){
    gsom_model[["data"]] = data
  }
  
  class(gsom_model) = "gsom"
  
  return(gsom_model)
  
}
