#######################################
#GSOM - Growing Self Organizing Maps
#train.r
#26/10/16 - Alex Hunziker
#######################################

# The Functions in this File are required in order to train the gsom model.
# gsom.train() is the main function, which should be called by the user.
# The performance intensive loop has been outsourced to C for performance reasons.

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

  t1 <- Sys.time()
  gsom_object <- grow.gsom(gsom_object, df, iterations, spreadFactor, alpha, beta, gridsize = gridsize, nhood=nhood, grow=grow, initrad = initrad)
  t2 <- Sys.time()
  print(t2-t1)

	gsom_object$nodes$codes <- t(apply(gsom_object$nodes$codes, 1, function(x){(x*sd+mean)}))
  
  norm_param <- data.frame(mean = mean, sd = sd)
  gsom_object[["norm_param"]] <- norm_param
  
  if(keepdata==TRUE){
    gsom_object[["data"]] = data
  }
  
  class(gsom_object) = "gsom"
  
  return(gsom_object)
  
}
