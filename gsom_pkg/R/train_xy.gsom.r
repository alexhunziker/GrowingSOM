#######################################
#GSOM - Growing Self Organizing Maps
#train_xy
#02/11/16 - Alex Hunziker
#######################################

# The Functions in this File are required in order to train the gsom model.
# gsom.train() is the main function, which should be called by the user.
# The performance intensive loop has been outsourced to C for performance reasons.

train_xy.gsom <- function(data, y, spreadFactor=0.8, keepdata=FALSE, iterations=50, alpha=0.5, gridsize = FALSE, nhood= "rect", initrad = NULL, ...){
  
  # Normalizing the training or testdata (min/max) in order to balance the impact
  # of the different properties of the dataframe
  minx <- apply(data, 2, function(x){min(x)})
  maxx <- apply(data, 2, function(x){max(x)})
  df <- t(apply(data, 1, function(x){(x-minx)/ifelse(maxx==minx,1,(maxx-minx))}))
  
  if(is.vector(y)) cy = 1
  else cy = ncol(y)
  y <- as.matrix(y, ncol = cy)
  miny <- apply(y, 2, function(x){min(x)})
  maxy <- apply(y, 2, function(x){max(x)})
  if(cy > 1) y <- t(apply(y, 1, function(x){(x-miny)/ifelse(maxy==miny,1,(maxy-miny))}))
  else y <- matrix(apply(y, 1, function(x){(x-miny)/ifelse(maxy==miny,1,(maxy-miny))}), ncol=1)
  
  if(gridsize==FALSE) grow=1
  else grow=2
  
  if(gridsize == FALSE) gridsize=2
  else if(!is.numeric(gridsize)){
    error("Grid size must be nummeric (for classical kohonen map) or FALSE (for Growing SOM).")
  }
  
  t1 <- Sys.time()
  gsom_object <- grow_xy.gsom(y, df, iterations, spreadFactor, alpha, gridsize, nhood, grow, initrad = initrad)
  t2 <- Sys.time()
  print(t2-t1)

	gsom_object$nodes$codes <- t(apply(gsom_object$nodes$codes, 1, function(x){(x*ifelse(maxx==minx,1,(maxx-minx))+minx)}))
	if(cy > 1) gsom_object$nodes$predict <- t(apply(gsom_object$nodes$predict, 1, function(x){(x*ifelse(maxy==miny,1,(maxy-miny))+miny)}))
  else gsom_object$nodes$predict <- matrix(apply(gsom_object$nodes$predict, 1, function(x){(x-miny)/ifelse(maxy==miny,1,(maxy-miny))}), ncol=1)
  
	#convert data types if only one variable to be predicted.
	if(cy==1) {
	  gsom_object$nodes$predict <- data.frame(gsom_object$nodes$predict)
	  data.frame(apply(gsom_object$nodes$predict, 1, function(x){(x*ifelse(maxy==miny,1,(maxy-miny))+miny)}))
	}
	else gsom_object$nodes$predict <- t(apply(gsom_object$nodes$predict, 1, function(x){(x*ifelse(maxy==miny,1,(maxy-miny))+miny)}))
  colnames(gsom_object$nodes$predict) = colnames(y)
	
  norm_param <- data.frame(min = minx, max = maxx)
  norm_param_y <- data.frame(miny = miny, maxy = maxy)
  gsom_object[["norm_param"]] <- norm_param
  gsom_object[["norm_param_y"]] <- norm_param_y
  
  if(keepdata==TRUE){
    gsom_object[["data"]] = data
  }
  
  class(gsom_object) = "gsom"
  
  return(gsom_object)
  
}
