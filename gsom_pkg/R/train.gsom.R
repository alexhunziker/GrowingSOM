train.gsom <-
function(data, spreadFactor=0.5, keepdata=FALSE, iterations=50, alpha, ...){
  
  # Normalizing the training or testdata (min/max) in order to balance the impact
  # of the different properties of the dataframe
  min <- apply(data, 2, function(x){min(x)})
  max <- apply(data, 2, function(x){max(x)})
  df <- t(apply(data, 1, function(x){(x-min)/ifelse(max==min,1,(max-min))}))
  
  t1 <- Sys.time()
  gsom_model <- grow.gsom(gsom_model, df, iterations, spreadFactor)
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
