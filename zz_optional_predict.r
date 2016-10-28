############
#predict.r
############

# More pseudocode than anything else

predict <- function(df, gsom_model, retaindata=FALSE){
  
  if(!gsom_model[["dep_val"]]) error("Wrong input, use trained gsom model for dependent values.")
  
  # Normalizing the training or testdata (min/max) in order to balance the impact
  # of the different properties of the dataframe
  min <- gsom_model$norm_param$min
  max <- gsom_model$norm_param$max
  df <- t(apply(data, 1, function(x){(x-min)/(max-min)}))
  
  codes <- rep(0, times=nrow(df))
  ndist <- rep(0, times=nrow(df))
  
  outc = .C("map_data",
            plendf = as.integer(nrow(df)),
            lennd = as.integer(nrow(gsom_model$nodes$weights)),
            dim = as.integer(ncol(gsom_model$nodes$weights)),
            df = as.double(df),
            weights =as.double(gsom_model$nodes$weights), 
            codes = as.double(codes),
            ndist = as.double(ndist)
  )
  
  dist <- outc$ndist
  code <- outc$codes
  
  predict <- For each element of df, based on the code, the depval. 
  
  gsom_mapped = list();
  gsom_mapped[["nodes"]] = gsom_model$nodes
  gsom_mapped[["mapped"]] = as.data.frame(codes, dist, row.names=c(codes, dist))
  gsom_mapped[["predict"]] = predict
  gsom_mapped[["norm_param"]] = gsom_model$scale
  if(retaindata) gsom_mapped[["data"]] == df;
  
  class(gsom_mapped) = "gsom_mapped"
  
  return(gsom_mapped)
  
}