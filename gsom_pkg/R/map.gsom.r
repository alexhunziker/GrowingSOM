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

map.gsom <- function(gsom_model, df, retaindata=FALSE){
  
  # Normalizing the training or testdata (min/max) in order to balance the impact
  # of the different properties of the dataframe
  min <- gsom_model$norm_param$min
  max <- gsom_model$norm_param$max
  df <- t(apply(df, 1, function(x){(x-min)/(max-min)}))
  
  codes <- rep(0, times=nrow(df))
  ndist <- rep(0, times=nrow(df))
  freq <- rep(0, times=nrow(gsom_model$nodes$weight))
  
  outc = .C("map_data",
            plendf = as.integer(nrow(df)),
            lennd = as.integer(nrow(gsom_model$nodes$weight)),
            dim = as.integer(ncol(gsom_model$nodes$weight)),
            df = as.double(df),
            weights =as.double(as.matrix(gsom_model$nodes$weight)), 
            codes = as.double(codes),
            ndist = as.double(ndist),
            freq = as.double(freq)
  )
  
  dist <- outc$ndist
  code <- outc$codes
  print(code)

  gsom_mapped = list();
  gsom_mapped[["nodes"]] = gsom_model$nodes
  gsom_mapped[["nodes"]]$error = NULL
  gsom_mapped[["nodes"]]$freq = outc$freq
  gsom_mapped[["mapped"]] = data.frame(codes=code, dist=dist)
  gsom_mapped[["predict"]] = predict
  gsom_mapped[["norm_param"]] = gsom_model$scale
  if(retaindata) gsom_mapped[["data"]] == df;
  
  class(gsom_mapped) = "gsom_mapped"
  
  return(gsom_mapped)
  
}