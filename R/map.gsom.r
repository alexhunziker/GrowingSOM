#######################################
#GSOM - Growing Self Organizing Maps
#map.r
#11/10/16 - Alex Hunziker
#######################################

# This Function maps new data onto a trained gsom_object without adjusting 
# the gsom_object itself.
# Requires: trained gsom_object and testdata (DataFrame)
# Returns: mapped_data, which includes the nodes with position of nodes, frequency and average errors
#   as well as the error and winning node for each node of the testdata

map.gsom <- function(gsom_object, df, retaindata=FALSE, ...){
  
  # Normalizing the training or testdata (mean/sd) in order to balance the impact
  # of the different properties of the dataframe
  mean <- gsom_object$norm_param$mean
  sd <- gsom_object$norm_param$sd
  dfs <- t(apply(df, 1, function(x){(x-mean)/ifelse(sd==0,1,sd)}))
  gsom_object$nodes$codes <- t(apply(gsom_object$nodes$codes, 1, function(x){(x-mean)/ifelse(sd==0,1,sd)}))
  
  bmn <- rep(0, times=nrow(df))
  ndist <- rep(0, times=nrow(df))
  freq <- rep(0, times=nrow(gsom_object$nodes$codes))
  
  outc = .C("map_data",
            plendf = as.integer(nrow(df)),
            lennd = as.integer(nrow(gsom_object$nodes$codes)),
            dim = as.integer(ncol(gsom_object$nodes$codes)),
            df = as.double(dfs),
            codes =as.double(as.matrix(gsom_object$nodes$codes)), 
            bmn = as.double(bmn),
            ndist = as.double(ndist),
            freq = as.double(freq)
  )
  
  dist <- outc$ndist
  bmn <- outc$bmn

  gsom_mapped = list();
  gsom_mapped[["nodes"]] = gsom_object$nodes
  gsom_mapped[["nodes"]]$distance = NULL
  gsom_mapped[["nodes"]]$freq = outc$freq
  gsom_mapped[["mapped"]] = data.frame(bmn=bmn, dist=dist)
  gsom_mapped[["norm_param"]] = gsom_object$norm_param
  if(retaindata) gsom_mapped[["data"]] = df

  gsom_mapped$nodes$codes <- t(apply(gsom_mapped$nodes$codes, 1, function(x){(x*sd+mean)}))
  
  class(gsom_mapped) = "gsom"
  
  return(gsom_mapped)
  
}
