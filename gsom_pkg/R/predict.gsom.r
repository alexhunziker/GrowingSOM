############
#predict.r
############

# More pseudocode than anything else

predict.gsom <- function(gsom_object, df, retaindata=FALSE){
  
  if(is.null(gsom_object$nodes$predict)) stop("Wrong input, use trained gsom model for dependent values.")
  
  # Normalizing the training or testdata (min/max) in order to balance the impact
  # of the different properties of the dataframe
  min <- gsom_object$norm_param$min
  max <- gsom_object$norm_param$max
  df <- as.matrix(df)
  df <- t(apply(df, 1, function(x){(x-min)/ifelse(max==min,1,(max-min))}))
  gsom_object$nodes$codes <- t(apply(gsom_object$nodes$codes, 1, function(x){(x-min)/ifelse(max==min,1,(max-min))}))
  miny <- gsom_object$norm_param_y$miny
  maxy <- gsom_object$norm_param_y$maxy

  bmn <- rep(0, times=nrow(df))
  ndist <- rep(0, times=nrow(df))
  freq <- rep(0, times=nrow(gsom_object$nodes$codes))
  
  outc = .C("map_data",
            plendf = as.integer(nrow(df)),
            lennd = as.integer(nrow(gsom_object$nodes$codes)),
            dim = as.integer(ncol(gsom_object$nodes$codes)),
            df = as.double(df),
            codes =as.double(as.matrix(gsom_object$nodes$codes)), 
            bmn = as.double(bmn),
            ndist = as.double(ndist),
            freq = as.double(freq)
  )
  
  dist <- outc$ndist
  bmn <- outc$bmn
  bmn = matrix(bmn, ncol= 1)
  
  predict = data.frame(matrix(ncol=ncol(gsom_object$nodes$predict), nrow=nrow(df)))
  for(i in (1:nrow(df))){
    predict[i,] = gsom_object$nodes$predict[bmn[i,1],]
  }
  preditc = 0;
  
  cy = ncol(gsom_object$nodes$predict)

  gsom_object$nodes$codes <- t(apply(gsom_object$nodes$codes, 1, function(x){(x*ifelse(max==min,1,(max-min))+min)}))
  
  gsom_mapped = list();
  gsom_mapped[["nodes"]] = gsom_object$nodes
  gsom_mapped[["nodes"]]$error = NULL
  gsom_mapped[["nodes"]]$freq = outc$freq
  gsom_mapped[["prediction"]] = data.frame(bmn=bmn, dist=dist, predict=predict)
  gsom_mapped[["norm_param"]] = gsom_object$scale
  if(retaindata) gsom_mapped[["data"]] == df;
  
  class(gsom_mapped) = "gsom"
  
  return(gsom_mapped)
  
}
