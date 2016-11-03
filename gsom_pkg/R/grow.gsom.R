grow.gsom <-
function(gsom_model, df, repet, spreadFactor){
  
  # Set some parameters
  lentr <- 10000
  lentn <- 10000
  
  lrinit <- 0.9
  alpha <- 0.9 #Learning Rate Depreciation factor.
  radius <- 3 #Initial Radius. Missing feature.
  
  df <- as.matrix(df)
  
  weights <- matrix(0, nrow=lentn, ncol=ncol(df))
  weights[1:4,] <- runif(4*ncol(df))
  
  distnd <- rep(0, lentn) #Error per node
  freq <- rep(0, lentn) #Frequ of nodes
  
  gt = -ncol(df) * log(spreadFactor) * 0.01*nrow(df)
  
  npos <- matrix(0, nrow=lentn, ncol=2)
  npos[1:4,] <- c(0, 1, 1, 0, 1, 0, 1, 0)
  
  if(repet > lentr) error("Max nr of iterations exceeded.")
  
  currtrain <- matrix(0, nrow=lentr, ncol=5)
  
  outc = .C("som_train_loop",
            df = as.double(df),
            weights = as.double(weights),
            distnd = as.double(distnd),
            prep = as.integer(repet), #repetitions
            plendf = as.integer(nrow(df)),
            plennd = as.integer(4),
            plrinit = as.double(lrinit),
            freq = as.double(freq),
            alpha = as.double(alpha), #for lr depreciation
            pdim = as.integer(ncol(df)),
            gt = as.double(gt), # Growth Rate
            npos = as.double(npos),
            pradius = as.double(radius), 
            plentn = as.integer(lentn), #Max num of nodes
            plentd = as.integer(nrow(df)),
            currtrain = as.double(currtrain),
            plentr = as.integer(lentr) #Max num of iterations
  )
  
  training <- matrix(outc$currtrain, ncol=5)
  training <- training[1:repet,]
  colnames(training) <- c("iteration", "training_stage", "meandist", "num_of_nodes", "nodegrowth")
  
  weights <- matrix(outc$weights, ncol=ncol(df))
  weights <- weights[1:outc$plennd,]
  colnames(weights) <- colnames(df)
  
  npos <- matrix(outc$npos, ncol=2)
  npos <- matrix(npos[1:outc$plennd,], ncol=2)
  colnames(npos) <- c("x", "y")
  
  freq <- outc$freq[1:outc$plennd]
  distnd <- outc$distnd[1:outc$plennd]
  
  nodes <- list(
    position = data.frame(npos),
    weight = data.frame(weights),
    error = distnd,
    freq = freq
  )
  
  training = data.frame(training)
  GT = outc$gt
  
  gsom_model <- list(nodes = nodes, training = training, GT = GT, norm_param=0)
  class(gsom_model) <- "gsom"
  
  return(gsom_model)
  
}
