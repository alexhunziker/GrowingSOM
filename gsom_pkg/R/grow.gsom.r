#######################################
#GSOM - Growing Self Organizing Maps
#train.r
#26/10/16 - Alex Hunziker
#######################################

# The Functions in this File are required in order to train the gsom model.
# gsom.train() is the main function, which should be called by the user.
# The performance intensive loop has been outsourced to C for performance reasons.

#Mainly calls the C loop and processes returned data
grow.gsom <- function(gsom_model, df, repet, spreadFactor, alpha, gridsize, nhood, grow){
  
  # Set some parameters
  lentr <- 10000
  lentn <- 10000
  
  lrinit <- 0.9
  alpha <- 0.9 #Learning Rate Depreciation factor.
  radius <- 3 #Initial Radius. Missing feature.
  if(grow==2) radius = sqrt(gridsize)
  
  initnodes <- gridsize*gridsize
  
  df <- as.matrix(df)
  
  weights <- matrix(0, nrow=lentn, ncol=ncol(df))
  weights[1:initnodes,] <- runif(initnodes*ncol(df))
  
  distnd <- rep(0, lentn) #Error per node
  freq <- rep(0, lentn) #Frequ of nodes
  
  gt = -ncol(df) * log(spreadFactor) * 0.01*nrow(df)
  
  npos <- matrix(0, nrow=lentn, ncol=2)
  #npos[1:initnodes,] <- c(0, 1, 1, 0, 1, 0, 1, 0)
  if(nhood=="rect"){
    for(i in 1:gridsize){
      for(j in 1:gridsize) npos[gridsize*(i-1)+j,] = c(i, j)
    }
  }else{
    for(i in 1:gridsize){
      if(i/2 - rounded(i/2,0) != 0) q=0.5
      else q=0
      for(j in 1:gridsize) npos[gridsize*(i-1)*j,] = c(i+q, j)
    }
  }
  
  if(repet > lentr) error("Max nr of iterations exceeded.")
  
  currtrain <- matrix(0, nrow=lentr, ncol=5)
  
  outc = .C("som_train_loop",
            df = as.double(df),
            weights = as.double(weights),
            distnd = as.double(distnd),
            prep = as.integer(repet), #repetitions
            plendf = as.integer(nrow(df)),
            plennd = as.integer(initnodes),
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
            plentr = as.integer(lentr), #Max num of iterations
            hex = as.integer(0),
            grow = as.integer(grow)
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
