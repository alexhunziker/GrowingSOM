#######################################
#GSOM - Growing Self Organizing Maps
#train.r
#26/10/16 - Alex Hunziker
#######################################

# The Functions in this File are required in order to train the gsom model.
# gsom.train() is the main function, which should be called by the user.
# The performance intensive loop has been outsourced to C for performance reasons.

gsom.train <- function(data, spreadFactor=0.5, keepdata=FALSE, iterations=50, alpha, ...){
  
  # Normalizing the training or testdata (min/max) in order to balance the impact
  # of the different properties of the dataframe
  min <- apply(data, 2, function(x){min(x)})
  max <- apply(data, 2, function(x){max(x)})
  df <- t(apply(data, 1, function(x){(x-min)/(max-min)}))
  
  gsom_model <- gsom.grow(gsom_model, df, iterations, spreadFactor)
  
  norm_param <- data.frame(min = min, max = max)
  gsom_model[["norm_param"]] <- norm_param
  
  if(keepdata==TRUE){
    gsom_model[["data"]] = data
  }
  
  return(gsom_model)
  
}


#Mainly calls the C loop and processes returned data
gsom.grow <- function(gsom_model, df, repet, spreadFactor){
  
  # Set some parameters
  lentr <- 10000
  lentn <- 10000
  
  lrinit <- 0.9
  alpha <- 0.99 #Learning Rate Depreciation factor.
  radius <- 3 #Initial Radius. Missing feature.
  
  df <- as.matrix(df)
  
  weights <- matrix(0, nrow=lentn, ncol=ncol(df))
  weights[1:4,] <- runif(4*ncol(df))
  
  distnd <- rep(0, lentn) #Error per node
  freq <- rep(0, lentn) #Frequ of nodes
  
  gt = -ncol(df) * log(spreadFactor)
  
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
            pradius = as.integer(radius), 
            plentn = as.integer(lentn), #Max num of nodes
            plentd = as.integer(nrow(df)),
            currtrain = as.double(currtrain),
            plentr = as.integer(lentr) #Max num of iterations
  )
  
  training <- matrix(outc$currtrain, ncol=5)
  training <- training[1:repet,]
  colnames(training) <- c("iteration", "training_stage", "meandist", "num_of_nodes", "nodegrowth")
  
  weights <- matrix(outc$weights, ncol=4)
  weights <- weights[1:outc$plennd,]
  
  npos <- matrix(outc$npos, ncol=2)
  npos <- npos[1:outc$plennd,]
  colnames(npos) <- c("x", "y")
  
  freq <- outc$freq[1:outc$plennd]
  distnd <- outc$distnd
  
  nodes <- list(
    position = data.frame(npos),
    weight = data.frame(weights),
    error = distnd,
    freq = freq
  )
  
  training = data.frame(training)
  GT = outc$gt
  
  gsom_model <- list(nodes = nodes, training = training, GT = GT, norm_param=0)
  
  return(gsom_model)
  
}