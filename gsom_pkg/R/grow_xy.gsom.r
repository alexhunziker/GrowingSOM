#######################################
#GSOM - Growing Self Organizing Maps
#grow_xy
#02/11/16 - Alex Hunziker
#######################################

# The Functions in this File are required in order to train the gsom model.
# gsom.train() is the main function, which should be called by the user.
# The performance intensive loop has been outsourced to C for performance reasons.

#Mainly calls the C loop and processes returned data
grow_xy.gsom <- function(y, df, repet, spreadFactor, alpha, gridsize, nhood, grow){

  # Set some parameters
  lentr <- 10000
  lentn <- 10000
  
  lrinit <- 0.9
  alpha <- 0.9 #Learning Rate Depreciation factor.
  radius <- 3 #Initial Radius. Missing feature.
  if(grow==2) radius = sqrt(gridsize)
  
  initnodes <- gridsize*gridsize
  
  df <- as.matrix(df)
  
  codes <- matrix(0, nrow=lentn, ncol=ncol(df))
  codes[1:initnodes,] <- runif(initnodes*ncol(df))
                                 
  predict <- matrix(0, nrow=lentn, ncol=ncol(y))
  codes[1:initnodes,] <- runif(initnodes*ncol(df))
  
  
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
    error("Missing Feature.")
    for(i in 1:gridsize){
      if(i/2 - rounded(i/2,0) != 0) q=0.5
      else q=0
      for(j in 1:gridsize) npos[gridsize*(i-1)*j,] = c(i+q, (j+q)*1.5)
    }
  }
  
  if(repet > lentr) error("Max nr of iterations exceeded.")
  
  currtrain <- matrix(0, nrow=lentr, ncol=5)
  
  as.integer(nrow(y))
  as.integer(ncol(y))
  as.double(predict)
  
  outc = .C("som_train_loop_xy",
            df = as.double(df),
            codes = as.double(codes),
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
            grow = as.integer(grow),
            y = as.double(y),
            leny = as.integer(nrow(y)),
            pydim = as.integer(ncol(y)),
            predict = as.double(predict)
  )
  
  training <- matrix(outc$currtrain, ncol=5)
  training <- training[1:repet,]
  colnames(training) <- c("iteration", "training_stage", "meandist", "num_of_nodes", "nodegrowth")
  
  codes <- matrix(outc$codes, ncol=ncol(df))
  codes <- codes[1:outc$plennd,]
  colnames(codes) <- colnames(df)
  
  predict <- matrix(outc$predict, ncol=ncol(y))
  predict <- predict[1:outc$plennd,]
  colnames(predict) <- colnames(y)
  
  npos <- matrix(outc$npos, ncol=2)
  npos <- matrix(npos[1:outc$plennd,], ncol=2)
  colnames(npos) <- c("x", "y")
  
  freq <- outc$freq[1:outc$plennd]
  distnd <- outc$distnd[1:outc$plennd]
  
  nodes <- list(
    position = data.frame(npos),
    codes = data.frame(codes),
    predict = data.frame(predict),
    error = distnd,
    freq = freq
  )
  
  training = data.frame(training)
  GT = outc$gt
  
  gsom_model <- list(nodes = nodes, training = training, GT = GT, norm_param=0)
  class(gsom_model) <- "gsom"
  
  return(gsom_model)
  
}
