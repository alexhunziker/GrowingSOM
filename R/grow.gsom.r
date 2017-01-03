#######################################
#GSOM - Growing Self Organizing Maps
#train.r
#26/10/16 - Alex Hunziker
#######################################

# The Functions in this File are required in order to train the gsom model.
# gsom.train() is the main function, which should be called by the user.
# The performance intensive loop has been outsourced to C for performance reasons.

#Mainly calls the C loop and processes returned data
grow.gsom <- function(gsom_model, df, repet, spreadFactor, alpha, beta, gridsize, nhood, grow, initrad){
  
  # Set some parameters
  lentr <- 10000
  lentn <- 10000
  
  lrinit <- 0.9
  
 initnodes <- gridsize*gridsize

	if(is.null(initrad)){
		if(grow==1) radius = 3  # Arbitrary default value
		else radius = sqrt(gridsize)
	} else{
	  radius = initrad
	}
  
  df <- as.matrix(df)
  
  codes <- matrix(0, nrow=lentn, ncol=ncol(df))
  codes[1:initnodes,] <- runif(initnodes*ncol(df))
  
  distnd <- rep(0, lentn) #Error per node
  freq <- rep(0, lentn) #Frequ of nodes
  
  gt = -ncol(df) * log(spreadFactor) * 0.01*nrow(df)
  
  npos <- matrix(0, nrow=lentn, ncol=2)
  #npos[1:initnodes,] <- c(0, 1, 1, 0, 1, 0, 1, 0)
  
  
  if(nhood=="rect"){
    hex = 0
    for(i in 1:gridsize){
      for(j in 1:gridsize) npos[gridsize*(i-1)+j,] = c(i, j)
    }
  }else{
    
    hex = 1
    counter = 1
    
    for(i in 1:gridsize){
      for(j in 1:floor(gridsize/2)){
        npos[counter,] = c(i, 1.5*j)
        npos[counter +1,] = c((i-0.5), (1.5*j+0.75))
        counter = counter + 2
      }
      if(gridsize %% 2 != 0){
        npos[counter,] = c(i, 1.5*(j+1))
        counter = counter + 1
      }
    }
  }
  
  if(repet > lentr) stop("Max nr of iterations exceeded.")
  
  currtrain <- matrix(0, nrow=lentr, ncol=5)
  
  outc = .C("som_train_loop",
            df = as.double(df),
            codes = as.double(codes),
            distnd = as.double(distnd),
            prep = as.integer(repet), #repetitions
            plendf = as.integer(nrow(df)),
            plennd = as.integer(initnodes),
            plrinit = as.double(lrinit),
            freq = as.double(freq),
            alpha = as.double(alpha), #for lr depreciation
            beta = as.double(beta), #for depreciation of neighbour learning
            pdim = as.integer(ncol(df)),
            gt = as.double(gt), # Growth Rate
            npos = as.double(npos),
            pradius = as.double(radius), 
            plentn = as.integer(lentn), #Max num of nodes
            plentd = as.integer(nrow(df)),
            currtrain = as.double(currtrain),
            plentr = as.integer(lentr), #Max num of iterations
            hex = as.integer(hex),
            grow = as.integer(grow),
            PACKAGE = "GrowingSOM"
  )
  
  training <- matrix(outc$currtrain, ncol=5)
  training <- training[1:repet,]
  colnames(training) <- c("iteration", "training_stage", "meandist", "num_of_nodes", "nodegrowth")
  
  codes <- matrix(outc$codes, ncol=ncol(df))
  codes <- codes[1:outc$plennd,]
  colnames(codes) <- colnames(df)
  
  npos <- matrix(outc$npos, ncol=2)
  npos <- matrix(npos[1:outc$plennd,], ncol=2)
  colnames(npos) <- c("x", "y")
  
  freq <- outc$freq[1:outc$plennd]
  distnd <- outc$distnd[1:outc$plennd]
  
  nodes <- list(
    position = data.frame(npos),
    codes = data.frame(codes),
    distance = distnd,
    freq = freq
  )
  
  training = data.frame(training)
  GT = outc$gt
  
  gsom_model <- list(nodes = nodes, training = training, GT = GT, norm_param=0, nhood = nhood)
  class(gsom_model) <- "gsom"
  
  return(gsom_model)
  
}
