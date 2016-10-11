#GSOM - Growing Self Organizing Maps
#train.r
#11/10/16 - Alex Hunziker
#######################################

# The Functions in this File are required in order to train the gsom model.
# gsom.train() is the main function, which should be called by the user.

gsom.train <- function(data, spreadFactor=0.5, keepdata=FALSE, rep=50){
  
  df <- gsom.normalize(data)
  gsom_model <- gsom.init(df, spreadFactor)
  gsom_model <- gsom.grow(gsom_model, df, rep)
  
  stop("Missing feature.")
  
  return(gsom_model)
  
}

# Normalizing the training or testdata (min/max) in order to balance the impact
# of the different properties of the dataframe
gsom.normalize <- function(df){
  min <- apply(df, 2, function(x){min(x)})
  max <- apply(df, 2, function(x){max(x)})
  df <- t(apply(df, 1, function(x){(x-min)/(max-min)}))
  return(df)
}

# Gnenerate the basic data structure for the GSOM model. 
# Creates for initial nodes that are initialized with random values.
# Returns trainable gsom model
gsom.init <- function(df, spreadFactor){
  #Initialize Network
  nodes <- list(
    position = data.frame(
      x=c(1, 1, 0, 0),
      y=c(1, 0, 1, 0)
    ),
    weight = matrix(runif((4*ncol(df))), 4),
    error = c(0, 0, 0, 0),
    freq = c(0, 0, 0, 0)
  )
  training = data.frame(
    iteration = integer(0),
    training_stage = integer(0),
    meandist = numeric(0),
    nodecount = integer(0),
    nodegrowth = integer(0)
  )
  herr = 0
  GT = -ncol(df) * log(spreadFactor)
  gsom_net <- list(nodes = nodes, training = training, herr = herr, GT = GT)
  #Calculate Growth Threshhold
  return(gsom_net)
}

# Main loop during stage 1 (Growing phase) of the training process
# Loops through training set and assigns every event a winner and keeps track of error and initializes growth of the grid
# Requires gsom_model, df. Tuning Parameters: Max no of repetitions, alpha (Parameter determining the learning rate (between 0 and 1))
# Returns gsom_model (fro smoothing phase)
gsom.grow <- function(gsom_model, df, rep=50, alpha = 0.5, ...){
  winner <- "Dummy"
  total_iterations <- nrow(df)*rep
  gsom_model$nodes$error <- rep(0, times=length(gsom_model$nodes$error))
  
  #Present Input
  #Replace with c code for speed gains
  for(i in 1:rep){
    
    nodegrow <- 0
    errorsum <- 0
    
    gsom_model$nodes$freq <- rep(0, times=length(gsom_model$nodes$freq))
    #gsom_model$nodes$error <- rep(0, times=length(gsom_model$nodes$error))
    t1 <- Sys.time()
    for(j in 1:nrow(df)){
      #Set Learning Rate
      #Recalculate Learning Rate
      #Contradiction to Paper. Learning rate is as defined in kohonen package
      learningRate <- (total_iterations-(i*j))/total_iterations * 1-(3.8/nrow(gsom_model$nodes$weight)) * alpha
      
      #CalcErrorValues
      errors <- sqrt(rowSums(sweep(gsom_model$nodes$weight, MARGIN = 2, df[j,], FUN="-")^2, dims=1))
      
      #Determine Winner
      minError=min(errors)
      errorsum = errorsum + minError
      winner <- which(grepl(minError, errors))
      #fixes problem when two errors are the same for two nodes
      if(length(winner)>1) winner <- winner[1]
      
      #Increase Error Value of Winner
      gsom_model$nodes$error[winner] <- gsom_model$nodes$error[winner] + minError
      gsom_model$nodes$freq[winner] <- gsom_model$nodes$freq[winner] + 1
      
      #Adjust Weights of Neighbourhood
      #Note: Just the direct neighbours are considered. The paper does not give information about the neighbourhood size used.
      winner.x <- gsom_model$nodes$position[winner,'x']
      winner.y <- gsom_model$nodes$position[winner,'y']
      self <- winner
      top <- which(gsom_model$nodes$position$x == winner.x & gsom_model$nodes$position$y == winner.y+1)
      bottom <- which(gsom_model$nodes$position$x == winner.x & gsom_model$nodes$position$y == winner.y-1)
      right <- which(gsom_model$nodes$position$x == winner.x+1 & gsom_model$nodes$position$y == winner.y)
      left <- which(gsom_model$nodes$position$x == winner.x-1 & gsom_model$nodes$position$y == winner.y)
      adjust <- c(self, top, bottom, left, right)
      for(k in 1:length(adjust)){
        gsom_model$nodes$weight[adjust[k],] <- gsom_model$nodes$weight[adjust[k],]+ (df[j,] - gsom_model$nodes$weight[adjust[k],]) * learningRate
      }
      
      #Update Maximum Error
      gsom_model$herr <- max(gsom_model$nodes$error[winner], gsom_model$herr)
      
      #Growth Condition
      if(gsom_model$herr > gsom_model$GT){
        gsom_model$herr <- 0
        if(length(adjust)>4){
          #Not a boundry node -> Spread it.
          gsom_model$nodes$error[self] = gsom_model$GT / 2
          #Note: This is ugly, since growth condidion is only checked, once the following nodes are winners again.
          #Also we are checking the topology once again...
          top <- which(gsom_model$nodes$position$x == winner.x & gsom_model$nodes$position$y == winner.y+1)
          bottom <- which(gsom_model$nodes$position$x == winner.x & gsom_model$nodes$position$y == winner.y-1)
          right <- which(gsom_model$nodes$position$x == winner.x+1 & gsom_model$nodes$position$y == winner.y)
          left <- which(gsom_model$nodes$position$x == winner.x-1 & gsom_model$nodes$position$y == winner.y)
          adjust <- c(self, top, bottom, left, right)
          for(m in 1:4){
            #The paper suggests values for gamma between 0 and 1
            gam <- 0.5
            gsom_model$nodes$error[adjust[m]] = gsom_model$nodes$error[adjust[m]] * (1+gam)
          }
        } else {
          #Boundry node -> Grow net.
          nodegrow = nodegrow + 1
          if(length(top)==0) gsom_model <- gsom.newnode(gsom_model, winner.x, winner.y+1)
          if(length(bottom)==0) gsom_model <- gsom.newnode(gsom_model, winner.x, winner.y-1)
          if(length(right)==0) gsom_model <- gsom.newnode(gsom_model, winner.x+1, winner.y)
          if(length(left)==0) gsom_model <- gsom.newnode(gsom_model, winner.x-1, winner.y)
          #Not clearly mentioned in the paper
          gsom_model$nodes$error[self]=0
        }
        #Distribute Error
        
      }
      
      rm(self, top, bottom, left, right)
      #print_crude(gsom_model)
      #Sys.sleep(0.1)
      
    }
    meandist <- errorsum/nrow(df)
    curr_train = c(iteration=i, training_stage=1, meandist=meandist, nodecount=nrow(gsom_model$nodes$position), nodegrow=nodegrow)
    print(curr_train)
    gsom_model$training[nrow(gsom_model$training)+1,] <- curr_train
    gsom_model <- gsom.emptyremove(gsom_model)
    t2 <- Sys.time()
    print(t2-t1)
    print_crude(gsom_model)
  }
  #Run until Growth is reduced to minimum level:
  #Initialize new node weights
  #Initialize learning rate
  return(gsom_model)
}

gsom.smooth <- function(...){
  #Reduce learning rate
  #fix starting neighbourhood
  #Find Winner
  #Adapt Weights
}

# Removes nodes, that were never marked as winners. These nodes are irrelevant as no
# Data is best represented by these nodes
# Unfortunatly the paper is not clear on when to remove empty nodes. It is assumed for this
# implementation, that the function shall be called after every iteration through the training df.
# Requires and returns gsom_models
gsom.emptyremove <- function(gsom_model){
  rmlist <- which(gsom_model$nodes$freq == 0)
  if(length(rmlist)!=0) {
    gsom_model$nodes$position <- gsom_model$nodes$position[-rmlist,]
    gsom_model$nodes$weight <- gsom_model$nodes$weight[-rmlist,]
    gsom_model$nodes$error <- gsom_model$nodes$error[-rmlist]
    gsom_model$nodes$freq <- gsom_model$nodes$freq[-rmlist]
  }
  return(gsom_model)
}

# Creates new nodes and initializes them.
# Requires gsom_model plus the coordinates where to grow the new node.
# Returns gsom_model. 
gsom.newnode <- function(gsom_model, x, y){
  #grow node
  gsom_model$nodes$position <- rbind(gsom_model$nodes$position, c(x, y))
  gsom_model$nodes$error <- c(gsom_model$nodes$error, 0)
  gsom_model$nodes$freq <- c(gsom_model$nodes$freq, 0)
  
  #Get Topology
  self <- which(gsom_model$nodes$position$x == x & gsom_model$nodes$position$y == y)
  top <- which(gsom_model$nodes$position$x == x & gsom_model$nodes$position$y == y+1)
  bottom <- which(gsom_model$nodes$position$x == x & gsom_model$nodes$position$y == y-1)
  left <- which(gsom_model$nodes$position$x == x+1 & gsom_model$nodes$position$y == y)
  right <- which(gsom_model$nodes$position$x == x-1 & gsom_model$nodes$position$y == y)
  neighbours <- c(top, bottom, left, right)
  
  #Weight calculation
  if(length(neighbours)>1){
    #Case B
    nweight <- (gsom_model$nodes$weight[neighbours[1],] + gsom_model$nodes$weight[neighbours[2],]) / 2
  }else{
    if(length(neighbours)==0)
      stop("No neighbours")
    posnx <- gsom_model$nodes$position[neighbours[1],'x']
    posny <- gsom_model$nodes$position[neighbours[1],'y']
    s2 <- which(gsom_model$nodes$position$x == posnx & gsom_model$nodes$position$y == posny)
    t2 <- which(gsom_model$nodes$position$x == posnx & gsom_model$nodes$position$y == posny+1)
    b2 <- which(gsom_model$nodes$position$x == posnx & gsom_model$nodes$position$y == posny-1)
    l2 <- which(gsom_model$nodes$position$x == posnx+1 & gsom_model$nodes$position$y == posny)
    r2 <- which(gsom_model$nodes$position$x == posnx-1 & gsom_model$nodes$position$y == posny)
    n2 <- c(b2, t2, r2, l2)
    if(which(n2==self))
      n2 <- n2[-which(n2==self)]
    ca <- which(gsom_model$nodes$position$x == ((posnx - x) + posnx) & gsom_model$nodes$position$y == ((posny - y) + posny))
    if(length(ca)!=0){
      #Case A
      nweight <- rep(-1, ncol(gsom_model$nodes$weight))
      for(l in 1:ncol(gsom_model$nodes$weight)){
        if(gsom_model$nodes$weight[ca,l]>gsom_model$nodes$weight[s2,l])
          nweight[l] <- gsom_model$nodes$weight[s2,l]-(gsom_model$nodes$weight[ca,l] - gsom_model$nodes$weight[s2,l])
        else
          nweight[l] <- gsom_model$nodes$weight[s2,l]+(gsom_model$nodes$weight[s2,l] - gsom_model$nodes$weight[ca,l])
      }
    } else if(length(n2) == 0){
      #Case D
      t1 <- apply(gsom_model$nodes$weight, 2, min)
      t2 <- apply(gsom_model$nodes$weight, 2, max)
      nweight <- (t1+t2)/2
    } else {
      #Case C
      nweight <- rep(-1, ncol(gsom_model$nodes$weight))
      for(l in 1:ncol(gsom_model$nodes$weight)){
        if(gsom_model$nodes$weight[n2[1],l]>gsom_model$nodes$weight[s2,l])
          nweight[l] <- gsom_model$nodes$weight[s2,l]-(gsom_model$nodes$weight[n2[1],l] - gsom_model$nodes$weight[s2,l])
        else
          nweight[l] <- gsom_model$nodes$weight[s2,l]+(gsom_model$nodes$weight[s2,l] - gsom_model$nodes$weight[n2[1],l])
      }
    }
  }
  
  #Enter new weights
  gsom_model$nodes$weight <- rbind(gsom_model$nodes$weight, nweight)
  
  return(gsom_model)
}