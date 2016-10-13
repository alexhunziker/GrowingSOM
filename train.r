#GSOM - Growing Self Organizing Maps
#train.r
#11/10/16 - Alex Hunziker
#######################################

# The Functions in this File are required in order to train the gsom model.
# gsom.train() is the main function, which should be called by the user.

gsom.train <- function(data, spreadFactor=0.5, keepdata=FALSE, iterations=50, alpha, ...){
  
  # Normalizing the training or testdata (min/max) in order to balance the impact
  # of the different properties of the dataframe
  min <- apply(data, 2, function(x){min(x)})
  max <- apply(data, 2, function(x){max(x)})
  df <- t(apply(data, 1, function(x){(x-min)/(max-min)}))
  
  gsom_model <- gsom.init(df, spreadFactor)
  gsom_model <- gsom.grow(gsom_model, df, iterations)
  
  norm_param <- data.frame(min = min, max = max)
  gsom_model[["norm_param"]] <- norm_param
  
  if(keepdata==TRUE){
    gsom_model[["data"]] = data
  }
  
  return(gsom_model)
  
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
  gsom_net <- list(nodes = nodes, training = training, herr = herr, GT = GT, norm_param=0)
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
    learningRate <- 0.8
    
    gsom_model$nodes$freq <- rep(0, times=length(gsom_model$nodes$freq))
    #gsom_model$nodes$error <- rep(0, times=length(gsom_model$nodes$error))
    t1 <- Sys.time()
    for(j in 1:nrow(df)){
      #Set Learning Rate
      #Recalculate Learning Rate
      #Contradiction to Paper. Learning rate is as defined in kohonen package
      #learningRate <- (total_iterations-(i*j))/total_iterations * 1-(3.8/nrow(gsom_model$nodes$weight)) * alpha
      learningRate <- alpha * (1-(3.8/nrow(gsom_model$nodes$weight))) * learningRate
      
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
      adjust <- gsom.get_neighours(gsom_model, winner.x, winner.y)
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
          gsom_model$nodes$error[adjust["self"]] = gsom_model$GT / 2
          
          #Note: This is ugly, since growth condidion is only checked, once the following nodes are winners again.
          #Also we are checking the topology once again...
          adjust <- gsom.get_neighours(gsom_model, winner.x, winner.y)
          
          for(m in 1:4){
            #The paper suggests values for gamma between 0 and 1
            gam <- 0.5
            gsom_model$nodes$error[adjust[m]] = gsom_model$nodes$error[adjust[m]] * (1+gam)
          }
          
        } else {
          
          #Boundry node -> Grow net.
          nodegrow = nodegrow + 1
          if(is.na(adjust["top"])) gsom_model <- gsom.newnode(gsom_model, winner.x, winner.y+1)
          if(is.na(adjust["bottom"])) gsom_model <- gsom.newnode(gsom_model, winner.x, winner.y-1)
          if(is.na(adjust["right"])) gsom_model <- gsom.newnode(gsom_model, winner.x+1, winner.y)
          if(is.na(adjust["left"])) gsom_model <- gsom.newnode(gsom_model, winner.x-1, winner.y)
          
          #Not clearly mentioned in the paper
          gsom_model$nodes$error[adjust["self"]]=0
          learningRate <- 0.8
        
        }
        
      }

    }
    
    print(errorsum)
    meandist <- errorsum/nrow(df)
    curr_train = c(iteration=i, training_stage=1, meandist=meandist, nodecount=nrow(gsom_model$nodes$position), nodegrow=nodegrow)
    print(curr_train)
    gsom_model$training[nrow(gsom_model$training)+1,] <- curr_train
    gsom_model <- gsom.emptyremove(gsom_model)
    t2 <- Sys.time()
    
    #Arbitrary!
    if(gsom_model$training$nodecount[i] <= gsom_model$training$nodecount[i-4] && i > 4) break
  }
  
  gsom.plot(gsom_model, type="property")
  gsom_model <- gsom.smooth(gsom_model, df, rep, i, alpha)
  
  return(gsom_model)
}

gsom.smooth <- function(gsom_model, df, rep, n, alpha){
  
  winner <- numeric()
  total_iterations <- nrow(df)*rep
  nodegrow <- 0

  #Present Input
  #Replace with c code for speed gains
  for(i in 1:rep){
    
    errorsum <- 0
    gsom_model$nodes$freq <- rep(0, times=length(gsom_model$nodes$freq))
    gsom_model$nodes$error <- rep(0, times=length(gsom_model$nodes$error))
    
    t1 <- Sys.time()
    
    for(j in 1:nrow(df)){
      #Set Learning Rate
      #Recalculate Learning Rate
      #Contradiction to Paper. Learning rate is as defined in kohonen package
      learningRate <- (total_iterations-(i*j)) / total_iterations * alpha / 3
      
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
      adjust <- gsom.get_neighours(gsom_model, winner.x, winner.y)
      # for(k in adjust){
      #  gsom_model$nodes$weight[k,] <- gsom_model$nodes$weight[k,]+ (df[j,] - gsom_model$nodes$weight[k,]) * learningRate
      # }
      #tempstor <- gsom_model$nodes$weight[adjust["self"],]
      gsom_model$nodes$weight[adjust["self"],] <- gsom_model$nodes$weight[adjust["self"],]+ 
        (df[j,] - gsom_model$nodes$weight[adjust["self"],]) * learningRate
      #temp2 <- gsom_model$nodes$weight[adjust["self"],]
      #adjust = adjust[2:length(adjust)]
      #adjust <- na.omit(adjust)
      #for(k in adjust){
      #  gsom_model$nodes$weight[k,] <- gsom_model$nodes$weight[k,] + (df[j,] - gsom_model$nodes$weight[k,]) * learningRate * 0.5
      #}
      
    }
    
    print(errorsum)
    meandist <- errorsum/nrow(df)
    curr_train = c(iteration=i+n, training_stage=2, meandist=meandist, nodecount=nrow(gsom_model$nodes$position), nodegrow=nodegrow)
    print(curr_train)
    gsom_model$training[nrow(gsom_model$training)+1,] <- curr_train
    t2 <- Sys.time()
    print(t2-t1)
  }
  
  gsom_model <- gsom.emptyremove(gsom_model)
  gsom.plot(gsom_model, type="property")

  return(gsom_model)
}

# Gets all the direct neighbours of a point with coord(x,y) and returns a vector with the rownumber of
# the point and its neighbours
gsom.get_neighours <- function(gsom_model, x, y){
  self <- which(gsom_model$nodes$position$x == x & gsom_model$nodes$position$y == y)
  top <- which(gsom_model$nodes$position$x == x & gsom_model$nodes$position$y == y+1)
  bottom <- which(gsom_model$nodes$position$x == x & gsom_model$nodes$position$y == y-1)
  right <- which(gsom_model$nodes$position$x == x+1 & gsom_model$nodes$position$y == y)
  left <- which(gsom_model$nodes$position$x == x-1 & gsom_model$nodes$position$y == y)
  neighbours <- c(self=self, top=top, bottom=bottom, left=left, right=right)
  return(neighbours)
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
  neighbours <- gsom.get_neighours(gsom_model, x, y)
  self <- neighbours[1]
  neighbours <- na.omit(neighbours[2:5])
  
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

# Function for debugging only. 
# Can be used to plot the gsom frequency at any stage druring the training process.
# All collumns of the dataframe are expected to contain numeric values.
print_crude <- function(gsom_model_t){
  plot(gsom_model_t$nodes$position$x, 
       gsom_model_t$nodes$position$y, 
       type="n", main="Training Plot (Frequency)", xlab="", ylab="", xaxt='n', yaxt='n',
       pch=16, cex=3, col=gray(gsom_model_t$nodes$freq/max(gsom_model_t$nodes$freq)))
  text(gsom_model_t$nodes$position$x, 
       gsom_model_t$nodes$position$y, 
       gsom_model_t$nodes$freq)
}