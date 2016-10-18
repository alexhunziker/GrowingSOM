#load("Q:/Abteilungsprojekte/eng/SWWData/Alex/Validation/SBR_raw.RData")
#df <- traindata
df <- iris[,1:4]

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

# Normalizing the training or testdata (min/max) in order to balance the impact
# of the different properties of the dataframe
gsom.normalize <- function(df){
  min <- apply(df, 2, function(x){min(x)})
  max <- apply(df, 2, function(x){max(x)})
  df <- t(apply(df, 1, function(x){(x-min)/(max-min)}))
  return(df)
}
#df <- gsom.normalize(df[1:800,3:10])
df <- gsom.normalize(df)

# Gnenerate the basic data structure for the GSOM model. 
# Creates for initial nodes that are initialized with random values.
# Returns trainable gsom model
gsom.init <- function(spreadFactor, df){
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
gsom_model <- gsom.init(0.3, df)

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

# Main loop during stage 1 (Growing phase) of the training process
# Loops through training set and assigns every event a winner and keeps track of error and initializes growth of the grid
# Requires gsom_model, df. Tuning Parameters: Max no of repetitions, alpha (Parameter determining the learning rate (between 0 and 1))
# Returns gsom_model (fro smoothing phase)
gsom.grow <- function(gsom_model, df, rep=50, alpha = 0.5, ...){
  winner <- "Dummy"
  total_iterations <- nrow(df)*rep
  gsom_model$nodes$error <- rep(0, times=length(gsom_model$nodes$error))
  
  #Present Input
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
#gsom_model_t = gsom.grow(gsom_model, df[1:800,], rep=100)
gsom_model_t = gsom.grow(gsom_model, df, rep=500)

gsom.smooth <- function(...){
  #Reduce learning rate
  #fix starting neighbourhood
  #Find Winner
  #Adapt Weights
}

gsom.geterrordist <- function(winner, newnode, ...){
  newerror = 0
  for(k in 1:ncol(winner$weight)){
    newerror = newerror + (newnode[,k]-winner$weight[,k])^2
  }
  newerror = sqrt(newerror)
  errordist = winner.error + newerror
  return(errordist)
}

gsom.updateherr <- function(errordist_i, herr){
  if(herr < errordist_i) return(errordist_i)
  else return(herr)
}

gsom.isborderneuron <- function(gsom_net){
  if(gsom_net[x]) TRUE
}

gsom.grownet <- function(...){
  if(gsom_net[gsom_net$x==(toextednd$x+1), gsom_net$y==gtoextend$y])
    rbind(gsom_net, c(toextend$x+1, toextend$y, NULL, 0))
  if(gsom_net[gsom_net$x==(toextednd$x-1), gsom_net$y==gtoextend$y])
    rbind(gsom_net, c(toextend$x-1, toextend$y, NULL, 0))
  if(gsom_net[gsom_net$x==toextednd$x, gsom_net$y==(gtoextend$y+1)])
    rbind(gsom_net, c(toextend$x, toextend$y+1, NULL, 0))
  if(gsom_net[gsom_net$x==toextednd$x, gsom_net$y==(gtoextend$y-1)])
    rbind(gsom_net, c(toextend$x, toextend$y-1, NULL, 0))
}

gsom.train <- function(df, sf, len=50){
  df <- gsom.normalize(df)
  model <- gsom.init(sf, df)
  gsom.grow(model, df, len)
}

gsom.propertyplot_crude <- function(gsom_model){
  plot(x=gsom_model_t$training$iteration, y=gsom_model_t$training$meandist, col=gsom_model_t$training$training_stage+1, type="l",
       main="Training Progress", xlab="Number of iterations", ylab="Mean Distance to Unit")
  plot(gsom_model_t$nodes$position$x, 
       gsom_model_t$nodes$position$y, 
       type="p", main=paste("Density"), xlab="", ylab="", xaxt='n', yaxt='n',
       pch=16, cex=3, col=grey((gsom_model_t$nodes$freq/max(gsom_model_t$nodes$freq))^2))
  for(n in 1:ncol(gsom_model_t$nodes$weight)){
    plot(gsom_model_t$nodes$position$x, 
       gsom_model_t$nodes$position$y, 
       type="p", main=paste("Property Nr:", n), xlab="", ylab="", xaxt='n', yaxt='n',
       pch=16, cex=3, col=gray((gsom_model_t$nodes$weight[,n]/max(gsom_model_t$nodes$weight[,n]))^2))
  }
}