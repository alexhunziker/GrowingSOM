#######################################
#GSOM - Growing Self Organizing Maps
#map.r
#11/10/16 - Alex Hunziker
#######################################

# This function will plot the standard gsom plots, namely:
#   -Count Frequency of nodes (count)
#   -Distance within a node (dist)
#   -Distance to neighbouring nodes (dist_neighbours)
#   -Learning Process per Iteration (Average Mean Error) (training)
#   -Plot of Properties (property)
# Requires: a trained gsom model or mapped gsom data. Furthermore it requires
#   the type of plot desired. For the property plot the dimension(s) to be 
#   plotted can be indicated.
# Returns: nothing

gsom.plot <- function(gsom_object, type="count", dim=0, main=""){
  
  if(!exists(gsom_object)) stop("GSOM object (trained model or mapped data) has to be provided.")
  
  if(type == "count"){
    
    plot(gsom_object$nodes$position$x, 
         gsom_object$nodes$position$y, 
         type="p", main=paste("Density"), xlab="", ylab="", xaxt='n', yaxt='n',
         pch=16, cex=3, col=grey((gsom_object$nodes$freq/max(gsom_object$nodes$freq))^2)
    )
    
  }else if(type = "dist"){
    
    plot(gsom_object$nodes$position$x, 
         gsom_object$nodes$position$y, 
         type="p", main=paste("Density"), xlab="", ylab="", xaxt='n', yaxt='n',
         pch=16, cex=3, col=grey((gsom_object$nodes$error/max(gsom_object$nodes$error))^2)
    )
    
  } else if(type = "dist_neighbours") {
    
    stop("Missing Feature. Sorry...")
    plot(gsom_object$nodes$position$x, 
         gsom_object$nodes$position$y, 
         type="p", main=paste("Density"), xlab="", ylab="", xaxt='n', yaxt='n',
         pch=16, cex=3, col=grey((gsom_object$nodes$error/max(gsom_object$nodes$error))^2)
    )
    
  } else if(type = "training") {
    
    if(!exists(gsom_object$training)) stop("Trained gsom model expected, but obtained different data structure.")
    
    stop("Missing Feature")
    
  } else if(type = "property") {
    
    if(any(dim > ncol(gsom_object$nodes$weight))) stop("Invalid value for parameter dim.")
    if(dim == 0) dim <- c(1:ncol(gsom_object$nodes$weight))
    
    for(i in dim){ #should eventually be changed to colnames. For works there as well
      if(main == "") main <- paste("Property Nr:", i)
      
      plot(gsom_object$nodes$position$x, 
           gsom_object$nodes$position$y, 
           type="p", main=main, xlab="", ylab="", xaxt='n', yaxt='n',
           pch=16, cex=3, col=gray((gsom_object$nodes$weight[,i]/max(gsom_object$nodes$weight[,i]))^2)
      )
      
    }
    
  } else {
    
    stop("Invalid value for parameter type.")
    
  }
  
}