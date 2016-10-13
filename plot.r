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
  
  if(!exists("gsom_object")) stop("GSOM object (trained model or mapped data) has to be provided.")
  
  if(type == "count"){
    
    plot(gsom_object$nodes$position$x, 
         gsom_object$nodes$position$y, 
         type="p", main=paste("Density"), xlab="", ylab="", xaxt='n', yaxt='n',
         pch=16, cex=3, col=color.scale(gsom_object$nodes$freq/max(gsom_object$nodes$freq),c(0,1,1),c(1,1,0),0)
    )
    
  }else if(type == "dist"){
    
    plot(gsom_object$nodes$position$x, 
         gsom_object$nodes$position$y, 
         type="p", main=paste("Density"), xlab="", ylab="", xaxt='n', yaxt='n',
         pch=16, cex=3, col=grey((gsom_object$nodes$error/max(gsom_object$nodes$error))^2)
    )
    
  } else if(type == "dist_neighbours") {
    
    stop("Missing Feature. Sorry...")
    plot(gsom_object$nodes$position$x, 
         gsom_object$nodes$position$y, 
         type="p", main=paste("Density"), xlab="", ylab="", xaxt='n', yaxt='n',
         pch=16, cex=3, col=grey((gsom_object$nodes$error/max(gsom_object$nodes$error))^2)
    )
    
  } else if(type == "training") {
    
    if(is.null(gsom_object[["training"]])) stop("Trained gsom model expected, but obtained different data structure.")
    
    if(main == "") main <- "Training Progress"
    
    plot(x=gsom_object$training$iteration[gsom_object$training$training_stage==1], 
         y=gsom_object$training$meandist[gsom_object$training$training_stage==1], col=2, type="l",
         main=main, xlab="Number of iterations", ylab="Mean Distance to Unit", xlim = c(0, length(gsom_object$training$iteration)),
         ylim = c(min(gsom_object$training$meandist), max(gsom_object$training$meandist)))
    points(x=gsom_object$training$iteration[gsom_object$training$training_stage==2], 
         y=gsom_object$training$meandist[gsom_object$training$training_stage==2], col=3, type="l")
    
  } else if(type == "property") {
    
    if(any(dim > ncol(gsom_object$nodes$weight))) stop("Invalid value for parameter dim.")
    if(dim == 0) dim <- c(1:ncol(gsom_object$nodes$weight))
    
    if(main == "") gennames = TRUE
    for(i in dim){ #should eventually be changed to colnames. For works there as well
      if(exists("gennames")) main <- paste("Property Nr:", i)
      
      plot(gsom_object$nodes$position$x, 
           gsom_object$nodes$position$y, 
           type="p", main=main, xlab="", ylab="", xaxt='n', yaxt='n',
           pch=16, cex=3, col=color.scale((gsom_object$nodes$weight[,i]/max(gsom_object$nodes$weight[,i])),c(0.3,0.9),c(0,0.95),c(0.3,0.95))
      )
      
    }
    
  } else {
    
    stop("Invalid value for parameter type.")
    
  }
  
}