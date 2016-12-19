#######################################
#GSOM - Growing Self Organizing Maps
#map.r
#11/10/16 - Alex Hunziker
#######################################

# This function will plot the standard gsom plots, namely:
#   -Count Frequency of nodes (count)
#   -Distance within a node (distance)
#   -Distance to neighbouring nodes (dist_neighbours)
#   -Learning Process per Iteration (Average Mean Error) (training)
#   -Plot of Properties (property)
# Requires: a trained gsom model or mapped gsom data. Furthermore it requires
#   the type of plot desired. For the property plot the dimension(s) to be 
#   plotted can be indicated.
# Returns: nothing

plot.gsom <- function(gsom_object, type="count", colors=NULL, dim=0, main="", ...){
  
  if(!exists("gsom_object")) stop("GSOM object (trained model or mapped data) has to be provided.")
  
  if(type == "count"){
    
    if(is.null(colors)) colors = list(c(0.9,0),c(0.9,0),c(0.9,1))
    
    hlim = max(max(gsom_object$nodes$position[, 1])-min(gsom_object$nodes$position[, 1]), 
               max(gsom_object$nodes$position[, 2])-min(gsom_object$nodes$position[, 2]))/2 + 0.5
    hx = (max(gsom_object$nodes$position[, 1])+min(gsom_object$nodes$position[, 1]))/2
    hy = (max(gsom_object$nodes$position[, 2])+min(gsom_object$nodes$position[, 2]))/2
    par(mar=c(5.1,4.1,4.1,5.5))
    plot(gsom_object$nodes$position$x, 
         gsom_object$nodes$position$y, 
         type="n", main=paste("Observations per node"), xlab="", ylab="", xaxt='n', yaxt='n',
         xlim = c(hx-hlim, hx+hlim), ylim=c(hy-hlim, hy+hlim), pch=16, cex=3, ...)
    symbols(gsom_object$nodes$position[, 1], gsom_object$nodes$position[, 2],
            circles = rep(0.4, nrow(gsom_object$nodes$position)), inches = FALSE,
            add = TRUE, bg=plotrix::color.scale(gsom_object$nodes$freq,colors[[1]],colors[[2]],colors[[3]]))
    plot_scale(zlim=c(min(gsom_object$nodes$freq),max(gsom_object$nodes$freq)),
               col=plotrix::color.scale(min(gsom_object$nodes$freq):max(gsom_object$nodes$freq),colors[[1]],colors[[2]],colors[[3]]))
    par(mar=c(5.1,4.1,4.1,2.1))
    
  }else if(type == "distance"){
    
    if(is.null(colors)) colors = list(c(0.9,0),c(0.9,0),c(0.9,1))
    
    hlim = max(max(gsom_object$nodes$position[, 1])-min(gsom_object$nodes$position[, 1]), 
               max(gsom_object$nodes$position[, 2])-min(gsom_object$nodes$position[, 2]))/2 + 0.5
    hx = (max(gsom_object$nodes$position[, 1])+min(gsom_object$nodes$position[, 1]))/2
    hy = (max(gsom_object$nodes$position[, 2])+min(gsom_object$nodes$position[, 2]))/2
    par(mar=c(5.1,4.1,4.1,5.5))
    plot(gsom_object$nodes$position$x, 
         gsom_object$nodes$position$y, 
         type="n", main=paste("Avgerage Distance From BMN"), xlab="", ylab="", xaxt='n', yaxt='n',
         xlim = c(hx-hlim, hx+hlim), ylim=c(hy-hlim, hy+hlim), pch=16, cex=3, 
         col=plotrix::color.scale(gsom_object$nodes$distance,colors[[1]],colors[[2]],colors[[3]]), ...
    )
    symbols(gsom_object$nodes$position[, 1], gsom_object$nodes$position[, 2],
            circles = rep(0.4, nrow(gsom_object$nodes$position)), inches = FALSE,
            add = TRUE, bg=plotrix::color.scale(gsom_object$nodes$distance,colors[[1]],colors[[2]],colors[[3]]))
    minattr <- min(gsom_object$nodes$distance)
    maxattr <- max(gsom_object$nodes$distance)
    scale <- seq(minattr, maxattr, by=(maxattr-minattr)/100)
    plot_scale(zlim=c(min(gsom_object$nodes$distance),max(gsom_object$nodes$distance)),
               col=plotrix::color.scale(scale,colors[[1]],colors[[2]],colors[[3]]))
    par(mar=c(5.1,4.1,4.1,2.1))
    
  } else if(type == "dist_neighbours") {
    
    stop("Plot cannot be generated. Reason: Missing Feature.")
    
  } else if(type == "training") {
    
    if(is.null(gsom_object[["training"]])) stop("Trained gsom model expected, but obtained different data structure.")
    
    if(main == "") main <- "Training Progress"
    
    plot(x=gsom_object$training$iteration[gsom_object$training$training_stage==1], 
         y=gsom_object$training$meandist[gsom_object$training$training_stage==1], col=2, type="l",
         main=main, xlab="Number of iterations", ylab="Mean Distance to Unit", xlim = c(0, length(gsom_object$training$iteration)),
         ylim = c(min(gsom_object$training$meandist), max(gsom_object$training$meandist)), ...)
    points(x=gsom_object$training$iteration[gsom_object$training$training_stage==2], 
         y=gsom_object$training$meandist[gsom_object$training$training_stage==2], col=3, type="l")
    legend(length(gsom_object$training$meandist)*0.65, max(gsom_object$training$meandist), 
           legend=c("Growing Phase", "Smoothing Phase"), col=c(2, 3), lty=1, cex=0.8, lwd=2)
    
  } else if(type == "property") {
    
    if(is.null(colors)) colors = list(c(0.15,0.95,0.7),c(0.4,0.95,0.1),c(0.65,0.95,0.15))
    
    par(mar=c(5.1,4.1,4.1,5.5))
    if(any(dim > ncol(gsom_object$nodes$codes))) stop("Invalid value for parameter dim.")
    if(dim == 0) dim <- c(1:ncol(gsom_object$nodes$codes))
    
    if(main == "") gennames = TRUE
    for(i in dim){ #should eventually be changed to colnames. For works there as well
      
      if(exists("gennames")) main <- paste("Property:", colnames(gsom_object$nodes$codes)[i])
      
      minattr <- gsom_object$norm_param[i,1]
      maxattr <- gsom_object$norm_param[i,2]
      scale <- seq(minattr, maxattr, by=(maxattr-minattr)/100)
      
      hlim = max(max(gsom_object$nodes$position[, 1])-min(gsom_object$nodes$position[, 1]), 
                 max(gsom_object$nodes$position[, 2])-min(gsom_object$nodes$position[, 2]))/2 + 0.5
      hx = (max(gsom_object$nodes$position[, 1])+min(gsom_object$nodes$position[, 1]))/2
      hy = (max(gsom_object$nodes$position[, 2])+min(gsom_object$nodes$position[, 2]))/2
      
      par(mar=c(5.1,4.1,4.1,5.5))
      
      plot(gsom_object$nodes$position$x, 
           gsom_object$nodes$position$y, 
           type="n", main=main, xlab="", ylab="", xaxt='n', yaxt='n',
           xlim = c(hx-hlim, hx+hlim), ylim=c(hy-hlim, hy+hlim), pch=16, cex=3, ...)
      symbols(gsom_object$nodes$position[, 1], gsom_object$nodes$position[, 2],
              circles = rep(0.4, nrow(gsom_object$nodes$position)), inches = FALSE,
              add = TRUE, bg=plotrix::color.scale(gsom_object$nodes$codes[,i],colors[[1]],colors[[2]],colors[[3]]))
      plot_scale(zlim=c(minattr,maxattr),
                 col=plotrix::color.scale(scale,colors[[1]],colors[[2]],colors[[3]]))
      
    }
    
    par(mar=c(5.1,4.1,4.1,2.1))
    
  } else if(type == "predict") {
    
    if(is.null(colors)) colors = list(c(0.5,1), c(1,0.5), c(0.4,0.8))
    
    par(mar=c(5.1,4.1,4.1,5.5))
    if(any(dim > ncol(gsom_object$nodes$predict))) stop("Invalid value for parameter dim.")
    if(dim == 0) dim <- c(1:ncol(gsom_object$nodes$predict))
    
    if(main == "") gennames = TRUE
    for(i in dim){ #should eventually be changed to colnames. For works there as well
      
      if(exists("gennames")) main <- paste("Prediction:", colnames(gsom_object$nodes$predict[i]))
      
      minattr <- gsom_object$norm_param_y[i,1]
      maxattr <- gsom_object$norm_param_y[i,2]
      scale <- seq(minattr, maxattr, by=(maxattr-minattr)/100)
      
      hlim = max(max(gsom_object$nodes$position[, 1])-min(gsom_object$nodes$position[, 1]), 
                 max(gsom_object$nodes$position[, 2])-min(gsom_object$nodes$position[, 2]))/2 + 0.5
      hx = (max(gsom_object$nodes$position[, 1])+min(gsom_object$nodes$position[, 1]))/2
      hy = (max(gsom_object$nodes$position[, 2])+min(gsom_object$nodes$position[, 2]))/2
      
      par(mar=c(5.1,4.1,4.1,5.5))
      
      plot(gsom_object$nodes$position$x, 
           gsom_object$nodes$position$y, 
           type="n", main=main, xlab="", ylab="", xaxt='n', yaxt='n',
           xlim = c(hx-hlim, hx+hlim), ylim=c(hy-hlim, hy+hlim), pch=16, cex=3, ...)
      symbols(gsom_object$nodes$position[, 1], gsom_object$nodes$position[, 2],
              circles = rep(0.4, nrow(gsom_object$nodes$position)), inches = FALSE,
              add = TRUE, bg=plotrix::color.scale(gsom_object$nodes$predict[,i],colors[[1]],colors[[2]],colors[[3]]))
      plot_scale(zlim=c(minattr,maxattr),
                 col=plotrix::color.scale(scale,colors[[1]],colors[[2]],colors[[3]]))
      
    }
    
    par(mar=c(5.1,4.1,4.1,2.1))
    
  } else {
    
    stop("Invalid value for parameter type.")
    
  }
  
}