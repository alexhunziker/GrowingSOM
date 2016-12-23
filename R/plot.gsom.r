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

plot.gsom <- function(x, type="count", colors=NULL, dim=0, main="", ...){
  
  if(!exists("x")) stop("GSOM object (trained model or mapped data) has to be provided.")
  
  if(type == "count"){
    
    if(is.null(colors)) colors = list(c(0.9,0),c(0.9,0),c(0.9,1))
    
    hlim = max(max(x$nodes$position[, 1])-min(x$nodes$position[, 1]), 
               max(x$nodes$position[, 2])-min(x$nodes$position[, 2]))/2 + 0.5
    hx = (max(x$nodes$position[, 1])+min(x$nodes$position[, 1]))/2
    hy = (max(x$nodes$position[, 2])+min(x$nodes$position[, 2]))/2
    par(mar=c(5.1,4.1,4.1,5.5))
    plot(x$nodes$position$x, 
         x$nodes$position$y, 
         type="n", main=paste("Observations per node"), xlab="", ylab="", xaxt='n', yaxt='n',
         xlim = c(hx-hlim, hx+hlim), ylim=c(hy-hlim, hy+hlim), pch=16, cex=3, ...)
    symbols(x$nodes$position[, 1], x$nodes$position[, 2],
            circles = rep(0.4, nrow(x$nodes$position)), inches = FALSE,
            add = TRUE, bg=plotrix::color.scale(x$nodes$freq,colors[[1]],colors[[2]],colors[[3]]))
    plot_scale(zlim=c(min(x$nodes$freq),max(x$nodes$freq)),
               col=plotrix::color.scale(min(x$nodes$freq):max(x$nodes$freq),colors[[1]],colors[[2]],colors[[3]]))
    par(mar=c(5.1,4.1,4.1,2.1))
    
  }else if(type == "distance"){
    
    if(is.null(colors)) colors = list(c(0.9,0),c(0.9,0),c(0.9,1))
    
    hlim = max(max(x$nodes$position[, 1])-min(x$nodes$position[, 1]), 
               max(x$nodes$position[, 2])-min(x$nodes$position[, 2]))/2 + 0.5
    hx = (max(x$nodes$position[, 1])+min(x$nodes$position[, 1]))/2
    hy = (max(x$nodes$position[, 2])+min(x$nodes$position[, 2]))/2
    par(mar=c(5.1,4.1,4.1,5.5))
    plot(x$nodes$position$x, 
         x$nodes$position$y, 
         type="n", main=paste("Avgerage Distance From BMN"), xlab="", ylab="", xaxt='n', yaxt='n',
         xlim = c(hx-hlim, hx+hlim), ylim=c(hy-hlim, hy+hlim), pch=16, cex=3, 
         col=plotrix::color.scale(x$nodes$distance,colors[[1]],colors[[2]],colors[[3]]), ...
    )
    symbols(x$nodes$position[, 1], x$nodes$position[, 2],
            circles = rep(0.4, nrow(x$nodes$position)), inches = FALSE,
            add = TRUE, bg=plotrix::color.scale(x$nodes$distance,colors[[1]],colors[[2]],colors[[3]]))
    minattr <- min(x$nodes$distance)
    maxattr <- max(x$nodes$distance)
    scale <- seq(minattr, maxattr, by=(maxattr-minattr)/100)
    plot_scale(zlim=c(min(x$nodes$distance),max(x$nodes$distance)),
               col=plotrix::color.scale(scale,colors[[1]],colors[[2]],colors[[3]]))
    par(mar=c(5.1,4.1,4.1,2.1))
    
  } else if(type == "dist_neighbours") {
    
    stop("Plot cannot be generated. Reason: Missing Feature.")
    
  } else if(type == "training") {
    
    if(is.null(x[["training"]])) stop("Trained gsom model expected, but obtained different data structure.")
    
    if(main == "") main <- "Training Progress"
    
    plot(x=x$training$iteration[x$training$training_stage==1], 
         y=x$training$meandist[x$training$training_stage==1], col=2, type="l",
         main=main, xlab="Number of iterations", ylab="Mean Distance to Unit", xlim = c(0, length(x$training$iteration)),
         ylim = c(min(x$training$meandist), max(x$training$meandist)), ...)
    points(x=x$training$iteration[x$training$training_stage==2], 
         y=x$training$meandist[x$training$training_stage==2], col=3, type="l")
    legend(length(x$training$meandist)*0.65, max(x$training$meandist), 
           legend=c("Growing Phase", "Smoothing Phase"), col=c(2, 3), lty=1, cex=0.8, lwd=2)
    
  } else if(type == "property") {
    
    if(is.null(colors)) colors = list(c(0.15,0.95,0.7),c(0.4,0.95,0.1),c(0.65,0.95,0.15))
    
    par(mar=c(5.1,4.1,4.1,5.5))
    if(any(dim > ncol(x$nodes$codes))) stop("Invalid value for parameter dim.")
    if(dim == 0) dim <- c(1:ncol(x$nodes$codes))
    
    if(main == "") gennames = TRUE
    for(i in dim){ #should eventually be changed to colnames. For works there as well
      
      if(exists("gennames")) main <- paste("Property:", colnames(x$nodes$codes)[i])
      
      minattr <- min(x$nodes$codes[,i])
      maxattr <- max(x$nodes$codes[,i])
      #minattr <- x$norm_param[i,1]
      #maxattr <- x$norm_param[i,2]
      scale <- seq(minattr, maxattr, by=(maxattr-minattr)/100)
      
      hlim = max(max(x$nodes$position[, 1])-min(x$nodes$position[, 1]), 
                 max(x$nodes$position[, 2])-min(x$nodes$position[, 2]))/2 + 0.5
      hx = (max(x$nodes$position[, 1])+min(x$nodes$position[, 1]))/2
      hy = (max(x$nodes$position[, 2])+min(x$nodes$position[, 2]))/2
      
      par(mar=c(5.1,4.1,4.1,5.5))
      
      plot(x$nodes$position$x, 
           x$nodes$position$y, 
           type="n", main=main, xlab="", ylab="", xaxt='n', yaxt='n',
           xlim = c(hx-hlim, hx+hlim), ylim=c(hy-hlim, hy+hlim), pch=16, cex=3, ...)
      symbols(x$nodes$position[, 1], x$nodes$position[, 2],
              circles = rep(0.4, nrow(x$nodes$position)), inches = FALSE,
              add = TRUE, bg=plotrix::color.scale(x$nodes$codes[,i],colors[[1]],colors[[2]],colors[[3]]))
      plot_scale(zlim=c(minattr,maxattr),
                 col=plotrix::color.scale(scale,colors[[1]],colors[[2]],colors[[3]]))
      
    }
    
    par(mar=c(5.1,4.1,4.1,2.1))
    
  } else if(type == "predict") {
    
    if(is.null(colors)) colors = list(c(0.5,1), c(1,0.5), c(0.4,0.8))
    
    par(mar=c(5.1,4.1,4.1,5.5))
    if(any(dim > ncol(x$nodes$predict))) stop("Invalid value for parameter dim.")
    if(dim == 0) dim <- c(1:ncol(x$nodes$predict))
    
    if(main == "") gennames = TRUE
    for(i in dim){ #should eventually be changed to colnames. For works there as well
      
      if(exists("gennames")) main <- paste("Prediction:", colnames(x$nodes$predict[i]))
      
      minattr <- min(x$nodes$codes[,i])
      maxattr <- max(x$nodes$codes[,i])
      scale <- seq(minattr, maxattr, by=(maxattr-minattr)/100)
      
      hlim = max(max(x$nodes$position[, 1])-min(x$nodes$position[, 1]), 
                 max(x$nodes$position[, 2])-min(x$nodes$position[, 2]))/2 + 0.5
      hx = (max(x$nodes$position[, 1])+min(x$nodes$position[, 1]))/2
      hy = (max(x$nodes$position[, 2])+min(x$nodes$position[, 2]))/2
      
      par(mar=c(5.1,4.1,4.1,5.5))
      
      plot(x$nodes$position$x, 
           x$nodes$position$y, 
           type="n", main=main, xlab="", ylab="", xaxt='n', yaxt='n',
           xlim = c(hx-hlim, hx+hlim), ylim=c(hy-hlim, hy+hlim), pch=16, cex=3, ...)
      symbols(x$nodes$position[, 1], x$nodes$position[, 2],
              circles = rep(0.4, nrow(x$nodes$position)), inches = FALSE,
              add = TRUE, bg=plotrix::color.scale(x$nodes$predict[,i],colors[[1]],colors[[2]],colors[[3]]))
      plot_scale(zlim=c(minattr,maxattr),
                 col=plotrix::color.scale(scale,colors[[1]],colors[[2]],colors[[3]]))
      
    }
    
    par(mar=c(5.1,4.1,4.1,2.1))
    
  } else {
    
    stop("Invalid value for parameter type.")
    
  }
  
}