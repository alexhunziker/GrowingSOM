\name{plot.gsom}
\alias{plot.gsom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Plot gsom models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.gsom(gsom_object, type = "count", dim = 0, main = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gsom_object}{
%%     ~~Describe \code{gsom_object} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{dim}{
%%     ~~Describe \code{dim} here~~
}
  \item{main}{
%%     ~~Describe \code{main} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (gsom_object, type = "count", dim = 0, main = "") 
{
    if (!exists("gsom_object")) 
        stop("GSOM object (trained model or mapped data) has to be provided.")
    if (type == "count") {
        par(mar = c(5.1, 4.1, 4.1, 5.5))
        plot(gsom_object$nodes$position$x, gsom_object$nodes$position$y, 
            type = "n", main = paste("Density (Observations per unit)"), 
            xlab = "", ylab = "", xaxt = "n", yaxt = "n", pch = 16, 
            cex = 3, col = plotrix::color.scale(gsom_object$nodes$freq, 
                c(0.9, 0), c(0.9, 0), c(0.9, 1)))
        symbols(gsom_object$nodes$position[, 1], gsom_object$nodes$position[, 
            2], circles = rep(0.4, nrow(gsom_object$nodes$position)), 
            inches = FALSE, add = TRUE, bg = plotrix::color.scale(gsom_object$nodes$freq, 
                c(0.9, 0), c(0.9, 0), c(0.9, 1)))
        image.plot(legend.only = TRUE, zlim = c(min(gsom_object$nodes$freq), 
            max(gsom_object$nodes$freq)), col = plotrix::color.scale(min(gsom_object$nodes$freq):max(gsom_object$nodes$freq), 
            c(0.9, 0), c(0.9, 0), c(0.9, 1)))
        par(mar = c(5.1, 4.1, 4.1, 2.1))
    }
    else if (type == "dist") {
        par(mar = c(5.1, 4.1, 4.1, 5.5))
        plot(gsom_object$nodes$position$x, gsom_object$nodes$position$y, 
            type = "n", main = paste("Avgerage Euclidan Distance From BMU"), 
            xlab = "", ylab = "", xaxt = "n", yaxt = "n", pch = 16, 
            cex = 3, col = plotrix::color.scale(gsom_object$nodes$error, 
                c(0.9, 0), c(0.9, 0), c(0.9, 1)))
        symbols(gsom_object$nodes$position[, 1], gsom_object$nodes$position[, 
            2], circles = rep(0.4, nrow(gsom_object$nodes$position)), 
            inches = FALSE, add = TRUE, bg = plotrix::color.scale(gsom_object$nodes$error, 
                c(0.9, 0), c(0.9, 0), c(0.9, 1)))
        minattr <- min(gsom_object$nodes$error)
        maxattr <- max(gsom_object$nodes$error)
        scale <- seq(minattr, maxattr, by = (maxattr - minattr)/100)
        image.plot(legend.only = TRUE, zlim = c(min(gsom_object$nodes$error), 
            max(gsom_object$nodes$error)), col = plotrix::color.scale(scale, 
            c(0.9, 0), c(0.9, 0), c(0.9, 1)))
        par(mar = c(5.1, 4.1, 4.1, 2.1))
    }
    else if (type == "dist_neighbours") {
        stop("Missing Feature. Sorry...")
        plot(gsom_object$nodes$position$x, gsom_object$nodes$position$y, 
            type = "n", main = paste("Density"), xlab = "", ylab = "", 
            xaxt = "n", yaxt = "n", pch = 16, cex = 3, col = grey((gsom_object$nodes$error/max(gsom_object$nodes$error))^2))
        symbols(gsom_object$nodes$position[, 1], gsom_object$nodes$position[, 
            2], circles = rep(0.4, nrow(gsom_object$nodes$position)), 
            inches = FALSE, add = TRUE, bg = grey((gsom_object$nodes$error/max(gsom_object$nodes$error))^2))
    }
    else if (type == "training") {
        if (is.null(gsom_object[["training"]])) 
            stop("Trained gsom model expected, but obtained different data structure.")
        if (main == "") 
            main <- "Training Progress"
        plot(x = gsom_object$training$iteration[gsom_object$training$training_stage == 
            1], y = gsom_object$training$meandist[gsom_object$training$training_stage == 
            1], col = 2, type = "l", main = main, xlab = "Number of iterations", 
            ylab = "Mean Distance to Unit", xlim = c(0, length(gsom_object$training$iteration)), 
            ylim = c(min(gsom_object$training$meandist), max(gsom_object$training$meandist)))
        points(x = gsom_object$training$iteration[gsom_object$training$training_stage == 
            2], y = gsom_object$training$meandist[gsom_object$training$training_stage == 
            2], col = 3, type = "l")
        legend(length(gsom_object$training$meandist) * 0.65, 
            max(gsom_object$training$meandist), legend = c("Growing Phase", 
                "Smoothing Phase"), col = c(2, 3), lty = 1, cex = 0.8, 
            lwd = 2)
    }
    else if (type == "property") {
        par(mar = c(5.1, 4.1, 4.1, 5.5))
        if (any(dim > ncol(gsom_object$nodes$weight))) 
            stop("Invalid value for parameter dim.")
        if (dim == 0) 
            dim <- c(1:ncol(gsom_object$nodes$weight))
        if (main == "") 
            gennames = TRUE
        for (i in dim) {
            if (exists("gennames")) 
                main <- paste("Property Nr:", i)
            minattr <- gsom_object$norm_param[i, 1]
            maxattr <- gsom_object$norm_param[i, 2]
            scale <- seq(minattr, maxattr, by = (maxattr - minattr)/100)
            plot(gsom_object$nodes$position$x, gsom_object$nodes$position$y, 
                type = "n", main = main, xlab = "", ylab = "", 
                xaxt = "n", yaxt = "n", pch = 16, cex = 3, col = plotrix::color.scale(gsom_object$nodes$weight[, 
                  i], c(0.3, 0.9), c(0, 0.95), c(0.3, 0.95)))
            symbols(gsom_object$nodes$position[, 1], gsom_object$nodes$position[, 
                2], circles = rep(0.4, nrow(gsom_object$nodes$position)), 
                inches = FALSE, add = TRUE, bg = plotrix::color.scale(gsom_object$nodes$weight[, 
                  i], c(0.3, 0.9), c(0, 0.95), c(0.3, 0.95)))
            image.plot(legend.only = TRUE, zlim = c(minattr, 
                maxattr), col = plotrix::color.scale(scale, c(0.3, 
                0.9), c(0, 0.95), c(0.3, 0.95)))
        }
        par(mar = c(5.1, 4.1, 4.1, 2.1))
    }
    else {
        stop("Invalid value for parameter type.")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
