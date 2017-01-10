########################################################
# This file originates from the package 'fields'
# and has been slightly modified
########################################################


# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
plot_scale <- function(..., add = FALSE,
                         breaks= NULL, nlevel = 64, col = NULL,  
                         horizontal = FALSE, legend.shrink = 0.9, legend.width = 1.2, 
                         legend.lab = NULL, legend.line= 2,                    
                         graphics.reset = FALSE, bigplot = NULL, smallplot = NULL, lab.breaks = NULL,
                         axis.args = NULL, legend.args = NULL, legend.cex=1.0, midpoint = FALSE) {
  # Thanks to S. Koehler and  S. Woodhead
  # for comments on making this a better function
  #
  # save current graphics settings
  old.par <- par(no.readonly = TRUE)
  # set defaults for color scale 
  # note this works differently than the image function.
  nlevel<- length( col)
  #  figure out zlim from passed arguments
  #  also set the breaks for colors if they have not been passed, 
  info <- imagePlotInfo(..., breaks=breaks, nlevel=nlevel)
  # breaks have been computed if not passed in the call
  breaks<- info$breaks
  legend.mar <- 5.1
  # figure out how to divide up the plotting real estate 
  temp <- imageplot.setup(add = add, legend.shrink = legend.shrink, 
                          legend.width = legend.width, legend.mar = legend.mar, 
                          horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
  # bigplot has plotting region coordinates for image
  # smallplot has plotting coordinates for legend strip
  smallplot <- temp$smallplot
  bigplot <- temp$bigplot
  ##
  ## check dimensions of smallplot
  if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
    par(old.par)
    stop("plot region too small to add legend\n")
  }
  # Following code draws the legend using the image function
  # and a one column image.
  # What might be confusing is the values of the "image" are the same 
  # as the locations on the legend axis.
  # Moreover the image values are in the middle of each breakpoint category
  # thanks to Tobias Nanu Frechen and Matthew Flickinger 
  # for sorting out some problems with the breaks position in the legend.
  ix <- 1:2
  iy<- breaks
  nBreaks<- length( breaks)
  midpoints<- (breaks[1:(nBreaks-1)] +  breaks[2:nBreaks] )/2
  iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints)) 
  
  # next par call sets up a new plotting region just for the legend strip
  # at the smallplot coordinates
  par(new = TRUE, pty = "m", plt = smallplot, err = -1)
  # draw color scales the two  cases are horizontal/vertical 
  # add a label if this is passed.
  if (!horizontal) {
    image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
          ylab = "", col = col, breaks=breaks)
  }
  else {
    image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
          ylab = "", col = col, breaks=breaks)
  }
  # create the argument list to draw the axis
  #  this avoids 4 separate calls to axis and allows passing extra
  # arguments.
  if (!is.null(lab.breaks)) {
    # axis with labels at break points
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
                        at = breaks, labels = lab.breaks), axis.args)
  }
  else {
    # If lab.breaks is not specified ( with or without breaks), pretty
    # tick mark locations and labels are computed internally,
    # or as specified in axis.args at the function call
    axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                        mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
                   axis.args)
  }
  #
  # now add the axis to the legend strip.
  # notice how all the information is in the list axis.args
  do.call("axis", axis.args)
  # add a box around legend strip
  box()
  #
  # add a label to the axis if information has been  supplied
  # using the mtext function. The arguments to mtext are
  # passed as a list like the drill for axis (see above)
  #
  if (!is.null(legend.lab)) {
    legend.args <- list(text = legend.lab, side = ifelse(horizontal, 
                                                         1, 4), line = legend.line, cex=legend.cex)
    #                    just guessing at a good default for line argument!
  }
  # add the label using mtext function
  if (!is.null(legend.args)) {
    do.call(mtext, legend.args)
  }
  #
  # clean up graphics device settings
  # reset to larger plot region with right user coordinates.
  mfg.save <- par()$mfg
  par(old.par)
  par(mfg = mfg.save, new = FALSE)
  invisible()
}



"imagePlotInfo" <- function(..., breaks = NULL, nlevel) {
  #NOTE:
  # image.plot.info 
  # has been renamed as imagePlotInfo to avoid confusion with
  # an S3 method
  temp <- list(...)
  #
  xlim <- NA
  ylim <- NA
  zlim <- NA
  poly.grid <- FALSE
  #
  # go through various cases of what these can be
  #
  ##### x,y,z list is first argument
  if (is.list(temp[[1]])) {
    xlim <- range(temp[[1]]$x, na.rm = TRUE)
    ylim <- range(temp[[1]]$y, na.rm = TRUE)
    zlim <- range(temp[[1]]$z, na.rm = TRUE)
    if (is.matrix(temp[[1]]$x) & is.matrix(temp[[1]]$y) & 
        is.matrix(temp[[1]]$z)) {
      poly.grid <- TRUE
    }
  }
  ##### check for polygrid first three arguments should be matrices
  #####
  if (length(temp) >= 3) {
    if (is.matrix(temp[[1]]) & is.matrix(temp[[2]]) & is.matrix(temp[[3]])) {
      poly.grid <- TRUE
    }
  }
  #####  z is passed without an  x and y  (and not a poly.grid!)
  #####
  if (is.matrix(temp[[1]]) & !poly.grid) {
    xlim <- c(0, 1)
    ylim <- c(0, 1)
    zlim <- range(temp[[1]], na.rm = TRUE)
  }
  ##### if x,y,z have all been passed find their ranges.
  ##### holds if poly.grid or not
  #####
  if (length(temp) >= 3) {
    if (is.matrix(temp[[3]])) {
      xlim <- range(temp[[1]], na.rm = TRUE)
      ylim <- range(temp[[2]], na.rm = TRUE)
      zlim <- range(temp[[3]], na.rm = TRUE)           
    }
  }
  # if constant z values perturb the range (1e-8) by epsilon to 
  # avoid other problems in drawing legend later on
  if( !is.na( zlim[1] ) ){
    if( zlim[1] == zlim[2]){
      if( zlim[1]==0){
        zlim[1]<- -1e-8
        zlim[2]<- 1e-8}
      else{		 
        delta<- .01*abs(zlim[1])
        zlim[1]<- zlim[1] - delta
        zlim[2]<- zlim[2] + delta
      }
    }
  }
  #### parse x,y,z if they are  named arguments
  # determine if  this is polygon grid (x and y are matrices)
  if (is.matrix(temp$x) & is.matrix(temp$y) & is.matrix(temp$z)) {
    poly.grid <- TRUE
  }
  # set limits from the usual $x $y $z format of image object    
  xthere <- match("x", names(temp))
  ythere <- match("y", names(temp))
  zthere <- match("z", names(temp))
  if (!is.na(zthere)) 
    zlim <- range(temp$z, na.rm = TRUE)
  if (!is.na(xthere)) 
    xlim <- range(temp$x, na.rm = TRUE)
  if (!is.na(ythere)) 
    ylim <- range(temp$y, na.rm = TRUE)
  
  # overwrite limits with passed values
  if (!is.null(temp$zlim)) 
    zlim <- temp$zlim
  if (!is.null(temp$xlim)) 
    xlim <- temp$xlim
  if (!is.null(temp$ylim)) 
    ylim <- temp$ylim
  # At this point xlim, ylim and zlim should be correct 
  # using all the different possibilities and defaults for these values
  #        
  #  Now set up the breaks
  if( is.null(breaks)){
    midpoints<- seq( zlim[1], zlim[2],,nlevel)
    delta<- (midpoints[2]- midpoints[1])/2
    # nlevel +1 breaks with the min and max as midpoints 
    # of the first and last bins.
    
    breaks <- c( midpoints[1]- delta, midpoints + delta)
  }        
  list(xlim = xlim, ylim = ylim, zlim = zlim, poly.grid = poly.grid,
       breaks=breaks)
}

# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
# NOTE:
# image.plot.plt<- function(...){
# this function has been renamed as imageplot.setup to avoid confusion with
# an S3 method
#   imageplot.setup(...)}

"imageplot.setup" <- function(x, add = FALSE, legend.shrink = 0.9, 
                              legend.width = 1, horizontal = FALSE, legend.mar = NULL, 
                              bigplot = NULL, smallplot = NULL, ...) {
  old.par <- par(no.readonly = TRUE)
  if (is.null(smallplot)) 
    stick <- TRUE
  else stick <- FALSE
  if (is.null(legend.mar)) {
    legend.mar <- ifelse(horizontal, 3.1, 5.1)
  }
  # compute how big a text character is
  char.size <- ifelse(horizontal, par()$cin[2]/par()$din[2], 
                      par()$cin[1]/par()$din[1])
  # This is how much space to work with based on setting the margins in the
  # high level par command to leave between strip and big plot
  offset <- char.size * ifelse(horizontal, par()$mar[1], par()$mar[4])
  # this is the width of the legned strip itself.
  legend.width <- char.size * legend.width
  # this is room for legend axis labels
  legend.mar <- legend.mar * char.size
  # smallplot is the plotting region for the legend.
  if (is.null(smallplot)) {
    smallplot <- old.par$plt
    if (horizontal) {
      smallplot[3] <- legend.mar
      smallplot[4] <- legend.width + smallplot[3]
      pr <- (smallplot[2] - smallplot[1]) * ((1 - legend.shrink)/2)
      smallplot[1] <- smallplot[1] + pr
      smallplot[2] <- smallplot[2] - pr
    }
    else {
      smallplot[2] <- 1 - legend.mar
      smallplot[1] <- smallplot[2] - legend.width
      pr <- (smallplot[4] - smallplot[3]) * ((1 - legend.shrink)/2)
      smallplot[4] <- smallplot[4] - pr
      smallplot[3] <- smallplot[3] + pr
    }
  }
  if (is.null(bigplot)) {
    bigplot <- old.par$plt
    if (!horizontal) {
      bigplot[2] <- min(bigplot[2], smallplot[1] - offset)
    }
    else {
      bottom.space <- old.par$mar[1] * char.size
      bigplot[3] <- smallplot[4] + offset
    }
  }
  if (stick & (!horizontal)) {
    dp <- smallplot[2] - smallplot[1]
    smallplot[1] <- min(bigplot[2] + offset, smallplot[1])
    smallplot[2] <- smallplot[1] + dp
  }
  return(list(smallplot = smallplot, bigplot = bigplot))
}
