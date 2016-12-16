\name{plot.gsom}
\alias{plot.gsom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Plot gsom objects
}
\description{
  Plot growing self-organising maps, obtained from functions of the GrowingSOM package. Several types of plots are supported.
}
\usage{
  plot.gsom(gsom_object, type = "count", colors=NULL, dim = 0, main = "")
}
\arguments{
  \item{gsom_object}{
    a gsom object
  }
  \item{type}{
    defines the type of the plot. Possible values are:\cr
    - training: shows the mean distance to the code of the best matching unit during the training.\cr
    - count: Graphical Representation of how many observations are put into each unit.\cr
    - property: properties of each unit can be calculated and shown in colour code. It can be used to visualise the similarity of one particular object to all units in the map, to show the mean similarity of all units and the objects mapped to them, etcetera\cr
    - distance: shows the mean distance of objects mapped to a unit to the codebook vector of that unit. The smaller the distances, the better the objects are represented by unit.\cr
    - predict: properties of the properties to be modelled can be shown similarly to the "property" plot.
  }
  \item{dim}{
    can be used if type is "property" to specifly which property should be plotted
  }
  \item{colors}{
    colors to be used to generate the plots. If set to NULL, the default colors will be used. Format: list(red, green, blue)\cr
    \code{Examples:\cr 
        colors = list(c(0.1, 0.9), c(0.1, 0.9), c(0.1, 0.9))   For gray scale\cr 
        colors = list(c(1.0, 0.0), c(0.0, 1.0), c(0.0, 0.0))   Red to green}
  }
  \item{main}{
    changes the title of the plots
  }
  \item{...}{
    Collects further graphical parameters for the plot
  }
}

\author{
  Alex Hunziker
}
%\note{
%  a plot for the distance of direct neighbours to the respective node should be introduced in the future.
%}

\examples{
  # Get some Data
  data(iris)
  iris_d = iris[,1:4]
  
  # Train a GSOM map
  gsom_map = train.gsom(iris_d_tr, spreadFactor = 0.99)
  
  # Counts
  plot(gsom_map)
  
  # Properties
  plot(gsom_map, type = "properties")
  
  # Training
  plot(gsom_map, type = "training")
  
  # Distances
  plot(gsom_map, type = "distance")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~GSOM }
\keyword{ ~SOM }
\keyword{ ~Self }
\keyword{ ~Organizing }
\keyword{ ~Map }
\keyword{ ~Growing }
\keyword{ ~Plot }