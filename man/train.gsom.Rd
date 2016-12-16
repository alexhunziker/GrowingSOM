\name{train.gsom}
\alias{train.gsom}

\title{
  Train a Growing Self-Organizing Map
}
\description{
  Computes a growing self-organizing map for mapping high-dimendional data to 2D. 
}
\usage{
  train.gsom(data, spreadFactor=0.9, keepdata=FALSE, 
      iterations=50, alpha=0.9, beta = 0.5, gridsize = FALSE, 
      nhood= "rect", initrad)
}

\arguments{
  \item{data}{
    a matrix or data.frame, with each row representing an observation and each column a dimension.
  }
  \item{spreadFactor}{
    the spread factor determines the rate with which new units are added to the map. Values close to 0 lead to few growth and therefore less nodes thatn values close to 1. The default value is 0.9. 
  }
  \item{keepdata}{
    if set to TRUE, a copy of the traindata will be stored in the gsom object.
  }
  \item{iterations}{
    number of times that the dataframe will be presented to the network. (Growing and Smoothing Phase combined)
  }
  \item{alpha}{
    discount factor for the learning rate during the growing phase of the training. Values should be between 0 and 1.
  }
  \item{beta}{
    propagation rate. Determines the rate at which the error of a node, that cannot grow any nodes, is passed on to its neighbours. Suggested values may range between 0 and 1.
  }
  \item{gridsize}{
    default value is FALSE. If a nummeric value is entered, the grid-size of the network will be preditermined as a square with length of gridsize. No growth of the network will take place in this case.
  }
  \item{nhood}{
    define how the grid will be built, and how the neighbourhood will look like consequently. Allowed values are "rect" (rectangular) and "hex" (hexagonal).
  }
  \item{initrad}{
    if the gridsize is predetermined, the initial radius can be chosen here. If left blank, the square root of gridsize will be taken. \cr
    A larger initrad can increase the quality of the clustering. However, the script can get very slow when a too large value is chosen. (number of necessary computations rises exponentially)
  }
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
  Euclidean distance is used in order to calculate the distances.
}
\value{
  an S3 object of class "gsom" with components:
  \item{nodes$position}{
    the location of the nodes on the map.
  }
  \item{nodes$codes}{
    codes that were established during the training for each node and dimension of the data.
  }
  \item{nodes$distance}{
    average distance of observations from their best matching units, per unit.
  }
  \item{nodes$freq}{
    how many times each node was the best matching unit during the last iteration.
  }
  \item{training}{
    reports on the progress of the training after each iteration over the data. "training_stage" indicates whether the algotithm was in the growing phase (1) or the smoothing phase (2), "meandist" records the average distance to the best matching unit, "num_of_nodes" stores the size of the map, and "nodegrowth" keeps track of the nodes grown during each iteration.
  }
  \item{GT}{
    the accumulated error per node that was required to grow a new unit during the training phase of the network.
  }
  \item{norm_param}{
    parameters which were used for each dimension to normalize the data. Needed for plotting the correct scale and to map future data.
  }
  \item{data}{
    data that was used to train the model, not normalized.
  }
}
\references{
 Damminda Alahakoon, Saman K. Halgamuge (2000): Dynamic Self-Organizing Maps with Controlled Growth for Knowledge Discovery. IEEE TRANSACTIONS ON NEURAL NETWORKS, VOL. 11.
}
\author{
  Alex Hunziker
}
\note{
  In contrast to the paper the algorighm is based on the following adjustments have been made:
  
  1) Changing the Formula to calculate the growth-rate from the SpreadingFactor, in order to take into consideration the number of observations
  
  2) The learning rate is reduced according to the formula used for regualr Kohonen networks during phase tow, in order to prevent too fast depreciation.
  
  The decrease of neighbourhood size only takes place in phase two and happens on a linear basis from the specified start value to 0.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  train_xy.gsom, map.gsom, plot.gsom
}
\examples{

# load data
data(iris)
iris <- iris[,1:4]

# Train Gsom Model
gsom_iris <- train.gsom(iris)

# Some more Parameters
gsom_iris <- train.gsom(iris, spreadFactor=0.95, keepdata=TRUE, iterations=30, alpha=0.5, gridsize = FALSE, nhood= "rect")

# Fixed Grid size
gsom_iris <- train.gsom(iris, iterations=30, gridsize = 10)

}

\keyword{ Growing }
\keyword{ GSOM }
\keyword{ Self-Organizing }
\keyword{ Map }
\keyword{ SOM }
\keyword{ training }
\keyword{ initialize }