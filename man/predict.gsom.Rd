\name{predict.gsom}
\alias{predict.gsom}
\title{
  Predict function for gsom
}
\description{
  Map objects to a supervised GrowingSOM map, and return for each object the property associated with the corresponding winning node. 
}
\usage{
\method{predict}{gsom}(object, df, retaindata = FALSE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    a trained supervised gsom object. (see train_xy.gsom)
  }
  \item{df}{
    Matrix or Dataframe of observations for which the Y values should be predicted.\cr
    Columns of this input argument must be corresponding to the ones used as data when training the map.
  }
  \item{retaindata}{
    logical. If set to TRUE a copy of the mapped data (df) is stored as a copy in the returned gsom object.
  }
  \item{...}{
    not used.
  }
}
\value{
  gsom object containing the following information:
  \item{nodes$position}{
    Coordinates of the units of the map (copied from the input gsom map)
  }
  \item{nodes$codes}{
    a matrix containing the codes for each variable of each node (copied from the input gsom map)
  }
  \item{nodes$predict}{
    a matrix containing the codes for each property that are to be modelled (copied from the input gsom map)
  }
  \item{nodes$freq}{
    the number of observations that were mapped to each node.
  }
  \item{prediction$bmn}{
    best matching node for each observation provided.
  }
  \item{prediction$dist}{
    distance measurement of observations to their best matching units.
  }
  \item{prediction$value}{
    prediction for the modelled properties for each observation.
  }
}

\author{
  Alex Hunziker
}
\note{
  Distance Measurement: The following Formula is used for determining the distance from the best matching unit: \cr
  \code{sqrt((sum(df-code))^2 / Dimenstions)}\cr
  This is a different scale than is used by the kohonen package.
}

\seealso{
  train_xy.gsom, map.gsom
}
\examples{
  # load data
  data("auto_mpg")
  s = sample(1:392, 300)
  train_set = auto_mpg[s,1:8]
  test_set = auto_mpg[-s,1:8]
  
  # Train Gsom Model (hexagonal grid)
  gsom_map <- train_xy.gsom(train_set[,2:8], train_set[,1], spreadFactor = 0.9, nhood="hex")
  
  print(gsom_map)
  plot(gsom_map, type = "predict")
  
  # Predict mpg for the test set
  gsom_predictions = predict.gsom(gsom_map, test_set[,2:8])
}

\keyword{ GSOM }
\keyword{ GrowingSOM }
\keyword{ Predict }
