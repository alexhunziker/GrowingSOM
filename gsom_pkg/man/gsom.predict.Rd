\name{predict.gsom}
\alias{predict.gsom}
\title{
  Predict function for gsom
}
\description{
  Map objects to a trained GSOM map, and return for each object the property associated with the corresponding winning unit. 
}
\usage{
  predict.gsom(gsom_object, df, retaindata = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
    Matrix or Dataframe of observations for which the Y values should be predicted.
  }
  \item{gsom_object}{
    a trained supervised gsom object. (see train_xy.gsom)
  }
  \item{retaindata}{
    logical. If set to TRUE a copy of the mapped data (df) is stored as a copy in the returned gsom object.
  }
}
\value{
  gsom object containing the following information:
  \item{nodes$position}{
    Coordinates of the units of the map (copied from the input gsom map)
  }
  \item{nodes$codes}{
    a matrix containing the codes for each variable of each unit (copied from the input gsom map)
  }
  \item{nodes$predict}{
    a matrix containing the codes for each property that are to be modelled (copied from the input gsom map)
  }
  \item{nodes$freq}{
    the number of observations that were mapped to each unit.
  }
  \item{prediction$bmn}{
    best matching unit for each observation provided.
  }
  \item{prediction$dist}{
    distance measurement of observations to their best matching units.
  }
  \item{prediction$predict}{
    prediction for the modelled properties for each observation.
  }
}
\references{
%% ~put references to the literature/web site here ~
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
  # Get some Data
  data(iris)
  
  # Feature in col 5 should be modelled
  iris_p = as.numeric(iris[,5])
  iris_d = iris[,1:4]
  
  # Devide into test and training sample
  testsample = sample(1:150, 50)
  iris_p_tr = iris_p[-testsample]
  iris_p_te = iris_p[testsample]
  iris_d_tr = iris_d[-testsample,]
  iris_d_te = iris_d[testsample,]
  
  # Train a GSOM map
  gsom_map = train_xy.gsom(iris_d_tr, spreadFactor = 0.99)
  
  # Predict
  gsom_predict = predict.gsom(gsom_map, iris_d_te)
  
  # The predicted values can now be compared with iris_p_te
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~GSOM }
\keyword{ ~SOM }
\keyword{ ~Self }
\keyword{ ~Organizing }
\keyword{ ~Map }
\keyword{ ~Growing }
\keyword{ ~Predict }
