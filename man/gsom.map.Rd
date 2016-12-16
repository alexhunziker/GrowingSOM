\name{map.gsom}
\alias{map.gsom}
\title{
  Map data to a supervised or unsupervised gsom object
}
\description{
  Map a data matrix onto a trained gsom object.
}
\usage{
## S3 method for class 'gsom'
map.gsom(gsom_object, df, retaindata = FALSE)
}

\arguments{
  \item{gsom_object}{
    a trained gsom map (either unsupervised, as obtained by train.gsom, or supervised, as returned by train_xy.gsom).
  }
  \item{df}{
    Data matrix, with rows corresponding to objects, and columns to the dimensions the gsom object was trained with
  }
  \item{retaindata}{
    if set to TRUE a copy of the mapped data (unscaled) will be added to the returned object.
  }
}

\value{
  Returns a S3 object of type "gsom" containing:
  \item{nodes$position}{
    the location of the nodes on the map.
  }
  \item{nodes$codes}{
    codes that were established during the training for each node and dimension of the data.
  }
  \item{nodes$freq}{
    how many times each node was the best matching node for the mapped matrix.
  }
  \item{mapped$bmn}{
    contains the best matching node for each of the data that was mapped.
  }
  \item{mapped$dist}{
    distance from best matching node for each row of the mapped data-matrix.
  }
  \item{data}{
    Unscaled copy of the data that was mapped. 
  }
}

\author{
  Alex Hunziker
}

\seealso{
  train.gsom, predict.gsom
}
\examples{
  #Get some data
  data(iris)
  s = sample(1:150, 100)
  train_set = iris[s,1:4]
  test_set = iris[-s,1:4]

  # Create a GSOM Model
  gsom_iris <- train.gsom(train_set, spreadfactor=0.75)

  # Mapping
  mapped_iris <- map(gsom_iris, test_set)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~GSOM }
\keyword{ ~SOM }
\keyword{ ~Self }
\keyword{ ~Organizing }
\keyword{ ~Map }
\keyword{ ~Growing }
\keyword{ ~Map }
