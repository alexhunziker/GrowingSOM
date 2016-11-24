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
map.gsom(gsom_model, df, retaindata = FALSE, ...)
}

\arguments{
  \item{gsom_model}{
    a trained gsom model (either unsupervised, as obtained by train.gsom, or supervised, as returned by train_xy.gsom).
  }
  \item{df}{
    Data matrix, with rows corresponding to objects, and columns to the dimensions the gsom object was trained with
  }
  \item{retaindata}{
    if set to TRUE a copy of the mapped data (unscaled) will be added to the returned object.
  }
  \item{...}{
    not used.
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
  \item{nodes$frequency}{
    how many times each node was the best matching unit for the mapped matrix.
  }
  \item{mapped$codes}{
    contains the best matching unit for each of the data that was mapped.
  }
  \item{mapped$dist}{
    distance from best matching unit for each row of the mapped data-matrix.
  }
  \item{...}{
    Other useful parameters to be added
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
rand.iris <- iris[sample(nrow(iris)),]
train <- rand.iris[1:100, 1:4]
test <- rand.iris[101:150, 1:4]

# Create a GSOM Model
gsom.iris <- train.gsom(train, spreadfactor=0.95)

#Actual Mapping
mapped.iris <- map(gsom.iris, test)
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
