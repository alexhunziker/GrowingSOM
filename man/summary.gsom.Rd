\name{summary.gsom}
\alias{summary.gsom}
\alias{print.gsom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Summary and print methods for gsom objects
}
\description{
Summary and print methods for gsom objects. The print method shows the dimensions and the topology of the map; if information on the training data is included, the summary method additionally prints information on the size of the data and the mean distance of an object to its closest node, which is an indication of the quality of the mapping.
}
\usage{
## S3 method for class 'gsom'
summary(gsom_object, ...)

## S3 method for class 'gsom'
print(gsom_object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gsom_object}{
  an object of class gsom.
}
\item{...}{
  not used.
}
}

\author{
  Alex Hunziker
}

\seealso{
  train.gsom, train_xy.gsom
}

\examples{
# load some data
data(iris)
iris <- iris[1:4]

#make a gsom map
gsom.iris <- train.gsom(iris)

#summarize
summary(gsom.iris)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~GSOM }
\keyword{ ~SOM }
\keyword{ ~Self }
\keyword{ ~Organizing }
\keyword{ ~Map }
\keyword{ ~Growing }
\keyword{ ~Summary }
\keyword{ ~Print }

