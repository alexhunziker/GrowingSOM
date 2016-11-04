\name{summary.gsom}
\alias{summary.gsom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Summary and print methods for gsom objects
}
\description{
Summary and print methods for gsom objects. The print method shows the dimensions and the topology of the map; if information on the training data is included, the summary method additionally prints information on the size of the data and the mean distance of an object to its closest codebookvector, which is an indication of the quality of the mapping.
}
\usage{
## S3 method for class 'gsom'
summary(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gsom_model}{
  a gsom object.
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
data(iris)
iris <- iris[1:4]

gsom.iris <- train.gsom(iris)

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

