\name{summary.gsom}
\alias{summary.gsom}
\alias{print.gsom}

\title{
 Summary and print methods for gsom objects
}
\description{
Summary and print methods for gsom objects. The print method shows the dimensions and the topology of the map; if information on the training data is included, the summary method additionally prints information on the size of the data and the mean distance of an object to its closest node, which is an indication of the quality of the mapping.
}
\usage{
\method{summary}{gsom}(object, \dots)

\method{print}{gsom}(x, \dots)
}

\arguments{
  \item{x}{
  an object of class gsom.
}
  \item{object}{
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

\keyword{ GSOM }
\keyword{ GrowingSOM }
\keyword{ Summary }
\keyword{ Print }

