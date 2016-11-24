\name{gsom.summary}
\alias{gsom.summary}
\alias{gsom.print}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Summarize a gsom model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
  summary.gsom(gsom_object, ...)
  
  print.gsom(gsom_object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gsom_object}{
    any valid gsom object
  }
  \item{...}{
    
  }
}

\value{
  Prints Information about the gsom_object
}

\author{
  Alex Hunziker
}

\seealso{
  train.gsom, train_xy.gsom, map.gsom, predict.gsom
}
\examples{
  # Get some Data
  data(iris)
  iris_d = iris[,1:4]
  
  # Train a GSOM map
  gsom_map = train.gsom(iris_d_tr, spreadFactor = 0.99)
  
  # Summary / Print
  print(gsom_map)
  summary(gsom_map)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~GSOM }
\keyword{ ~SOM }
\keyword{ ~Self }
\keyword{ ~Organizing }
\keyword{ ~Map }
\keyword{ ~Growing }
\keyword{ ~Print }
\keyword{ ~Summary }
