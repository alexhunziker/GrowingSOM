\name{GrowingSOM}
\alias{GSOM}
\alias{GrowingSelfOrganizingMap}
\docType{package}
\title{
\packageTitle{GrowingSOM}
}
\description{
  \packageDescription{GrowingSOM}
}
\details{

  The GSOM was developed to address the issue of identifying a suitable map size in the SOM.\cr\cr
    It starts with a minimal number of nodes and grows new nodes on the boundary based on a heuristic. By using the value called Spread Factor (SF), the growth of the GSOM can be controlled.\cr\cr
    This R-Package supports both a supervised and an unsupervised version of Growing Self-Organizing maps.\cr
    Maps can be trained, new data can be mapped to existing maps and in the case of supervised GSOM maps, the dependant values can be predicted.
  \packageIndices{GrowingSOM}
}
\author{
  Package Author: \packageAuthor{GrowingSOM}

  Maintainer: \packageMaintainer{GrowingSOM}
}
\references{
  This package is conceptually mostly based on the version of Growing Self-Organizing Maps as described in this paper: \cr
  \cr
  Damminda Alahakoon, Saman K. Halgamuge (2000): Dynamic Self-Organizing Maps with Controlled Growth for Knowledge Discovery. IEEE TRANSACTIONS ON NEURAL NETWORKS, VOL. 11.
}
\keyword{ package }
\keyword{ growing }
\keyword{ self-organizing }
\keyword{ Map }
\keyword{ GrowingSOM }
\keyword{ GSOM }
%%\seealso{
%%    Non-growing Self-Organizing maps:\cr
%%    %\pkg{roxygen2}
%%      %%\code{\link[<kohonen>:<kohonen>-package]{<kohonen>-package}}
%%}
\examples{

  # Get some data
  data("iris")
  s = sample(1:150, 100)
  train_set = iris[s,1:4]
  test_set = iris[-s,1:4]
  
  # Train GrowingSOM
  gsom_map <- train.gsom(train_set, spreadFactor=0.8)
  
  # Some Plots
  plot(gsom_map, type="count")
  par(mfrow=c(2,2))
  plot(gsom_map, type="property")
  par(mfrow=c(1,1))
  
  # Map new data
  gsom_mapped = map.gsom(gsom_map, test_set)
  plot(gsom_mapped)
  
}
