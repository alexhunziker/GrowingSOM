\name{auto_mpg}
\alias{auto_mpg}
\title{
  Auto MPG Data Set
}

\description{
"The data concerns city-cycle fuel consumption in miles per gallon, to be predicted in terms of 3 multivalued discrete and 5 continuous attributes." (Quinlan, 1993)
}

\usage{
  data("auto_mpg")
}

\format{
Number of Instances: 398\cr
Number of Attributes: 8\cr
\cr
1. mpg: continuous\cr
2. cylinders: multi-valued discrete\cr
3. displacement: continuous\cr
4. horsepower: continuous\cr
5. weight: continuous\cr
6. acceleration: continuous\cr
7. model year: multi-valued discrete\cr
8. origin: multi-valued discrete\cr
9. car name: string (unique for each instance)
}

\details{
  This dataset was taken from the StatLib library which is maintained at Carnegie Mellon University. The dataset was used in the 1983 American Statistical Association Exposition. Missing values have been removed.\cr
}

\references{
  Quinlan,R. (1993). Combining Instance-Based and Model-Based Learning. In Proceedings on the Tenth International Conference of Machine Learning, 236-243, University of Massachusetts, Amherst. Morgan Kaufmann. 
}

\examples{
  # load data
  data("auto_mpg")
  
  # show head
  head(auto_mpg)
}

\keyword{ autos_mpg }
\keyword{ data }

