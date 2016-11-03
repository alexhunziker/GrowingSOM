\name{gsom.summary}
\alias{gsom.summary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Summarize a gsom model
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gsom.summary(gsom_model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gsom_model}{
%%     ~~Describe \code{gsom_model} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (gsom_model) 
{
    observations <- sum(gsom_model$nodes$freq)
    dimenstions <- ncol(gsom_model$nodes$weight)
    nodes <- nrow(gsom_model$nodes$position)
    iterations <- nrow(gsom_model$training)
    distance <- gsom_model$training$meandist[iterations]
    if (is.null(gsom_model[["data"]])) 
        sorstat = "is not"
    else sorstat = "is"
    print(paste("Growing SOM map with", nodes, "nodes."))
    print(paste("Training data used:", observations, "Observations with", 
        dimenstions, "Dimensions"))
    print(paste("Training data", sorstat, "stored in the model."))
    print(paste("Mean Distance to the closest unit in the map is:", 
        distance, "(after", iterations, "iterations)"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
