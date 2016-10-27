\name{gsom.map}
\alias{gsom.map}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gsom.map(newdf, gsom_model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{newdf}{
%%     ~~Describe \code{newdf} here~~
}
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
function (newdf, gsom_model) 
{
    nodes <- list(position = gsom_model$nodes$position, weight = gsom_model$nodes$weight, 
        avgerror = rep(0, nrow(gsom_model$nodes$position)), freq = rep(0, 
            nrow(gsom_model$nodes$position)))
    observation <- data.frame(winner = numeric(), error = numeric())
    mapped_data <- list(nodes = nodes, observation = observation)
    newdf <- t(apply(newdf, 1, function(x) {
        (x - gsom_model$norm_param$min)/(gsom_model$norm_param$max - 
            gsom_model$norm_param$min)
    }))
    for (j in 1:nrow(newdf)) {
        errors <- sqrt(rowSums(sweep(gsom_model$nodes$weight, 
            MARGIN = 2, newdf[j, ], FUN = "-")^2, dims = 1))
        minError = min(errors)
        winner <- which(grepl(minError, errors))
        if (length(winner) > 1) 
            winner <- winner[1]
        mapped_data$nodes$freq[winner] <- mapped_data$nodes$freq[winner] + 
            1
        mapped_data$nodes$avgerror[winner] <- mapped_data$nodes$avgerror[winner] + 
            minError
        mapped_data$observation[j, ] <- c(winner = winner, error = minError)
    }
    mapped_data$nodes$avgerror <- mapped_data$nodes$avgerror/mapped_data$nodes$freq
    return(mapped_data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
