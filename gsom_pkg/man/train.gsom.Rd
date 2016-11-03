\name{train.gsom}
\alias{train.gsom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Train
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
train.gsom(data, spreadFactor = 0.5, keepdata = FALSE, iterations = 50, alpha, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{spreadFactor}{
%%     ~~Describe \code{spreadFactor} here~~
}
  \item{keepdata}{
%%     ~~Describe \code{keepdata} here~~
}
  \item{iterations}{
%%     ~~Describe \code{iterations} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (data, spreadFactor = 0.5, keepdata = FALSE, iterations = 50, 
    alpha, ...) 
{
    min <- apply(data, 2, function(x) {
        min(x)
    })
    max <- apply(data, 2, function(x) {
        max(x)
    })
    df <- t(apply(data, 1, function(x) {
        (x - min)/ifelse(max == min, 1, (max - min))
    }))
    t1 <- Sys.time()
    gsom_model <- grow.gsom(gsom_model, df, iterations, spreadFactor)
    t2 <- Sys.time()
    print(t2 - t1)
    norm_param <- data.frame(min = min, max = max)
    gsom_model[["norm_param"]] <- norm_param
    if (keepdata == TRUE) {
        gsom_model[["data"]] = data
    }
    class(gsom_model) = "gsom"
    return(gsom_model)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
