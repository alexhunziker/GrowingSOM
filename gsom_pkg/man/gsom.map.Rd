\name{gsom.map}
\alias{gsom.map}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Mapping Function for gsom
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
gsom.map(df, gsom_model, retaindata = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
}
  \item{gsom_model}{
%%     ~~Describe \code{gsom_model} here~~
}
  \item{retaindata}{
%%     ~~Describe \code{retaindata} here~~
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
function (df, gsom_model, retaindata = FALSE) 
{
    min <- gsom_model$norm_param$min
    max <- gsom_model$norm_param$max
    df <- t(apply(df, 1, function(x) {
        (x - min)/(max - min)
    }))
    codes <- rep(0, times = nrow(df))
    ndist <- rep(0, times = nrow(df))
    freq <- rep(0, times = nrow(gsom_model$nodes$weight))
    outc = .C("map_data", plendf = as.integer(nrow(df)), lennd = as.integer(nrow(gsom_model$nodes$weight)), 
        dim = as.integer(ncol(gsom_model$nodes$weight)), df = as.double(df), 
        weights = as.double(as.matrix(gsom_model$nodes$weight)), 
        codes = as.double(codes), ndist = as.double(ndist), freq = as.double(freq))
    dist <- outc$ndist
    code <- outc$codes
    print(code)
    gsom_mapped = list()
    gsom_mapped[["nodes"]] = gsom_model$nodes
    gsom_mapped[["nodes"]]$error = NULL
    gsom_mapped[["nodes"]]$freq = outc$freq
    gsom_mapped[["mapped"]] = data.frame(codes = code, dist = dist)
    gsom_mapped[["predict"]] = predict
    gsom_mapped[["norm_param"]] = gsom_model$scale
    if (retaindata) 
        gsom_mapped[["data"]] == df
    class(gsom_mapped) = "gsom_mapped"
    return(gsom_mapped)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
