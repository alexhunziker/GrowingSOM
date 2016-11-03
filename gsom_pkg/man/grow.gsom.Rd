\name{grow.gsom}
\alias{grow.gsom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Grow the gsom model
}
\description{
  Internal Function that guides the growth and training of the gsom model. Not to be called by the user.
  Any gsom model should be calculated by calling gsom().
}
\usage{
grow.gsom(gsom_model, df, repet, spreadFactor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gsom_model}{
%%     ~~Describe \code{gsom_model} here~~
}
  \item{df}{
%%     ~~Describe \code{df} here~~
}
  \item{repet}{
%%     ~~Describe \code{repet} here~~
}
  \item{spreadFactor}{
%%     ~~Describe \code{spreadFactor} here~~
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
function (gsom_model, df, repet, spreadFactor) 
{
    lentr <- 10000
    lentn <- 10000
    lrinit <- 0.9
    alpha <- 0.9
    radius <- 3
    df <- as.matrix(df)
    weights <- matrix(0, nrow = lentn, ncol = ncol(df))
    weights[1:4, ] <- runif(4 * ncol(df))
    distnd <- rep(0, lentn)
    freq <- rep(0, lentn)
    gt = -ncol(df) * log(spreadFactor) * 0.01 * nrow(df)
    npos <- matrix(0, nrow = lentn, ncol = 2)
    npos[1:4, ] <- c(0, 1, 1, 0, 1, 0, 1, 0)
    if (repet > lentr) 
        error("Max nr of iterations exceeded.")
    currtrain <- matrix(0, nrow = lentr, ncol = 5)
    outc = .C("som_train_loop", df = as.double(df), weights = as.double(weights), 
        distnd = as.double(distnd), prep = as.integer(repet), 
        plendf = as.integer(nrow(df)), plennd = as.integer(4), 
        plrinit = as.double(lrinit), freq = as.double(freq), 
        alpha = as.double(alpha), pdim = as.integer(ncol(df)), 
        gt = as.double(gt), npos = as.double(npos), pradius = as.double(radius), 
        plentn = as.integer(lentn), plentd = as.integer(nrow(df)), 
        currtrain = as.double(currtrain), plentr = as.integer(lentr))
    training <- matrix(outc$currtrain, ncol = 5)
    training <- training[1:repet, ]
    colnames(training) <- c("iteration", "training_stage", "meandist", 
        "num_of_nodes", "nodegrowth")
    weights <- matrix(outc$weights, ncol = ncol(df))
    weights <- weights[1:outc$plennd, ]
    colnames(weights) <- colnames(df)
    npos <- matrix(outc$npos, ncol = 2)
    npos <- matrix(npos[1:outc$plennd, ], ncol = 2)
    colnames(npos) <- c("x", "y")
    freq <- outc$freq[1:outc$plennd]
    distnd <- outc$distnd[1:outc$plennd]
    nodes <- list(position = data.frame(npos), weight = data.frame(weights), 
        error = distnd, freq = freq)
    training = data.frame(training)
    GT = outc$gt
    gsom_model <- list(nodes = nodes, training = training, GT = GT, 
        norm_param = 0)
    class(gsom_model) <- "gsom"
    return(gsom_model)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
