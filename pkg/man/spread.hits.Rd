\name{spread.hits}
\alias{spread.hits}
\alias{permute.hits}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spread hits over a network
}
\description{
Spread hits over a network
}
\usage{
spread.hits(g, h = 10, lambda = 1, distmethod = "shortest.paths", start.node = NULL, hitColor = "red", D = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
  \item{distmethod}{
%%     ~~Describe \code{distmethod} here~~
}
  \item{start.node}{
%%     ~~Describe \code{start.node} here~~
}
  \item{hitColor}{
%%     ~~Describe \code{hitColor} here~~
}
  \item{D}{
%%     ~~Describe \code{D} here~~
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

1+1

## The function is currently defined as
function (g, h = 10, lambda = 1, distmethod = "shortest.paths", 
    start.node = NULL, hitColor = "red", D = NULL) 
{
    if (class(g) != "igraph") 
        stop("Input 'g' is not a graph: Stop feeding me crap!")
    if (is.null(V(g)$hits)) 
        V(g)$hits <- 0
    if (is.null(V(g)$color)) 
        V(g)$color <- "grey"
    if (is.null(start.node)) 
        start.node <- sample(as.character(V(g)[hits == 0]), 1)
    if (is.null(D)) 
        D <- dist.graph(g, method = distmethod)
    prob <- lambda * exp(-lambda * D[start.node, ])
    prob[start.node] <- 0
    if (sum(prob) > 0 && sum(prob > 0) < h) 
        prob[prob == 0] <- min(prob[prob != 0])
    if (sum(prob) > 0) 
        prob <- prob/sum(prob)
    else prob <- rep(1, length(prob))/length(prob)
    sampled <- c(start.node, sample(as.character(V(g)), size = h - 
        1, prob = prob))
    sampled <- as.numeric(sampled)
    V(g)[sampled]$hits <- 1
    V(g)[sampled]$color <- hitColor
    return(g)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
