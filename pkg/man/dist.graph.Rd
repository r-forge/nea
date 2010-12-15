\name{dist.graph}
\alias{dist.graph}
\alias{graph.diffusion}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Distances on graph
}
\description{
Distances on graphs
}
\usage{
dist.graph(g, method = "shortest.paths", correctInf = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{correctInf}{
%%     ~~Describe \code{correctInf} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
bla
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
bla bla bla
}
\references{
http://www.markowetzlab.org/software
}
\author{
Florian Markowetz
}
\note{
bla bla bla
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
function (g, method = "shortest.paths", correctInf = TRUE, ...) 
{
    D <- switch(method, shortest.paths = shortest.paths(g, ...), 
        diffusion = graph.diffusion(g, ...)$dist)
    if (correctInf) 
        D[D == Inf] <- max(D[D != Inf]) + 1
    dimnames(D) <- list(V(g), V(g))
    return(D)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
