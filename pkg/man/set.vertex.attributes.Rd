\name{set.vertex.attributes}
\alias{set.vertex.attributes}
\alias{get.all.vertex.attributes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get/set more than one vertex attribute in an igraph object
}
\description{
get/set more than one vertex attribute in an igraph object
}
\usage{
set.vertex.attributes(graph, D)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{graph}{
%%     ~~Describe \code{graph} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
1+1

## The function is currently defined as
function (graph, D) 
{
    if (!is.igraph(graph)) 
        stop("Not a graph object")
    if (!is.data.frame(D)) 
        stop("Not a dataframe")
    for (i in colnames(D)) graph <- set.vertex.attribute(graph, 
        i, value = D[, i])
    return(graph)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
