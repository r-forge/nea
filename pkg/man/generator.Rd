\name{generator}
\alias{generator}
\alias{sumInDistBins}
\alias{panel.smoothscatter}
\alias{panel.cor}
\alias{my.mget}
\alias{is.wholenumber}
\alias{GOnodeAnno}
\alias{read.delim.BioGRID}
\alias{my.as.igraph.vs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Small helper functions
}
\description{
Small helper functions
}
\usage{
generator(g)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
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
function (g) 
{
    if (is.null(E(g)$weight)) 
        E(g)$weight = 1
    A <- get.adjacency(g, attr = "weight")
    D <- diag(apply(A, 1, sum))
    H <- A - D
    return(H)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
