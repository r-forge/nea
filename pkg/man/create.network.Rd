\name{create.network}
\alias{create.network}
\alias{create.grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create networks
}
\description{
Create networks
}
\usage{
create.network(n = 100, type = "barabasi", edge.weights = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{edge.weights}{
%%     ~~Describe \code{edge.weights} here~~
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
bla bla bla
}
\references{
http://www.markowetzlab.org/software
}
\author{
Florian Markowetz
}
\note{
bla
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
y=1
x=2
x+y

## The function is currently defined as
function (n = 100, type = "barabasi", edge.weights = NULL, ...) 
{
    if (!type \%in\% c("barabasi", "erdos.renyi", "grid")) 
        stop("Input 'type' unknown: Stop feeding me crap!")
    g <- switch(type, barabasi = barabasi.game(n, directed = FALSE, 
        ...), erdos.renyi = erdos.renyi.game(n, p.or.m = (2 * 
        n - 2 * floor(sqrt(n))), type = "gnm", directed = FALSE, 
        ...), grid = create.grid(n))
    if (is.null(edge.weights) || edge.weights == 0) 
        E(g)$weight <- sample(1:10, length(E(g)), replace = TRUE)
    else E(g)$weight <- edge.weights
    return(g)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }

