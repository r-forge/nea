\name{Knet}
\alias{Knet}
\alias{Knet.fct}
\alias{plotKnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
K-function on networks
}
\description{
K-function on networks
}
\usage{
Knet(g, nperm = 100, distmethod = "shortest.paths", prob = c(0, 0.05, 0.5, 0.95, 1), phenotype = "hits", verbose = TRUE, parallel = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{nperm}{
%%     ~~Describe \code{nperm} here~~
}
  \item{distmethod}{
%%     ~~Describe \code{distmethod} here~~
}
  \item{prob}{
%%     ~~Describe \code{prob} here~~
}
  \item{phenotype}{
%%     ~~Describe \code{phenotype} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{parallel}{
%%     ~~Describe \code{parallel} here~~
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
function (g, nperm = 100, distmethod = "shortest.paths", prob = c(0, 
    0.05, 0.5, 0.95, 1), phenotype = "hits", verbose = TRUE, 
    parallel = NULL, ...) 
{
    if (verbose) 
        cat("Computing distances on graph\n")
    D <- dist.graph(g, method = distmethod)
    K <- Knet.fct(g, distmethod = distmethod, phenotype = phenotype, 
        D = D, ...)
    K.obs <- K$K
    AUK.obs <- K$AUK
    if (!is.null(nperm) && nperm > 1) {
        if (verbose) 
            cat(paste("Knet with", nperm, "permutations\n"))
        if (verbose) 
            if (nperm < 50) 
                cat("Don't be shy! Try some more permutations next time ...\n")
        if (is.null(parallel)) {
            K.perm <- sapply(1:nperm, function(i) Knet.fct(permute.hits(g, 
                phenotype), phenotype = phenotype, D = D, ...)$K)
        }
        else {
            if (!is.wholenumber(parallel)) 
                stop("Argument 'parallel' is not a whole number")
            if (verbose) 
                cat("Setting up cluster of size", parallel)
            cl <- makeCluster(parallel, type = "SOCK")
            clusterEvalQ(cl, source("/Users/markow01/Dropbox/Papers in preparation/Network Enrichment Analysis/functions_2010-12-10.R"))
            clusterExport(cl, c("g", "D", "phenotype"))
            if (verbose) 
                cat("Running permutations on cluster\n")
            K.perm <- parSapply(cl, 1:nperm, function(i) Knet.fct(permute.hits(g, 
                phenotype), phenotype = phenotype, D = D)$K)
            stopCluster(cl)
        }
        K.quan <- apply(K.perm, 1, function(x) quantile(x, prob = prob))
        AUK.perm <- apply(K.perm, 2, function(x) sum(x)/length(x))
        pval <- max(sum(AUK.perm > AUK.obs), 1)/nperm
    }
    else {
        if (verbose) 
            cat(paste("Knet without permutations\n"))
        K.perm <- NA
        K.quan <- NA
        AUK.perm <- NA
        pval <- NA
    }
    res <- list(graph = K$g, K.obs = K.obs, AUK.obs = AUK.obs, 
        K.perm = K.perm, AUK.perm = AUK.perm, K.quan = K.quan, 
        nodeK = K$nodeK, nodeAUK = K$nodeAUK, pval = pval)
    class(res) <- "Knet"
    invisible(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
