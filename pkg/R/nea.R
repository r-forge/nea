####################################################################################
##
## Generate networks according to different models
##
create.grid <- function(n=100){
	nr <- floor(sqrt(n))
	## empty graph:
	g <- graph.empty(n,directed=FALSE) 
	## add edges from i to i+1 (within layer)
	#g <- add.edges(g,c(0,rep(1:(n-2),each=2),(n-1)))
	g <- add.edges(g,rbind(0:(n-2),1:(n-1))[,c(rep(TRUE,nr-1),FALSE)])
	## add edges between layers
	g <- add.edges(g,as.numeric(rbind(0:(n-nr-1),nr:(n-1))))
	return(g)
} # end function create.grid

## main function:
create.network <- function(n=100,type="barabasi", edge.weights=NULL,...){
	if(!type%in%c("barabasi","erdos.renyi","grid")) stop("Input 'type' unknown: Stop feeding me crap!")
	##
	g <- switch(type,
		barabasi    = barabasi.game(n,directed=FALSE, ...),
		erdos.renyi = erdos.renyi.game(n,p.or.m=(2*n-2*floor(sqrt(n))),type="gnm",directed=FALSE,...), ## ~nr of edges of grid
		grid        = create.grid(n)	
		)
	##
	## add edge weights
	if (is.null(edge.weights)||edge.weights==0) E(g)$weight <- sample(1:10,length(E(g)),replace=TRUE) else E(g)$weight <- edge.weights

	##
	return(g)
} # end function create.network

####################################################################################
##
## Distance measures on graphs
##

## produces generator matrix for diffusion kernel
## if all edge weights = 1, this is the same as -graph.laplacian(g)
generator <- function(g){	
	if(is.null(E(g)$weight)) E(g)$weight = 1
	A <- get.adjacency(g,attr="weight")
	D <- diag(apply(A,1,sum))
	H <- A-D
	return(H)
} ## end function generator


## Diffusion Kernels on Graphs and Other Discrete Structures, Kondor and Lafferty, 2002
graph.diffusion <- function(g,beta=1,correctNeg=TRUE){
		H <- generator(g)
		x <- eigen(H) ## H ~ x$vectors %*% diag(x$values) %*% t(x$vectors)
		K <- x$vectors %*% diag(exp(beta * x$values)) %*% t(x$vectors)
		D2 <- outer(diag(K),diag(K),"+") - 2*K
		if (any(D2<0) && correctNeg){
			warning("Negative 'distances' set to Zero!")
			D2[D2<0] <- 0
			}
		D <- sqrt(D2)
		dimnames(D) <- list(V(g),V(g))
		return(list(kernel=K,dist=D))
} ## end function graph.diffusion

## main distance function that interfaces to all others
dist.graph <- function(g,method="shortest.paths",correctInf=TRUE,...){
	##	
	D <- switch(method,
		##floyd.warshall = floyd.warshall.all.pairs.sp(igraph.to.graphNEL(g)),
		shortest.paths = shortest.paths(g, ...),
		diffusion      = graph.diffusion(g, ...)$dist
	)
	##
	if (correctInf) D[D==Inf] <- max(D[D!=Inf])+1  ## D=Inf happens if there are unconnected nodes
	##
	dimnames(D) <- list(V(g),V(g))
	return(D)
	##
} ## end function dist.graph

####################################################################################
##
## spread hits over a given network
##	
spread.hits <- function(g,
						h=10,
						lambda=1,
						distmethod="shortest.paths",
						start.node=NULL,
						hitColor="red",
						D = NULL
						){
	##
	if (class(g)!="igraph")  stop("Input 'g' is not a graph: Stop feeding me crap!")
	if (is.null(V(g)$hits))  V(g)$hits <- 0
	if (is.null(V(g)$color)) V(g)$color <- "grey"
	if (is.null(start.node)) start.node <- sample(as.character(V(g)[hits==0]),1)
	##
	if (is.null(D)) D <- dist.graph(g,method=distmethod)
	prob <- lambda*exp(-lambda * D[start.node,]) ## prob decays exponentially with distance to start.node
	prob[start.node] <- 0                        ## don't sample start node (again)
	if (sum(prob)>0 && sum(prob>0)<h) prob[prob==0] <- min(prob[prob!=0])
	if (sum(prob)>0) prob <- prob/sum(prob) else prob <- rep(1,length(prob))/length(prob)   ## normalize 
	sampled <- c(start.node,sample(as.character(V(g)),size=h-1,prob=prob)) ## start.node plus h-1 others
	sampled <- as.numeric(sampled)
	##
	V(g)[sampled]$hits <- 1
	V(g)[sampled]$color <- hitColor
	##
	
	##
	return(g)
	## 	
} ## end function spread.hits
	
####################################################################################
##
## Enrichment statistic
##

## discrete and weighted hits
Knet.fct <- function(g,distmethod="shortest.paths",D=NULL,node.contributions=FALSE,phenotype="hits",nsteps=100,nTopNodes=10,colorTopNodes="red4"){
	
	##
	pheno  <- get.vertex.attribute(g,phenotype) ## phenotypes as specified by user
	if (is.null(pheno)) stop("no phenotype in graph!")
	if (!is.numeric(pheno)){ 
		warning("phenotype converted to numeric")
		pheno <- as.numeric(pheno)
	}
	
	## setup:
	hits   <- (pheno != 0) & !is.na(pheno) ## 0's and NA's are no hits
	phits  <- pheno[hits]                  ## non-zero phenotypes
	nhits  <- sum(pheno,na.rm=TRUE)        ## complete phenotype
	nnodes <- vcount(g)                    ## number of nodes
	
	## compute graph distance:
	if (is.null(D)) D <- dist.graph(g,method=distmethod)
	if (vcount(g)!=ncol(D) || ncol(D) !=nrow(D)) stop("Distance matrix D must be square with nrow(D)=vcount(g)")
	D.hits <- D[hits,hits]
	
	## compute K:
	breaks <- seq(from=0,to=max(D),length.out=min(nsteps,max(D)+1))
	S <- sumInDistBins(D.hits,phits,breaks)		
	Z <- nhits * (nhits/nnodes) ## normalizing constant
	K <- 1/Z * S                ## K-stat on graphs
	
	## normalize K to lie in [0,1]
	theo.max <- (nnodes * (nhits^2 - nhits))/nhits^2
	if (round(max(K),2) != round(theo.max,2)) warning("theoretical max not observed max")
	K       <- K/theo.max
	AUK     <- sum(K)/length(K)
	
	## compute individual node contributions
	nodeK   <- apply(D.hits,1,function(x) sumInDistBins(x,phits,breaks,normalize=TRUE))
	nodeAUK <- apply(nodeK,2,function(x)sum(x)/length(x))
	
	##
	if (is.null(V(g)$shape)) V(g)$shape="circle"
	if (is.null(V(g)$color)) V(g)$color="white"
	topNodes <- as.numeric(names(sort(nodeAUK,decreasing=TRUE)[1:min(nhits,nTopNodes)]))+1
	V(g)$shape[topNodes] <- "rectangle"
	V(g)$color[topNodes] <- colorTopNodes
	
	##	
	return(list(graph=g,K=K,AUK=AUK,nodeK=nodeK,nodeAUK=nodeAUK))
	##
} ## end function Knet

sumInDistBins <- function(D,phits,breaks,normalize=FALSE){
	times <- ifelse(is.null(nrow(D)),1,nrow(D)) 
	S0 <- tapply(rep(phits,times),cut(D,breaks),FUN=sum) ## hits in each distance bin
	S <- cumsum(ifelse(is.na(S0),0,S0))
	if (normalize) S <- S/max(S)		
	return(S)
} ## end function sumInDistBins

##
##
##
permute.hits <- function(g,phenotype="hits") set.vertex.attribute(g,name=phenotype,value=sample(get.vertex.attribute(g,phenotype)))

##
##
##
Knet <- function(	g,
					nperm=100,
					distmethod="shortest.paths",
					prob=c(0,0.05,0.5,0.95,1),
					phenotype="hits",
					verbose=TRUE,
					parallel=NULL,
					...
				){
	##
	if (verbose) cat("Computing distances on graph\n")
	D <- dist.graph(g,method=distmethod)
	K <- Knet.fct(g,distmethod=distmethod,phenotype=phenotype,D=D,...)
	K.obs <- K$K
	AUK.obs <- K$AUK
	##
	if (!is.null(nperm) && nperm>1){
		if (verbose) cat(paste("Knet with",nperm,"permutations\n"))
		if (verbose) if (nperm < 50) cat("Don't be shy! Try some more permutations next time ...\n")
		if (is.null(parallel)){
			if (verbose) cat("Running",nperm,"permutations\n")
			K.perm <- sapply(1:nperm,function(i) Knet.fct(permute.hits(g,phenotype),phenotype=phenotype,D=D,...)$K)
		} else {
			require(snow)
			if (!is.wholenumber(parallel)) stop("Argument 'parallel' is not a whole number\n")
			if (verbose) cat("Setting up cluster of size",parallel,"\n")
			##browser()
			cl <- makeCluster(parallel,type="SOCK",verbose=FALSE)
			clusterEvalQ(cl,library(nea))
			##clusterExport(cl, c("g","D","phenotype")) 
			if (verbose) cat("Running",nperm,"permutations on cluster\n")
			K.perm <- parSapply(cl,1:nperm,function(i) Knet.fct(permute.hits(g,phenotype),phenotype=phenotype,D=D)$K)
			stopCluster(cl)
			}
		K.quan <- apply(K.perm,1,function(x)quantile(x,prob=prob))
		AUK.perm <- apply(K.perm,2,function(x)sum(x)/length(x))
		pval <- max(sum(AUK.perm > AUK.obs),1)/nperm
	} else {
		if (verbose) cat(paste("Knet without permutations\n"))
		K.perm   <- NA
		K.quan   <- NA
		AUK.perm <- NA
		pval     <- NA	
	}
	##
	res <- list(graph=K$g,K.obs=K.obs,AUK.obs=AUK.obs,K.perm=K.perm,AUK.perm=AUK.perm,K.quan=K.quan,nodeK=K$nodeK,nodeAUK=K$nodeAUK,pval=pval)
	class(res) <- "Knet"
	invisible(res)
	##
} ## end function Knet.permute

plot.Knet <- function(x,...){ ## plot.Knet <- function(x,sequential=FALSE){
	res <- x
	##
	if (!all(c("K.perm","K.quan","AUK.perm")%in%names(res))){
		plot(res$K.obs,main="Knet function",type="l",lwd=3,col="red",xlab="distance on graph",ylab="Knet")
	}else{
		##	
		K.obs    <- res$K.obs
		K.perm   <- res$K.perm
		K.quan   <- res$K.quan
		AUK.obs  <- res$AUK.obs
		AUK.perm <- res$AUK.perm
		##
		if(!sequential) par(mfrow=c(1,2))
		## plot 1
		if (ncol(K.quan)==5) lty<-c("dotted","solid","solid","solid","dotted") else lty<-"solid"
		matplot(t(K.quan),type="l",lty=lty,col="grey",lwd=c(2,1,3,1,2),xlab="Distance in network",ylab="K-function",main="Network enrichment analysis")
		polygon(x=c(1:ncol(K.quan),ncol(K.quan):1),y=c(K.quan["5%",],rev(K.quan["95%",])),col="grey",density=20)
		lines(K.obs,col="red",lwd="3")
		## plot 2
		hist(AUK.perm,xlim=c(min(c(.4,min(AUK.perm))),1),col="grey",border="grey",xlab="Area under K-curve",main="AUK: observed v. permuted")
		abline(v=AUK.obs,lwd=2,col="red")
	}
} ## end plot.Knet


####################################################################################
read.delim.BioGRID <- function(path.to.BioGRID.tab.delim.file){
	
	## TO DO: use graph.data.frame
	
	## read graph from file
	biogrid <- read.delim(path.to.BioGRID.tab.delim.file, colClasses ="character", header=TRUE)
	mypast <- function(x)c(paste(x[1],x[3],x[5],x[10],sep="&"),paste(x[2],x[4],x[6],x[11],sep="&"))
	fromTo <- t(apply(biogrid,1,mypast))
	g <- graph.edgelist(fromTo,directed=FALSE)
	
	## set graph attributes
	g <- set.graph.attribute(g,"sourceDatabase","BioGRID")
	g <- set.graph.attribute(g,"sourceFile",path.to.BioGRID.tab.delim.file)
	
	## set node attributes
	node.anno <- sapply(V(g)$name,strsplit,split="&")
	V(g)$geneName  <- sapply(node.anno,function(x)strsplit(strsplit(x[2],"\\|")[[1]][1],":")[[1]][2])
	V(g)$RefSeq    <- NA
	V(g)$entrezID  <- sapply(node.anno,function(x)strsplit(strsplit(x[1],"\\|")[[1]][1],":")[[1]][2])
	V(g)$databaseID <- sapply(node.anno,function(x)strsplit(strsplit(x[1],"\\|")[[1]][2],":")[[1]][2])
	V(g)$name      <- V(g)$geneName
	
	## set edge attributes TODO
	
	##
	return(g)
	} ## end function read.BioGRID.graph


panel.cor <- function(x, y, digits=2, prefix="", cex.cor, scale.r=TRUE,...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    if (scale.r) cex.r <- cex.cor * r else cex.r = 2
    text(0.5, 0.5, txt, cex = cex.r)
}

panel.smoothscatter <- function (x, y, col = par("col"), bg = NA, pch = par("pch"),  cex = 1, col.smooth = "red", span = 2/3, iter = 3, ...) {
    smoothScatter(x,y,add=TRUE)
    #points(x, y, pch = pch, col = col, bg = bg, cex = cex)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
        lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
            col = col.smooth, ...)
}

##
## put GO anno on vertexes
##
GOnodeAnno <- function(nodeNames,llim=10,ulim=300,envir){
	goterms <- ls(envir)
	fct <- function(x,llim,ulim){
		tmp <- nodeNames %in% get(x,envir)
		if (sum(tmp) >=llim & sum(tmp) <= ulim) return(tmp) else return(NULL)
		}
	tmp <- lapply(goterms,fct,llim,ulim)
	names(tmp) <- goterms
	nu <- !sapply(tmp,is.null)
	if (!any(nu)) stop("no overlap within given limits")
	res <- matrix(unlist(tmp[nu]),ncol=sum(nu),byrow=FALSE)
	dimnames(res) <- list(nodeNames,names(tmp[nu]))
	return(res*1)
	}

get.all.vertex.attributes <- function(graph){
	if (!is.igraph(graph)) stop("Not a graph object")
	A <- as.data.frame(graph[[9]][[3]])
	if (is.null(V(g)$name)) rownames(A) <- 0:(vcount(g)-1) else rownames(A) <- V(g)$name
	return(A)
	}

set.vertex.attributes <- function(graph,D){
	if (!is.igraph(graph)) stop("Not a graph object")
	if (!is.data.frame(D)) stop("Not a dataframe")
	for (i in colnames(D)) graph <- set.vertex.attribute(graph,i,value=D[,i])
	return(graph)
	}

my.mget <- function(x,envir){
	xdim <- dim(x)
	y <- as.vector(as.matrix(x))
	tmp <- mget(y,envir,ifnotfound=NA)
	tmp[is.na(tmp)] <- names(tmp)[is.na(tmp)]
	tmp <- sapply(tmp,function(x)x[1]) ## not unlist, but pick first one 
	if (is.null(xdim)) z <- tmp else z<-matrix(tmp,ncol=xdim[2],nrow=xdim[1])
	return(z)
	}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

#