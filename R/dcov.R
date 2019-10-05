wdcor.test <-
function(x, y, G.list=NULL,R=NULL) {
    ## check for valid number of replicates R
    index=1.0
    
    if(is.null(G.list))
    G.list<-c(0)
    
    
    # distance covariance test for multivariate independence
    
    x <- as.matrix(x)
    y <- as.matrix(y)
    G.list<-as.vector(G.list)
    
    n <- nrow(x)
    m <- nrow(y)
    
    if (n != m) stop("Sample sizes must agree")
    if (! (all(is.finite(c(x, y)))))
    stop("Data contains missing or infinite values")
    
    stat <- dcorr <- reps1<-reps2 <- 0
    dcov <- rep(0, 4)
    if (R > 0) reps1 <- rep(0, R)
    if (R > 0) reps2 <- rep(0, R)
    
    pval <- rep(1,2)
    
    dims <- c(n, ncol(x), ncol(y), R,length(G.list))
    
    # dcov = [dCov,dCor,dVar(x),dVar(y)]
    a <- .C("wdCOVtest",
    x = as.double(t(x)),
    y = as.double(t(y)),
    byrow = as.integer(TRUE),
    dims = as.integer(dims),
    gammalist = as.double(G.list),
    index = as.double(index),
    reps1 = as.double(reps1),
    reps2 = as.double(reps2),
    DCOV = as.double(dcov),
    pval = as.double(pval),
    PACKAGE = "wdcor")
    # test statistic is n times the square of dCov statistic
    dcorr <- a$DCOV
    V <- dcorr[[2]]
    names(V) <- "wdcor"
    dataname <- paste("replicates ", R, sep="")
    method <- ifelse(R > 0, "wdCor test of independence",
    "Specify the number of replicates R>0 to perform the test of independence")
    pval <- ifelse (R < 1, NA, a$pval[2])
    e <- list(
    method = method,
    statistic = V,
    p.value = pval,
    replicates = reps2,
    n = n,
    data.name = dataname)
    class(e) <- "htest"
    return(e)
}

.wdcor <-
function(x, y, G.list=NULL) {
    # distance covariance statistic for independence
    # dcov = [dCov,dCor,dVar(x),dVar(y)]   (vector)
    # this function provides the fast method for computing dCov
    # it is called by the dcov and dcor functions
    if(is.null(G.list))
    G.list<-c(0)
    
    index=1.0
    
    x <- as.matrix(x);
    y <- as.matrix(y);
    G.list<-as.vector(G.list);
    
    n <- nrow(x)
    m <- nrow(y)
    if (n != m) stop("Sample sizes must agree")
    if (! (all(is.finite(c(x, y)))))
    stop("Data contains missing or infinite values")
    dims <- c(n, ncol(x), ncol(y),length(G.list))
    idx <- 1:dims[1]
    DCOV <- numeric(4)
    a <- .C("wdCOV",
    x = as.double(t(x)),
    y = as.double(t(y)),
    byrow = as.integer(TRUE),
    dims = as.integer(dims),
    gammalist = as.double(G.list),
    index = as.double(index),
    idx = as.double(idx),
    DCOV = as.double(DCOV),
    PACKAGE = "wdcor")
    return(a$DCOV)
}



wdcor <-
function(x, y, G.list=NULL) {
    # distance correlation statistic for independence
    if(is.null(G.list))
    G.list<-c(0)
    
    index=1.0
    
    x <- as.matrix(x)
    y <- as.matrix(y)
    G.list<-as.vector(G.list)
    
    return(.wdcor(x, y, G.list)[2])
}

