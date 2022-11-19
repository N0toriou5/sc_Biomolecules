library(DESeq2)


#' vst - Variance Stabilizing Transformation (based on DESeq2)
#'
#' This function applies Variance-Stabilizing Transformation on RNASeq counts.
#' Plus, it adds white gaussian noise to the data, in order to break ties for further differential expression analysis.
#'
#' @param rawcounts An integer matrix, with genes as rows and samples as columns
#' @return A numeric matrix
#' @export
vst <- function(rawcounts) {
    set.seed(1)
    if(!any(class(rawcounts)=="matrix")){
        stop("x must be an integer matrix object", call.=F)
    }
    coldata<-as.data.frame(colnames(rawcounts))
    coldata[,1]<-as.factor(coldata[,1])
    colnames(coldata)<-"samplename"
    x<-DESeqDataSetFromMatrix(countData=rawcounts,
                              colData=coldata,
                              design=~samplename
    )
    #
    x<-varianceStabilizingTransformation(x,blind=TRUE)
    x<-assay(x)
    #
    tmp <- apply(x, 2, function(x) {
        x <- sort(unique(x))
        x <- cbind(x[1:(length(x)-1)], x[2:length(x)])
        x <- cbind(x[, 1], sqrt(f.rvar.na(x)))
        return(x)
    })
    tmp <- cbind(unlist(lapply(tmp, function(x) x[, 1]), use.names=F), unlist(lapply(tmp, function(x) x[, 2]), use.names=F))
    tmp1 <- smooth.spline(tmp[, 1], tmp[, 2], spar=.5)
    tmp[tmp[, 1]>tmp1$x[which.min(tmp1$y)], 2] <- 0
    tmp1 <- smooth.spline(tmp[, 1], tmp[, 2], spar=.5)
    return(x+rnorm(length(x))*predict(tmp1, x)$y)
}


#' Variance of rows for arrays with NA values
#'
#' This function computes the variance by rows ignoring NA values
#'
#' @param x Numeric matrix
#' @return Numeric vector
#' @export
f.rvar.na <- function(x) {
    ave <- as.vector(f.rmean.na(x))
    pos <- which(is.na(x))
    largo <- f.rlength.na(x)
    x[pos] <- rep(ave,ncol(x))[pos]
    (x-ave)^2 %*% rep(1,ncol(x))/(largo-1)
}

#' Mean of rows for arrays with NA values
#'
#' This function compute the mean by rows ignoring NA values
#'
#' @param x Numeric matrix
#' @return Numeric vector
#' @export
f.rmean.na <- function(x) {
    largo <- f.rlength.na(x)
    x[is.na(x)] <- 0
    res <- x %*% rep(1,ncol(x)) / largo
    names(res) <- rownames(x)
    res
}
#' Length of rows for arrays with NA values
#'
#' This function report the length of rows of a matrix ignoring the NA values
#'
#' @param x Matrix
#' @return Integer indicating the number of non-NA rows
#' @export
f.rlength.na <- function(x) {
    r <- x/x
    r[x==0] <- 1
    r[!is.finite(r)] <- 0
    r %*% rep(1,ncol(r))
}