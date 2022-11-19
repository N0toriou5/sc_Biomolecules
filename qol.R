# Quality of Life functions
windowsFonts("xkcd" = windowsFont("xkcd-Regular"))
options(java.parameters = "-Xmx8000m")


slice<-function(x){
  return(x[1:5,1:5])
}

# Correlation to P-value
r2p<-function(r,n){
  t<-(r*sqrt(n-2)) / sqrt(1-(r^2))
  p<-2*pt(t,df=n-2,lower=FALSE)
  return(p)
}

# Plot a matrix composed of color entries
colomatplot<-function(topMatrix){
  plot(0,col="white",xlim=c(0,ncol(topMatrix)),ylim=c(0,nrow(topMatrix)),xaxt="n",yaxt="n",type="n",frame.plot=TRUE,xlab="",ylab="",xaxs="i",yaxs="i")
  for(i in (nrow(topMatrix):1)){
    for (j in 1:ncol(topMatrix)){
      rect(
        xleft=j-0.9,
        ybottom=i-1,
        xright=j+0.1,
        ytop=i,
        col=topMatrix[nrow(topMatrix)-i+1,j],
        border=topMatrix[nrow(topMatrix)-i+1,j]
      )
    }
  }
}




# Modified function to separate labels on scatter plot
library(wordcloud)
textplot3<-function(x, y, words, cex = 1, pch = 16, pointcolor = "#FFFFFF00",line.width=1, 
                    new = FALSE, show.lines = TRUE, line.col="#00000033",rstep=0.5,
                    padding=" ",...) 
{
  if (new) {
    plot(x, y, type = "n", ...)
  }
  words<-paste0(padding,words,padding)
  
  lay <- wordlayout(x, y, words, cex, rstep=rstep, ...)
  if (show.lines) {
    for (i in seq_len(length(x))) {
      xl <- lay[i, 1]
      yl <- lay[i, 2]
      w <- lay[i, 3]
      h <- lay[i, 4]
      if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] > 
          yl + h) {
        points(x[i], y[i], pch = pch, col = pointcolor, 
               cex = 0.5)
        nx <- xl + 0.5 * w
        ny <- yl + 0.5 * h
        lines(c(x[i], nx), c(y[i], ny), col = line.col,lwd=line.width)
      }
    }
  }
  text(lay[, 1] + 0.5 * lay[, 3], lay[, 2] + 0.5 * lay[, 4], 
       words, cex = cex, ...)
}



# Only title plot
titplot<-function(title,cex=2,srt=0){
  opar<-par()$mar
  par(mar = c(0,0,0,0))
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x=0.5,y=0.5,title, cex = cex, col = "black",srt=srt)
  par(mar=opar)
}


#' z2p
#'
#' This function gives a gaussian p-value corresponding to the provided Z-score
#'
#' @param z a Z score
#' @return a p-value
#' @examples
#' z<-1.96
#' z2p(z)
#' @export
z2p <- function(z) {
    pnorm(abs(z), lower.tail = FALSE) * 2
}

#' p2z
#'
#' This function gives a gaussian Z-score corresponding to the provided p-value
#' Careful: sign is not provided
#'
#' @param p a p-value
#' @return z a Z score
#' @examples
#' p<-0.05
#' p2z(p)
#' @export
p2z <- function(p) {
    qnorm(p/2, lower.tail = FALSE)
}

#' Stouffer integration of Z scores
#'
#' This function gives a gaussian Z-score corresponding to the provided p-value
#' Careful: sign is not provided
#'
#' @param x a vector of Z scores
#' @return Z an integrated Z score
#' @examples
#' zs<-c(1,3,5,2,3)
#' stouffer(zs)
#' @export
stouffer <- function(x) {
    Z <- sum(x)/sqrt(length(x))
    return(Z)
}


#' kmgformat - Nice Formatting of Numbers
#'
#' This function will convert thousand numbers to K, millions to M, billions
#' to G, trillions to T, quadrillions to P
#'
#' @param input A vector of values
#' @param roundParam How many decimal digits you want
#' @return A character vector of formatted numebr names
#' @examples
#' # Thousands
#' set.seed(1)
#' a<-runif(1000,0,1e4)
#' plot(a,yaxt='n')
#' kmg<-kmgformat(pretty(a))
#' axis(2,at=pretty(a),labels=kmg)
#'
#' # Millions to Billions
#' set.seed(1)
#' a<-runif(1000,0,1e9)
#' plot(a,yaxt='n',pch=20,col=val2col(a))
#' kmg<-kmgformat(pretty(a))
#' axis(2,at=pretty(a),labels=kmg)
#' @export
kmgformat <- function(input, roundParam = 1) {
    signs <- sign(input)
    signs[signs == 1] <- ""
    signs[signs == -1] <- "-"
    absinput <- abs(input)
    output <- c()
    for (i in absinput) {
        if (i < 1000) {
            output <- c(output, i)
        } else if (i < 1e+06) {
            i <- round(i/1000, roundParam)
            i <- paste0(i, "K")
            output <- c(output, i)
        } else if (i < 1e+09) {
            i <- round(i/1e+06, roundParam)
            i <- paste0(i, "M")
            output <- c(output, i)
        } else if (i < 1e+12) {
            i <- round(i/1e+09, roundParam)
            i <- paste0(i, "G")
            output <- c(output, i)
        } else if (i < 1e+15) {
            i <- round(i/1e+12, roundParam)
            i <- paste0(i, "T")
            output <- c(output, i)
        } else if (i < 1e+18) {
            i <- round(i/1e+15, roundParam)
            i <- paste0(i, "P")
            output <- c(output, i)
        } else {
            output <- c(output, i)
        }
    }
    output <- paste0(signs, output)
    return(output)
}

val2col<-function (z, col1 = "navy", col2 = "white", col3 = "red3", 
                   nbreaks = 1000, center = TRUE, rank = FALSE) 
{
  isMatrix <- FALSE
  if (is.matrix(z)) {
    isMatrix <- TRUE
    oriz <- z
  }
  if (is.character(z)) {
    z <- as.numeric(as.factor(z))
  }
  if (rank) {
    z <- rank(z)
  }
  if (center) {
    extreme = round(max(abs(z)))
    breaks <- seq(-extreme, extreme, length = nbreaks)
    z <- z - mean(z)
  }
  else {
    breaks <- seq(min(z), max(z), length = nbreaks)
  }
  ncol <- length(breaks) - 1
  col <- gplots::colorpanel(ncol, col1, col2, col3)
  CUT <- cut(z, breaks = breaks)
  colorlevels <- col[match(CUT, levels(CUT))]
  names(colorlevels) <- names(z)
  if (isMatrix) {
    colormatrix <- matrix(colorlevels, ncol = ncol(oriz), 
                          nrow = nrow(oriz))
    dimnames(colormatrix) <- dimnames(oriz)
    return(colormatrix)
  }
  return(colorlevels)
}

scatter2<-function (x, y, method = "pearson", threshold = 0.01, showLine = TRUE, 
                    bgcol = FALSE, pch = 20, extendXlim = FALSE, ...) {
  common <- intersect(names(x), names(y))
  x <- x[common]
  y <- y[common]
  if (!extendXlim) {
    plot(x, y, pch = pch, ...)
  }
  else {
    plot(x, y, pch = pch, xlim = 1.2 * c(min(x), max(x)), 
         ...)
  }
  cc <- cor.test(x, y, method = method)
  ccp <- signif(cc$p.value, 3)
  cccor <- signif(cc$estimate, 3)
  #mtext(paste0("CC=", cccor, " (p=", ccp, ")"), cex = 0.6)
  if (bgcol) {
    if (cccor >= 0) {
      bgcol <- "#FF000033"
    }
    else {
      bgcol <- "#0000FF33"
    }
    if (ccp > threshold) {
      bgcol <- "#FFFFFF00"
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2], 
         par("usr")[4], col = bgcol)
  }
  grid(col = "gray10")
  if (showLine) {
    lm1 <- lm(y ~ x)
    abline(lm1$coef)
  }
}


