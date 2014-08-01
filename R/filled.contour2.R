#'  @title Plot color filled contours
#'  @description
#'  \code{filled.contour2} is modification by Ian Taylor of the filled.contour function
#'  to remove the key and facilitate overplotting with contour(). Bug fixed by I. Lopez-Coto
#'  @param see \link{filled.contour}

#'  @author  Ian Taylor and I. Lopez-Coto, 2013 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export


filled.contour2 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
  {
    require(field)

    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    on.exit(par(par.orig))
    w <- (3 + mar.orig[2]) * par("csi") * 2.54
    par(las = las)
    mar <- mar.orig
    plot.new()
    par(mar=mar)
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
      stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
      storage.mode(z) <- "double"
    .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                            col = col)
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) 
      box()
    if (missing(plot.title)) 
      title(...)
    else plot.title
    invisible()
  }

#'  @title Plot color filled contours
#'  @description
#'  \code{filledContour2} is a wrapper to filled.contour2 with raster objects
#'  @param see \link{filled.contour}

#'  @author I. Lopez-Coto, 2013 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export
#'  
#'  
filledContour2 <- function (x, y = 1, maxpixels = 1e+05, ...) 
{
  if (nlayers(x) > 1) {
    y <- min(max(1, y), nlayers(x))
    x <- raster(x, y)
  }
  x <- sampleRegular(x, maxpixels, asRaster = TRUE, useGDAL = TRUE)
  X <- xFromCol(x, 1:ncol(x))
  Y <- yFromRow(x, nrow(x):1)
  Z <- t(matrix(getValues(x), ncol = x@ncols, byrow = TRUE)[nrow(x):1, 
                                                            ])
  filled.contour2(x = X, y = Y, z = Z, ...)
}
