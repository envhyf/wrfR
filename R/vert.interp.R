#'  @title Interpolate to fixed z coordinates
#'  @description
#'  \code{vert.interp} interpolate the heights to fixed z coordinatesn
#'  @param subset 3D array with the variable to interpolate (lon, lat, vert)
#'  @param zabg 3D array with heights above ground for each grid cell (lon, lat, vert) 
#'  @param zout vector with the desired heights
#'  @return 3D array with the interpolated field

#'  @author I. Lopez-Coto, 2014 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export

vert.interp <- function(subset,zabg,zout)
{
  
  tempout <- array(0,dim=dim(zabg))
  for (n in 1:dim(subset)[1])
  { 
    for (m in 1:dim(subset)[2])
    { 
      tempout[n,m,] <- approx(zabg[n,m,],subset[n,m,1:dim(zabg)[3]], xout=zout,rule = 2)$y
    }
  }
  return(tempout)
}