#'  @title Rasterize and Reproject 2D fields
#'  @description
#'  \code{wrf2d.raster} reprojects and rasterizes WRF 2D fields to lat lon
#'  @param subset a matrix with the input 2D field 
#'  @param lon_WRF a matrix with the WRF longitudes 
#'  @param lat_WRF a matrix with the WRF latitudes
#'  @param proj.latlon proj4 string witht the desired projection
#'  @param WRFproj proj4 string witht the current projection
#'  @param reproject \code{logical} whether or not to reproject 
#'  @return a reprojected raster object
#'  @author I. Lopez-Coto, 2014 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export
#  @examples
#  \dontrun{
#  WRFgrid <- getGRID(paste(pathWRF,"wrfinput_",domain[d],sep=''))
#  for (i in length(names(WRFgrid))){
#  assign(names(WRFgrid)[[i]],WRFgrid[[i]])
#  } 
#  }

wrf2d.raster <- function(subset,lon_WRF,lat_WRF,proj.latlon,WRFproj,reproject=TRUE)
{
  require(rgdal)
  require(raster)
  sp.data <- data.frame(lon=as.vector(lon_WRF), lat=as.vector(lat_WRF), conc=as.vector(subset))
  
  coordinates(sp.data) <- ~ lon+lat
  proj4string(sp.data) = CRS(proj.latlon)
  
  sp.data.proj <- spTransform(sp.data,CRSobj=CRS(WRFproj))
  
  ## Create raster object
  r <- raster(ncols=ncol(lat_WRF),nrows=nrow(lat_WRF))
  r[] <- 0
  extent(r) <- extent(sp.data.proj)
  projection(r) <- projection(sp.data.proj)   
  r.p <- rasterize(coordinates(sp.data.proj),r,sp.data.proj$conc)
  if (reproject){
    r.p <- projectRaster(r.p,crs=proj.latlon)
  }
  
  return(r.p)
}