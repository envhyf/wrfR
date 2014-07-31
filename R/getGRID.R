#'  @title Get the WRF grid
#'  @description
#'  \code{getGRID} is a function to get the WRF grid and related projection parameters
#'  @param filename input filename
#'  @return a list containing the lat and lon arrays and projection parameters

#'  @author I. Lopez-Coto, 2013 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export
#'  @examples
#'  \dontrun{
#'  WRFgrid <- getGRID(paste(pathWRF,"wrfinput_",domain[d],sep=''))
#'  for (i in length(names(WRFgrid))){
#'  assign(names(WRFgrid)[[i]],WRFgrid[[i]])
#'  } 
#'  }

getGRID <- function(filename)
{
  require('RNetCDF')
  nc <- open.nc(filename)
  lat_WRF <- var.get.nc(nc,"XLAT",start=c(1,1,1), count=c(NA,NA,1))
  lon_WRF <- var.get.nc(nc,"XLONG",start=c(1,1,1), count=c(NA,NA,1))
  cosalpha <- var.get.nc(nc,"COSALPHA",start=c(1,1,1), count=c(NA,NA,1)) #[,,1]
  sinalpha <- var.get.nc(nc,"SINALPHA",start=c(1,1,1), count=c(NA,NA,1)) #[,,1]
  HGT_WRF <- var.get.nc(nc,"HGT",start=c(1,1,1), count=c(NA,NA,1))  
  lat_1 <- att.get.nc(nc,'NC_GLOBAL','TRUELAT1')
  lat_2 <- att.get.nc(nc,'NC_GLOBAL','TRUELAT2')
  lat_0 <- att.get.nc(nc,'NC_GLOBAL','CEN_LAT')
  lon_0 <- att.get.nc(nc,'NC_GLOBAL','CEN_LON')
  map.proj <- att.get.nc(nc,'NC_GLOBAL','MAP_PROJ')
  close.nc(nc)
  
  if (map.proj == 1){  #  "Lambert Conformal"
    WRFproj <- paste('+proj=lcc +lat_1=',lat_1,' +lat_2=',lat_2,' +lat_0=',lat_0,' +lon_0=',lon_0,' +a=6370000 +es=0.0',sep='')
  }  
  if (map.proj == 2){  #  "Polar Stereographic" 
    WRFproj <- paste('+proj=stere',' +lat_ts=',lat_2,' +lat_0=',lat_0,' +lon_0=',lon_0,' +a=6370000 +es=0.0',sep='')
  } 
  
  output <- list(lat_WRF=lat_WRF,lon_WRF=lon_WRF,HGT_WRF=HGT_WRF,WRFproj=WRFproj,cosalpha=cosalpha,sinalpha=sinalpha)
  return(output)
}