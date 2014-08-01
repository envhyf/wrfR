#'  @title Extract WRF data at one station location 
#'  @description
#'  \code{getCONC4st} extracts 4D WRF variables at ONE station location from a WRF file. The purpose of this function is to extract a time serie as fast as possible
#'  @param filename a string with the WRF filename 
#'  @param lon_WRF 2D WRF longitudes matrix as produced by \link{getGRID}
#'  @param lat_WRF 2D WRF latitudes matrix as produced by \link{getGRID}
#'  @param lon station longitude
#'  @param lat station latitude 
#'  @param agl station height above ground level 
#'  @param var vector of strings containing the name(s) of 4D WRF variable(s) to extract 
#'  @param vlev number of vertical levels (default = 10)
#'  @return list containing objects, mostly rasters, for the variables \cr
#'  c('U3D','V3D','W3D','WS3D','WD3D','P3D','T3D','TH3D','TD3D',var4d,var3d)   

#'  @author I. Lopez-Coto, 2014 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export



getCONC4st <- function(filename,lon_WRF,lat_WRF,lon,lat,agl,var=c('CO2_BIO','CO2_ANT'), vlev=10)
  
  ## Funtion to extract the concentration of GHGs from WRF for the provided station/tower
  # Inputs:
  
  # filename -> string with the filename to process
  # lon_WRF -> 2D Array with WRF longitudes
  # lat_WRF -> 2D Array with WRF latitudes
  # lon -> station longitude
  # lat -> station latitude 
  # agl -> station height above ground level
  # var -> Vector of tracers to extract (c('CO2_BIO','CO2_ANT'))
# vlev -> vertical levels to process

# Output: ...


{
  require(raster)
  
  # Get the station location
  dist <- pointDistance(c(lon,lat),cbind(as.vector(lon_WRF),as.vector(lat_WRF)),longlat=TRUE,r=6370*1e3) # Uses WRF earth radius
  dist <- matrix(dist,ncol=ncol(lon_WRF),nrow=nrow(lon_WRF))
  stloc <- which(dist==min(dist),arr.ind = TRUE)
  
  nc <-open.nc(filename)
  
  times_temp <- var.get.nc(nc,"Times")
  ## Needed to get the real heights 
  PH <- var.get.nc(nc,"PH",start=c(stloc[1],stloc[2],1,1), count=c(1,1,vlev,NA))
  PHB <- var.get.nc(nc,"PHB",start=c(stloc[1],stloc[2],1,1), count=c(1,1,vlev,NA))
  P <- var.get.nc(nc,"P",start=c(stloc[1],stloc[2],1,1), count=c(1,1,vlev,NA))
  PB <- var.get.nc(nc,"PB",start=c(stloc[1],stloc[2],1,1), count=c(1,1,vlev,NA))
  HGT <- var.get.nc(nc,"HGT",start=c(stloc[1],stloc[2],1), count=c(1,1,1))
  ## Extract the tracers
  
  for (i in 1:length(var))
  {
    assign(paste(var[i],'_temp',sep=''),var.get.nc(nc,var[i],start=c(stloc[1],stloc[2],1,1), count=c(1,1,vlev,NA)))
  }
  close.nc(nc)
  
  nlvl<-dim(P)[1]
  z <- array(0,dim=c(nlvl,dim(P)[2]))
  
  
  for (k in 1:(nlvl-1)){
    z[k,]<-((((PH[k,]+PH[k+1,])/2) + ((PHB[k,]+PHB[k+1,]) / 2) )/ 9.81)
  }
  
  zout <- max(agl+HGT,z[1,])
  
  output<-NULL
  temp<-array(0,dim(z)[2])
  
  for (i in 1:length(var))
  {
    
    for (j in 1:dim(z)[2])
    {
      temp[j] <- approx(z[,j],get(paste(var[i],'_temp',sep=''))[,j], xout = zout)$y 
    }
    
    assign(var[i],temp) # save the interpolated tracer    
    temp[]<-0  # clean the variable 
    
    output <- cbind(output,get(var[i])) 
  }
  
  output<-data.frame(output)
  
  dates <- strptime(times_temp, tz = "GMT", format = "%Y-%m-%d_%H:%M:%S") 
  
  output<-cbind(dates,output)
  dimnames(output)[[2]]<-c('dates',var)
  
  return(output) 
  
}
