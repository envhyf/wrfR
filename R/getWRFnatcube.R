#'  @title Convert WRF data to 3D cube (WRF_lon, WRF_lat, z) 
#'  @description
#'  \code{getWRFnatcube} extracts WRF default (meteorological variables) + additional variables to 3D cube (WRF_lon, WRF_lat, z)  from WRF files. It conserves the original lat lon
#'  \strong{WARNING:} this function may be slow and require a lot of memory depending on the input file size.
#'  @param nc pointer to the WRF netcdf file as produced by \link{open.nc}
#'  
#'  @param time_ind \code{integer} with the time index within the file being processed
#'  @param vlev number of vertical levels (default = 10)
#'  @param HGT_WRF 2D terrain elevation array as produced by \link{getGRID}
#'  @param lon_WRF 2D WRF longitudes matrix as produced by \link{getGRID}
#'  @param lat_WRF 2D WRF latitudes matrix as produced by \link{getGRID}
#'  @param WRFproj proj4 string witht the current projection as produced by \link{getGRID}
#'  @param var4d vector of strings containing the name(s) of 4D WRF variable(s) to extract 
#'  @param var3d vector of strings containing the name(s) of 3D WRF variable(s) to extract 
#'  @return list containing matrix objects for the variables \cr
#'  c('U3D','V3D','W3D','WS3D','WD3D','P3D','T3D','TH3D','TD3D',var4d,var3d)   

#'  @author I. Lopez-Coto, 2014 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export


getWRFnatcube<-function(nc,time_ind,vlev=59,HGT_WRF,zout,var4d=c('QKE','tracer_10'),var3d=c('PBLH'))
{
  require(RNetCDF)
  
  ## Constants
  
  P0 = 1000  #hPa
  kappa = 0.2854
  
  A = 2.53e8 #kPa
  B = 5.42e3 #K (Kelvin)
  e = 0.622 #(approximated from R'/Rv)
  R = 287.05 #J/kg dry air
  
  proj.latlon<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
  
  ### Read the WRF file
  
  TWRF <- var.get.nc(nc,"T",start=c(1,1,1,time_ind), count=c(NA,NA,vlev,1))
  QVAPOR <- var.get.nc(nc,"QVAPOR",start=c(1,1,1,time_ind), count=c(NA,NA,vlev,1))      
  
  U.stag <- var.get.nc(nc,"U",start=c(1,1,1,time_ind), count=c(NA,NA,vlev,1))
  U <-apply(U.stag,c(2,3),diff)/2+U.stag[1:(dim(U.stag)[1])-1,,]
  
  V.stag <- var.get.nc(nc,"V",start=c(1,1,1,time_ind), count=c(NA,NA,vlev,1))
  V <-aperm(apply(V.stag,c(1,3),diff)/2,c(2,1,3))+V.stag[,1:(dim(V.stag)[2])-1,]
  
  W.stag <- var.get.nc(nc,"W",start=c(1,1,1,time_ind), count=c(NA,NA,vlev+1,1))
  W <-aperm(apply(W.stag,c(1,2),diff)/2,c(2,3,1))+W.stag[,,1:(dim(W.stag)[3]-1)]
  
  Uearth <- 0*U
  Vearth <- 0*V
  
  for (nlev in 1:vlev)
  {
    Uearth[,,nlev] <- U[,,nlev]*cosalpha + V[,,nlev]*sinalpha
    Vearth[,,nlev] <- V[,,nlev]*cosalpha - U[,,nlev]*sinalpha
  }
  
  PH.stag <- var.get.nc(nc,"PH",start=c(1,1,1,time_ind), count=c(NA,NA,vlev+1,1))
  PH <-aperm(apply(PH.stag,c(1,2),diff)/2,c(2,3,1))+PH.stag[,,1:(dim(PH.stag)[3])-1]
  PHB.stag <- var.get.nc(nc,"PHB",start=c(1,1,1,time_ind), count=c(NA,NA,vlev+1,1))
  PHB <-aperm(apply(PHB.stag,c(1,2),diff)/2,c(2,3,1))+PHB.stag[,,1:(dim(PHB.stag)[3])-1]     
  
  P <- var.get.nc(nc,"P",start=c(1,1,1,time_ind), count=c(NA,NA,vlev,1))
  PB <- var.get.nc(nc,"PB",start=c(1,1,1,time_ind), count=c(NA,NA,vlev,1))  
  
  p <- (P+PB)/100
  
  z <- array(0,dim=c(dim(P)[1],dim(P)[2],dim(P)[3]-1))
  zabg<-z
  Temp <- array(0,dim=c(dim(P)[1],dim(P)[2],dim(P)[3]-1))
  Td <- Temp
  nlvl<-dim(P)[3]
  
  for (k in 1:(nlvl-1)){
    z[,,k]<-((((PH[,,k]+PH[,,k+1])/2) + ((PHB[,,k]+PHB[,,k+1]) / 2) )/ 9.81)
    zabg[,,k]<-z[,,k]-HGT_WRF        
    Temp[,,k] = (TWRF[,,k]+300)*((P[,,k]+PB[,,k])/(P0*100))^kappa-273.15
    Td[,,k] = B/log(A*e/(QVAPOR[,,k]*(P[,,k]+PB[,,k])/1000))-273.15
  }
  
  ## Process the data
  
  U3D <- vert.interp(Uearth,zabg,zout)    
  V3D <- vert.interp(Vearth,zabg,zout)
  W3D <- vert.interp(W,zabg,zout) 
  P3D <- vert.interp(p,zabg,zout)
  T3D <- vert.interp(Temp,zabg,zout)
  TD3D <- vert.interp(Td,zabg,zout)        
  TH3D <- vert.interp(TWRF,zabg,zout)+300-273.15 
  
  WS3D <- sqrt(U3D^2 + V3D^2)
  WD3D <- atan2(U3D,V3D)*180/(3.1459)+180
  
 # WD3D <- WS3D
  
#  for (vl in 1:dim(U3D)[3]){
#   WD3D[[vl]] <- atan2(U3D[[1]],V3D[[1]])*180/(3.1459)+180
#  }
  
  st.WRF.temp <-list(U3D,V3D,W3D,WS3D,WD3D,P3D,T3D,TH3D,TD3D)
  varnames <- c('U3D','V3D','W3D','WS3D','WD3D','P3D','T3D','TH3D','TD3D')
  ## Extract 4D variables and interpolate (var4d)
  
  if (length(var4d)>0)
  {
    for (v in 1:length(var4d))
    {
      assign(paste(var4d[v],'_temp',sep=''),var.get.nc(nc,var4d[v],start=c(1,1,1,time_ind), count=c(NA,NA,vlev-1,1)))
      assign(paste(var4d[v],'3D',sep=''),vert.interp(get(paste(var4d[v],'_temp',sep='')),zabg,zout))
      #d4temp <- cbind(d4temp,temp2)
      st.WRF.temp[[length(st.WRF.temp)+1]] <- get(paste(var4d[v],'3D',sep=''))
      varnames<-c(varnames,var4d[v])
    }
  }       
  
  ## Extract 3D variables (var3d)
  
  if (length(var3d)>0)
  {
    d3temp <- NULL
    for (v in 1:length(var3d))
    {
      assign(paste(var3d[v],'_temp',sep=''),var.get.nc(nc,var3d[v],start=c(1,1,time_ind), count=c(NA,NA,1)))
      assign(paste(var3d[v],'2D',sep=''),get(paste(var3d[v],'_temp',sep='')))
      st.WRF.temp[[length(st.WRF.temp)+1]] <- get(paste(var3d[v],'2D',sep=''))
      varnames<-c(varnames,var3d[v])
    }
  }
  
  names(st.WRF.temp) <- varnames
  
  return(st.WRF.temp)
  
}  ## end function