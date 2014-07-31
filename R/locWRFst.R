#'  @title Localize stations indexes
#'  @description
#'  \code{locWRFst} localizes the lat and lon indexes for each station within the WRF grid
#'  @param fileWRF a string containing the WRF filename
#'  @param st.data a data.frame  [  st.data <- data.frame(st.data$Lat,st.data$Lon,st.data$AGL,st.data$Name)  ]
#'  @return A dataframe containing stations information \cr
#'  dimnames(st.loc)[[1]] = stations names # \cr 
#'  dimnames(st.loc)[[2]] = c('lat','lon','elev','lat_indx','lon_indx')  
#'  @author I. Lopez-Coto, 2013 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export

locWRFst <- function(fileWRF,st.data)
  
{
  
require(raster)
  
  ## Get the WRF grid
  
  WRFgrid <- getGRID(fileWRF)
  for (i in 1:length(names(WRFgrid))){
    assign(names(WRFgrid)[[i]],WRFgrid[[i]]) 
  }
  
  ## Localize stations indexes 
  
  st.loc <- matrix(0,nrow=dim(st.data)[1],ncol=6)
  st.names <- matrix(0,nrow=dim(st.data)[1],ncol=1)
  
  st.loc[,2:4] <- cbind(st.data$Lat,st.data$Lon,st.data$AGL)
  st.loc[,1] <- st.data$Name
  
  
  for (st in 1:dim(st.data)[1])
  {
    
    # Get the station location
    dist <- pointDistance(c(st.data$Lon[st],st.data$Lat[st]),cbind(as.vector(lon_WRF),as.vector(lat_WRF)),longlat=TRUE,r=6370*1e3) # Uses WRF earth radius
    dist <- matrix(dist,ncol=ncol(lon_WRF),nrow=nrow(lon_WRF))
    stloc <- which(dist==min(dist),arr.ind = TRUE)
    
    lat_indx<-stloc[2]
    lon_indx<-stloc[1]
    st.loc[st,5:6]<-c(lat_indx,lon_indx)
  }
  
  dimnames(st.loc)[[2]]<-c('name','lat','lon','elev','lat_indx','lon_indx')
  
  st.loc<-data.frame(st.loc,stringsAsFactors=F)
  class(st.loc[,2:6])<-'numeric'
  
  return(st.loc)
}
