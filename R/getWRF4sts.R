#'  @title Convert and extract WRF data at stations locations
#'  @description
#'  \code{getWRF4sts} extracts WRF default (meteorological variables) + additional variables at stations locations from WRF files. 
#'  @param filenames a vector of WRF filenames (It can be length(filenames) == 1)
#'  @param st.loc dataframe as produced by locWRFst() 
#'  @param vlev number of vertical levels (default = 10)
#'  @param var4d vector of strings containing the name(s) of 4D WRF variable(s) to extract 
#'  @param var3d vector of strings containing the name(s) of 3D WRF variable(s) to extract 
#'  @param masl \code{logical} whether or not the elevation is above sea level
#'  @return list containing the dataframes for each station \cr
#'  c('name','lat','lon','date','u','v','w','ws','wd','p','T','theta','Td',var4d,var3d)   

#'  @author I. Lopez-Coto, 2014 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export

getWRF4sts <- function(filenames,st.loc,var4d=c('QKE'),var3d=c('SWDOWN','PBLH'),vlev=10,masl=FALSE)
{
  require(RNetCDF)
  
  ## Constants
  
  P0 = 1000  #hPa
  kappa = 0.2854
  
  A = 2.53e8 #kPa
  B = 5.42e3 #K (Kelvin)
  e = 0.622 #(approximated from R'/Rv)
  R = 287.05 #J/kg dry air
  
  ### Create variables for output allocation 
  
  varsnum <- length(var3d) + length(var4d)
  
  st.WRF.temp <- data.frame(vector(),vector(),vector(),as.POSIXlt(vector(),tz = "GMT"),vector(),vector(),vector(),vector(),vector(),vector(),vector(),vector(),vector(),matrix(vector(), 0,   varsnum))
  
  dimnames(st.WRF.temp)[[2]] <- c('name','lat','lon','date','u','v','w','ws','wd','p','T','theta','Td',var4d,var3d)
  
  st.WRF.list <- list()
  
  for (i in 1:dim(st.loc)[1]){   
    st.WRF.list[[i]]<-st.WRF.temp #,times=dim(st.loc)[1]))
  }
  
  tind<-1 ## Initialize counter for storing each time stamp 
  
  
  ## Get the WRF grid
  
  WRFgrid <- getGRID(filenames[1])
  for (i in 1:length(names(WRFgrid))){
    assign(names(WRFgrid)[[i]],WRFgrid[[i]]) 
  }
  
  for (f in 1:(length(filenames)))
    
  {
    
    fileWRF<-filenames[f]
    
    nc <-open.nc(fileWRF)
    times <- var.get.nc(nc,"Times")
    date.in.WRF <- strptime(times, tz = "GMT", format = "%Y-%m-%d_%H:%M:%S")
    
    for (time_ind in 1:length(date.in.WRF))
    {
      
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
      
      #WS<-sqrt(Uearth^2 + Vearth^2)
      #WD<-atan2(Uearth,Vearth)*180/(3.1459)+180
      
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
      
      
      # Extract the stations data and interpolate to the correct height (m)
      
      
      for (i in 1:dim(st.loc)[1]){
        
        ## First these fields
        
        if (masl){
          zout <- max(as.numeric(as.vector(st.loc$elev[i])),z[st.loc$lon_indx[i],st.loc$lat_indx[i],1])
        } else {
          zout <- max(as.numeric(as.vector(st.loc$elev[i]))+HGT_WRF[st.loc$lon_indx[i],st.loc$lat_indx[i]],z[st.loc$lon_indx[i],st.loc$lat_indx[i],1])
        }
        
        U.st <- approx(z[st.loc$lon_indx[i],st.loc$lat_indx[i],],Uearth[st.loc$lon_indx[i],st.loc$lat_indx[i],1:dim(z)[3]], xout=zout)$y
        V.st <- approx(z[st.loc$lon_indx[i],st.loc$lat_indx[i],],Vearth[st.loc$lon_indx[i],st.loc$lat_indx[i],1:dim(z)[3]], xout=zout)$y
        w.st <- approx(z[st.loc$lon_indx[i],st.loc$lat_indx[i],],W[st.loc$lon_indx[i],st.loc$lat_indx[i],1:dim(z)[3]], xout=zout)$y
        p.st <- approx(z[st.loc$lon_indx[i],st.loc$lat_indx[i],],p[st.loc$lon_indx[i],st.loc$lat_indx[i],1:dim(z)[3]], xout=zout)$y
        t.st <- approx(z[st.loc$lon_indx[i],st.loc$lat_indx[i],],Temp[st.loc$lon_indx[i],st.loc$lat_indx[i],1:dim(z)[3]], xout=zout)$y
        td.st <- approx(z[st.loc$lon_indx[i],st.loc$lat_indx[i],],Td[st.loc$lon_indx[i],st.loc$lat_indx[i],1:dim(z)[3]], xout=zout)$y
        th.st <- approx(z[st.loc$lon_indx[i],st.loc$lat_indx[i],],TWRF[st.loc$lon_indx[i],st.loc$lat_indx[i],1:dim(z)[3]]+300-273.15, xout=zout)$y
        
        
        ws<-sqrt(U.st^2 + V.st^2)
        wd<-atan2(U.st,V.st)*180/(3.1459)+180
        
        ## Extract 4D variables and interpolate (var4d)
        
        if (length(var4d)>0)
        {
          d4temp <- NULL
          for (v in 1:length(var4d))
          {
            assign(paste(var4d[v],'_temp',sep=''),var.get.nc(nc,var4d[v],start=c(st.loc$lon_indx[i],st.loc$lat_indx[i],1,time_ind), count=c(1,1,vlev-1,1)))
            temp2 <- approx(z[st.loc$lon_indx[i],st.loc$lat_indx[i],], get(paste(var4d[v],'_temp',sep='')), xout=zout)$y
            d4temp <- cbind(d4temp,temp2)
          }
        }       
        
        ## Extract 3D variables (var3d)
        
        if (length(var3d)>0)
        {
          d3temp <- NULL
          for (v in 1:length(var3d))
          {
            assign(paste(var3d[v],'_temp',sep=''),var.get.nc(nc,var3d[v],start=c(st.loc$lon_indx[i],st.loc$lat_indx[i],time_ind), count=c(1,1,1)))
            d3temp <- cbind(d3temp,get(paste(var3d[v],'_temp',sep='')))
          }
        }
        
        ## Store it as a data.frame
        
        st.WRF.temp <- data.frame(st.loc$name[i],st.loc$lat[i],st.loc$lon[i],date.in.WRF[time_ind],U.st,V.st,w.st,ws,wd,p.st,t.st,th.st,td.st,stringsAsFactors=F)
        
        if (length(var4d)>0){st.WRF.temp<-cbind(st.WRF.temp,d4temp)}
        if (length(var3d)>0){st.WRF.temp<-cbind(st.WRF.temp,d3temp)}
        
        st.WRF.list[[i]][tind,] <- st.WRF.temp
        
      }
      
      tind <- tind + 1  ## Increase the global time counter 
      
      print(paste(dim(st.WRF.temp)[[2]]-4, ' variables and ', dim(st.loc)[1] ,' station(s) processed for ', date.in.WRF[time_ind],sep=''))
      
    }  ## End time_ind
    
    close.nc(nc)
    
  } ## End filenames (f)
  
  ## Remove duplicate dates (just in case)
  
  for (i in 1:dim(st.loc)[1]){
    st.WRF.list[[i]]<-st.WRF.list[[i]][!(duplicated(st.WRF.list[[i]]$date)),]
  }
  
  names(st.WRF.list) <- st.loc$name
  
  return(st.WRF.list)
} ## End function
