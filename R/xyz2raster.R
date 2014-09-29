#'  @title Rasterize xyz data to a regular lat lon grid
#'  @description
#'  \code{xyz2raster}This function performs an ordinary kriging of the (log) concentration to "predict" the spatial distribution of the input value
#'  @param data data.frame(lon=lon,lat=lat)  ## All the stations coordinates
#'  @param subset data.frame(lon=lon,lat=lat,conc=conc) ##  Subset of stations information
#'  @param log \code{logical} whether or not to perform the fitting in log-scale 
#'  @param proj.utm proj4 string \cr
#'  '+proj=utm +zone=31 +units=km +datum=WGS84 +no_defs' # ETEX  \cr
#'  '+proj=utm +zone=18 +units=km +datum=NAD83 +no_defs' # ANATEX 
#'  @return raster object with the "predicted" field

#'  @author I. Lopez-Coto, 2014 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export


xyz2raster  <- function(data,subset,p.res=100,proj.utm='+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=km +nadgrids=@null +no_defs',log=TRUE){

  
library(sp)
library(rgdal)
library(gstat)
library(raster)
 
# We need to perform the kriging in length units (m, km ...), no in degrees. 
# So we'll need to reproject from lat lon to other projection. 
  
proj.latlon<- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

# Filter duplicate data
subset<-subset[!(duplicated(subset$lon) & duplicated(subset$lat)),]
duplicated(subset$name)

if (log){
  sp.data <- data.frame(lon=subset$lon, lat=subset$lat, conc=log(subset$conc))
} else {
  sp.data <- data.frame(lon=subset$lon, lat=subset$lat, conc=(subset$conc))  
}
#sp.data$conc[!(is.finite(sp.data$conc))]<-0

coordinates(sp.data) <- ~ lon+lat

proj4string(sp.data) = CRS(proj.latlon)

sp.data.proj <- spTransform(sp.data,CRSobj=CRS(proj.utm))

# Create new grid 

g1 <- expand.grid(x = seq(range(data$lon)[1]-3,range(data$lon)[2]+3, length=p.res), y = seq(range(data$lat)[1]-3,range(data$lat)[2]+3, length=p.res))
g1.proj <- spTransform(SpatialPoints(g1,proj=CRS(proj.latlon)),CRSobj=CRS(proj.utm))
r <- raster(ncol=p.res,nrow=p.res) 
extent(r)<-extent(g1.proj)
projection(r)<-projection(g1.proj)
r[]<-0

## make gstat object:
g <- gstat(id="conc", formula=conc ~ 1, data=sp.data.proj)

# create directional variograms at 0, 45, 90, 135 degrees from north (y-axis)
v <- variogram(g, alpha=c(0,45,90,135))

# fit the model
v.fit <- fit.variogram(v, model=vgm(model='Lin' , anis=c(0, 0.5)))

## update the gstat object:
g <- gstat(g, id="conc", model=v.fit )

## perform ordinary kriging prediction:
p <- predict(g, model=v.fit, newdata=g1.proj)

# Reproject back to lat lon
p.ll <- spTransform(p,CRSobj=CRS(proj.latlon))
rr.p <- projectRaster(r,crs=proj.latlon)

## Convert to raster
r.p<-rasterize(p.ll,rr.p)
projection(r.p)

r.p.temp<-raster(ncol=100,nrow=100,crs=projection(r.p),ext=extent(g1))

r.p2<-resample(crop(r.p[[2]],extent(g1)),r.p.temp)

if (log){
  r.out <- exp(r.p2)
} else {
  r.out <- r.p2
}

return(r.out)

}


