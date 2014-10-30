#'  @title Retrieve ISD stations information 
#'  @description
#'  \code{ISDstations} download the station list from the ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ish-history.csv
#'  @return data.frame with ISD stations list
#'  @author I. Lopez-Coto, 2013 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export

ISDstations <- function(){stations <- read.csv('ftp://ftp.ncdc.noaa.gov/pub/data/noaa/isd-history.csv')}

#'  @title Retrieve the ISD data file 
#'  @description
#'  \code{getISD} download the station data from the ftp://ftp.ncdc.noaa.gov/pub/data/noaa/...
#'  @param pathout 
#'  @param year 4 digits year
#'  @param USAF usaf station code (6 digits but it accepts less)
#'  @param WBAN wban station code (5 digits but it accepts less)
#'  @return character vector with one element per raw ISD report
#'  @author I. Lopez-Coto, 2013 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export

getISD <- function(year=2014,USAF = 84490,WBAN=13025,stations=NULL)
{
  require(RCurl)
 meta<-NULL 
if (!is.null(stations)){info<-stations[stations$USAF==USAF & stations$WBAN == WBAN,]; print(info)}

if (nchar(USAF)<6){USAF <- paste(rep('0',6-nchar(USAF)),USAF,sep='')}
if (nchar(WBAN)<5){WBAN <- paste(rep('0',5-nchar(WBAN)),WBAN,sep='')}

url <- paste('http://www1.ncdc.noaa.gov/pub/data/noaa/',year,'/',sep='')
file <- paste(USAF,'-',WBAN,'-',year,'.gz',sep='')

if (url.exists(paste(url,file,sep=''))){
con <- gzcon(url(paste(url,file,sep=''), "r"))
meta <- readLines(con)
close(con)
}

# if (file.exists(file)) {
# command<-paste('rm ',file,sep='')
# system(command)
# }
# command<-paste('wget ',url,file,sep='')
# system(command)
# command<-paste('mv ',file,' ',pathout,sep='')
# system(command)
# command<-paste('gunzip ',pathout,file,sep='')
# system(command)

#return(paste(pathout,file,sep=''))
return(meta)

}


#'  @title Read and reformat the ISD data file 
#'  @description
#'  \code{readISD} read and reformat the ISD data file 
#'  @param stdata character vector with one element per ISD raw report
#'  @param QC \code{logical}
#'  @return data.frame with ISD station data
#'  @author I. Lopez-Coto, 2013 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export
#'  @examples
#'  \dontrun{
#'  stdata<-getISD()
#'  data <- readISD(stdata,QC=FALSE)
#'  } 

readISD <- function(stdata,QC=TRUE){
 
  format.ISD<-read.csv(system.file('data/format_ISD.csv', package = "wrfR"))
  ## ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ISD/ish-format-document.pdf
  ## ftp://ftp.ncdc.noaa.gov/pub/data/noaa/ISD/software/isd_display.pl
  
#  if (file.exists(filename))
#  {
#    stdata<-read.csv(filename)
#    stdata<-as.character(stdata[,])  
    
    df <- data.frame(matrix(0,nrow=length(stdata),ncol=length(format.ISD$field)), stringsAsFactors=F)
    dimnames(df)[[2]]<-format.ISD$field
    
    for (j in 1:(length(format.ISD$field)-1))
    {
      
      df[,j]<-substring(stdata,format.ISD$start[j]+1,format.ISD$start[j]+format.ISD$count[j])
      
    }
    
    if (QC){
    #keep the 'Passed all quality control checks' data  (flag 1 or 5) 
    
    indx<-which(df$sea_levp_flag == 1 | df$sea_levp_flag == 5)
    df<-df[indx,]
    indx<-which(df$wind_dir_flag == 1 | df$wind_dir_flag == 5)
    df<-df[indx,]
    indx<-which(df$wind_speed_flag == 1 | df$wind_speed_flag == 5)
    df<-df[indx,]
    indx<-which(df$air_temp_flag == 1 | df$air_temp_flag == 5)
    df<-df[indx,]
    indx<-which(df$dew_point_flag == 1 | df$dew_point_flag == 5)
    df<-df[indx,]
    }
    
    ## Apply the scaling factor
    
    df$wind_dir<-as.numeric(df$wind_dir)
    df$wind_speed<-as.numeric(df$wind_speed)/10
    df$air_temp<-as.numeric(df$air_temp)/10
    df$dew_point<-as.numeric(df$dew_point)/10
    df$sea_lev_press<-as.numeric(df$sea_lev_press)/10
    df$lat<-as.numeric(df$lat)/1000
    df$long<-as.numeric(df$long)/1000
    df$elev<-as.numeric(df$elev)  
    
    date<-strptime(paste(df$date,df$gmt,sep=''),tz = "GMT", format = "%Y%m%d%H%M")
    
    df<-data.frame(df$usaf_id,df$wban,df$call_letters,df$lat,df$long,df$elev,date,df$wind_speed,df$wind_dir,df$sea_lev_press,df$air_temp,df$dew_point,stringsAsFactors=F)
    
    dimnames(df)[[2]][1:6]<-substring(dimnames(df)[[2]][1:6],4,30)
    dimnames(df)[[2]][8:12]<-substring(dimnames(df)[[2]][8:12],4,30)
    
#  } else {
#    df<-NULL 
#  }
  return(df)
}

