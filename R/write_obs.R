
#'  @title Write observations to littleR format
#'  @description
#'  \code{write_obs} write one observations report to littleR format. This function is based on the example upa.f 
#'   originary from MM5 and currently distributed with obsgrid
#'  @param xlat Latitude
#'  @param xlon Longitude
#'  @param string1 Station description (Code and Name)
#'  @param string2 
#'  @param string3 Observation type ('FM-35 TEMP')
#'  @param string4 
#'  @param ter Terrain elevation
#'  @param kx Number of vertical levels
#'  @param mdate Date 
#'  @param slp Sea level pressure (Pa)
#'  @param p Pressures vector (Pa) 
#'  @param z Heights vector (m) 
#'  @param t Temperatures vector (K) 
#'  @param td Dew point vector (K) 
#'  @param spd Wind speed vector (m/s)
#'  @param dir Wind direction vector (deg)
#'  @param filename filename string
#'  @return write one observations report to littleR format\cr
#'  @author I. Lopez-Coto, 2014 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export
#'  @examples
#'  \dontrun{
#'  
#' kx<-10
#' 
#' p  <- c(1000.,   850.,   700.,   500.,   400.,   300., 250.,   200.,   150.,   100.)
#' z  <- c(100.,  1500.,  3000.,  5500.,  7000.,  9000., 10500., 12000., 13500., 16000.)
#' t  <- c( 14.,     6.,    -4.,   -21.,   -32.,   -45.,  -52.,   -57.,   -57.,   -57.)
#' td <- c( 13.,     3.,    -9.,   -28.,   -41.,   -55.,  -62.,   -67.,   -67.,   -67.)
#' spd<- c(  1.,     3.,     5.,     7.,     9.,    11., 13.,    15.,   17.,     19.)
#' dir<- c(  0.,    30.,    60.,    90.,   120.,   150.,  180.,   210.,  240.,    270.)
#' 
#' p=p*100.
#' t=t+273.15
#' td=td+273.15  
#' 
#' slp <- c(101325.)
#' ter <- 1.
#' xlat <- 22.
#' xlon <- 115.
#' mdate <- '1995073018'
#' 
#' string1 <- '99001  Maybe more site info'
#' string2 <- 'SOUNDINGS FROM ????????? SOURCE'
#' string3 <- 'FM-35 TEMP'
#' string4 <- '    '
#' 
#' filename <- 'test.dat'
#' 
#' 
#' write_obs(xlat,xlon, string1,string2,string3,string4,ter,kx,mdate, 
#'           slp, p, z, t,td,spd,dir,filename)  
#'  }


write_obs <- function(xlat,xlon, string1=NULL , string2=NULL , 
                      string3=NULL , string4=NULL , ter, kx, mdate , 
                      slp, p, z, t,td,spd,dir,filename)
{
  
  if (file.exists(filename)){
    append <- TRUE
  } else {
    append <- FALSE    
  }
  
  date_char<-paste('      ',mdate,'0000',sep='')
  
  if (is.null(string1)) {string1<- ' '}
  if (is.null(string2)) {string2<- ' '}
  if (is.null(string3)) {string3<- ' '}
  if (is.null(string4)) {string4<- ' '}
  
  
  ### Write header report
  
  if (kx > 1){
    is.sounding <- 'T'
  } else {
    is.sounding <- 'F'
  }
  
  header_data <- data.frame(xlat,xlon, string1 , string2 , 
                            string3 , string4 , ter, kx, 0,0,iseq_num=0,0, 
                            'T','F',is.sounding, 
                            -888888, -888888, date_char , 
                            slp,0,-888888.,0, -888888.,0, -888888.,0, -888888.,0, 
                            -888888.,0, 
                            -888888.,0, -888888.,0, -888888.,0, -888888.,0, 
                            -888888.,0, 
                            -888888.,0, -888888.,0)
  
  
  write.mat(m=header_data, c("%20.5f","%20.5f","%40s","%40s","%40s","%40s","%20.5f",rep("%10i",5),
                             rep("%10s",3),"%10i","%10i","%20s",rep(c("%13.5f" , "%7i"),13), "\n"), file = filename,append=append)
  
  ### Write measurements
  
  for (k in 1:kx) {
    meas_data <- t(matrix(c( p[k], 0, z[k],0, t[k],0, td[k],0, spd[k],0, dir[k],0,-888888.,0, -888888.,0,-888888.,0, -888888.,0)))
    write.mat(m=meas_data, c(rep(c("%13.5f" , "%7i"),10), "\n"), file = filename,append=TRUE)
  }
  
  ### Write end measurements
  
  meas_data <- t(matrix(c( -777777.,0, -777777.,0,kx,0, -888888.,0, -888888.,0, -888888.,0,  -888888.,0, -888888.,0, -888888.,0,  -888888.,0)))
  write.mat(m=meas_data, c(rep(c("%13.5f" , "%7i"),10), "\n"), file = filename,append=TRUE) 
  
  ### Write end report
  end_data <- t(matrix(c(kx,0,0)))
  write.mat(m=end_data, rep(c("%7i",  "\n"), c(3, 1)), file = filename,append=TRUE) 
  
}


#'  @title Write fortran format 
#'  @description
#'  \code{write.mat} write fortran formatted string
#'  @param m matrix or data.frame
#'  @param code a vector of codes
#'  @param file file name string
#'  @param ...
#'  @return write formatted output\cr
#'  @author Gabor Grothendieck, 2010 (http://r.789695.n4.nabble.com/write-fortran-td1587119.html)
#'  @export
#'  
write.mat <- function(m, codes, sep = "", ...) {
  s <- do.call(sprintf, unname(c(paste(codes, collapse = ""),
                                 as.data.frame(m))))
  if (length(list(...)) > 0) cat(s, sep = sep, ...) else s
}
