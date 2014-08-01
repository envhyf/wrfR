#'  @title Pressure to Sea Level Pressure
#'  @description
#'  \code{p2slp} convert from pressure to sea level pressure 
#'  @param p Pressures vector (hPa)
#'  @param elev Station height above sea level (m)
#'  @param Temp Temperature vector (C) (Use to be the 12 hours mean)
#'  @param n \code{integer} number of elements to be used on the Temp running average
#'  @return pressure corrected to sea level
#'  @author I. Lopez-Coto, 2013 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export

p2slp  <- function(p,elev,Temp,n)
{
  Tcap = filter(Temp,rep(1/n,n), sides=1) + 273.15
  pseal = p*exp(elev/(29.3*Tcap))
  return(pseal)
} 


#'  @title Sea Level Pressure to Pressure
#'  @description
#'  \code{slp2p} convert from pressure to sea level pressure 
#'  @param pseal Pressures vector (hPa)
#'  @param elev Station height above sea level (m)
#'  @param Temp Temperature vector (C) (Use to be the 12 hours mean)
#'  @param n \code{integer} number of elements to be used on the Temp running average
#'  @return pressure corrected to sea level
#'  @author I. Lopez-Coto, 2013 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export

slp2p  <- function(pseal,elev,Temp,n)
{
  Tcap = filter(Temp,rep(1/n,n), sides=1) + 273.15
  p = pseal*exp(-elev/(29.3*Tcap))
  return(p)
} 
