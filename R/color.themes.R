#'  @title Color Themes
#'  @description
#'  \code{color.themes} provides different color schemes
#'  @param theme \code{string} or \code{integer} indicating the desired color scheme c("radioactive","galaxyblue","bw",'pinkpanther',"iceblue",'apocalypse_now','lime','luna','bossanova','trideca')
#'  @param n \code{integer} number of colors
#'  @param rev \code{logical} reverse order
#'  @return a list with a vector of colors (colsch), and 3 independent colors (col1, col2, st.marker.col). In addition, load the library \link{RColorBrewer}
#'  @author I. Lopez-Coto, 2014 (israel.lopez@@dfa.uhu.es / inl@@nist.gov)
#'  @export

color.themes <- function(theme,n,rev){

  library(RColorBrewer)

if (theme == "radioactive" | theme == 1)
{
  colsch <- colorRampPalette(c("green", "black"))(n)
  if (rev){colsch <- rev(colsch)}
  col1<- 'red'
  col2<- 'white'
  st.marker.col <- 'white'
}

if (theme == "galaxyblue"  | theme == 2)
{
  colsch <- colorRampPalette(c("#55FFFF", "grey10"))(n)
  if (rev){colsch <- rev(colsch)}
  col1<- 'orange'
  col2<- 'red'
  st.marker.col <- 'white'
}

if (theme == "bw"  | theme == 3)
{
  colsch <- colorRampPalette(c("white", "black"))(n)
  if (rev){colsch <- rev(colsch)}
  col1<- 'red'
  col2<- 'white'
  st.marker.col <- 'white'
}

if (theme == "pinkpanther"  | theme == 4)
{
  colsch <- colorRampPalette(c("pink", "black"))(n)
  if (rev){colsch <- rev(colsch)}
  col1<- 'pink'
  col2<- 'white'
  st.marker.col <- 'white'
}

if (theme == "iceblue"  | theme == 5)
{
  colsch <- colorRampPalette(c('white', rgb(0,0.184,0.46,1)))(n)
  if (rev){colsch <- rev(colsch)}
  col1<- 'red'
  col2<- 'white'
  st.marker.col <- 'white'
}

if (theme == 'apocalypse_now'  | theme == 6) #"mini_B60D")
{
  colsch <- colorRampPalette(c('white','#D9D6CF', '#CCC9C3', '#BFBDB6', '#B3B0AA', '#FFB60D', '#FFB60D','#F28100', '#F24900'))(n) 
  if (rev){colsch <- rev(colsch)}
  col1<- 'red'
  col2<- 'grey30'
  st.marker.col <- 'grey30'
}

if (theme == 'lime'  | theme == 7) 
{
  colsch <- colorRampPalette(c('#FEE82C', '#549116', '#135C1D', '#01383E', '#000126'))(n)
  if (rev){colsch <- rev(colsch)}
  col1<- 'white'
  col2<- 'white'
  st.marker.col <- 'grey30'
}

if (theme == 'luna'  | theme == 8) 
{
  colsch <- colorRampPalette(c('#F8EEE1', '#C40F0F', '#110D0D', '#C4C2BC', '#9A9894'))(n)
  if (rev){colsch <- rev(colsch)}
  col1<- 'white'
  col2<- 'white'
  st.marker.col <- 'grey30'
}

if (theme == 'bossanova'  | theme == 9) 
{
  colsch <- colorRampPalette(c('#A9CDB7', '#AAC491', '#97985A', '#E9839C', '#D04768'))(n)
  if (rev){colsch <- rev(colsch)}
  col1<- 'white'
  col2<- 'white'
  st.marker.col <- 'grey40'
}

if (theme == 'trideca'  | theme == 10) 
{
  colsch <- colorRampPalette(c('#DDD6D5', '#F9C71F', '#F8A318', '#D83123', '#2F5083'))(n)
  if (rev){colsch <- rev(colsch)}
  col1<- '#2F5083'
  col2<- '#2F5083'
  st.marker.col <-  '#2F5083' #  'grey30'
}


return(list(colsch=colsch,col1=col1,col2=col2,st.marker.col=st.marker.col))
}