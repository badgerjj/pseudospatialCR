# Packages 


using<-function(...) {
  libs<-unlist(list(...))
  req<-unlist(lapply(libs,require,character.only=TRUE))
  need<-libs[req==FALSE]
  n<-length(need)
  if(n>0){
    libsmsg<-if(n>2) paste(paste(need[1:(n-1)],collapse=", "),",",sep="") else need[1]
    print(libsmsg)
    if(n>1){
      libsmsg<-paste(libsmsg," and ", need[n],sep="")
    }
    libsmsg<-paste("The following packages could not be found: ",libsmsg,"\n\r\n\rInstall missing packages?",collapse="")
    if(winDialog(type = c("yesno"), libsmsg)=="YES"){       
      install.packages(need)
      lapply(need,require,character.only=TRUE)
    }
  }
}

using("coda",
      "crawl",
      "doFuture",
      "doRNG",
      "foreach",
      "fuzzyjoin",
      "ggplot2",
      "ggspatial",
      "ks",
      "lubridate",
      "mapview",
      "mgcv",
      "pathroutr",
      "posterior",
      "raster",
      "readxl",
      "sf",
      "sfheaders",
      "spdep",
      "stars",
      "tidyverse",
      "terra",
      "units")

if (!require('progressr')) install.packages('progressr'); library('progressr'); handlers(global = TRUE)
registerDoFuture()


argosDiag2Cov_manual = function(Major, Minor, Orientation){
  a=Major
  b=Minor
  if(any(b<=.Machine$double.eps,na.rm=TRUE)) stop("There are very small (or 0) values for the minor ellipse lengths! These may need to be removed.")
  theta=Orientation
  if(any(theta < 0 | theta > 180, na.rm = TRUE)) stop("Argos diagnostic data orientation outside of [0,180]!")
  if(any(a < 0, na.rm = TRUE)) stop("Argos diagnostic data major axis < 0!")
  if(any(b < 0, na.rm = TRUE)) stop("Argos diagnostic data minor axis < 0!")
  theta = pi*(theta/180)
  k=sqrt(2)
  v1 = (a/k)^2*sin(theta)^2 + (b/k)^2*cos(theta)^2
  v2 = (a/k)^2*cos(theta)^2 + (b/k)^2*sin(theta)^2
  c12 = ((a^2 - b^2)/k^2)*cos(theta)*sin(theta)
  rho = c12/(sqrt(v1)*sqrt(v2))
  check = (v1*v2-c12^2) > 0 
  if(any(rho > 1 | rho < -1, na.rm=TRUE)) stop("Faulty Argos error correlation calculated from 'argosDiag2Cov' function")
  return(data.frame(ln.sd.x=log(sqrt(v1)), ln.sd.y=log(sqrt(v2)), error.corr=rho, diag.check=check))
}

not_within <- function(x,y) !st_within(x,y)


crw2ctmm <- function(data, time.name, silent=TRUE){
  tmp <- data  %>% st_transform(4326) %>% st_coordinates() %>% 
    cbind(data,.) %>% dplyr::select(.data[[time.name]], X, Y, rep) %>% st_drop_geometry() %>% 
    mutate(timestamp=as.character(.data[[time.name]])) %>% 
    rename(
      `location-long` = X,
      `location-lat` = Y,
      `tag-local-identifier` = rep
    )
  if(silent){
    out <- suppressMessages(
      ctmm::as.telemetry(
        tmp, 
        projection=st_crs(data)$proj4string)
    )
  } else{
    out <- ctmm::as.telemetry(
      tmp, 
      projection=st_crs(data)$proj4string
    )
  }
  return(out)
}

get_segments <- function(x, gap=7, time_unit="days"){
  dt <- diff(x) %>% `units<-`(time_unit)
  time_diff <- c(0, dt)
  seg_id <- (time_diff >= gap) %>% {cumsum(.)+1}
  return(seg_id)
}

location_rate <- function(x, time_unit="day", stat=mean, ...){
  out <- table(lubridate::round_date(x$timestamp, time_unit))
  return(stat(out,...))
}

hr_cut <- function(x, cut=0.95){
  if(sum(x)!=1) stop("UD values must be normalized to find home range cut!")
  revx <- rev(sort(x))
  out <- revx[min(which(cumsum(revx)>=cut))]
}

