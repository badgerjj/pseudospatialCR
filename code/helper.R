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
    if(utils::winDialog(type = c("yesno"), libsmsg)=="YES"){       
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
      "MASS",
      "mgcv",
      "mvtnorm",
      "pathroutr",
      "plotKML",
      "posterior",
      "postpack",
      "raster",
      "readxl",
      "rgdal",
      "R2jags",
      "sf",
      "sfheaders",
      "shinystan",
      "sp",
      "spdep",
      "stars",
      "tidyverse",
      "terra",
      "truncnorm",
      "units",
      "viridis")

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

crwMLE_batch <- function(tdata){
  p <- progressor(nrow(tdata))
  fits <- foreach(i=1:nrow(tdata), .packages="sf") %dorng%{
    temp_dat <- tdata$data[[i]] %>% dplyr::arrange(timestamp)
    
    # Establish fixed parameter values
    if(temp_dat$loc_type[1] == "Least squares"){
      err.model <- list(
        x =  ~ 0 + ln.sd.x + lc0 + lcA + lcB#,
        # y =  ~ 0 + ln.sd.y + lc0 + lcA + lcB,
        # rho =  ~ error.corr
      )
      fixPar <- c(1,NA,NA,NA,NA,NA)
      constr <- list(
        lower=c(rep(0,3), -Inf, -Inf), 
        upper=c(rep(Inf,3), Inf, Inf)
      )
      theta <- c(rep(log(1.2),3),9,0)
    } else {
      err.model <- list(
        x =  ~ 0 + ln.sd.x,
        y =  ~ 0 + ln.sd.y,
        rho =  ~ error.corr)
      fixPar <- c(1,1, NA,NA)
      constr <- list(lower=-Inf, upper=Inf)
      theta <- c(9,0)
    }
    
    # Fit ctcrw model
    suppressMessages(
      out <- crawl::crwMLE(
        mov.model = ~1, err.model = err.model, data = temp_dat, Time.name="timestamp",
        fixPar = fixPar, constr = constr, theta = theta,
        control = list(maxit=2000), initialSANN = list(maxit=1500, temp=10),
        attempts=10, method = "L-BFGS-B")
    )
    p()
    out 
  }
  return(fits)
}



crwPredict_batch <- function(fit_list, predTime, barrier, vis_graph){
  p <- progressor(length(fit_list))
  pred <- foreach(i=1:length(fit_list), .packages=c("sf","dplyr"))%dorng%{
    pred <- crawl::crwPredict(fit_list[[i]], predTime=predTime, return.type="flat") %>% 
      crawl::crw_as_sf(ftype="POINT", locType="p") %>% 
      # mutate(
      #   tmp_seg1 = zoo::na.locf(seg_id),
      #   tmp_seg2 = zoo::na.locf(seg_id, fromLast=TRUE),
      #   seg_id = ifelse(tmp_seg1!=tmp_seg2, NA, tmp_seg1)
      # ) %>% select(-tmp_seg1, -tmp_seg2) %>% 
      pathroutr::prt_trim(iso20)
    fix <- pathroutr::prt_reroute(pred, barrier, vis_graph, blend=FALSE)
    pred <- pathroutr::prt_update_points(fix, pred) 
    p()
    pred
  }
  return(pred)
}


crwSample_batch <- function(size, fit_list, predTime, barrier, vis_graph){
  p <- progressor(length(fit_list))
  slist <- foreach(i=1:length(fit_list), .packages=c("sf","dplyr"))%dorng%{
    simObj <- crawl::crwSimulator(fit_list[[i]], parIS = 0, predTime="15 min")
    out <- foreach(j=1:size)%do%{
      samp <- crawl::crwPostIS(simObj, fullPost = FALSE) %>% 
        crawl::crw_as_sf(ftype="POINT", locType="p") %>% 
        pathroutr::prt_trim(iso20)
      samp_fix <- pathroutr::prt_reroute(samp, iso20, vis_graph, blend=FALSE)
      samp <- pathroutr::prt_update_points(samp_fix, samp) %>% dplyr::mutate(rep=j)
      samp
    }
    p()
    out
  }
  return(slist)
}


