dir<-"~/pseudospatialCR"                       
source(file.path(dir, "code","helper.R"))

### ~~~~~~~~~~~~~~~~
### Read in data
### ~~~~~~~~~~~~~~~~
load(file.path(dir, "data", "mov.RData")) 

#project to mercator
mov <- st_as_sf(mov, coords=c("long","lat")) %>% sf::st_set_crs(4326) 
mov_box <- sf::st_bbox(mov)
mov <- st_transform(mov,3857) %>% rename(timestamp=date)

# Get Bathymetry and shoreline data
load(file.path(dir, "data", "habitat_layers.RData")) 


ggplot()+
  ggspatial::layer_spatial(bathy)+
  geom_sf(data=mov, alpha=0.2, size=0.2, color="lawngreen")+
  geom_sf(data=iso20)+
  theme_classic()+ theme(legend.position="none")+
  ggtitle("Telemetry data points for single \nfalse killer whale deployment in the MHI")


# Visability graph for rerouting paths around barrier (iso20)
iso20_b <- st_buffer(iso20, 10000)
aug_pts <- raster::rasterToContour(bathy, levels=c(-100, -200, -300)) %>% st_as_sf() %>% 
  st_sample(300, type='regular') %>% st_cast("POINT")
aug_pts <- aug_pts[st_intersects(aug_pts, iso20_b, sparse=F)]  
vis_graph <- pathroutr::prt_visgraph(iso20, centroids=TRUE, aug_points=aug_pts)


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Create error structure for all types of tags (LS, KF, FastGPS)
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mov <- cbind(mov,
          argosDiag2Cov_manual(mov$semi_major,
                               mov$semi_minor,
                               mov$ellipse_orient)
        ) %>% dplyr::mutate(
  ln.sd.x = case_when(
    argos_lc %in% c("DP","G") ~ log(50),
    is.na(ln.sd.x) & argos_lc=="L3" ~ log(250),
    is.na(ln.sd.x) & argos_lc=="L2" ~ log(500),
    is.na(ln.sd.x) & argos_lc %in% c("L1","L0","LA","LB") ~ log(1500),
    TRUE ~ ln.sd.x
  ),
  ln.sd.y = case_when(
    argos_lc %in% c("DP","G") ~ log(50),
    is.na(ln.sd.y) & argos_lc=="L3" ~ log(250),
    is.na(ln.sd.y) & argos_lc=="L2" ~ log(500),
    is.na(ln.sd.y) & argos_lc %in% c("L1","L0","LA","LB") ~ log(1500),
    TRUE ~ ln.sd.y
  ),
  error.corr = ifelse(is.na(error.corr), 0, error.corr),
  lc0 = ifelse(argos_lc %in% c("L0","LA","LB"), 1, 0),
  lcA = ifelse(argos_lc %in% c("LA","LB"), 1, 0),
  lcB = ifelse(argos_lc=="LB", 1, 0)
)


### ~~~~~~~~~~~~~~~
### Fit crawl model
### ~~~~~~~~~~~~~~~

temp_dat <- mov %>% dplyr::arrange(timestamp)
if(temp_dat$loc_type[1] == "Least squares"){
  err.model <- list(
    x =  ~ 0 + ln.sd.x + lc0 + lcA + lcB#,
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
  fit <- crawl::crwMLE(
    mov.model = ~1, err.model = err.model, data = temp_dat, Time.name="timestamp",
    fixPar = fixPar, constr = constr, theta = theta,
    control = list(maxit=2000), initialSANN = list(maxit=1500, temp=10),
    attempts=10, method = "L-BFGS-B")
)
fit 


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Predict tracks and reroute around land
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pred <- crawl::crwPredict(fit, predTime="1 hour", return.type="flat") %>% 
  crawl::crw_as_sf(ftype="POINT", locType="p") %>% 
    pathroutr::prt_trim(iso20)
fix <- pathroutr::prt_reroute(pred, iso20, vis_graph, blend=FALSE)
pred <- pathroutr::prt_update_points(fix, pred) 


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Sample from track posterior and reroute around land
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

simObj <- crawl::crwSimulator(fit, parIS = 0, predTime="15 min")
sims <- foreach(j=1:8)%do%{
  samp <- crawl::crwPostIS(simObj, fullPost = FALSE) %>% 
    crawl::crw_as_sf(ftype="POINT", locType="p") %>% 
    pathroutr::prt_trim(iso20)
  samp_fix <- pathroutr::prt_reroute(samp, iso20, vis_graph, blend=FALSE)
  samp <- pathroutr::prt_update_points(samp_fix, samp) %>% dplyr::mutate(rep=j)
  samp
}


### ~~~~~~~~~~~~
### Save data 
### ~~~~~~~~~~~~

save(mov, fit, sims, pred, vis_graph, file=file.path(dir, "out_step1.RData"))
