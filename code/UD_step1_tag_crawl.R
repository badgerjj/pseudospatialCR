dir<-"~/Dropbox/work/fkw"
source(file.path(dir, "published_scripts/crawl_helper.R"))
source(file.path(dir, "published_scripts/load_packages.R"))

### ~~~~~~~~~~~~~~~~
### Read in data
### ~~~~~~~~~~~~~~~~
datadir<-file.path(dir,"data")
mov <- readr::read_csv(file.path(datadir,"tag data", "PcTag001-083_DAF_for_PIFSC_2023Jan11.csv")) 

#project to mercator
mov <- st_as_sf(mov, coords=c("long","lat")) %>% sf::st_set_crs(4326) 
mov_box <- sf::st_bbox(mov)
mov <- st_transform(mov,3857) %>% rename(timestamp=date)

# Get Bathymetry and shoreline data
land <- sf::read_sf(file.path(datadir, "environmental", dsn ="MHI coastline_USGS/Coastline.shp")) %>% 
  sf::st_transform(3857) %>% 
  sf::st_simplify(preserveTopology=TRUE, dTolerance=250)


ext<-st_bbox(land)%>%
  st_as_sfc%>%
  st_buffer(200000)%>%
  st_bbox()


load(file.path(datadir, "environmental","bathy_raster3-MHI.rdata"))
bathy<-raster::projectRaster(bathy_raster3, crs="EPSG:3857")%>%
  crop(ext)

iso20 <- raster::rasterToContour(bathy, levels=-20) %>% st_as_sf() %>% 
  st_cast("POLYGON") %>% st_make_valid()
iso20 <- st_union(st_geometry(land), iso20) %>% st_union()

ggplot()+
  ggspatial::layer_spatial(bathy)+
  ggspatial::layer_spatial(mov, alpha=0.2, size=0.2, color="lawngreen")+
  ggspatial::layer_spatial(iso20)+
  theme_classic()+ theme(legend.position="none")


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
          argos2diag
)

# argosDiag2Cov = function(Major, Minor, Orientation){
#   a=Major
#   b=Minor
#   if(any(b<=.Machine$double.eps,na.rm=TRUE)) stop("There are very small (or 0) values for the minor ellipse lengths! These may need to be removed.")
#   theta=Orientation
#   if(any(theta < 0 | theta > 180, na.rm = TRUE)) stop("Argos diagnostic data orientation outside of [0,180]!")
#   if(any(a < 0, na.rm = TRUE)) stop("Argos diagnostic data major axis < 0!")
#   if(any(b < 0, na.rm = TRUE)) stop("Argos diagnostic data minor axis < 0!")
#   theta = pi*(theta/180)
#   k=sqrt(2)
#   v1 = (a/k)^2*sin(theta)^2 + (b/k)^2*cos(theta)^2
#   v2 = (a/k)^2*cos(theta)^2 + (b/k)^2*sin(theta)^2
#   c12 = ((a^2 - b^2)/k^2)*cos(theta)*sin(theta)
#   rho = c12/(sqrt(v1)*sqrt(v2))
#   check = (v1*v2-c12^2) > 0 
#   if(any(rho > 1 | rho < -1, na.rm=TRUE)) stop("Faulty Argos error correlation calculated from 'argosDiag2Cov' function")
#   return(data.frame(ln.sd.x=log(sqrt(v1)), ln.sd.y=log(sqrt(v2)), error.corr=rho, diag.check=check))
  

mov <- mov %>% dplyr::mutate(
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

### ~~~~~~~~~~~~~~~~~~~~~~~~
### Nest data by individuals
### ~~~~~~~~~~~~~~~~~~~~~~~~

mov <- mov %>% 
  dplyr::group_by(tag_id, ptt, population, cluster, tag_type, pseu_clust) %>% 
  tidyr::nest() %>% ungroup()%>%
  filter(pseu_clust=="Include",
         population=="MHI")

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Find big gaps and locations rate summary
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mov <- mov %>% mutate(
  data = map(data, ~{
    .x$seg_id_1 = get_segments(.x$timestamp,1)
    .x$seg_id_2 = get_segments(.x$timestamp,2)
    .x
    }),
  num_seg_1 = map_dbl(data, ~{max(.x$seg_id_1)}),
  num_seg_2 = map_dbl(data, ~{max(.x$seg_id_2)}),
  loc_rate = map_dbl(data, ~{location_rate(.x, "hour")}),
  num_locs = map_int(data, nrow),
  num_days = map_int(data, ~{
    lubridate::round_date(.x$timestamp, "days") %>% 
      unique(.) %>% length(.)
    })
)


war### ~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Fit crawl model in parallel
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~

# crwMLE_batch is in 'crawl_helper.R'
plan("multisession", workers=6)
mov$fits <- crwMLE_batch(mov)
plan("sequential")

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Predict tracks and reroute around land
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plan("multisession", workers=6)
mov$pred <- crwPredict_batch(mov$fits, predTime="1 hour", barrier=iso20, vis_graph=vis_graph)
plan("sequential")  


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Sample from track posterior and reroute around land
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plan("multisession", workers=6)
mov$mi_samp <- crwSample_batch(8, fit_list=mov$fits, predTime="1 hour", barrier=iso20, vis_graph=vis_graph)
plan("sequential")  


### ~~~~~~~~~
### Save data 
### ~~~~~~~~~

save(mov, iso20, land, file="processed_data/fkw_telem_data.RData")
