dir<-"~/pseudospatialCR"
source(file.path(dir,"code","helper.R"))


### ~~~~~~~~~~~~
### load ud maps
### ~~~~~~~~~~~~
load(file.path(dir,"out_step1.RData"))

### ~~~~~~~~~~~~~~~~~
### load habitat data
### ~~~~~~~~~~~~~~~~~
load(file.path(dir, "data", "habitat_layers.RData")) 

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Establish neighborhood structure for modeling
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

grid <- habitat_sf %>% dplyr::select(geometry) %>% st_centroid() %>% st_coordinates()


mov <- mov %>% mutate(       
  seg_id_2 = get_segments(timestamp,2))

seg_tbl <- mov %>% group_by(seg_id_2) %>% st_drop_geometry() %>% 
    summarize(start=min(timestamp), end=max(timestamp))
  
ud_<- foreach(j = c(1:length(sims)), .combine='rbind', 
              .packages = c("dplyr","sf"))%do%{
                mi_tmp <- sims[[j]] %>% dplyr::filter(locType=="p") %>% 
                  fuzzyjoin::fuzzy_left_join(
                    seg_tbl, 
                    by=c(timestamp="start",timestamp="end"),
                    match_fun = list(`>=`, `<=`)
                  ) 
              
                  
                  xy <- st_coordinates(mi_tmp) 
                  # n.ess <- mean(effectiveSize(mcmc(xy)))/nrow(xy)
                  f.ess <- posterior::ess_basic(xy)/nrow(xy)
                  xy <- xy[!is.na(mi_tmp$seg_id_2),]
                  n.ess <- nrow(xy)*f.ess
                  H = var(xy)*n.ess^(-1/3)
                  ud_ <- ks::kde(xy, H=H, eval.points=grid, binned=F)$estimate
                  ud_ <- ud_/sum(ud_); ud_
                }
ud<- list(ud = colMeans(ud_), var=apply(ud_, 2, var))

save(mov, sims, pred, fit, vis_graph, grid, ud, 
     file=file.path(dir,"out_step2.RData"))
