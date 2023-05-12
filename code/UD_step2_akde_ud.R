dir<-"~/Dropbox/work/fkw/github"
source(file.path(dir,"helper.R"))


### ~~~~~~~~~~~~
### load ud maps
### ~~~~~~~~~~~~
load(file.path(dir,"processed_data","crawl_out.RData"))

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Establish neighborhood structure for modeling
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

grid <- habitat_sf %>% dplyr::select(geometry) %>% st_centroid() %>% st_coordinates()


seg_tbl <- mov %>% group_by(seg_id_2) %>% st_drop_geometry() %>% 
    summarize(start=min(timestamp), end=max(timestamp))
  
ud_<- foreach(j = c(1:length(samp)), .combine='rbind', 
              .packages = c("dplyr","sf"))%do%{
                mi_tmp <- samp[[j]] %>% dplyr::filter(locType=="p") %>% 
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
# plan("sequential")

save(mov, habitat_sf, grid, land, iso20, ud, file=file.path(dir, "processed_data","kde_out.RData"))
