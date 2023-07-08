dir<-"~/pseudospatialCR"
source(file.path(dir,"code","helper.R"))

#dir.create(path=file.path(dir, "processed_data","cluster_uds"))

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### load ud maps & habitat layers
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

load(file.path(dir, "data","kde_data_full_old.RData"))     # all individuals 
load(file.path(dir, "data", "habitat_layers.RData"))

N <- nrow(mov)
 
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Plot each individual KDE UD first
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


buff200 <- st_union(land) %>% st_buffer(200000) 
buff200_l <- st_cast(buff200, "LINESTRING")

for(i  in 1:N){
  tot_hrs <- diff(as.numeric(range(mov$mi_samp[[i]][[1]]$timestamp)))/(3600)
  use <- mov$ud[[i]] %>% as_tibble() %>% 
    mutate(
      geometry = st_geometry(habitat_sf),
    ) %>% st_as_sf() %>% 
    st_intersection(buff200) %>% 
    mutate(ud = ud/sum(ud)) %>% filter(ud>=hr_cut(ud,0.99)) %>% 
    mutate(ud = tot_hrs*ud)
  
   p_use <- ggplot() +
     ggspatial::layer_spatial(use, aes(fill = ud), lwd=0) + 
     ggspatial::layer_spatial(land, fill = "gray20", color="gray50", size=0.1) +
     ggspatial::layer_spatial(buff200_l) +
     scale_fill_viridis_c("Hours", option="C") +
     scale_alpha_continuous(guide = 'none') + 
     theme_void() + 
     ggtitle(paste("ID: ",mov$tag_id[[i]]," / Total hrs: ", tot_hrs," / Cluster: ",mov$cluster[[i]]))
  
   file.nm<-paste0("kde_",mov$tag_id[[i]], ".png", sep="")
   ggsave(p_use, file=file.path(dir, "processed_data", "ud_plots", file.nm), width=6.5, height=4, dpi="retina")
  cat(i, " ")
}


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Average KDEs using meta-analysis
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


use_df <- mov %>% dplyr::select(tag_id, old_cluster) %>% 
  mutate(
    ud = foreach(i=1:N)%do%{
      tot_hrs <- diff(as.numeric(range(mov$mi_samp[[i]][[1]]$timestamp)))/(3600)
      min_ud <- min(mov$ud[[i]]$ud)
      use <- mov$ud[[i]] %>% as_tibble() %>% 
        mutate(
          geometry = st_geometry(habitat_sf),
        ) %>% st_as_sf() %>% 
        mutate(
          ud = ud/sum(ud),
        ) %>% st_centroid() %>% sfheaders::sf_to_df(fill=TRUE)
      use
    }
  )


### ~~~~~~~~~~~~~~~~~~~~~~~~~~
### Social group (cluster) UDs
### ~~~~~~~~~~~~~~~~~~~~~~~~~~


size_n<-apply(table(mov[c("tag_id","old_cluster")]), 2, sum)

nms <- paste0("cluster_ud_", 1:max(use_df$old_cluster))

for(i in 1:max(use_df$old_cluster, na.rm=T)){
  dtmp <- use_df %>% filter(old_cluster==i) %>% 
    unnest(cols=ud) %>% mutate(tag_id = factor(tag_id)) %>% 
    droplevels() %>% group_by(point_id) %>% nest() %>% ungroup()
  
  dtmp <- dtmp %>% mutate(
    ud = map_dbl(data, ~{mean(.x$ud)}),
    cv = map_dbl(data, ~{sd(.x$ud)/sqrt(n())})/ud) %>% 
    dplyr::select(-point_id, -data) %>% mutate(geometry=st_geometry(habitat_sf)) %>% 
    st_as_sf() 
  
  assign(nms[i], dtmp)
  file.nm<-paste0(nms[i],".RData")
  save(list=nms[i], file=file.path(dir, "processed_data","cluster_uds",file.nm))
  
  dtmp <- filter(dtmp, ud>=hr_cut(ud,0.95))
  
  p_use <- ggplot() +
    ggspatial::layer_spatial(dtmp, aes(fill = ud), lwd=0) + 
    ggspatial::layer_spatial(land, fill = "gray20", color="gray50", size=0.1) +
    ggspatial::layer_spatial(buff200_l) +
    scale_fill_viridis_c("PMF", option="C") +
    scale_alpha_continuous(guide = 'none') + 
    theme_void() + ggtitle(paste0("95% contour / Cluster: ",i, ", n = ", size_n[i]))
  file.nm<-paste0("kde_cluster_", i, ".png")
  ggsave(p_use, file=file.path(dir,"processed_data","cluster_uds",file.nm), width=6.5, height=4, dpi="retina")
}


