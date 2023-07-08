dir<-"~/pseudospatialCR"
source(file.path(dir, "code","helper.R"))

n.clust<-3
years<-c(1999:2015)
ts_length<-length(years)

for(k in 1:n.clust){
  load(file.path(dir,"data","processed_data","cluster_uds", paste0("cluster_ud",k,".RData")))
}

ref<-st_crs(ud1)

mhi<-st_read(file.path(datadir,"environmental","MHI coastline_USGS/Coastline.shp"))%>%
  `st_crs<-`(4326) %>% drop_units()%>%
  st_transform(crs=ref)%>% 
  sf::st_simplify(preserveTopology=TRUE, dTolerance=250)

buff200 <- st_union(mhi) %>% st_buffer(200000) 
buff200_l <- st_cast(buff200, "LINESTRING")
buff100 <- st_buffer(mhi, 100000) %>% st_union() %>% st_convex_hull()%>%st_transform(crs=ref)

surveys.all<-readRDS(file.path(datadir,"effort", "surveys_1999_2021_PWF_CRC_PIFSC_WDF.RData"))

#checking
unique(st_geometry_type(surveys.all))
hist(st_coordinates(surveys.all)[,1])
hist(st_coordinates(surveys.all)[,2])


colors<-c("CRC" = "cyan4",
          "PIFSC" = "aquamarine", 
          "PIFSC Ships" = "grey80",
          "PWF" = "darkturquoise",
          "WDF" = "cornflowerblue",
          "Other" = "black")

trks<- ggplot()+
  geom_sf(data=surveys.all[surveys.all$source=="PIFSC Ships",],
          aes(color="PIFSC Ships"), alpha=0.05, size=0.05)+
  geom_sf(data=surveys.all[ surveys.all$source=="CRC",], 
          aes(color="CRC"), size=0.05)+
  geom_sf(data=surveys.all[surveys.all$source=="PWF",], 
          aes(color="PWF"), size=0.05)+
  geom_sf(data=surveys.all[surveys.all$source=="PIFSC",], 
          aes(color="PIFSC"), size=0.05)+
  geom_sf(data=surveys.all[surveys.all$source=="Other",], 
          aes(color="Other"), alpha=0.05, size=0.05)+
  geom_sf(data=surveys.all[surveys.all$source=="WDF",], 
          aes(color="WDF"), size=0.05)+
  scale_color_manual(values = colors, 
                     name = "Source")+
  geom_sf(data=mhi, color="grey")+
  labs(title="", x="Longitude", y="Latitude")+
  guides(colour = guide_legend(override.aes = list(size=3)))+
  theme_classic()
trks 
ggsave(filename = "effort.png", dpi = 500)

# 
# s<-surveys.all[surveys.all$year=="2011",]
# trks<- ggplot()+
#   geom_sf(data=s[s$source=="PIFSC Ships",],
#           aes(color="PIFSC Ships"), alpha=0.05, size=0.05)+
#   geom_sf(data=s[ s$source=="CRC",], 
#           aes(color="CRC"), size=0.05)+
#   geom_sf(data=s[s$source=="PWF",], 
#           aes(color="PWF"), size=0.05)+
#   geom_sf(data=s[s$source=="PIFSC",], 
#           aes(color="PIFSC"), size=0.05)+
#   geom_sf(data=s[s$source=="Other",], 
#           aes(color="Other"), alpha=0.05, size=0.05)+
#   geom_sf(data=s[s$source=="WDF",], 
#           aes(color="WDF"), size=0.05)+
#   scale_color_manual(values = colors, 
#                      name = "Source")+
#   geom_sf(data=mhi, color="grey")+
#   labs(title="", x="Longitude", y="Latitude")+
#   guides(colour = guide_legend(override.aes = list(size=3)))+
#   theme_classic()
# trks 
# ggsave(filename = "effort.png", dpi = 500)
# 


effort<-list()

#kernel density variables
extent_vec <- st_bbox(ud1)[c(1,3,2,4)]
cellsize=2000

n_y <- ceiling((extent_vec[4]-extent_vec[3])/cellsize)
n_x <- ceiling((extent_vec[2]-extent_vec[1])/cellsize)

extent_vec[2] <- extent_vec[1]+(n_x*cellsize)-cellsize
extent_vec[4] <- extent_vec[3]+(n_y*cellsize)-cellsize

for(t in 1:ts_length){
  
  effort.t<-surveys.all[surveys.all$year==years[t],]
  unique(effort.t$source) #checking
  
  #make kernel density 
  
  coords <- st_coordinates(effort.t)
  matrix <- kde2d(coords[,1],coords[,2],
                  h=150000, 
                  n = c(n_x,n_y), 
                  lims = extent_vec)
  contour(matrix)
  points(coords)
  
  kd<-raster(matrix)
  crs(kd)<-ref
  
  
  kde<-kd%>%st_as_stars(dx=cellsize, dy=cellsize)
  kde$layer[kde$layer<quantile(kde$layer,0.05,na.rm=T)]<-0
  st_crs(kde)<-ref
  effort[[t]]<-kde%>%st_crop(.,buff100)
  
  
}


#plot kds
t.plot<-13
kd_plot<- ggplot()+
  geom_stars(data=effort[[t.plot]], show.legend=FALSE)+
  scale_fill_viridis(na.value="white") +
  geom_sf(data=st_geometry(surveys.all[surveys.all$year==years[t.plot],]), color="black", size=0.025)+
  geom_sf(data=mhi, color="grey")+
  geom_sf(data=buff100, fill=NA)+
  theme_classic()+
  labs(title=paste("Effort kernel, ",years[t.plot], sep=""), x="", y="")
kd_plot 


saveRDS(effort, file.path(datadir,"processed_data", "effort_1999_2021_PWF_CRC_PIFSC_WDF.rds"))


#crop at 1000m from shore or something ( look at devins code)

ud1_<-st_rasterize(ud1["ud"], dx=cellsize, dy=cellsize)%>%st_crop(.,buff100)
ud1_$ud<-ud1_$ud/sum(ud1_$ud, na.rm=T)
ud1_$ud[ud1_$ud<quantile(ud1_$ud,0.05, na.rm=T)]<-0


ud2_<-st_rasterize(ud2["ud"], dx=cellsize, dy=cellsize)%>%st_crop(.,buff100)
ud2_$ud<-ud2_$ud/sum(ud2_$ud, na.rm=T)
ud2_$ud[ud2_$ud<quantile(ud2_$ud,0.05, na.rm=T)]<-0


ud3_<-st_rasterize(ud3["ud"], dx=cellsize, dy=cellsize)%>%st_crop(.,buff100)
ud3_$ud<-ud3_$ud/sum(ud3_$ud, na.rm=T)
ud3_$ud[ud3_$ud<quantile(ud3_$ud,0.05, na.rm=T)]<-0


ud4_<-st_rasterize(ud4["ud"], dx=cellsize, dy=cellsize)%>%st_crop(.,buff100)
ud4_$ud<-ud4_$ud/sum(ud4_$ud, na.rm=T)
ud4_$ud[ud4_$ud<quantile(ud4_$ud,0.05, na.rm=T)]<-0


#plot uds
ud_plot3<- ggplot()+
  ggspatial::layer_spatial(ud3_)+
  scale_fill_viridis(na.value="white",option="magma") +
  ggspatial::layer_spatial(mhi, fill="gray20", color="gray50")+
  geom_sf(data=buff100, fill=NA)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(title="Normalized Cluster 3 space use kernel (n = 12) ", fill="", x="Longitude", y="Latitude")
ud_plot3 

ud_plot2<- ggplot()+
  ggspatial::layer_spatial(ud2_)+
  scale_fill_viridis(na.value="white",option="magma") +
  ggspatial::layer_spatial(mhi, fill="gray20", color="gray50")+
  geom_sf(data=buff100, fill=NA)+
  theme_classic()+
  theme(legend.position = "none",axis.text.y = element_blank(), axis.text.x = element_blank())+
  labs(title="Normalized Cluster 2 space use kernel (n = 4) ", y="")
ud_plot2

ud_plot1<- ggplot()+
  ggspatial::layer_spatial(ud1_)+
  scale_fill_viridis(na.value="white",option="magma") +
  ggspatial::layer_spatial(mhi, fill="gray20", color="gray50")+
  geom_sf(data=buff100, fill=NA)+
  theme_classic()+
  theme(legend.position = "none",axis.text.x = element_blank())+
  labs(title="Normalized Cluster 1 space use kernel (n = 25) ", x="",y="Latitude")
ud_plot1

ud_plot4<- ggplot()+
  ggspatial::layer_spatial(ud4_)+
  scale_fill_viridis(na.value="white",option="magma") +
  ggspatial::layer_spatial(mhi, fill="gray20", color="gray50")+
  geom_sf(data=buff100, fill=NA)+
  theme_classic()+
  theme(legend.position = "none",axis.text.y = element_blank())+
  labs(title="Normalized Cluster 4 space use kernel (n = 4) ", x="Longtidue",y="")
ud_plot4

library(patchwork)
(ud_plot1+ud_plot2)/(ud_plot3+ud_plot4)

ggsave(filename = "clust_uds_4.png", dpi = 500, width = 9, height = 9, units = "in")



overlap<-matrix(NA, nrow=ts_length, ncol=n.clust)

x1<-c(ud1_$ud)
x2<-c(ud2_$ud)
x3<-c(ud3_$ud)
x4<-c(ud4_$ud)

dxdy<-cellsize*cellsize

for(t in 1:ts_length){
  
  effort.t<-effort[[t]]
  y<-c(effort.t$layer)
  
  tmp<-as.data.frame(cbind(y, x1,x2,x3,x4))
  tmp<-na.omit(tmp)
  
  overlap[t,1]<-sum(sqrt(tmp$y)*sqrt(tmp$x1))*dxdy   
  overlap[t,2]<-sum(sqrt(tmp$y)*sqrt(tmp$x2))*dxdy   
  overlap[t,3]<-sum(sqrt(tmp$y)*sqrt(tmp$x3))*dxdy  
  overlap[t,4]<-sum(sqrt(tmp$y)*sqrt(tmp$x4))*dxdy  
  
  
  
}

overlap

dat<-as.data.frame(overlap)
colnames(dat)<-c("Clust1","Clust2","Clust3","Clust4")
dat$year<-1999:2021
dat<-reshape2::melt(dat, id.vars="year")
colnames(dat)<-c("year","cluster","overlap")

o<-ggplot(dat, aes(x=year, y=overlap))+
  geom_line(aes(colour=cluster))+
  ylab("Effort Overlap")+
  theme_classic()
o


#for(c in 1:n.clust){
# overlap[,c]<-(overlap[,c]-mean(overlap[,c]))/sd(overlap[,c])
#}

write.table(dat,file.path(datadir,"processed_data","overlap_4clust_raw.txt"))

overlap.s<-matrix(NA, nrow=dim(overlap)[1], ncol=dim(overlap)[2])
for(t in 1:ts_length){
  overlap.s[t,]<-(overlap[t,]-mean(overlap[t,]))/sd(overlap[t,])
}

write.table(overlap.s,file.path(datadir,"processed_data","overlap_4clust.txt"))

## Overlap varibale CRC full time series

s<-surveys.all[surveys.all$source=="CRC",]

effort.CRC<-list()

#kernel density variables
extent_vec <- st_bbox(ud1)[c(1,3,2,4)]
cellsize=2000

n_y <- ceiling((extent_vec[4]-extent_vec[3])/cellsize)
n_x <- ceiling((extent_vec[2]-extent_vec[1])/cellsize)

extent_vec[2] <- extent_vec[1]+(n_x*cellsize)-cellsize
extent_vec[4] <- extent_vec[3]+(n_y*cellsize)-cellsize

years<-as.integer(sort(unique(s$year)))
ts_length<-length(years)

for(t in 1:ts_length){
  
  effort.t<-s[s$year==years[t],]
  unique(effort.t$source) #checking
  
  #make kernel density 
  
  coords <- st_coordinates(effort.t)
  matrix <- kde2d(coords[,1],coords[,2],
                  h=150000, 
                  n = c(n_x,n_y), 
                  lims = extent_vec)
  contour(matrix)
  points(coords)
  
  kd<-raster(matrix)
  crs(kd)<-ref
  
  
  kde<-kd%>%st_as_stars(dx=cellsize, dy=cellsize)
  kde$layer[kde$layer<quantile(kde$layer,0.05,na.rm=T)]<-0
  st_crs(kde)<-ref
  effort.CRC[[t]]<-kde%>%st_crop(.,buff100)
  
  
}


overlap.CRC<-matrix(NA, nrow=ts_length, ncol=n.clust)


dxdy<-cellsize*cellsize

for(t in 1:ts_length){
  
  effort.t<-effort.CRC[[t]]
  y<-c(effort.t$layer)
  
  tmp<-as.data.frame(cbind(y, x1,x2,x3,x4))
  tmp<-na.omit(tmp)
  
  overlap.CRC[t,1]<-sum(sqrt(tmp$y)*sqrt(tmp$x1))*dxdy   
  overlap.CRC[t,2]<-sum(sqrt(tmp$y)*sqrt(tmp$x2))*dxdy   
  overlap.CRC[t,3]<-sum(sqrt(tmp$y)*sqrt(tmp$x3))*dxdy  
  overlap.CRC[t,4]<-sum(sqrt(tmp$y)*sqrt(tmp$x4))*dxdy  
  
  
  
}

overlap.CRC

dat<-as.data.frame(overlap.CRC)
colnames(dat)<-c("Clust1","Clust2","Clust3","Clust4")
dat$year<-as.integer(sort(unique(s$year)))
dat<-reshape2::melt(dat, id.vars="year")
colnames(dat)<-c("year","cluster","overlap")

o<-ggplot(dat, aes(x=year, y=overlap))+
  geom_line(aes(colour=cluster))+
  ylab("Effort Overlap, CRC")+
  theme_classic()
o

write.table(dat,file.path(datadir,"processed_data","overlap_4clust_raw_CRC.txt"))

overlap.CRC.s<-matrix(NA, nrow=dim(overlap.CRC)[1], ncol=dim(overlap.CRC)[2])
for(t in 1:ts_length){
  overlap.CRC.s[t,]<-(overlap.CRC[t,]-mean(overlap.CRC[t,]))/sd(overlap.CRC[t,])
}

write.table(overlap.CRC.s,file.path(datadir,"processed_data","overlap_4clust_CRC.txt"))




##############################################################################

### Overlap for old cluster assignments

n.clust<-3
years<-c(2000:2015)
ts_length<-length(years)

for(k in 1:n.clust){
  load(file.path(datadir,"processed_data","old_cluster_uds", paste("old_cluster_ud_",k,".RData", sep="")))
}


ud1<-old_cluster_ud_1
ud2<-old_cluster_ud_2
ud3<-old_cluster_ud_3


rm(old_cluster_ud_1,old_cluster_ud_2,old_cluster_ud_3)
ref<-st_crs(ud1)

mhi<-st_read(file.path(datadir,"MHI coastline_USGS/Coastline.shp"))%>%
  st_transform(crs=ref)%>% 
  sf::st_simplify(preserveTopology=TRUE, dTolerance=250)

bathy <- terra::rast(file.path(datadir, "mhi_mbsyn_bathytopo_1km_v21.nc"))%>% 
  wrap()%>%rast()
bathy <- terra::project(bathy, y="EPSG:3857")

buff200 <- st_union(mhi) %>% st_buffer(200000) 
buff200_l <- st_cast(buff200, "LINESTRING")
buff100 <- st_buffer(mhi, 100000) %>% st_union() %>% st_convex_hull()%>%st_transform(crs=ref)

surveys<-readRDS(file.path(datadir, "effort", "surveys_1999_2021_PWF_CRC_PIFSC_WDF.RData"))%>%
  filter(year < 2016 & year > 1999)%>%
  filter(source == "CRC")

#checking
unique(st_geometry_type(surveys))
hist(st_coordinates(surveys)[,1])
hist(st_coordinates(surveys)[,2])


trks<- ggplot()+
  geom_sf(data=surveys, color="wheat4", size=0.025)+
  geom_sf(data=mhi, color="grey")+
  geom_sf(data=buff100, fill=NA)+
  labs(x="Longitude", y="Latitude")+
  theme_classic()+
  ggsn::north(data=surveys, location="bottomleft")+
  ggsn::scalebar(data=surveys, location = "bottomleft", dist = 50, 
                 dist_unit="nm", transform=FALSE, model="WGS84", st.size=3,
                 st.dist = 0.05, anchor=c(x=-17885000, y=2100000))
trks 

ggsave(filename = "effort_CRC_2000_2015.png", dpi = 400)


effort<-list()

#kernel density variables
extent_vec <- st_bbox(ud1)[c(1,3,2,4)]
cellsize=2000

n_y <- ceiling((extent_vec[4]-extent_vec[3])/cellsize)
n_x <- ceiling((extent_vec[2]-extent_vec[1])/cellsize)

extent_vec[2] <- extent_vec[1]+(n_x*cellsize)-cellsize
extent_vec[4] <- extent_vec[3]+(n_y*cellsize)-cellsize

for(t in 1:ts_length){
  
  effort.t<-surveys[surveys$year==years[t],]
  unique(effort.t$source) #checking
  
  #make kernel density 
  
  coords <- st_coordinates(effort.t)
  matrix <- kde2d(coords[,1],coords[,2],
                  h=150000, 
                  n = c(n_x,n_y), 
                  lims = extent_vec)
  contour(matrix)
  points(coords)
  
  kd<-raster(matrix)
  crs(kd)<-ref
  
  
  kde<-kd%>%st_as_stars(dx=cellsize, dy=cellsize)
  kde$layer[kde$layer<quantile(kde$layer,0.05,na.rm=T)]<-0
  st_crs(kde)<-ref
  effort[[t]]<-kde%>%st_crop(.,buff100)
  
  
}


#plot kds
t.plot<-16
kd_plot<- ggplot()+
  geom_stars(data=effort[[t.plot]], show.legend=FALSE)+
  scale_fill_viridis(na.value="white") +
  geom_sf(data=st_geometry(surveys[surveys$year==years[t.plot],]), color="black", size=0.025)+
  geom_sf(data=mhi, color="grey")+
  geom_sf(data=buff100, fill=NA)+
  theme_classic()+
  labs(title=paste("Effort kernel, ",years[t.plot], sep=""), x="", y="")
kd_plot 



#crop at 1000m from shore or something ( look at devins code)

ud1_<-st_rasterize(ud1["ud"], dx=cellsize, dy=cellsize)%>%st_crop(.,buff100)
ud1_$ud<-ud1_$ud/sum(ud1_$ud, na.rm=T)
ud1_$ud[ud1_$ud<quantile(ud1_$ud,0.05, na.rm=T)]<-0


ud2_<-st_rasterize(ud2["ud"], dx=cellsize, dy=cellsize)%>%st_crop(.,buff100)
ud2_$ud<-ud2_$ud/sum(ud2_$ud, na.rm=T)
ud2_$ud[ud2_$ud<quantile(ud2_$ud,0.05, na.rm=T)]<-0


ud3_<-st_rasterize(ud3["ud"], dx=cellsize, dy=cellsize)%>%st_crop(.,buff100)
ud3_$ud<-ud3_$ud/sum(ud3_$ud, na.rm=T)
ud3_$ud[ud3_$ud<quantile(ud3_$ud,0.05, na.rm=T)]<-0


#plot uds
ud_plot3<- ggplot()+
  ggspatial::layer_spatial(ud3_)+
  scale_fill_viridis(na.value="white",option="magma") +
  ggspatial::layer_spatial(mhi, fill="gray20", color="gray50")+
  geom_sf(data=buff100, fill=NA)+
  theme_classic()+
  labs(title="Normalized Cluster 3 space use kernel (n = 9) ", fill="", x="Longitude", y="Latitude")
ud_plot3 

ud_plot2<- ggplot()+
  ggspatial::layer_spatial(ud2_)+
  scale_fill_viridis(na.value="white",option="magma") +
  ggspatial::layer_spatial(mhi, fill="gray20", color="gray50")+
  geom_sf(data=buff100, fill=NA)+
  theme_classic()+
  theme(legend.position = "none",axis.text.x = element_blank())+
  labs(title="Normalized Cluster 2 space use kernel (n = 4) ", y="Latitude")
ud_plot2

ud_plot1<- ggplot()+
  ggspatial::layer_spatial(ud1_)+
  scale_fill_viridis(na.value="white",option="magma") +
  ggspatial::layer_spatial(mhi, fill="gray20", color="gray50")+
  geom_sf(data=buff100, fill=NA)+
  theme_classic()+
  theme(legend.position = "none",axis.text.x = element_blank())+
  labs(title="Normalized Cluster 1 space use kernel (n = 31) ", x="",y="Latitude")
ud_plot1

library(patchwork)
ud_plot1/ud_plot2/ud_plot3

ggsave(filename = "clust_uds.png", dpi = 500, width = 7, height = 9, units = "in")


overlap<-matrix(NA, nrow=ts_length, ncol=n.clust)

x1<-c(ud1_$ud)
x2<-c(ud2_$ud)
x3<-c(ud3_$ud)


dxdy<-cellsize*cellsize

for(t in 1:ts_length){
  
  effort.t<-effort[[t]]
  y<-c(effort.t$layer)
  
  # overlap[t,1]<-sum(y[x1>0],na.rm = T)*dxdy   # PHR: prob. effort is in cluster 1's UD
  # overlap[t,2]<-sum(y[x2>0],na.rm = T)*dxdy   
  # overlap[t,3]<-sum(y[x3>0],na.rm = T)*dxdy   
  # overlap[t,4]<-sum(y[x4>0],na.rm = T)*dxdy   
  # overlap[t,5]<-sum(y[x5>0],na.rm = T)*dxdy   
  tmp<-as.data.frame(cbind(y, x1,x2,x3))
  tmp<-na.omit(tmp)
  
  overlap[t,1]<-sum(sqrt(tmp$y)*sqrt(tmp$x1))*dxdy   
  overlap[t,2]<-sum(sqrt(tmp$y)*sqrt(tmp$x2))*dxdy   
  overlap[t,3]<-sum(sqrt(tmp$y)*sqrt(tmp$x3))*dxdy  
}

overlap

dat<-as.data.frame(overlap)
colnames(dat)<-c("Clust1","Clust2","Clust3")
dat$year<-years
dat<-reshape2::melt(dat, id.vars="year")
colnames(dat)<-c("year","cluster","overlap")

o<-ggplot(dat, aes(x=year, y=overlap))+
  geom_line(aes(colour=cluster))+
  labs(y="Overlap (not standardized)", x="Year")+
  scale_color_discrete(name = "Social Cluster", labels = c("Cluster 1", "Cluster 2", "Cluster 3"))+
  theme_classic()
o

dat$overlap.s<-(dat$overlap-mean(dat$overlap))/sd(dat$overlap)

o<-ggplot(dat, aes(x=year, y=overlap.s))+
  geom_line(aes(colour=cluster))+
  labs(y="Overlap (standardized)", x="Year")+
  scale_color_manual(values = c("blue4", "cyan3","darkorange"),
                     name = "Social Cluster", labels = c("Cluster 1", "Cluster 2", "Cluster 3"))+
  theme_classic()
o

ggsave(filename = "overlap_3clust.png", dpi = 400)


for(c in 1:n.clust){
  overlap[,c]<-(overlap[,c]-mean(overlap[,c]))/sd(overlap[,c])
}




write.table(overlap,file.path(datadir,"processed_data","overlap_3clust.txt"))

