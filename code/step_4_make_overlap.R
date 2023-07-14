dir<-"~/pseudospatialCR"
source(file.path(dir, "code","helper.R"))

n.clust<-3
years<-c(2000:2015)
ts_length<-length(years)

for(k in 1:n.clust){
  load(file.path(dir,"data","processed_data", paste0("ud_",k,".RData")))
}

ref<-st_crs(ud_1)

load(file.path(dir, "data","habitat_layers.RData"))
load(file=file.path(dir, "data","surveys_CRC_2000_2015.RData"))
surveys<-surveys.CRC%>%
  dplyr::mutate(year = as.integer(year))%>%
  dplyr::select(-source)

#checking
unique(st_geometry_type(surveys))
hist(st_coordinates(surveys)[,1])       #make sure no crazy outliers 
hist(st_coordinates(surveys)[,2])

trks<- ggplot()+
  geom_sf(data=surveys, color="wheat4", size=0.025)+
  geom_sf(data=land, color="grey")+
  geom_sf(data=buff100, fill=NA)+
  labs(x="Longitude", y="Latitude")+
  theme_classic()+
  ggsn::north(data=surveys, location="bottomleft")+
  ggsn::scalebar(data=surveys, location = "bottomleft", dist = 50, 
                 dist_unit="nm", transform=FALSE, model="WGS84", st.size=3,
                 st.dist = 0.05, anchor=c(x=-17885000, y=2100000))
trks 


effort<-list()

#kernel density variables
extent_vec <- st_bbox(ud_1)[c(1,3,2,4)]
cellsize=2000

n_y <- ceiling((extent_vec[4]-extent_vec[3])/cellsize)
n_x <- ceiling((extent_vec[2]-extent_vec[1])/cellsize)

extent_vec[2] <- extent_vec[1]+(n_x*cellsize)-cellsize
extent_vec[4] <- extent_vec[3]+(n_y*cellsize)-cellsize

for(t in 1:ts_length){
  effort.t<-surveys[surveys$year==years[t],]
  
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
  effort[[t]]<-kde
}


#plot kde
t.plot<-13
kd_plot<- ggplot()+
  ggspatial::layer_spatial(data=effort[[t.plot]]%>%st_crop(.,buff100))+
  viridis::scale_fill_viridis(na.value="white") +
  geom_sf(data=surveys[surveys$year==years[t.plot],], color="black", size=0.025)+
  geom_sf(data=land, color="grey")+
  geom_sf(data=buff100, fill=NA)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(title=paste("Effort kernel, ",years[t.plot], sep=""), x="", y="")
kd_plot 


saveRDS(effort, file.path(dir,"data","processed_data", "effort_kdes.RData"))


#crop UDs at 1000m from shore & at upper 95% of probability density (to ease computation)
# and normalize so it sums to 1

ud_1_crop<-st_rasterize(ud_1["ud"], template = effort[[1]])%>%st_crop(.,buff100)
ud_1_crop$ud<-ud_1_crop$ud/sum(ud_1_crop$ud, na.rm=T)           
ud_1_crop$ud[ud_1_crop$ud<quantile(ud_1_crop$ud,0.05, na.rm=T)]<-0


ud_2_crop<-st_rasterize(ud_2["ud"], template = effort[[1]])%>%st_crop(.,buff100)
ud_2_crop$ud<-ud_2_crop$ud/sum(ud_2_crop$ud, na.rm=T)
ud_2_crop$ud[ud_2_crop$ud<quantile(ud_2_crop$ud,0.05, na.rm=T)]<-0


ud_3_crop<-st_rasterize(ud_3["ud"], template = effort[[1]])%>%st_crop(.,buff100)
ud_3_crop$ud<-ud_3_crop$ud/sum(ud_3_crop$ud, na.rm=T)
ud_3_crop$ud[ud_3_crop$ud<quantile(ud_3_crop$ud,0.05, na.rm=T)]<-0


#Compute overlap using Bhattacharya's affinity
overlap<-matrix(NA, nrow=ts_length, ncol=n.clust)

x1<-c(ud_1_crop$ud)
x2<-c(ud_2_crop$ud)
x3<-c(ud_3_crop$ud)

dxdy<-cellsize*cellsize


for(t in 1:ts_length){
  
  effort.t<-effort[[t]]%>%st_crop(.,buff100)
  y<-c(effort.t$layer)
  
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
  ylab("Effort Overlap")+
  theme_classic()
o



write.table(overlap, file.path(dir,"data","processed_data","overlap_3clust.txt"))

