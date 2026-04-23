rm(list=ls())
library(humboldt)
library(terra)
library(raster)
library(tidyverse)
library(tcltk)

### READ IN  YOUR DISTRIBUTION DATA: 

## read in your distribution data: 
Dist <- read_excel("Dist.csv")


# create a points dataframe: 
points<- Dist %>% select(longitude, latitude)

old_dist<- Dist %>% filter(year<=1950)
new_dist<- Dist %>% filter(year > 1950)

## get just the lat long of your years: 
old_dist<- Dist %>% filter(year<=1950) %>% select(-year)
new_dist<- Dist %>% filter(year > 1950) %>% select(-year)

## READ IN YOUR ENVIRONMENTAL DATA: 
## load rasters: 
Past<- stack("Env_T1.tif")
Current<- stack("Env_T2.tif")




#####                                     ####
########### HUMBOLDT FORMATTING ##############

## get these env data into humboldt format:
### 1) 
##convert one of the rasters to a point dataframe to sample.  Use any raster input.
env.points_Past<-rasterToPoints(Past[[1]], fun=NULL, spatial=TRUE)
env.points_Current<- rasterToPoints(Current[[1]], fun=NULL, spatial=TRUE)

### 2) 
##Extract values to points from rasters
## Past: 
PAST_humb<-data.frame(raster::extract(Past, env.points_Past))
PAST_humb<-cbind(env.points_Past@coords,PAST_humb)

## Current: 
CURRENT_humb<- data.frame(raster::extract(Current, env.points_Current))
CURRENT_humb<- cbind(env.points_Current@coords, CURRENT_humb)

## make sure the variable names are the same: 
colnames(CURRENT_humb)==colnames(PAST_humb)

## get rid of NAs: 
CURRENT<- na.omit(CURRENT_humb)
PAST<-na.omit(PAST_humb)

## get the correct ordinal values for land use
CURRENT<- full_join(CURRENT, Classes)
CURRENT<- CURRENT %>% select(-landUse, -landcover)
PAST<- full_join(PAST, Classes)
PAST<- PAST %>% select(-landUse, -landcover)


### 3) 
### NOW get your points organized correctly: 
old_points<- old_dist
names(old_points)<-c('lon', 'lat')
old_points$taxon_name<- "porcupine_B1950"
old_points<- old_points %>% select(taxon_name, lon, lat)
old_points<- na.omit(old_points)

new_points<- new_dist 
names(new_points)<-c('lon', 'lat')
new_points$taxon_name<- "porcupine_A1950"
new_points<- new_points %>% select(taxon_name, lon, lat)
new_points<- na.omit(new_points)

#### 4) 
## run humboldt and DO RARIFY (SLOWER SLOWER) ## also print all supplemental plots: 
### ful extent: (NOT):
full<-humboldt.doitall(inname="full_extent", env1=PAST, 
                       env2=CURRENT, sp1=old_points, sp2=new_points, rarefy.dist=10, rarefy.units="km", 
                       env.reso=0.04166667, reduce.env=0, reductype="PCA", non.analogous.environments="YES", 
                       correct.env=T, env.trim=T,  env.trim.type="RADIUS", trim.buffer.sp1=500, trim.buffer.sp2=500, 
                       pcx=1, pcy=2, col.env=e.var, e.var=c(3:9), R=100, kern.smooth=1, e.reps=100, b.reps=100, nae="YES",
                       thresh.espace.z=0.0001, p.overlap=T, p.boxplot=T, p.scatter=T, run.silent=T, ncores=16)



##### shared AE ONLY (NDT)
Shared_AE<-humboldt.doitall(inname="shared_espace_ae", env1=PAST, 
                            env2=CURRENT, sp1=old_points, sp2=new_points, rarefy.dist=10, rarefy.units="km", 
                            env.reso=0.04166667, reduce.env=2, reductype="PCA", non.analogous.environments="NO", 
                            correct.env=T, env.trim=T,  env.trim.type="RADIUS", trim.buffer.sp1=500, trim.buffer.sp2=500, 
                            pcx=1, pcy=2, col.env=e.var, e.var=c(3:9), R=100, kern.smooth=1, e.reps=100, b.reps=100, nae="NO",
                            thresh.espace.z=0.0001, p.overlap=T, p.boxplot=T, p.scatter=T, run.silent=T, ncores=8)



## 2: INTERPRETING AND PLOTTING YOUR RESULTS: 

          ### find the points that are outside of your shared niche: 

cm1 <- as.data.frame(zz$scores.sp1)
cm2 <- as.data.frame(zz$scores.sp2)
cm1$env<- "env1"
cm2$env<-"env2"
cm12<- rbind(cm1, cm2)

## get a convex hull of each area: 

CH1 = chull(cm1$Axis1, cm1$Axis2)
CH2 = chull(cm2$Axis1, cm2$Axis2)

# find the intersection: 
library(geometry)
TEST<- intersectn(
  cm1[1:2],
  cm2[1:2],
  tol = 0,
  return.chs = FALSE,
  options = "Tv",
  fp = NULL,
  autoscale = FALSE
)

points(cm12$Axis1, cm12$Axis2)

## now deterine whether those points are in the convex hull: 
cm1$inH<- inhulln(TEST$ch, as.matrix(cm1[1:2]))
cm2$inH<- inhulln(TEST$ch, as.matrix(cm2[1:2]))
## now get only the overlap points: 
cm1_O<- cm1 %>% filter(inH==F)
cm2_O<- cm2 %>% filter(inH==F)


#### now get a plot of what's going on: 
test<- rast("land_Past.tif")
plot(test)
points(cm2_O[3:4])
points(cm1_O[3:4], col='red')
write.csv(cm1_O, "points_OutsideOfHumboldt_PAST_Shared.csv", row.names=F)
write.csv(cm2_O, "points_OutsideofHumboldt_CURRENT_Shared.csv", row.names=F)


          ### plot your results: 
PAST<- zz$scores.sp1
PRESENT<- zz$scores.sp2



### make your own plots in ggplot of humboldt results: 
library(MASS)
past_contour<- kde2d(PAST$Axis1, PAST$Axis2, n=50)
filled.contour(past_contour)

present_contour<- kde2d(PRESENT$Axis1, PRESENT$Axis2, n=50)
filled.contour(present_contour)


### okay, now do it in ggplot: 
ggplot() +
  stat_density_2d(data=PAST, mapping=aes(x = Axis1, y = Axis2, fill = ..level.., alpha=..level..), geom = "polygon")+
  stat_density_2d(data=PAST, mapping=aes(x = Axis1, y = Axis2, alpha=.3), color="dark blue", alpha=.3)+
  scale_fill_gradient(high="dark blue", low="light blue", name="past") + 
  new_scale_fill()+
  stat_density_2d(data=PRESENT, mapping=aes(x = Axis1, y = Axis2, fill = ..level.., alpha=..level..), geom = "polygon")+
  stat_density_2d(data=PRESENT, mapping=aes(x = Axis1, y = Axis2), color="dark red", alpha=.3)+
  scale_fill_gradient(high="dark red", low="pink", name="present")+
  theme_bw()+
  xlab("PC1")+
  ylab("PC2")+
  guides(alpha=FALSE)+
  ggtitle("Species density in E space")



### make a boxplot of PC occupancy: 

library(reshape2)
dat_PAST<- melt(PAST)
dat_PAST$group="Past"

dat_PRESENT<- melt(PRESENT)
dat_PRESENT$group="Present"

dat_BOTH<- rbind(dat_PRESENT, dat_PAST)
dat_BOTH <- dat_BOTH %>% filter(variable=="Axis1" | variable=="Axis2")

ggplot(dat_BOTH, aes(x=variable, y=value, fill=group))+ 
  geom_boxplot()+
  scale_fill_manual(values=c("dark blue", "dark red"))+
  xlab("PC axis")+ 
  ylab("score")+
  coord_flip()


#### do a density plot: 
PC1_dens<- ggplot(dat_BOTH %>% filter(variable=="Axis1"), aes(x=value, fill=group))+
  geom_density(alpha=.4)+
  xlab("PC1")+
  scale_fill_manual(values=c("dark blue", "dark red"))

ggplot(dat_BOTH %>% filter(variable=="Axis2"), aes(x=value, fill=group))+
  geom_density(alpha=.4)+
  xlab("PC2")+
  scale_fill_manual(values=c("dark blue", "dark red"))


### plot variable contribution: 
contrib=zz$pca.cal$co
eigen= zz$pca.cal$eig
library(ade4)
s.corcircle(contrib[, 1:2]/max(abs(contrib[, 1:2])), grid = F)
title(main = "Variable Contribution", sub = paste("PC1 = ", round(eigen[1]/sum(eigen) * 
                                                                    100, 2), "%", "PC2 = ", round(eigen[2]/sum(eigen) * 100, 2), "%"))


