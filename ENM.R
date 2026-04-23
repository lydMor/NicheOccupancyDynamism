rm(list=ls())

# variables for porcs:
# load your packages: 
library(dismo)
library(tidyverse)
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
gc()
library(rJava)
library(terra)
library(openxlsx)

#points:######
Dist<- read.xlsx("E.dorsatum_clean.xlsx")
Dist2<- Dist %>% filter(!(is.na(latitude)))
# create a points dataframe: 
points<- Dist2 %>% dplyr::select(longitude, latitude)
points_spthin <- spThin::thin.algorithm(points[,1:2], thin.par = 0.5, reps = 5)
# check how  many points you have left: 
dim(points_spthin[[1]])
## 10101 
# create the points dataframe: 
points_thin<- data.frame(points_spthin[1])
#reload env data because SpatRasters don't save as RDS files:#######
#These classes hold a C++ pointer to the data and they cannot be directly saved to a ".Rds" file or
#used in cluster computing. They cannot be recovered from a saved R session either. See wrap or
#writeRaster to work around that limitation."


## crop your raster to the porcupine extent: 
xrange <- c(min(Dist2$longitude) - 1, max(Dist2$longitude) + 1)
yrange <- c(min(Dist2$latitude) - 1, max(Dist2$latitude) + 1)
vars<- rast("~/Documents/porcupines/rast/AllEnvVars_2.5.tif")
land = rast("~/Documents/porcupines/rast/global_LULC_2015_correctResolution.tif")
puma = rast("~/Documents/porcupines/allModels/puma_BEST_R1spThin_kFold_porcRange.asc")
puma = subst(puma, NaN, NA)
unique(values(puma)[is.na(values(puma))])
names(puma) = "puma"
crs(puma) = crs(vars)
plot(puma)
fisher = rast("~/Documents/porcupines/allModels/fisher_BEST_R1spThin_kFold_porcRange.asc")
fisher = subst(fisher, NaN, NA)
unique(values(fisher)[is.na(values(fisher))])
names(fisher) = "fisher"
crs(fisher) = crs(vars)
unique(land)
landUse = subst(land, 1:7, c("Water",
                             "Forest",
                             "Grassland",
                             "Barren",
                             "Cropland",
                             "Urban",
                             "PermSnowIce"))
unique(landUse)
names(landUse) = "land_use"
Rasts1 = terra::subset(vars, subset=c("bio1",
                               "bio2",
                               "bio3",
                               "bio4",
                               "bio5",
                               "bio6",
                               "bio7",
                               "bio8",
                               "bio9",
                               "bio10",
                               "bio11",
                               "bio12",
                               "bio13",
                               "bio14",
                               "bio15",
                               "bio16",
                               "bio17",
                               "bio18",
                               "bio19",
                               "alt"))
Rasts1 = c(Rasts1, landUse)
Rasts1<- crop(Rasts1, c(xrange, yrange))
SpatRast1 = Rasts1
Rasts1 = raster::stack(SpatRast1)
r = ratify(Rasts1$land_use)
rat <- levels(r)[[1]]
rat$landcover <- c("Water",
                   "Forest",
                   "Grassland",
                   "Barren",
                   "Cropland",
                   "Urban",
                   "PermSnowIce")
rat$code <- c(1:7)
levels(r) <- rat
x <- deratify(r, "landcover")
is.factor(x)

Rasts1$land_use = x
is.factor(Rasts1$land_use)




#now get 2, 3


SpatRast2<- terra::subset(SpatRast1, subset=c('alt', 
                                 'bio8', 
                                 'bio4', 
                                 'bio12', 
                                 'bio7', 
                                 'bio9', 
                                 'bio6', 
                                 "land_use"))


SpatRast2_noLU<- terra::subset(SpatRast1, subset=c('alt', 
                                      'bio8', 
                                      'bio4', 
                                      'bio12', 
                                      'bio7', 
                                      'bio9', 
                                      'bio6'))

SpatRast3<- terra::subset(SpatRast1, subset=c('bio3', 'bio9', 'bio18', 'bio15', 'bio1', 'bio13', 'alt', 'land_use'))

is.factor(SpatRast2$land_use)
is.factor(SpatRast3$land_use)

SpatRast2P = c(SpatRast2, puma)
SpatRast2F = c(SpatRast2, fisher)
SpatRast2PF = c(SpatRast2P, fisher)
SpatRast2P_noLU = c(SpatRast2_noLU, puma)
SpatRast2F_noLU = c(SpatRast2_noLU, fisher)
SpatRast2PF_noLU = c(SpatRast2P_noLU, fisher)
SpatRast3P = c(SpatRast3, puma)
SpatRast3F = c(SpatRast3, fisher)
SpatRast3PF = c(SpatRast3P, fisher)

is.factor(SpatRast3PF$land_use)

#maxent doesn't work on a SpatRaster, only a *Raster
Rasts2 = raster::stack(SpatRast2)
Rasts2_noLU = raster::stack(SpatRast2_noLU)
Rasts_three = raster::stack(SpatRast3)
Rasts2P = raster::stack(SpatRast2P)
Rasts2F = raster::stack(SpatRast2F)
Rasts2PF = raster::stack(SpatRast2PF)
Rasts2P_noLU = raster::stack(SpatRast2P_noLU)
Rasts2F_noLU = raster::stack(SpatRast2F_noLU)
Rasts2PF_noLU = raster::stack(SpatRast2PF_noLU)
Rasts_threeP = raster::stack(SpatRast3P)
Rasts_threeF = raster::stack(SpatRast3F)
Rasts_threePF = raster::stack(SpatRast3PF)


bg_thin <-  randomPoints(Rasts1, nrow(points_thin), p = points) 
#Run the models::::########
## 2) REDUCED ENVIRONMENTAL VARIABLES: (Rasts2)
# a) spatially thinned points: 
fullModel2A_list = list()
eval2A_list = list()
for(i in 1:10){
  fold <- kfold(points, k=5) # witholding a 20% sample for testing
  points_test_thin <- points_thin[fold == 1, ]
  points_train_thin <- points_thin[fold != 1, ]
  model_tmp<- maxent(Rasts2, points_train_thin, removeDuplicates=T)
  fullModel2A_list[[i]] = model_tmp
  eval_tmp =evaluate(p=points_test_thin, a=bg_thin, model=model_tmp, x=Rasts2)
  eval2A_list[[i]] = eval_tmp
  print(i)
}
saveRDS(fullModel2A_list, "~/Documents/porcupines/allModels/porc_fullModel2A_list.RData")
saveRDS(eval2A_list, "~/Documents/porcupines/allModels/porc_Eval2A_list.RData")


#2aP 
fullModel2AP_list = list()
eval2AP_list = list()
for(i in 1:10){
  fold <- kfold(points, k=5) # witholding a 20% sample for testing
  points_test_thin <- points_thin[fold == 1, ]
  points_train_thin <- points_thin[fold != 1, ]
  model_tmp<- maxent(Rasts2P, points_train_thin, removeDuplicates=T)
  fullModel2AP_list[[i]] = model_tmp
  eval_tmp =evaluate(p=points_test_thin, a=bg_thin, model=model_tmp, x=Rasts2P)
  eval2AP_list[[i]] = eval_tmp
  print(i)
}
saveRDS(fullModel2AP_list, "~/Documents/porcupines/allModels/porc_fullModel2AP_list.RData")
saveRDS(eval2AP_list, "~/Documents/porcupines/allModels/porc_Eval2AP_list.RData")

#2AF = FISHER
fullModel2AF_list = list()
eval2AF_list = list()
for(i in 1:10){
  fold <- kfold(points, k=5) # witholding a 20% sample for testing
  points_test_thin <- points_thin[fold == 1, ]
  points_train_thin <- points_thin[fold != 1, ]
  model_tmp<- maxent(Rasts2F, points_train_thin, removeDuplicates=T)
  fullModel2AF_list[[i]] = model_tmp
  eval_tmp =evaluate(p=points_test_thin, a=bg_thin, model=model_tmp, x=Rasts2F)
  eval2AF_list[[i]] = eval_tmp
  print(i)
}

saveRDS(fullModel2AF_list, "~/Documents/porcupines/allModels/porc_fullModel2AF_list.RData")
saveRDS(eval2AF_list, "~/Documents/porcupines/allModels/porc_Eval2AF_list.RData")


#2APF = BOTH
fullModel2APF_list = list()
eval2APF_list = list()
for(i in 1:10){
  fold <- kfold(points, k=5) # witholding a 20% sample for testing
  points_test_thin <- points_thin[fold == 1, ]
  points_train_thin <- points_thin[fold != 1, ]
  model_tmp<- maxent(Rasts2PF, points_train_thin, removeDuplicates=T)
  fullModel2APF_list[[i]] = model_tmp
  eval_tmp =evaluate(p=points_test_thin, a=bg_thin, model=model_tmp, x=Rasts2PF)
  eval2APF_list[[i]] = eval_tmp
  print(i)
}
saveRDS(fullModel2APF_list, "~/Documents/porcupines/allModels/porc_fullModel2APF_list.RData")
saveRDS(eval2APF_list, "~/Documents/porcupines/allModels/porc_Eval2APF_list.RData")


## 2_noLU) REDUCED ENVIRONMENTAL VARIABLES, NO LAND USE: (Rasts2)
# a) spatially thinned points: 
fullModel2A_noLU_list = list()
eval2A_noLU_list = list()
for(i in 1:10){
  fold <- kfold(points, k=5) # witholding a 20% sample for testing
  points_test_thin <- points_thin[fold == 1, ]
  points_train_thin <- points_thin[fold != 1, ]
  model_tmp<- maxent(Rasts2_noLU, points_train_thin, removeDuplicates=T)
  fullModel2A_noLU_list[[i]] = model_tmp
  eval_tmp =evaluate(p=points_test_thin, a=bg_thin, model=model_tmp, x=Rasts2_noLU)
  eval2A_noLU_list[[i]] = eval_tmp
  print(i)
}

saveRDS(fullModel2A_noLU_list, "~/Documents/porcupines/allModels/porc_fullModel2A_noLU_list.RData")
saveRDS(eval2A_noLU_list, "~/Documents/porcupines/allModels/porc_Eval2A_noLU_list.RData")

#2AP NO LU = PUMAS NO LAND USE
fullModel2AP_noLU_list = list()
eval2AP_noLU_list = list()
for(i in 1:10){
  fold <- kfold(points, k=5) # witholding a 20% sample for testing
  points_test_thin <- points_thin[fold == 1, ]
  points_train_thin <- points_thin[fold != 1, ]
  model_tmp<- maxent(Rasts2P_noLU, points_train_thin, removeDuplicates=T)
  fullModel2AP_noLU_list[[i]] = model_tmp
  eval_tmp =evaluate(p=points_test_thin, a=bg_thin, model=model_tmp, x=Rasts2P_noLU)
  eval2AP_noLU_list[[i]] = eval_tmp
  print(i)
}
saveRDS(fullModel2AP_noLU_list, "~/Documents/porcupines/allModels/porc_fullModel2AP_noLU_list.RData")
saveRDS(eval2AP_noLU_list, "~/Documents/porcupines/allModels/porc_Eval2AP_noLU_list.RData")


#2AP NO LU = FISHERS NO LAND USE
fullModel2AF_noLU_list = list()
eval2AF_noLU_list = list()
for(i in 1:10){
  fold <- kfold(points, k=5) # witholding a 20% sample for testing
  points_test_thin <- points_thin[fold == 1, ]
  points_train_thin <- points_thin[fold != 1, ]
  model_tmp<- maxent(Rasts2F_noLU, points_train_thin, removeDuplicates=T)
  fullModel2AF_noLU_list[[i]] = model_tmp
  eval_tmp =evaluate(p=points_test_thin, a=bg_thin, model=model_tmp, x=Rasts2F_noLU)
  eval2AF_noLU_list[[i]] = eval_tmp
  print(i)
}
saveRDS(fullModel2AF_noLU_list, "~/Documents/porcupines/allModels/porc_fullModel2A_noLUF_list.RData")
saveRDS(eval2AF_noLU_list, "~/Documents/porcupines/allModels/porc_Eval2AF_noLU_list.RData")

#2APF NO LU = BOTH PREDS, NO LAND USE
fullModel2APF_noLU_list = list()
eval2APF_noLU_list = list()
for(i in 1:10){
  fold <- kfold(points, k=5) # witholding a 20% sample for testing
  points_test_thin <- points_thin[fold == 1, ]
  points_train_thin <- points_thin[fold != 1, ]
  model_tmp<- maxent(Rasts2PF_noLU, points_train_thin, removeDuplicates=T)
  fullModel2APF_noLU_list[[i]] = model_tmp
  eval_tmp =evaluate(p=points_test_thin, a=bg_thin, model=model_tmp, x=Rasts2PF_noLU)
  eval2APF_noLU_list[[i]] = eval_tmp
  print(i)
}
saveRDS(fullModel2APF_noLU_list, "~/Documents/porcupines/allModels/porc_fullModel2APF_noLU_list.RData")
saveRDS(eval2APF_noLU_list, "~/Documents/porcupines/allModels/porc_Eval2APF_noLU_list.RData")





## 3) REDUCED ENVIRONMENTAL VARIABLES: (Rasts_three)

#a) spatially thinned points
fullMode3A_list = list()
eval3A_list = list()
for(i in 1:10){
  fold <- kfold(points, k=5) # witholding a 20% sample for testing
  points_test_thin <- points_thin[fold == 1, ]
  points_train_thin <- points_thin[fold != 1, ]
  model_tmp<- maxent(Rasts_three, points_train_thin, removeDuplicates=T)
  fullMode3A_list[[i]] = model_tmp
  eval_tmp =evaluate(p=points_test_thin, a=bg_thin, model=model_tmp, x=Rasts_three)
  eval3A_list[[i]] = eval_tmp
  print(i)
}
saveRDS(fullMode3A_list, "~/Documents/porcupines/allModels/porc_fullMode3A_list.RData")
saveRDS(eval3A_list, "~/Documents/porcupines/allModels/porc_Eval3A_list.RData")

#3aP 
fullModel3AP_list = list()
eval3AP_list = list()
for(i in 1:10){
  fold <- kfold(points, k=5) # witholding a 20% sample for testing
  points_test_thin <- points_thin[fold == 1, ]
  points_train_thin <- points_thin[fold != 1, ]
  model_tmp<- maxent(Rasts_threeP, points_train_thin, removeDuplicates=T)
  fullModel3AP_list[[i]] = model_tmp
  eval_tmp =evaluate(p=points_test_thin, a=bg_thin, model=model_tmp, x=Rasts_threeP)
  eval3AP_list[[i]] = eval_tmp
  print(i)
}

saveRDS(fullModel3AP_list, "~/Documents/porcupines/allModels/porc_fullModel3AP_list.RData")
saveRDS(eval3AP_list, "~/Documents/porcupines/allModels/porc_Eval3AP_list.RData")


#3AF = FISHER
fullModel3AF_list = list()
eval3AF_list = list()
for(i in 1:10){
  fold <- kfold(points, k=5) # witholding a 20% sample for testing
  points_test_thin <- points_thin[fold == 1, ]
  points_train_thin <- points_thin[fold != 1, ]
  model_tmp<- maxent(Rasts_threeF, points_train_thin, removeDuplicates=T)
  fullModel3AF_list[[i]] = model_tmp
  eval_tmp =evaluate(p=points_test_thin, a=bg_thin, model=model_tmp, x=Rasts_threeF)
  eval3AF_list[[i]] = eval_tmp
  print(i)
}
saveRDS(fullModel3AF_list, "~/Documents/porcupines/allModels/porc_fullModel3AF_list.RData")
saveRDS(eval3AF_list, "~/Documents/porcupines/allModels/porc_Eval3AF_list.RData")


#3APF = BOTH
fullModel3APF_list = list()
eval3APF_list = list()
for(i in 1:10){
  fold <- kfold(points, k=5) # witholding a 20% sample for testing
  points_test_thin <- points_thin[fold == 1, ]
  points_train_thin <- points_thin[fold != 1, ]
  model_tmp<- maxent(Rasts_threePF, points_train_thin, removeDuplicates=T)
  fullModel3APF_list[[i]] = model_tmp
  eval_tmp =evaluate(p=points_test_thin, a=bg_thin, model=model_tmp, x=Rasts_threePF)
  eval3APF_list[[i]] = eval_tmp
  print(i)
}
saveRDS(fullModel3APF_list, "~/Documents/porcupines/allModels/porc_fullModel3APF_list.RData")
saveRDS(eval3APF_list, "~/Documents/porcupines/allModels/porc_Eval3APF_list.RData")




#save.image("~/Documents/porcupines/ENM_noMasking_porc2_iterativeKfold.RData")


fullModel2A_list[[1]]


#find best model######################
listoflists = list(eval2A_list,
                   eval2AP_list,
                   eval2AF_list,
                   eval2APF_list,
                   eval2A_noLU_list,
                   eval2AP_noLU_list,
                   eval2AF_noLU_list,
                   eval2APF_noLU_list,
                   eval3A_list,
                   eval3AP_list,
                   eval3AF_list,
                   eval3APF_list)
names(listoflists) = c("2A",
                       "2AP",
                       "2AF",
                       "2APF",
                       "2A_noLU",
                       "2AP_noLU",
                       "2AF_noLU",
                       "2APF_noLU",
                       "3A",
                       "3AP",
                       "3AF",
                       "3APF")
dfList = list()
x = expand.grid(1:length(listoflists),1:length(listoflists[[1]])) %>% arrange(Var1)

for(k in 1:nrow(x)){
  j = x$Var1[k]
  i = x$Var2[k]
  tmp = listoflists[[j]][[i]]
  df_tmp = data.frame(model = names(listoflists)[j],
                      iter = i,
                      numPresences = tmp@np, #num presences
                      numAbsences = tmp@na, #num absences
                      corCoef = tmp@cor[["cor"]], #Correlation coefficient
                      p_cor = tmp@pcor, #p-value for correlation coefficient
                      AUC = tmp@auc, #Area under the receiver operator (ROC) curve)
                      p_AUC = NA,
                      max_TPR_plus_TNR_at = tmp@t[which.max(tmp@TPR + tmp@TNR)] ) #the probability threshold at which our model maximizes the True Positive Rate and the True Negative Rate
  
  if(length(tmp@pauc) > 0 ){
    df_tmp$p_AUC = tmp@pauc #p-value for the AUC (for the Wilcoxon test W statistic
  }
  dfList[[k]] = df_tmp
}

library(dplyr)
df = do.call("bind_rows", dfList)
means = df %>% group_by(model) %>% summarise(meanAUC = mean(AUC)) %>% arrange(desc(meanAUC))
means[1,]
#best model is 2PF

#rerun best model(s) on full set of points########
bestModel <- maxent(Rasts2PF, points_thin, removeDuplicates=T)
Pred_bestModel<- predict(bestModel, Rasts2PF)
secondBest <- maxent(Rasts_threeF, points_thin, removeDuplicates=T)
Pred_secondBest <- predict(secondBest, Rasts_threeF)

plot(Pred_bestModel)
plot(Pred_secondBest)

saveRDS(bestModel, "~/Documents/porcupines/allModels/porc_bestModel_ENM.RData")
saveRDS(secondBest, "~/Documents/porcupines/allModels/porc_secondBestModel_ENM.RData")
writeRaster(Pred_bestModel, "~/Documents/porcupines/allModels/porc_BEST_R2PFspThin_kFold.asc", format = 'ascii', overwrite=T)
writeRaster(Pred_secondBest, "~/Documents/porcupines/allModels/porc_secondBest_R3FspThin_kFold.asc", format = 'ascii', overwrite=T)
#save.image("~/Documents/porcupines/ENM_noMasking_porc3_iterativeKfold.RData")

#future porc######
rm(list=ls())
load("~/Documents/porcupines/ENM_noMasking_porc3_iterativeKfold.RData")
Rasts_best = Rasts2PF
Rasts_secondBest = Rasts_threeF
rm(list= ls()[!(ls() %in% c('bestModel','Rasts_secondBest', 'Rasts_best' ,'yrange', 'xrange', 'Pred_secondBest', 'secondBest', 'Pred_bestModel'))])

#bioclim:
f = list.files("~/Documents/porcupines/rast/microc6", full.names = T)
bioClimListAll = lapply(f, function(x){
  t <- rast(x) 
})
names(bioClimListAll) = gsub(".tif", "", list.files("~/Documents/porcupines/rast/microc6", full.names = F))

#preds:
f = list.files("~/Documents/porcupines/FutureProjectionsNew/puma", pattern = "R1",full.names = T)
puma = lapply(f, function(x){
  t <- rast(x) 
})
names(puma) = gsub(".asc", "", stringr::word(list.files("~/Documents/porcupines/FutureProjectionsNew/puma", pattern = "R1",full.names = F), -2,-1, "_"))
for(i in 1:length(puma)){
  puma[[i]] = subst(puma[[i]], NaN, NA)
  print(is.nan(unique(values(puma[[i]])[is.na(values(puma[[i]]))])))
  names(puma[[i]]) = "puma"
  crs(puma[[i]]) = crs(bioClimListAll[[1]]$bio1)
}

f = list.files("~/Documents/porcupines/FutureProjectionsNew/fisher", pattern = "R1",full.names = T)
fisher = lapply(f, function(x){
  t <- rast(x) 
})
names(fisher) = gsub(".asc", "", stringr::word(list.files("~/Documents/porcupines/FutureProjectionsNew/fisher", pattern = "R1",full.names = F), -2,-1, "_"))
for(i in 1:length(fisher)){
  fisher[[i]] = subst(fisher[[i]], NaN, NA)
  print(is.nan(unique(values(fisher[[i]])[is.na(values(fisher[[i]]))])))
  names(fisher[[i]]) = "fisher"
  crs(fisher[[i]]) = crs(bioClimListAll[[1]]$bio1)
}


#land use:
a =rast("~/Documents/porcupines/rast/global_SSP5_RCP85_2080_correctResolution.tif")
b =rast("~/Documents/porcupines/rast/global_SSP5_RCP85_2100_correctResolution.tif")
c =rast("~/Documents/porcupines/rast/global_SSP1_RCP26_2080_correctResolution.tif")
d =rast("~/Documents/porcupines/rast/global_SSP1_RCP26_2100_correctResolution.tif")
res(a)
res(d)
LUList = list(a,b,c,d)
names(LUList) = c("SSP585_2080", "SSP585_2100", "SSP126_2080", "SSP126_2100")
for(i in 1:length(LUList)){
  names(LUList[[i]]) = "land_use"
}
rm(a,b,c,d,f)

alt = rast("~/Documents/porcupines/rast/RCP85_SSP5_2080_global.tif")$alt

scenarios = c("SSP585_2080", "SSP585_2100", "SSP126_2080", "SSP126_2100")

#get them into a new list, with each scenario:
scenarioListStacks = list()
scenarioList = list()
for(i in 1:length(scenarios)){
  bioclim =bioClimListAll[[which( names(bioClimListAll)==scenarios[i])]]
  landUse = LUList[[which( names(LUList)==scenarios[i])]]
  rat = data.frame(ID=1:7, land_use = c("Water",
                                          "Forest",
                                          "Grassland",
                                          "Barren",
                                          "Cropland",
                                          "Urban",
                                          "PermSnowIce"))
  landUse = as.factor(landUse)
  levels(landUse) = rat
  
  p = puma[[which( names(puma)==scenarios[i])]]
  f = fisher[[which( names(fisher)==scenarios[i])]]
  scen_stack = c(crop(bioclim, c(xrange, yrange)), 
                 crop(landUse, c(xrange, yrange)),
                 crop(alt, c(xrange, yrange)),
                 p, f)
  scenarioListStacks[[i]] = scen_stack
  names(scenarioListStacks)[i] = scenarios[i]
  
  rast = raster::stack(scen_stack)
  r = ratify(rast$land_use)
  rat <- levels(r)[[1]]
  rat
  rat$land_use <- c("Water",
                     "Forest",
                     "Grassland",
                     "Barren",
                     "Cropland",
                     "Urban",
                     "PermSnowIce")
  rat$code <- c(1:7)
  levels(r) <- rat
  x <- deratify(r, "land_use")
  if(!is.factor(x)){print("x is not a factor!")}
  rast$land_use = x
  if(!is.factor(rast$land_use)){print("land use is not a factor!")}
  
  scenarioList[[i]] = rast
  names(scenarioList)[i] = scenarios[i]
  print(i)
}
names(scenarioList[[i]])
#make sure it got pumas and fishers right....
puma[[which( names(puma)==scenarios[1])]]@ptr[["range_min"]]
puma[[which( names(puma)==scenarios[2])]]@ptr[["range_min"]]

rm(list= ls()[!(ls() %in% c('bestModel',
                            'Rasts_secondBest',
                            'Rasts_best' ,
                            'scenarioListStacks',
                            'scenarioList',
                            'Pred_secondBest',
                            'secondBest',
                            'Pred_bestModel',
                            'scenarios'))])

####BEST: Rasts2PF#
#futureList_best = list()
futureProjList_best = list()
#futureList_secondBest = list()
futureProjList_secondBest = list()
for(i in 1:length(scenarioList)){
  bigRast = scenarioList[[which( names(scenarioList)==scenarios[i])]]
  if(!is.factor(bigRast$land_use)){print("land use is not a factor!")}
  
  #best model:
  bestR = raster::subset(bigRast, names(Rasts_best))
  if(!all(sort(names(Rasts_best)) == sort(names(bestR)))){"variables don't match!"}
  #futureList_best[[i]] = bestR
  futureProj1<- predict(bestModel, bestR)
  futureProjList_best[[i]] = futureProj1
  writeRaster(futureProj1, paste0("~/Documents/porcupines/FutureProjectionsNew/porc_Future_best_", 
                                 scenarios[i],
                                 ".asc"), format = 'ascii', overwrite = T)
  
  #second best model:
  secondBestR = raster::subset(bigRast, names(Rasts_secondBest))
  if(!all(sort(names(Rasts_secondBest)) == sort(names(secondBestR)))){"variables don't match!"}
  #futureList_secondBest[[i]] = secondBestR
  futureProj2<- predict(secondBest, secondBestR)
  futureProjList_secondBest[[i]] = futureProj2
  writeRaster(futureProj2, paste0("~/Documents/porcupines/FutureProjectionsNew/porc_Future_secondBest_", 
                                 scenarios[i],
                                 ".asc"), format = 'ascii', overwrite = T)

  print(i)
  
  
}
#names(futureList_best) = scenarios
names(futureProjList_best) = scenarios
#names(futureList_secondBest) = scenarios
names(futureProjList_secondBest) = scenarios


#current:
plot(Pred_bestModel)
#SSP585_2080:
plot(futureProjList_best[[1]])
#SSP585_2100
plot(futureProjList_best[[2]])
#SSP126_2080
plot(futureProjList_best[[3]])
#SSP126_2100
plot(futureProjList_best[[4]])

plot(Pred_secondBest)
#SSP585_2080:
plot(futureProjList_secondBest[[1]])
#SSP585_2100
plot(futureProjList_secondBest[[2]])
#SSP126_2080
plot(futureProjList_secondBest[[3]])
#SSP126_2100
plot(futureProjList_secondBest[[4]])

#save.image("~/Documents/porcupines/future_porc.RData")


#PAST:: porc#########
rm(list=ls())
load("~/Documents/porcupines/ENM_noMasking_porc3_iterativeKfold.RData")
Rasts_best = Rasts2PF
Rasts_secondBest = Rasts_threeF
rm(list= ls()[!(ls() %in% c('bestModel','Rasts_secondBest', 'Rasts_best' ,'yrange', 'xrange', 'Pred_secondBest', 'secondBest', 'Pred_bestModel'))])

vars<- rast("~/Documents/porcupines/rast/AllEnvVars_2.5.tif")
alt = vars$alt
bioclim = subset(vars, names(vars)[grepl("bio", names(vars))])
land = rast("~/Documents/porcupines/rast/rastAnthromes_adjusted_withWater.tif")
names(land) = "land_use"
puma = rast("~/Documents/porcupines/PastProjections/puma_Past_R1spThin.asc")
is.nan(unique(values(puma)[is.na(values(puma))]))
puma = subst(puma, NaN, NA)
is.nan(unique(values(puma)[is.na(values(puma))]))
names(puma) = "puma"
crs(puma) = crs(vars)
fisher = rast("~/Documents/porcupines/PastProjections/fisher_Past_R1spThin.asc")
is.nan(unique(values(fisher)[is.na(values(fisher))]))
fisher = subst(fisher, NaN, NA)
is.nan(unique(values(fisher)[is.na(values(fisher))]))
names(fisher) = "fisher"
crs(fisher) = crs(vars)

pastStack = c(crop(bioclim, c(xrange, yrange)), 
              crop(land, c(xrange, yrange)),
              crop(alt, c(xrange,yrange)),
              puma, fisher)

PastRastAll = raster::stack(pastStack)
names(PastRastAll)

is.factor(PastRastAll$land_use)
unique(PastRastAll$land_use)
r = ratify(PastRastAll$land_use)
rat <- levels(r)[[1]]
rat
rat$landcover <- c("Water",
                   "Forest",
                   "Grassland",
                   "Barren",
                   "Cropland",
                   "Urban",
                   "PermSnowIce")
rat$code <- c(1:7)
levels(r) <- rat
plot(r)
x <- deratify(r, "landcover")
plot(x)
is.factor(x)
PastRastAll$land_use = x
is.factor(PastRastAll$land_use)

PastRastBest = raster::subset(PastRastAll, names(Rasts_best))
all(sort(names(PastRastBest)) == sort(names(Rasts_best)))
Past_best<- predict(bestModel, PastRastBest)

PastRastSecondBest = raster::subset(PastRastAll, names(Rasts_secondBest))
all(sort(names(PastRastSecondBest)) == sort(names(Rasts_secondBest)))
Past_secondbest<- predict(secondBest, PastRastSecondBest)
plot(Pred_bestModel)
plot(Past_best)
plot(Pred_secondBest)
plot(Past_secondbest)

writeRaster(Pred_bestModel, "~/Documents/porcupines/PastProjections/porc_Past_bestModel.asc", format = 'ascii', overwrite=T)
writeRaster(Pred_secondBest, "~/Documents/porcupines/PastProjections/porc_Past_secondBest.asc", format = 'ascii', overwrite=T)

#save.image("~/Documents/porcupines/past_porc.RData")

#
#TRASH########



#
#saveRDS(fullModel2A_list, "~/Documents/porcupines/allModels/porc_fullModel2A_list.RData")
#saveRDS(fullModel2AP_list, "~/Documents/porcupines/allModels/porc_fullModel2AP_list.RData")
#saveRDS(fullModel2AF_list, "~/Documents/porcupines/allModels/porc_fullModel2AF_list.RData")
#saveRDS(fullModel2APF_list, "~/Documents/porcupines/allModels/porc_fullModel2APF_list.RData")#

#saveRDS(fullModel2A_noLU_list, "~/Documents/porcupines/allModels/porc_fullModel2A_noLU_list.RData")
#saveRDS(fullModel2AP_noLU_list, "~/Documents/porcupines/allModels/porc_fullModel2AP_noLU_list.RData")
#saveRDS(fullModel2AF_noLU_list, "~/Documents/porcupines/allModels/porc_fullModel2A_noLUF_list.RData")
#saveRDS(fullModel2APF_noLU_list, "~/Documents/porcupines/allModels/porc_fullModel2APF_noLU_list.RData")#

#saveRDS(fullMode3A_list, "~/Documents/porcupines/allModels/porc_fullMode3A_list.RData")
#saveRDS(fullModel3AP_list, "~/Documents/porcupines/allModels/porc_fullModel3AP_list.RData")
#saveRDS(fullModel3AF_list, "~/Documents/porcupines/allModels/porc_fullModel3AF_list.RData")
#saveRDS(fullModel3APF_list, "~/Documents/porcupines/allModels/porc_fullModel3APF_list.RData")#

#saveRDS(eval2A_list, "~/Documents/porcupines/allModels/porc_Eval2A_list.RData")
#saveRDS(eval2AP_list, "~/Documents/porcupines/allModels/porc_Eval2AP_list.RData")
#saveRDS(eval2AF_list, "~/Documents/porcupines/allModels/porc_Eval2AF_list.RData")
#saveRDS(eval2APF_list, "~/Documents/porcupines/allModels/porc_Eval2APF_list.RData")#

#saveRDS(eval2A_noLU_list, "~/Documents/porcupines/allModels/porc_Eval2A_noLU_list.RData")
#saveRDS(eval2AP_noLU_list, "~/Documents/porcupines/allModels/porc_Eval2AP_noLU_list.RData")
#saveRDS(eval2AF_noLU_list, "~/Documents/porcupines/allModels/porc_Eval2AF_noLU_list.RData")
#saveRDS(eval2APF_noLU_list, "~/Documents/porcupines/allModels/porc_Eval2APF_noLU_list.RData")#

#saveRDS(eval3A_list, "~/Documents/porcupines/allModels/porc_Eval3A_list.RData")
#saveRDS(eval3AP_list, "~/Documents/porcupines/allModels/porc_Eval3AP_list.RData")
#saveRDS(eval3AF_list, "~/Documents/porcupines/allModels/porc_Eval3AF_list.RData")
#saveRDS(eval3APF_list, "~/Documents/porcupines/allModels/porc_Eval3APF_list.RData")

#future porc######
rm(list=ls())
load("~/Documents/porcupines/ENM_noMasking_porc3_iterativeKfold.RData")
rm(list= ls()[!(ls() %in% c('bestModel','Rasts2PF', 'Rasts_threeF' ,'yrange', 'xrange', 'Pred_secondBest', 'secondBest', 'Pred_bestModel'))])

#bioclim:
f = list.files("~/Documents/porcupines/rast/microc6", full.names = T)
bioClimListAll = lapply(f, function(x){
  t <- rast(x) 
})
names(bioClimListAll) = gsub(".tif", "", list.files("~/Documents/porcupines/rast/microc6", full.names = F))

bioClimList_R3 = list()
bioClimList_R2  = list()
for(i in 1:length(bioClimListAll)){
  names_tmp3 = names(Rasts_threeF)[grepl("bio", names(Rasts_threeF))]
  tmp3 = subset(bioClimListAll[[i]], names_tmp3)
  print(all(names(tmp3) == names_tmp3))
  bioClimList_R3 [[i]] = tmp3
  
  names_tmp2 = names(Rasts2PF)[grepl("bio", names(Rasts2PF))]
  tmp2 = subset(bioClimListAll[[i]], names_tmp2)
  print(all(names(tmp2) == names_tmp2))
  bioClimList_R2 [[i]] = tmp2
}
names(bioClimList_R3) = names(bioClimListAll)
names(bioClimList_R2) = names(bioClimListAll)
all(names(bioClimList_R3[[1]]) == c('bio3', 'bio9', 'bio18', 'bio15', 'bio1', 'bio13'))


#preds:
f = list.files("~/Documents/porcupines/FutureProjectionsNew/puma", pattern = "R1",full.names = T)
puma = lapply(f, function(x){
  t <- rast(x) 
})
names(puma) = gsub(".asc", "", stringr::word(list.files("~/Documents/porcupines/FutureProjectionsNew/puma", pattern = "R1",full.names = F), -2,-1, "_"))
for(i in 1:length(puma)){
  puma[[i]] = subst(puma[[i]], NaN, NA)
  print(is.nan(unique(values(puma[[i]])[is.na(values(puma[[i]]))])))
  names(puma[[i]]) = "puma"
  crs(puma[[i]]) = crs(bioClimListAll[[1]]$bio1)
}

f = list.files("~/Documents/porcupines/FutureProjectionsNew/fisher", pattern = "R1",full.names = T)
fisher = lapply(f, function(x){
  t <- rast(x) 
})
names(fisher) = gsub(".asc", "", stringr::word(list.files("~/Documents/porcupines/FutureProjectionsNew/fisher", pattern = "R1",full.names = F), -2,-1, "_"))
for(i in 1:length(fisher)){
  fisher[[i]] = subst(fisher[[i]], NaN, NA)
  print(is.nan(unique(values(fisher[[i]])[is.na(values(fisher[[i]]))])))
  names(fisher[[i]]) = "fisher"
  crs(fisher[[i]]) = crs(bioClimListAll[[1]]$bio1)
}


#land use:
a =rast("~/Documents/porcupines/rast/global_SSP5_RCP85_2080_correctResolution.tif")
b =rast("~/Documents/porcupines/rast/global_SSP5_RCP85_2100_correctResolution.tif")
c =rast("~/Documents/porcupines/rast/global_SSP1_RCP26_2080_correctResolution.tif")
d =rast("~/Documents/porcupines/rast/global_SSP1_RCP26_2100_correctResolution.tif")
res(a)
res(d)
LUList = list(a,b,c,d)
names(LUList) = c("SSP585_2080", "SSP585_2100", "SSP126_2080", "SSP126_2100")
for(i in 1:length(LUList)){
  names(LUList[[i]]) = "land_use"
}
rm(a,b,c,d,f)

alt = rast("~/Documents/porcupines/rast/RCP85_SSP5_2080_global.tif")$alt

scenarios = c("SSP585_2080", "SSP585_2100", "SSP126_2080", "SSP126_2100")
####BEST: Rasts2PF#
futureList_best = list()
futureProjList_best = list()
for(i in 1:length(scenarios)){
  bioclim =bioClimList_R2[[which( names(bioClimList_R2)==scenarios[i])]]
  landUse = LUList[[which( names(LUList)==scenarios[i])]]
  p1 = raster::stack(puma[[which( names(puma)==scenarios[i])]])
  f1 = raster::stack(fisher[[which( names(fisher)==scenarios[i])]])
  Future1 = c(crop(bioclim, c(xrange, yrange)), crop(alt, c(xrange, yrange)))
  Future1 = c(Future1, crop(landUse, c(xrange, yrange)))
  Future1 = stack(Future1, p1)
  Future1 = stack(Future1, f1)
  
  r = ratify(Future1$land_use)
  rat <- levels(r)[[1]]
  rat
  rat$landcover <- c("Water",
                     "Forest",
                     "Grassland",
                     "Barren",
                     "Cropland",
                     "Urban",
                     "PermSnowIce")
  rat$code <- c(1:7)
  levels(r) <- rat
  x <- deratify(r, "landcover")
  if(!is.factor(x)){print("x is not a factor!")}
  Future1$land_use = x
  if(!is.factor(Future1$land_use)){print("land use is not a factor!")}
  if(!all(sort(names(Rasts2PF)) == sort(names(Future1)))){"variables don't match!"}
  futureList_best[[i]] = Future1
  futureProj<- predict(bestModel, Future1)
  futureProjList_best[[i]] = futureProj
  writeRaster(futureProj, paste0("~/Documents/porcupines/FutureProjectionsNew/porc_Future_R2PFspThin_", 
                                 scenarios[i],
                                 ".asc"), format = 'ascii', overwrite = T)
  print(i)
}
names(futureList_best) = scenarios
names(futureProjList_best) = scenarios


#current:
plot(Pred_bestModel)
#SSP585_2080:
plot(futureProjList_best[[1]])
#SSP585_2100
plot(futureProjList_best[[2]])
#SSP126_2080
plot(futureProjList_best[[3]])
#SSP126_2100
plot(futureProjList_best[[4]])

#SECOND BEST##
futureList_secondBest = list()
futureProjList_secondBest = list()
for(i in 1:length(scenarios)){
  bioclim =bioClimList[[which( names(bioClimList)==scenarios[i])]]
  landUse = LUList[[which( names(LUList)==scenarios[i])]]
  p1 = raster::stack(puma[[which( names(puma)==scenarios[i])]])
  Future1 = c(crop(bioclim, c(xrange, yrange)), crop(alt, c(xrange, yrange)))
  Future1 = c(Future1, crop(landUse, c(xrange, yrange)))
  Future1 = stack(Future1, f1)
  
  r = ratify(Future1$land_use)
  rat <- levels(r)[[1]]
  rat
  rat$landcover <- c("Water",
                     "Forest",
                     "Grassland",
                     "Barren",
                     "Cropland",
                     "Urban",
                     "PermSnowIce")
  rat$code <- c(1:7)
  levels(r) <- rat
  x <- deratify(r, "landcover")
  if(!is.factor(x)){print("x is not a factor!")}
  Future1$land_use = x
  if(!is.factor(Future1$land_use)){print("land use is not a factor!")}
  if(!all(names(Rasts_threeP) == names(Future1))){"variables don't match!"}
  futureList_secondBest[[i]] = Future1
  futureProj<- predict(bestModel, Future1)
  futureProjList_secondBest[[i]] = futureProj
  writeRaster(futureProj, paste0("~/Documents/porcupines/FutureProjectionsNew/porc_Future_R3PspThin_", 
                                 scenarios[i],
                                 ".asc"), format = 'ascii', overwrite = T)
  print(i)
}
names(futureList_secondBest) = scenarios
names(futureProjList_secondBest) = scenarios

plot(Pred_secondBest)
#SSP585_2080:
plot(futureProjList_secondBest[[1]])
#SSP585_2100
plot(futureProjList_secondBest[[2]])
#SSP126_2080
plot(futureProjList_secondBest[[3]])
#SSP126_2100
plot(futureProjList_secondBest[[4]])

#save.image("~/Documents/porcupines/future_porc.RData")


bioClimList_R3 = list()
bioClimList_R2  = list()
for(i in 1:length(bioClimListAll)){
  names_tmp3 = names(Rasts_secondBest)[grepl("bio", names(Rasts_secondBest))]
  tmp3 = subset(bioClimListAll[[i]], names_tmp3)
  print(all(names(tmp3) == names_tmp3))
  bioClimList_R3 [[i]] = tmp3
  
  names_tmp2 = names(Rasts_best)[grepl("bio", names(Rasts_best))]
  tmp2 = subset(bioClimListAll[[i]], names_tmp2)
  print(all(names(tmp2) == names_tmp2))
  bioClimList_R2 [[i]] = tmp2
}
names(bioClimList_R3) = names(bioClimListAll)
names(bioClimList_R2) = names(bioClimListAll)
all(names(bioClimList_R3[[1]]) == c('bio3', 'bio9', 'bio18', 'bio15', 'bio1', 'bio13'))

