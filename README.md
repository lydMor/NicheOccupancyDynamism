# NicheOccupancyDynamism
Scripts and code associated with the generation of explicit tests for niche occupancy dynamism

The code provided here can be used to conduct temporal niche comparisons within the same species. In order to conduct appropriate comparisons with your own data, you should ensure that you have a reasonable inflection point (i.e. a time period that you suspect changes in niche occupancy may have occurred), which we differentiate using T1 (time 1) and T2 (time 2). Data requried for these analyses include: 

1) occurrence records for T1 and T2

2) Environmental data for T1 and T2

Here, you can use the provided files and associated R scripts to reproduce our original temporal niche comparisons. 

**Provided files:** 

TemporalNicheComp.R - R script detailing temporal niche comparison process. 

Dist_T1.csv - occurrence records for a focal species at T1 

Dist_T2.csv - occurrence records for a focal species at T2 

Env_T1.csv - environmental data associated with focal species occurrence records at T1

Env_T2.csv - environmental data associated with focal species occurence records at T2

**Additional data and code:**

On this page, you can also find data and code associated with the publication "Using temporal niche comparisons to capture niche occupancy dynamism: a case study in the North American Porcupine". We've provided the full porcupine distribution dataset, our environmental data, and the scripts used to conduct traditional ENMs and temporal forecasting. 

ENM.csv - script associated with generating an Ecological Niche Model using MAXENT in R. 

Porcupine.csv - full porcupine distribution data

Puma.csv - full Puma distribution data

Fisher.csv - full fisher distribution data

Env.tif - environmental rasters used to build porcupine, puma, and fisher ENMs


