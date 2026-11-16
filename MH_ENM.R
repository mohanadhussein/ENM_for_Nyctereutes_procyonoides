setwd("D:/İTÜ/7th semester/BIO 468E current topics in Eco")

install.packages("tidyterra")
library(tidyterra)
library(biomod2)
library(raster)
library(dismo)
library(sp)
library(dplyr)


#downloading data using API, don't include fossil data
Nyctereutes <- gbif(URLencode("Nyctereutes Temminck"), URLencode("Nyctereutes procyonoides"), ext = NULL, args = "basisOfRecord=HUMAN_OBSERVATION", geo = TRUE, sp = TRUE, removeZeros = TRUE, download = TRUE, ntries = 5, nrecs = 300, start = 1, end = 1000 )

#reparing data
Nyctereutes.df <- data.frame(Nyctereutes)
NyctereutesXY <- Nyctereutes.df[,c('lon','lat')]
Nyctereutesdfs <- slice_sample(NyctereutesXY, n = 300, replace = FALSE) #select 300 random rows only once

### Preparing data for biomod2
myRespVar <- 'Nyctereutes.procyonoides'
myResp <- SpatialPoints(coords = Nyctereutesdfs, proj4string = CRS(as.character(NA)), bbox = NULL)
myRespXY <- data.frame(Nyctereutesdfs)

#Why did you choose these parameters, write clearly 2 4 5 6 14 15 18 19
setwd("D:/İTÜ/7th semester/BIO 468E current topics in Eco/wc2.1_5m_bio")
bio2 <- raster("wc2.1_5m_bio_2.tif")
bio4 <- raster("wc2.1_5m_bio_4.tif") 
bio5 <- raster("wc2.1_5m_bio_5.tif") 
bio6 <- raster("wc2.1_5m_bio_6.tif")
bio14 <- raster("wc2.1_5m_bio_14.tif")
bio15 <- raster("wc2.1_5m_bio_15.tif") 
bio18 <- raster("wc2.1_5m_bio_18.tif")
bio19 <- raster("wc2.1_5m_bio_19.tif")

#explanatory variable global view
myExpl <- stack(bio2, bio4, bio5, bio6, bio14, bio15, bio18, bio19)
plot(myExpl)

#explanatory variable regional view
turkey_shape <- getData("GADM", country = "Turkey", level = 0)
myExpl_turkey <- crop(myExpl, turkey_shape) # Crop the explanatory variables to the extent of Turkey
myExpl_turkey_stack <- stack(myExpl_turkey) # convert into RasterStack required by BIOMOD
plot(myExpl_turkey_stack)

# format for BIOMOD
myBiomodData <- BIOMOD_FormatingData(resp.var=myResp,
                                     expl.var = myExpl,
                                     resp.xy= myRespXY,
                                     resp.name=myRespVar,
                                     PA.nb.rep = 1,
                                     PA.nb.absences = 1000,
                                     PA.strategy = 'SRE',
                                     filter.raster = TRUE,
                                     na.rm = TRUE)

#check cross validation using this table
crossval= bm_CrossValidation(
  myBiomodData,
  strategy = "random",
  nb.rep = 1,
  perc = 0.8,
  k = 0,
  balance = "presences",
  env.var = NULL,
  strat = "both",
  user.table = NULL,
  do.full.models = FALSE
)


myBiomodModelout <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = paste(myRespVar,"Modeling",sep=""),
                                    models = 'MAXNET',
                                    CV.strategy = "random",
                                    CV.nb.rep = 3,
                                    CV.perc = 0.75,
                                    prevalence = 0.5,
                                    var.import = 5,
                                    metric.eval = c('ROC', 'TSS'),
                                    scale.models = TRUE,
                                    CV.do.full.models = TRUE)

myBiomodEval <- BIOMOD_Evaluate(bm.mod = myBiomodModelout,
                                data = myBiomodData, # Use the model output and formatted data
                                eval.metric = c('ROC', 'TSS'), # Specify evaluation metrics
                                do.full.models = TRUE)

evaluation_df <- myBiomodModelout@models.failed

  
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelout,
                                  new.env = myExpl_turkey_stack,
                                  proj.name = 'current',
                                  models.chosen = 'all',
                                  metric.binary = c('ROC', 'TSS'),
                                  compress = TRUE,
                                  clamping.mask=F)

plot(myBiomodProj)

#future projection  2 4 5 6 14 15 18 19
abio <- raster::stack("wc2.1_5m_bioc_CMCC-ESM2_ssp585_2041-2060.tif")
a_bio2 <- subset(abio,2)
a_bio4 <- subset(abio,4)
a_bio5 <- subset(abio,5)
a_bio6 <- subset(abio,6)
a_bio14 <- subset(abio,14)
a_bio15 <- subset(abio,15)
a_bio18 <- subset(abio,18)
a_bio19 <- subset(abio,19)

myExplFuture <- stack(a_bio2, a_bio4, a_bio5, a_bio6, a_bio14,a_bio15,a_bio18,a_bio19)
names(myExplFuture) = names(myExpl)
plot(myExplFuture)

#future projection for turkey
turkey_shape <- getData("GADM", country = "Turkey", level = 0)
myExplFut_turkey <- crop(myExplFuture, turkey_shape) # Crop the explanatory variables to the extent of Turkey
myExplFut_turkey_stack <- stack(myExplFut_turkey) # convert into RasterStack required by BIOMOD
plot(myExplFut_turkey_stack)


myBiomodProjFut <- BIOMOD_Projection(bm.mod = myBiomodModelout,
                                     new.env = myExplFut_turkey_stack,
                                     proj.name = 'future',
                                     models.chosen = 'all',
                                     metric.binary = c('ROC', 'TSS'),
                                     compress = TRUE,
                                     clamping.mask= F)
plot(myBiomodProjFut)
 


myBiomodModelEval <- get_evaluations(myBiomodModelout)
myBiomodModelEval
