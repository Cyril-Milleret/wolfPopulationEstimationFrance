# R code to reproduce the analysis to estimate wolf population size 
# Milleret, C., Duchamp, D., Vandel, JM., Gimenez, O. 2025. Mise à jour des estimations démographiques et des effectifs de la population de loups en France lors de l’hiver 2024/25. OFB/CEFE-CNRS. XX pages. Disponible sur : <https:/ www.loupfrance.fr/>
rm(list=ls())
#LOAD PACKAGES
library(nimble)
library(nimbleSCR)
library(raster)
library(sf)
library(stars)
library(terra)
library(R.utils)
library(basicMCMCplots)
library(coda)
library(tidyverse)
library(RANN)

library(readxl)

## WORKING DIRECTORY & MODEL NAME
# set your working directory here
WD <- "C:/Personal_Cloud/OneDrive/Work/CNRS/Reports/GitHub/wolfPopulationEstimationFrance"
myVars <- list(
  # HABITAT SPECIFICATIONS
  HABITAT = list( habResolution = 10000,
                  habBuffer = 20000),
  WD = "C:/Personal_Cloud/OneDrive/Work/CNRS/SCR",
  
  # NGS DATA SPECIFICATIONS
  DATA = list( years = 2025,
               sex = c("F","M"),                   ## "Hann","Hunn","Ukjent"
               samplingMonths = list(11,12,1:3)), ## list(10:12,1:4), list(1:XXX), list(XX:XX,YY:YY)
  
  # DETECTORS SPECIFICATIONS
  DETECTORS = list( detSubResolution = 1000,
                    detResolution = 5000,
                    detDeadResolution = 15000),
  
  # DATA GENERATION
  DETECTIONS = list( maxDetDist = 40000,
                     resizeFactor = 0,
                     aug.factor = 4),
  
  ## MISCELLANEOUS
  plot.check = TRUE)

years <- myVars$DATA$years
nYears <- length(years)
YEARS <- lapply(years, function(x)c(x,x+1))

##source functions
sourceDirectory(file.path(WD,"Functions"), modifiedOnly = FALSE)

#SOURCE CPP FUNCTIONS FOR FAST MAPPING
# dir.function.cpp <- "C:/My_documents/rovquant/analyses/Rgit/RovQuant/Temp/CM/functions/cpp" ## The Temp folder
# for(i in list.files(dir.function.cpp)){sourceCpp(filePath(dir.function.cpp,i))}


## ==== 1. LOAD DATA ####
## ====  1.1 DNA DATA ####
## ====   1.1.1 ALIVE DATA ####
load(file.path(WD,"2025/Data/DNA.RData"))
write.csv(DNA,file = file.path(WD,"2025/Data/DNA.csv"))

## ====  1.2 GIS DATA ####
## ====   1.2.1 10*10 km SAMPLING GRID #####
grid1010 <- read_sf(file.path(WD,"2025/Data/","plan_echantillonnage_reg_2024-2025.shp"))
grid1010 <- grid1010[!is.na(grid1010$Meute), ]

st_crs(DNA) <- st_crs(grid1010)
plot(grid1010$geometry)

## ====   1.2.2 EUROPE ####
Europe <- read_sf(file.path(WD,"GISLayers","EuropeCountries"))
##subset to neighboring countries 
Europe <- Europe[Europe$NAME %in%c("France","Spain","Italy","Switzerland","Germany","Luxembourg","Belgium","Andorra","Monaco"),]
Europe <- st_transform(Europe,crs = st_crs(DNA))
## ====   1.2.3 FRANCE ####
France <- Europe[Europe$NAME %in%c("France"),]
France <- st_transform(France,crs = st_crs(DNA))
##REMOVE ALL ISLANDS 
parts = st_cast(France, "POLYGON", warn = FALSE)
# Step 2: compute areas
parts = parts |> 
  mutate(area = st_area(geometry))
# remove small islands
# keep only france
France = parts |> 
  filter(area == max(parts$area))

plot(France$geometry)
plot(DNA$geometry,add=T,col="red")

## ====   1.2.4 DEPARTEMENT #####
Departement <- read_sf(file.path(WD,"GISLayers","departements-20180101-shp","departements-20180101.shp"))
Departement <- st_transform(Departement,crs = st_crs(DNA))

plot(Departement$geometry,add=T)

## ====   1.2.5 REGION #####
Region <- read_sf(file.path(WD,"GISLayers","Region","regions_2015_metropole_region.shp"))
Region <- st_set_crs(Region, 27572)#required to fix the issue with the projection
Region <- st_transform(Region,crs = st_crs(DNA))

plot(Region$geometry,add=T,col="red")


## ====   1.2.6 COMMUNES #####
Communes <- read_sf(file.path(myVars$WD,"GIS","AdminBoundaries","Communes20220101","communes-20220101FR.shp"))
CommunesOriginal <- Communes <- st_transform(Communes, crs = st_crs(DNA))
Communes <- st_crop( Communes, st_buffer(France,dist = 20000)$geometry)
Communes <- st_simplify(Communes, preserveTopology = T, dTolerance = 200)
#CHECK IF THERE IS AN ISSUE WITH THE GEOMETRY
any(sf::st_is_empty(Communes))

## ====   1.2.7 DEM  #####
#dem was downloaded from elevatr::get_elev_raster
DEM <- raster(file.path(myVars$WD,"GIS","DEM","DEM.tif"))
# dem <- crop(DEM,st_buffer(France,dist = 80000))
#use focal to get an approx 10km resolution 
DEM <- raster::focal(DEM,w=matrix(1,nrow=11,ncol=11),fun=mean, na.rm=T, pad=T)
# DEM1 <- st_as_stars(DEM)
# carmel_mean15 = focal2(carmel, matrix(1, 15, 15), "mean")
plot(DEM)
plot(France$geometry,add=T)

## ====   1.2.8 LCIE WOLF DISTRIBUTION #####
#2023 DISTRIBUTION MAP 
LCIEDistribution2023 <- read_sf(file.path(WD,"GISLayers","LCIE_wolfDistribution","Wolf_LCIE_2023","Wolf_2017-2022_LCIE.shp"))
LCIEDistribution2023$PRESENCE1 <- 0
LCIEDistribution2023$PRESENCE1[LCIEDistribution2023$PRESENCE %in% "Sporadic"] <- 1
LCIEDistribution2023$PRESENCE1[LCIEDistribution2023$PRESENCE %in% "Permanent"] <- 3


#2018 DISTRIBUTION MAP 
LCIEDistribution2018 <- read_sf(file.path(WD,"GISLayers","LCIE_wolfDistribution","Wolf_LCIE_2018","2018_06_06_Wolf_IUCN_RedList.shp"))
str(LCIEDistribution2018)
LCIEDistribution2018$PRESENCE1 <- ifelse(LCIEDistribution2018$SPOIS %in% "Sporadic",1,3)
#2012 DISTRIBUTION MAP 
LCIEDistribution2012permanent <- read_sf(file.path(WD,"GISLayers","LCIE_wolfDistribution","Wolf_LCIE_2012","Clip_2012_12_01_Wolves_permanent.shp"))
LCIEDistribution2012permanent$SPOIS <- 3
LCIEDistribution2012sporadic <- read_sf(file.path(WD,"GISLayers","LCIE_wolfDistribution","Wolf_LCIE_2012","Clip_2012_12_01_Wolves_sporadic.shp"))
LCIEDistribution2012sporadic$SPOIS <- 1
LCIEDistribution2012 <- rbind(LCIEDistribution2012permanent, LCIEDistribution2012sporadic)
LCIEDistribution2012$PRESENCE1 <- LCIEDistribution2012$SPOIS

LCIEDistribution2012 <- st_transform(LCIEDistribution2012, crs=st_crs(DNA))
LCIEDistribution2018 <- st_transform(LCIEDistribution2018, crs=st_crs(DNA))
LCIEDistribution2023 <- st_transform(LCIEDistribution2023, crs=st_crs(DNA))


## ====   1.2.7 EFFORT #####
## ====     1.2.7.1 POTENTIAL Bauduin et al 2023 #####
#LOAD GRID
load(file.path(WD,"2025","Data","Effort.RData"))

## ====   1.2.10 SNOW #####
load(file.path(WD,"2025/Data/SNOW.RData"))
plot(SNOW)




## ====   1.2.12 ROAD DENSITY #####
load(file.path(WD,"GISLayers/roads.RData"))


## ====   1.2.13 FOREST #####
#Corinne land cover data should be downloaded https://land.copernicus.eu/en/products/corine-land-cover
#Layers are not provided given their size 
#CLC CLIPPED TO FRANCE IN QGIS
# CLC_FR <- st_read(file.path("C:/Personal_Cloud/OneDrive/Work/CNRS/SCR","GIS","CLC","u2018_clc2018_v2020_20u1_fgdb","CLC_FR.shp"))
#CLCForest <- CLC_FR[CLC_FR$Code_18 %in% c("311", "312", "313"),]#,"324
# ## ====   1.2.12 GRASSLAND #####
# CLCGrassland <- CLC_FR[CLC_FR$Code_18 %in% c("321", "322" ),]
## ====   1.2.15 Agriculture #####
#CLCAgri <- CLC_FR[as.numeric(CLC_FR$Code_18)  >199 &
#                    as.numeric(CLC_FR$Code_18)  <300,]
## ====   1.2.16 Human #####
# CLCHum <- CLC_FR[as.numeric(CLC_FR$Code_18)  >99 &
#                    as.numeric(CLC_FR$Code_18)  <200,]

## ====   1.2.14 AIRE DE PRESENCE ====
AirePresence <- read_sf(file.path(WD,"GISLayers","AirePresenceFR",
                                  "Massif_France_Mailles_Reg_Sce_2_2.shp"))

## ==== 2. Check NGS DATA ####
barplot(table(DNA$Year))
table(DNA$mois,DNA$Year)
NDet <- tapply(DNA$Id, DNA$Id, length)
mean(NDet)

## ====     2.1 PLOT AND SUMMARY TO CHECK EVERYTHING IS BEING USED ####
nrow(DNA)
plot(France$geometry)
plot(DNA$geometry,pch=16, cex=0.5, col="red",add=T)

#MAKE A TABLE WITH SUMMARY OF DETECTIONS PER MONTHS AND SEXE
matDets <- matrix(NA,nrow=4,ncol=6)
colnames(matDets) <- c("Nov","Dec","Jan","Feb","Mar","Total")
rownames(matDets) <- c("F","M","NI","Total")
tab <- table(DNA$SEXE,DNA$mois)
matDets[1:3,1:5]<- tab[c(2,3,1),c(4,5,1,2,3)]
matDets[,"Total"] <- rowSums(matDets,na.rm=T)
matDets["Total",] <- colSums(matDets,na.rm=T)

apply(table(DNA$Id,DNA$SEXE),2,function(x) sum(x>0))



## ==== 3. GENERATE HABITAT ====
## ====   3.1 DEFINE AREA SEARCHED ====
areaSearched <- st_buffer(DNA,dist = 600000)
areaSearched <- st_union(areaSearched)
#REMOVE ALL AREAS OUTSIDE OF FRANCE (NOT SEARCHED)
areaSearched <- st_intersection(areaSearched, France)
areaSearched <- st_union(areaSearched)
areaSearched <- st_as_sf(areaSearched)
plot(areaSearched)


## ====   3.2 DEFINE HABTIAT EXTENT ====
#BUFFER AREA SEARCHED
habitatExtent <- st_buffer(areaSearched,dist = myVars$HABITAT$habBuffer)
habitatExtent <- st_union(habitatExtent)
#REMOVE THE SEA FROM THE BUFFER
habitatExtent <- st_intersection(habitatExtent, Europe)
habitatExtent <- st_union(habitatExtent)
habitatExtent <- st_as_sf(habitatExtent)
plot(habitatExtent)
plot(France$geometry,add=T)

## ====   3.3 RASTERIZE HABITAT ====
#BASE THE HABITAT ON THE 10*10 KM GRID USED FOR THE MONITORING 
idx <- st_intersects(AirePresence, habitatExtent)
gridhabitat <- AirePresence[lengths(idx) > 0, ]
gridhabitat <- st_rasterize(gridhabitat, dx=10000, dy=10000 )

# #SOME SEARCHED GRID CELLS ARE OUTSIDE THE HABITAT (VERY LOW OVERLAP WITH HABITAT; DO NOT INCLUDE THEM)
# #keep habitat prop>50%
template.rSpa <- rast(gridhabitat[1])
template.rSpa <- raster(template.rSpa)

template.rSpa <- raster::rasterize(habitatExtent, template.rSpa, getCover=TRUE)        ## Add field=1 so its 1 for study area and NA for the rest.
template.rSpa[template.rSpa < 0.50] <- 0
template.rSpa[template.rSpa >= 0.50] <- 1
template.r <- template.rSpa

#store habitat 
habitat.r <- template.r#
habitat.r1 <- aggregate(habitat.r,fact=1)
habitat.r[habitat.r[]%in%0] <- NA
# GET THE HABITAT OBJECTS NECESSARY FOR NIMBLESCR
habitatxy <- data.frame(raster::coordinates(habitat.r)[!is.na(habitat.r)[],])#,df=T)#[is.na(habitat.r[][,1]),]
habitatxysf <- st_as_sf(habitatxy,coords = c("x", "y"), crs = st_crs(DNA))#double check the projection system
habitatSF.pol <- sf::st_as_sf(stars::st_as_stars(habitat.r), 
                              as_points = FALSE, merge = F)

#PLOT CHECK 
plot(habitat.r)
plot(habitatxysf$geometry,add=T)
plot(France$geometry,add=T)
plot(DNA$geometry,add=T,pch=16,col="red")
plot(grid1010$geometry,add=T,col=NA)

## ====   3.4 DEFINE DETECTOR GRID ====
detector.r <- rasterize(areaSearched, template.rSpa)
#BASE THE DETECTOR GRID ON SUBDETECTORS
subdetector.r <- disagg(rast(detector.r), fact=40)#250 resolution for the PAB
myDetectors <- MakeSearchGrid( subdetector.r=subdetector.r,
                               detResolution=myVars$DETECTORS$detResolution,
                               plot = F
)
myDetectors$detector.xy <- st_coordinates(myDetectors$main.detector.sf)
colnames(myDetectors$detector.xy) <- c("x","y")

#PLOT CHECK 
plot(habitat.r)
plot(France$geometry,add=T)
plot(myDetectors$main.detector.sf$geometry,add=T,col="red",pch=16,cex=0.1)


# ==== 4. GENERATE HABITAT-LEVEL COVARIATES ====
## ====   4.1 WOLF PRESENCE ====
## EXTRACT WOLF PERMANENT PRESENCE  IN EACH HABITAT CELL
habitatGrid  <- sf::st_as_sf(stars::st_as_stars(template.rSpa), 
                             as_points = FALSE, merge = F)
habitatGrid <- st_transform(habitatGrid,st_crs(LCIEDistribution2012))

#2012
sizehabitatCell <- st_area(habitatGrid)[1]
habitatGrid$id <- 1:nrow(habitatGrid)
intersection <- st_intersection(habitatGrid, LCIEDistribution2012) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  summarise(IUCN2012 = sum(iucn*PRESENCE1)/sizehabitatCell)
#summarise(IUCN = sum(SPOIS))

habitatGrid <- habitatGrid %>%
  left_join(intersection, by = "id")
habitatGrid$IUCN2012[is.na(habitatGrid$IUCN2012)] <- 0
plot(habitatGrid[,"IUCN2012"])

#2018
## Extract LCIE wolf presence in each habitat grid cell
intersection <- st_intersection(habitatGrid, LCIEDistribution2018) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  #summarise(iucn_2 = sum(SPOIS))
  summarise(IUCN2018 = sum(iucn*PRESENCE1)/sizehabitatCell)

tmp <- habitatGrid %>%
  left_join(intersection, by = "id")
tmp$IUCN2018[is.na(tmp$IUCN2018)] <- 0

#2023
## Extract LCIE wolf presence in each habitat grid cell
#combine aire de présence 2023 with LCIE
LCIECODE3 <- LCIEDistribution2023$CELLCODE[LCIEDistribution2023$PRESENCE1>2]
LCIECODE1 <- LCIEDistribution2023$CELLCODE[LCIEDistribution2023$PRESENCE1%in%1]
AirePresence$reg2023LCIE <- AirePresence$Reg_2023
AirePresence$reg2023LCIE[AirePresence$reg2023LCIE%in% 2] <- 3 
AirePresence$reg2023LCIE[AirePresence$CellCode%in% LCIECODE3 & AirePresence$reg2023LCIE%in% 0] <- 3
AirePresence$reg2023LCIE[AirePresence$CellCode%in% LCIECODE1 & AirePresence$reg2023LCIE%in% 0] <- 1



intersection <- st_intersection(habitatGrid, AirePresence) %>%
  mutate(iucn = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(id) %>%
  #summarise(iucn_2 = sum(SPOIS))
  summarise(IUCN2023 = sum(iucn*reg2023LCIE)/sizehabitatCell)

tmp1 <- habitatGrid %>%
  left_join(intersection, by = "id")
tmp1$IUCN2023[is.na(tmp1$IUCN2023)] <- 0

##sum all habitat grid
habitatGrid$IUCN <- habitatGrid$IUCN2012 + tmp$IUCN2018 + tmp1$IUCN2023
plot(habitatGrid[habitatGrid$layer>0,"IUCN"])

habitatGrid$IUCN <- scale(habitatGrid$IUCN)

## ====   4.4 PROPORTION OF FOREST ====
# habitatSF.pol1 <- st_transform(habitatSF.pol, st_crs(CLCForest))
# habitatSF.pol1$ID <- 1:nrow(habitatSF.pol1) 
# 
# cellForest <- st_intersection(habitatSF.pol1 , CLCForest) %>%
#   mutate(areaInter = st_area(.)) %>%
#   group_by(ID) %>%
#   summarise(areaforestCell = sum(areaInter)) %>%
#   mutate(propForestCell = areaforestCell / sizehabitatCell) %>%
#   as_tibble() %>%
#   select(ID, propForestCell)
# 
# #PLOT CHECK  
# forest.r <- habitat.r
# forest.r[] <- 0
# 
# forest.r[which(!is.na(habitat.r[]))[cellForest$ID]] <- cellForest$propForestCell
# plot(forest.r)
# 
# ## ====   4.5 PROPORTION OF GRASSLAND ====
# cellGrassland <- st_intersection(habitatSF.pol1 , CLCGrassland) %>%
#   mutate(areaInter = st_area(.)) %>%
#   group_by(ID) %>%
#   summarise(areaGrasslandCell = sum(areaInter)) %>%
#   mutate(propGrasslandCell = areaGrasslandCell / sizehabitatCell) %>%
#   as_tibble() %>%
#   select(ID, propGrasslandCell)
# #plot#check 
# grassland.r <- habitat.r
# grassland.r[] <- 0
# 
# grassland.r[which(!is.na(habitat.r[]))[cellGrassland$ID]] <- cellGrassland$propGrasslandCell
# plot(grassland.r)
# 
# ## ====   4.6 PROPORTION OF Agri ====
# cellAgri<- st_intersection(habitatSF.pol1 , CLCAgri) %>%
#   mutate(areaInter = st_area(.)) %>%
#   group_by(ID) %>%
#   summarise(areaAgriCell = sum(areaInter)) %>%
#   mutate(propAgriCell = areaAgriCell / sizehabitatCell) %>%
#   as_tibble() %>%
#   select(ID, propAgriCell)
# #plot#check 
# Agri.r <- habitat.r
# Agri.r[] <- 0
# 
# Agri.r[which(!is.na(habitat.r[]))[cellAgri$ID]] <- cellAgri$propAgriCell
# plot(Agri.r)
# 
# ## ====   4.7 PROPORTION OF human ====
# cellHum <- st_intersection(habitatSF.pol1 , CLCHum) %>%
#   mutate(areaInter = st_area(.)) %>%
#   group_by(ID) %>%
#   summarise(areaHumCell = sum(areaInter)) %>%
#   mutate(propHumCell = areaHumCell / sizehabitatCell) %>%
#   as_tibble() %>%
#   select(ID, propHumCell)
# #plot#check 
# Hum.r <- habitat.r
# Hum.r[] <- 0
# 
# Hum.r[which(!is.na(habitat.r[]))[cellHum$ID]] <- cellHum$propHumCell
# plot(Hum.r)
# save(grassland.r,forest.r,Hum.r,Agri.r, 
# file=file.path(WD,"GISLayers","LandCoverRaster.RData"))
load(file.path(WD,"GISLayers","LandCoverRaster.RData"))

## ====   4.4 BIND COVARIATES ====
habCovs <-  cbind(habitatGrid$IUCN[habitatGrid$layer>0],
                  Hum.r[!is.na(habitat.r[])],
                  grassland.r[!is.na(habitat.r[])],
                  forest.r[!is.na(habitat.r[])],
                  Agri.r[!is.na(habitat.r[])])
#SCALE COVARIATES 
habCovs <- scale(habCovs)
colnames(habCovs) <- c(#"DEM",
  "Wolf","Human","Grassland","Forest","Agriculture")#,"distCore")

##PlotCheckAndSave
pdf(file=file.path(WD,"2025","output","CovDens.pdf"))
tmp <- habitat.r
for(i in 1:dim(habCovs)[2]){
  tmp[!is.na(habitat.r[])] <- habCovs[,i]
  plot(tmp,main=colnames(habCovs)[i])
  plot(France$geometry,add=T)
}
dev.off()

## ==== 5. GENERATE DECTOR-LEVEL COVARIATES ====
myDetectors$main.detector.sf <- st_transform(myDetectors$main.detector.sf,st_crs(DNA))

## ====   5.1 EFFORT ====
 myDetectors$grid.poly <- st_transform(myDetectors$grid.poly,
                                       st_crs(CommunesBillGEACODNA))

intersection <- st_intersection(myDetectors$grid.poly,
                                CommunesBillGEACODNA) %>%
  mutate(area = st_area(.)) %>%
  st_drop_geometry() %>%
  group_by(layer) %>%
  #summarise(iucn_2 = sum(SPOIS))
  summarise(Nsorties = sum(area*nAllDNA)/sizehabitatCell)

tmpdet <- myDetectors$grid.poly %>%
  left_join(intersection, by = "layer")
# tmp$IUCN2018[is.na(tmp$IUCN2018)] <- 0
tmpdet[is.na(tmpdet)] <- 0
#check
plot(tmpdet["Nsorties"],border=NA)

## ====   5.2 REGION ====
Departement$id <- 1:nrow(Departement)
# Departement <- st_transform(Departement,
#                             st_crs(myDetectors$main.detector.sf))

tmpDis <- st_distance(myDetectors$main.detector.sf,
                      st_simplify(Departement,dTolerance = 1000,preserveTopology = T))
# tmpDis[500,]
DetRegion <- 0
for(i in 1:nrow(tmpDis)){
  DetRegion[i] <- which.min(tmpDis[i,])
}
DetRegion <- apply(tmpDis, 1, which.min)

table(DetRegion)
centr <- st_centroid(Departement)$geometry
plot(France$geometry)
text(st_coordinates(centr),Departement$code_insee)

Departement$code_insee[c(53,11,17,16,42,100,99,79,78,87,60)]
Departement$code_insee[c(88,89,95,97,4,52,18,37)]
Departement$code_insee[c(27,34,2,19,57,10,15,24,96,3,8,14,7,101,102,71,64,65,59)]

inseRegion1 <- c("01","39","71","25","21","70","90","68","69D","69M","58","89","42","75",
                 "45","77","02","59","08","51","10","52","88","54","57","67","62","80","76","95","92","93",
                 "94","91","78","28","27","76","60","55","50","14","61","95")

inseRegion2 <- c("07","26")

#inseRegion2 <- c("26","05","84","04","13","83","06")
inseRegion3 <- c("66","09","11","34","31","76","81","12","30","32","47",
                 "48","43","15","63","19","46","82","71","84","13","24",
                 "87","16","23","86","36","18","03","63","43","65","64","40","33","17",
                 "79","85","29","22","35","53","72","44","37","49","41","56")


DetRegion[DetRegion %in% which(Departement$code_insee %in%c(inseRegion1)) ] <- "A"
DetRegion[DetRegion %in% which(Departement$code_insee %in%c(inseRegion2)) ] <- "Z"
DetRegion[DetRegion %in% which(Departement$code_insee %in%c(inseRegion3)) ] <- "C"

# DetRegion[DetRegion %in% which(Departement$code_insee %in%c("84","13")) ] <- "84"

# DetRegion[DetRegion %in% which(Departement$code_insee %in%inseRegion2) ] <- "B"
# DetRegion[DetRegion %in% which(Departement$code_insee %in%inseRegion3) ] <- "C"
Departement[]
table(DetRegion)

trapIntercept <- as.numeric(as.factor(DetRegion))
par(mar=c(0,0,0,0))
plot(France$geometry)
col <- rainbow(max(as.numeric(as.factor(DetRegion)))+3)
plot(myDetectors$main.detector.sf$geometry,col=col[as.numeric(as.factor(DetRegion))+2],add=T)
plot(Departement$geometry,add=T)
text(st_coordinates(centr),Departement$code_insee)
# mapview::mapview(Departement)

## ====   5.3 SNOW ====
snowCov <- raster::extract(SNOW, st_transform(myDetectors$main.detector.sf,crs = st_crs(SNOW)))

## if NA returns the average value of the cells within 20000m 
snowXY <- coordinates(SNOW) 
snowXY <- st_as_sf(data.frame(snowXY), coords=c("x","y"), crs=st_crs(SNOW))

isna <- which(is.na(snowCov))#, 1, function(x)any(is.na(x))))
if(length(isna)>0){
  whichClose <- apply(st_distance(snowXY, st_transform(myDetectors$main.detector.sf[isna, ],crs=st_crs(SNOW))),2,function(x){order(x) [1:8]})
  for(i in 1:length(isna)){
    snowCov[isna[i]] <- mean(SNOW[whichClose[,i]],na.rm=T)
  }
}

##check 
tmpr <- myDetectors$maindetector.r
tmpr[!is.na(tmpr)] <- snowCov
plot(tmpr)



## ====   5.4 ROADS ====
roads <- focal(roads, matrix(1,3,3),mean,na.rm=T)
# plot(Road.r)
RoadCov <- terra::extract(roads, st_transform(myDetectors$main.detector.sf,crs = st_crs(roads)))
isna <- which(is.na(RoadCov))#, 1, function(x)any(is.na(x))))


roadXY <- coordinates(roads) 
roadXY <- st_as_sf(data.frame(roadXY), coords=c("x","y"), crs=st_crs(roads))

isna <- which(is.na(RoadCov))#, 1, function(x)any(is.na(x))))
if(length(isna)>0){
  whichClose <- apply(st_distance(roadXY, st_transform(myDetectors$main.detector.sf[isna, ],crs=st_crs(roads))),2,function(x){order(x) [1:8]})
  for(i in 1:length(isna)){
    RoadCov[isna[i]] <- mean(roads[whichClose[,i]],na.rm=T)
  }
}
tmpr[!is.na(tmpr[])]<- RoadCov
plot(tmpr)
## ====   5.2 BIND COVARIATES ====
trapCov <- scale(cbind(as.numeric(tmpdet$Nsorties), snowCov, log(RoadCov)))#myDetectors$grid.poly$nbVisits
colnames(trapCov) <- c("Effort","Snow","Roads")

cor(trapCov)

##PlotCheck
pdf(file=file.path(WD,"2025","output","CovTraps.pdf"))
tmp <- myDetectors$maindetector.r
for(i in 1:dim(trapCov)[2]){
  tmp[!is.na(myDetectors$maindetector.r[])] <- trapCov[,i]
  plot(tmp,main=colnames(trapCov)[i])
  plot(France$geometry,add=T)
  
}

tmp[!is.na(myDetectors$maindetector.r[])] <- trapIntercept
plot(tmp,main="Region")
plot(France$geometry,add=T)
dev.off()


## ==== 6. ASSIGN DETECTIONS TO DETECTORS ====
DNA$Year <- DNA$saisonyear
myData.alive <- AssignDetectors_v3sf( myData = DNA
                                      ,                
                                      myDetectors = myDetectors$main.detector.sf
                                      ,
                                      mysubDetectors = myDetectors$detector.sf
                                      ,
                                      radius = myVars$DETECTORS$detResolution)

## ==== 7. MAKE Y ==== 
y.ar <- MakeYsf( myData = myData.alive$myData.sp,
                 myDetectors = myDetectors$main.detector.sf,
                 method = "Binomial",
                 returnIdvector = TRUE)

## ==== 8. INDIVIDUAL COVARIATES ====
## ====   8.1 SEX ====
sex <- 0
i=107
for(i in 1:length(y.ar$Id.vector)){
  tmp <- DNA[DNA$Id %in% y.ar$Id.vector[i],]
  sex[i]  <- unique(tmp$SEXE)
}
# sex[sex%in% "XX"] <- "F"
# sex[sex%in% "XY"] <- "M"
#LET MODEL ASSIGN THE SEX OF THE INDIVIDUAL
sex[sex%in% "NI"] <- NA
sex[sex%in% "#N"] <- NA

sex[sex%in% "F"] <- "0"
sex[sex%in% "M"] <- "1"
sex <- as.numeric(sex)
#CHECK IT 
table(sex,useNA="always")

## ====   8.2 PREVIOUS DETECTION ====
#check if individuals was detected during the previous years 4 years
load(file.path(WD,"2025","Data","DNAAllPev.RData"))
#SOME CHECKS 
tab1 <- table(DNAAllPev$Id,DNAAllPev$saisonyear)
whichDets <- apply(tab1,1,function(x) sum(x)>1)#ids detected at least during two winters previous to 2021-2023
idGen <- names(whichDets)[whichDets]#unique(DNA1$GENOTYPE_TRANS)
PrevDets <- as.numeric(y.ar$Id.vector %in% idGen)
sum(PrevDets)/length(PrevDets)#PROPORTION OF IDS DETECTED THIS YEAR VS LAST YEAR

#mean number of detections 
mean(rowSums(y.ar$y.ar[y.ar$Id.vector %in% idGen,]))


## ==== 9. CHECK PATTERNS OF INDIVIDUAL DETECTIONS ====
distances <- CheckDistanceDetectionsV2sf( y = y.ar$y.ar, 
                                          detector.xy = myDetectors$main.detector.sf, 
                                          max.distance = 40000,
                                          method = "pairwise",
                                          plot.check = F)

#EXLCUDE DETECTIONS OF INDIVIDUALS DOING "LARGE MOVEMENT" DURING THE SAMPLING. 
par(mfrow=c(1,1),mar=c(1,1,1,1))
if(sum(distances$y.flagged) > 0){
  affected.ids <- which(apply(distances$y.flagged,1,sum)>0)
  for(i in affected.ids){
    plot(st_geometry(myDetectors$main.detector.sf),pch=16, cex=0.1)
    plot(st_geometry(France), add = T,border="red")
    
    tmp <- DNA[DNA$Id == y.ar$Id.vector[i], ]
    tmp <- tmp[order(tmp$date), ]
    tmp.xy <- st_coordinates(tmp)
    n.det <- nrow(tmp.xy)
    
    plot(st_geometry(tmp),add=T, col = "orange", pch = 16, cex = 1)
    
    arrows(x0 = tmp.xy[1:(n.det-1),1], y0 = tmp.xy[1:(n.det-1),2],
           x1 = tmp.xy[2:n.det,1], y1 = tmp.xy[2:n.det,2], length = 0.1, lwd = 1,col="orange")
    
    plot(st_geometry(myDetectors$main.detector.sf[which(y.ar$y.ar[i,] > 0), ]), pch = 16, col = "red",add=T)
    
    tmp2 <- myDetectors$main.detector.sf[which(y.ar$y.ar[i,] > 0 & distances$y.flagged[i,] == 1), ]
    plot(st_geometry(tmp2), col = "blue", pch = 13, cex = 1.5, lwd = 1,add=T)
  }#i
}#if
#if plot.check
y.ar$y.ar <- y.ar$y.ar* (1-distances$y.flagged)


## ==== 10. DATA AUGMENTATION ====
y.alive <- MakeAugmentation(y = y.ar$y.ar, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = 0)
SEX <- MakeAugmentation(y = sex, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = NA)
PrevDets <- MakeAugmentation(y = PrevDets, aug.factor = myVars$DETECTIONS$aug.factor, replace.value = NA)

## ==== 11. PREPARE NIMBLE SCR OBJECTS ====
## ====   11.1 SCALED DETECTORS ====
ScaledDetectors <- scaleCoordsToHabitatGrid(coordsData = myDetectors$detector.xy,
                                            coordsHabitatGridCenter = habitatxy,
                                            scaleToGrid =T )

## ====   11.2 GET WINDOW COORDINATES ====
ScaledLowUpCoords <- getWindowCoords(scaledHabGridCenter = ScaledDetectors$coordsHabitatGridCenterScaled,
                                     scaledObsGridCenter = ScaledDetectors$coordsDataScaled)


habitat.mx <- ScaledLowUpCoords$habitatGrid 
habitat.mx[habitat.mx[]>0]<-1

## ====   11.3 GET LOCAL OBJECTS ====
LocalDetectors <- getLocalObjects(habitatMask = habitat.mx,
                                  coords = ScaledDetectors$coordsDataScaled,
                                  dmax =  7,
                                  resizeFactor = 1,
                                  plot.check = F
)
## ====   11.4 GET SPARSE OBJECTS ====
ySparse <- getSparseY(y.alive)


## ==== 12. DEFINE MODEL CODE ====
modelCode <- nimbleCode({
  ##---- SPATIAL PROCESS 
  ## Prior for AC distribution parameter
  for(i in 1:nhabCov){
    habCoeffSlope[i] ~ dunif(-10,10)
  }
  
  habIntensity[1:numHabWindows] <- exp(habCovs[1:numHabWindows,1:nhabCov] %*% habCoeffSlope[1:nhabCov])
  
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  ## AC distribution
  for(i in 1:M){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],
      numGridRows =  numGridRows,
      numGridCols = numGridCols
    )
  }
  
  ##---- DEMOGRAPHIC PROCESS
  ## Prior for data augmentation
  psi ~ dunif(0,1)
  probMale ~ dunif(0,1)
  teta ~ dunif(0,1)
  ## Data augmentation
  for (i in 1:M){
    z[i] ~ dbern(psi)
    indCov[i,1] ~ dbern(teta)#prev dets
    indCov[i,2] ~ dbern(probMale)#sex
    
    
  }
  
  ##---- DETECTION PROCESS
  ## Priors for detection parameters
  for(c in 1:nCovTraps){
    betaTraps[c] ~ dunif(-5,5)
  }
  
  sigma ~ dunif(0, 50)
  
  for(c in 1:nCounties){
    p0[c] ~ dunif(0,0.50)
  }
  
  
  for(c in 1:2){
    indBetas[c] ~ dunif(-5,5)
  }
  
  ## Detection process
  for (i in 1:M){
    y[i, 1:lengthYCombined] ~  dbinomLocal_normalCovsis(size = trials[1:n.traps],
                                                        p0Traps = p0[1:nCounties],
                                                        sigma = sigma,
                                                        s = sxy[i,1:2],
                                                        trapCoords = trapCoords[1:n.traps,1:2],
                                                        localTrapsIndices = trapIndex[1:n.cells,1:maxNBDets],
                                                        localTrapsNum = nTraps[1:n.cells],
                                                        resizeFactor = ResizeFactor,
                                                        habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                        indicator = z[i],
                                                        lengthYCombined = lengthYCombined,
                                                        allowNoLocal = 0,
                                                        trapCovs =  trapCov[1:n.traps,1:nCovTraps],
                                                        trapCovsIntercept =  trapIntercept[1:n.traps],
                                                        trapBetas = betaTraps[1:nCovTraps],
                                                        indCov = indCov[i,1:2],
                                                        indBetas = indBetas[1:2]
                                                        
    )
    
  }
  
  ##---- DERIVED QUANTITIES
  ## Number of individuals in the population
  N <- sum(z[1:M])
  #density
  dens[1:numHabWindows] <- calculateDensity(s=sxy[1:M,1:2],
                                            habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],
                                            indicator = z[1:M],
                                            numWindows = numHabWindows,
                                            nIndividuals = M )
  #get N within the french border
  for(i in 1:numHabWindowsFR){
    densFr[i] <-  dens[HabWindowsFr[i]]
  }
  NFrance <- sum(densFr[1:numHabWindowsFR])
  
  
})
## ==== 13. INITIAL VALUES ====
## ====   13.1 Z ====
z <- ifelse(apply(y.alive, 1,sum)>0,1,NA)
z <- z
zinits <- ifelse(!is.na(z),NA,rbinom(sum(is.na(z)),1,0.5))

## ====   13.2 SEX ====
sexInits <- ifelse(!is.na(SEX),NA,rbinom(sum(is.na(SEX)),1,0.5))
PrevDetsInits <- ifelse(!is.na(PrevDets),NA,rbinom(sum(is.na(PrevDets)),1,0.5))

## ====   13.3 S ====
#IDENTIFY THE AUGMENTED IDS
idAugmented <- which(names(z) %in% "Augmented")

#sxy 
xy <- cbind(myData.alive$myData.sp$x, myData.alive$myData.sp$y)
colnames(xy) <- c("x","y")
ScaledCoords <- scaleCoordsToHabitatGrid(coordsData =  xy,
                                         coordsHabitatGridCenter = habitatxy,
                                         scaleToGrid =T )
myData.alivetmp <- myData.alive$myData.sp
myData.alivetmp$x <- ScaledCoords$coordsDataScaled[,1]
myData.alivetmp$y <- ScaledCoords$coordsDataScaled[,2]


#ADD A FAKE YEAR TO GET getSInits TO WORK
tmp1 <- myData.alivetmp[1:10,]
tmp1$Year <- as.numeric(tmp1$Year) +1

# GET INITIAL SXY VALUES 
sxy.init <- getSInits( AllDetections = rbind(myData.alivetmp,tmp1),
                       Id.vector = y.ar$Id.vector,
                       idAugmented = idAugmented,
                       lowerCoords = as.matrix(ScaledLowUpCoords$lowerHabCoords),
                       upperCoords = as.matrix(ScaledLowUpCoords$upperHabCoords),
                       habitatGrid = ScaledLowUpCoords$habitatGrid,
                       intensity = NULL,
                       sd = 4,
                       movementMethod = "dbernppACmovement_normal"
                       
)


## ====   13.4 SOTHER INITS VALUES ====
nimInits <- list( "sxy" = sxy.init[,,1],
                  "z" = zinits,
                  "sigma" = runif(1,1,1.1),
                  "habCoeffSlope" = runif(dim(habCovs)[2],-0.1,0.1),#[CM]#0,
                  "p0" = runif(max(trapIntercept),0.2,0.4),#array(runif(max(trapIntercept),0.2,0.4),c(max(trapIntercept),2)),#[CM]rep(0,dim(detCovs)[3]),
                  "betaTraps"  = runif(dim(trapCov)[2], 0.4, 0.5),
                  "psi" = runif(1,0.4,0.5),#,
                  "probMale"= runif(1,0.4,0.5),
                  "teta" = runif(1,0.4,0.5),
                  "indCov" = cbind(PrevDetsInits,sexInits),
                  "indBetas" = runif(2,0.2,0.3)
)

## ==== 14. NIMDATA ====
nimData <- list(y = ySparse$yCombined[,,1],
                trapCov = trapCov,
                habCovs = habCovs,
                z = z,
                indCov = cbind(PrevDets,SEX),
                trapIntercept = trapIntercept
)

# get the cells that in France to extract abundance. 
habitatBuffer <- habitat.r 
id <- raster::extract(habitatBuffer, areaSearched, cellnumbers=T)[[1]][,1]
habitatBuffer[id] <- 2
plot(habitatBuffer)

## SET RESOLUTION FOR EXTRACTION 
habDensity.r <- raster::disaggregate(habitatBuffer, fact=1)
habDensity.r[habDensity.r%in% 0] <- 1
cellInFrance <- which(raster::extract(habDensity.r,habitatxysf)>1)


## ==== 15. NIMCONSTANTS ====
nimConstants <- list( M = dim(ySparse$yCombined)[1],
                      nhabCov= dim(nimData$habCovs)[2],
                      numHabWindows = dim(ScaledLowUpCoords$upperHabCoords)[1],
                      nCounties = max(trapIntercept),
                      habitatGrid = ScaledLowUpCoords$habitatGrid,
                      numGridRows = dim(ScaledLowUpCoords$habitatGrid)[1],
                      numGridCols = dim(ScaledLowUpCoords$habitatGrid)[2],
                      lowerHabCoords = as.matrix(ScaledLowUpCoords$lowerHabCoords),
                      upperHabCoords = as.matrix(ScaledLowUpCoords$upperHabCoords),
                      y.maxDet = dim(LocalDetectors$habitatGrid)[1],
                      x.maxDet = dim(LocalDetectors$habitatGrid)[2],
                      ResizeFactor = LocalDetectors$resizeFactor,
                      n.cells = dim(LocalDetectors$localIndices)[1],
                      maxNBDets = LocalDetectors$numLocalIndicesMax,
                      trapIndex = LocalDetectors$localIndices,
                      nTraps = LocalDetectors$numLocalIndices,
                      habitatIDDet = LocalDetectors$habitatGrid,
                      lengthYCombined = ySparse$lengthYCombined,
                      trials = myData.alive$n.trials[[1]],
                      trapCoords = ScaledDetectors$coordsDataScaled,
                      n.traps=dim(ScaledDetectors$coordsDataScaled)[1],
                      nCovTraps=dim(trapCov)[2],
                      HabWindowsFr = cellInFrance,
                      numHabWindowsFR = length(cellInFrance)
)

## ==== 16. BUILD AND FIT MODEL  ====
## ==== 17. NIMPARAMETERS  ====
nimParams <- c("N","sigma","psi",
               "p0", "habCoeffSlope","betaTraps","teta","probMale","indBetas","NFrance")
# SAVE LESS ITERATIONS FOR AC LOCAITON AND INDIVIDUAL STATES
nimParams2 <- c("sxy","z")

# sum(nimData$trapCov)
# sum(nimData$habCovs)
# sum(nimData$y)
# sum(nimData$indCov,na.rm = T)
# sum(nimConstants$lowerHabCoords,na.rm = T)
# sum(nimConstants$habitatGrid,na.rm = T)
# sum(nimConstants$trapCoords,na.rm = T)
# sum(nimConstants$trapIndex,na.rm = T)
# sum(nimConstants$,na.rm = T)

## ==== 18. SAVE MODEL  ====
for(c in 1:4){
  save(nimData,
       nimConstants,
       nimParams,
       nimParams2,
       modelCode,
       nimInits,
       file = file.path(myVars$WD,"Output","2025", myVars$modelName,
                        paste(unlist(strsplit(myVars$modelName, '_'))[1],"Chain", c, ".RData", sep = "")))
}

## ==== 19. BUILD AND FIT MODEL  ====
## ====   19.1 BUILD AND CACULATE  ====
## ==== 19. BUILD AND FIT MODEL  ====
## ====   19.1 BUILD AND CACULATE  ====
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
model$calculate()
model$initializeInfo()
## ====   19.2 SOME CHECKS AND THE REST  ====
cmodel <- compileNimble(model)
MCMCconf <- configureMCMC(model = model,
                          monitors  = nimParams,
                          control = list(reflective = TRUE),
                          monitors2= nimParams2,
                          thin2 = 10,
                          thin = 1)
MCMC <- buildMCMC(MCMCconf)

cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)

## Run MCMC
MCMCRuntime <- system.time(samples <- runMCMC( mcmc = cMCMC,
                                               nburnin = 5000,
                                               niter = 37500,#run for longer# 500 iterations= 10 mins(M1, RAM= 64Go)
                                               nchains = 4,
                                               samplesAsCodaMCMC = TRUE))
### to access to full posterior, you can download them here: 
#https://sdrive.cnrs.fr/s/oCPC6mdeNtemXjZ
#Once loaded they can be loaded:
load("C:/Personal_Cloud/OneDrive/Work/CNRS/SCR/Output/2025/s54_[W,Hum,Gras,Agri,F]sigma[sex]p0[reg[dpt1],geaco,snow,rds.indcov[prev,sex]2025LargeHab/Posteriors.RData")
load("Posteriors.RData")

#process the nimbleoutput to get array of posteriors
myResults <- ProcessCodaOutput(samples$samples,params.omit = c("sxy","z"))
myResultsSZ <- ProcessCodaOutput(samples$samples2,params.omit = c("sxy","z"))


#check convergence 
#CHECK RHAT
myResults$Rhat
# visual check convergence of parameters 
basicMCMCplots::chainsPlot(samples$samples,var=c("N","betaTraps[1]" ,"betaTraps[2]","betaTraps[3]",  "habCoeffSlope[1]",
                                                 "habCoeffSlope[2]","habCoeffSlope[3]" ,"habCoeffSlope[4]" ,"indBetas[1]","indBetas[2]"))

basicMCMCplots::chainsPlot(samples$samples,var=c("p0[1]", "p0[2]", "p0[3]", "p0[4]",  "p0[5]",
                                                 "p0[6]", "p0[7]", "p0[8]", "p0[9]",  "probMale",
                                                 "psi", "sigma", "teta"))


# Extract abundance only within the french boundary
# Compile function to do it outside of nimble model from sxy and z
CcalculateDensity <- compileNimble(calculateDensity)

dens <- list()
NFrance <- 0
#loop through each iterations
for(i in 1:dim(myResultsSZ$sims.list$sxy)[1]){
  dens[[i]] <- CcalculateDensity(s = myResultsSZ$sims.list$sxy[i, ,], 
                                 habitatGrid = nimConstants$habitatGrid, 
                                 indicator = myResultsSZ$sims.list$z[i, ],
                                 numWindows = nimConstants$numHabWindows,
                                 nIndividuals = nimConstants$M)  
  NFrance[i] <-  sum(dens[[i]][nimConstants$HabWindowsFr])
}

## abundance estimates
mean(NFrance)# 1081.998
quantile(NFrance, probs=c(0.025,0.975))
#
##2.5% 97.5% 
##989  1187

## ====     20.5.5 P0 #### 
det.r <- myDetectors$maindetector.r

rdet1 <- det.rp0 <- rdet2 <- rdet3  <- raster(det.r) 
det.rp0[!is.na(rdet1[])]  <- nimData$trapIntercept
det.pol <- sf::st_as_sf(stars::st_as_stars(det.rp0), 
                        as_points = FALSE, merge = F)
det.pol <- det.pol %>%    group_by(layer) %>%summarize()

nregions <- dim(myResults$sims.list$p0)[2]
par(mfrow=c(1,2),mar=c(5,5,0.5,1))
plot(-1000,xlim=c(0, nregions+1), ylim=c(0,0.002),ylab="p0",xlab="Region",xaxt="n")
regNames <- 0

for(c in 1:nregions){
  regNames[c] <- DetRegion[which(trapIntercept %in% c)[1]]
  
  #axis(1, at=c(1:dim(habCovs)[2]), labels = colnames(habCovs) )
  tmp <- myResults$sims.list$p0[,]
  col <- viridis::viridis(nregions)# c("red","blue","green4","purple")
  plotBars(x=tmp[,c],
           at = c,
           quantile = c(0.0275, 0.975),
           quantile1 = c(0.25, 0.75),
           widthBar = 0.25,
           col = col[c],alpha= 0.5,
           alpha1= 0.8
  )
}

labels <- Departement$code_insee[as.numeric(regNames)]
labels[labels%in% NA] <- c("A","B","C")
axis(1, at=c(1:nregions), labels = labels,cex.axis=0.7 )

#PLOT REGIONS 
par(mar=c(0,1,1,0))

plot(st_geometry(habitatExtent))
plot(st_geometry(France),add=T)
plot(myDetectors$main.detector.sf[,]$geometry,
     col=col[trapIntercept],cex=0.8,add=T,pch=16)


## ====     20.5.2 SIGMA #### 
par(mar=c(5,5,5,5))
plot(-1000,xlim=c(0,3), ylim=c(0,5000),ylab="Sigma (m)",xlab="Sex",xaxt="n")
axis(1,at=c(1:2),labels = c("F","M"))
tmp <- myResults$sims.list$sigma*res(habitat.r)[1]
plotBars(x=tmp,
         at = 1,
         quantile = c(0.0275, 0.975),
         quantile1 = c(0.25, 0.75),
         widthBar = 0.1,
         col = "red",alpha= 0.5,
         alpha1= 0.8
)

## ====     20.5.3 DETECTOR COVARIATE #### 
par(mar=c(5,5,5,5))
plot(-1000,xlim=c(0,dim(trapCov)[2]+1), ylim=c(-1,1),ylab="BetaTraps",xlab="",xaxt="n")
abline(h=0)
axis(1,at=c(1:dim(trapCov)[2]),labels = colnames(trapCov))
tmp <- myResults$sims.list$betaTraps
col <- c("red","blue","green4")
for(c in 1:dim(trapCov)[2]){
  plotBars(x=tmp[,c],
           at = c,
           quantile = c(0.0275, 0.975),
           quantile1 = c(0.25, 0.75),
           widthBar = 0.1,
           col = col[c],alpha= 0.5,
           alpha1= 0.8
  )
}

## ====     20.5.4 DENSITY COVARIATE #### 
par(mar=c(5,5,5,5))
plot(-1000,xlim=c(0, dim(habCovs)[2]+1), ylim=c(-2,2),ylab="BetaHab",xlab="",xaxt="n")
abline(h=0)


axis(1, at=c(1:dim(habCovs)[2]), labels = colnames(habCovs) )
tmp <- myResults$sims.list$habCoeffSlope
col <- c("red","blue","green4","purple")
for(c in 1:dim(habCovs)[2]){
  plotBars(x=tmp[,c],
           at = c,
           quantile = c(0.0275, 0.975),
           quantile1 = c(0.25, 0.75),
           widthBar = 0.1,
           col = col[c],alpha= 0.5,
           alpha1= 0.8
  )
}

## ====     20.5.6 PROB DETECTED  #### 
par(mfrow=c(1,1),mar=c(5,5,0.5,1))
plot(-1000,xlim=c(0, 2), ylim=c(0,1),
     ylab="teta",
     xlab="",
     xaxt="n")
tmp <- myResults$sims.list$teta
col <- c("red","blue","green4")
plotBars(x=tmp,
         at = 1,
         quantile = c(0.0275, 0.975),
         quantile1 = c(0.25, 0.75),
         widthBar = 0.2,
         col = col[1],alpha= 0.5,
         alpha1= 0.8
)
## ====     20.5.7 EFFECT PREV DETE ON DETECTION  #### 
plot(-1000,xlim=c(0, 2), ylim=c(-1,1),
     ylab="IndCov(PrevCovs)",
     xlab="",
     xaxt="n")
abline(h=0)
tmp <- myResults$sims.list$indBetas
col <- c("red","blue","green4")
plotBars(x=tmp[,1],
         at = 1,
         quantile = c(0.0275, 0.975),
         quantile1 = c(0.25, 0.75),
         widthBar = 0.2,
         col = col[1],alpha= 0.5
)

## ====     20.5.8 PROB MALE #### 
par(mfrow=c(1,1),mar=c(5,5,0.5,1))
plot(-1000,xlim=c(0, 2), ylim=c(0,1),
     ylab="probMale",
     xlab="",
     xaxt="n")
tmp <- myResults$sims.list$probMale
col <- c("red","blue","green4")
plotBars(x=tmp,
         at = 1,
         quantile = c(0.0275, 0.975),
         quantile1 = c(0.25, 0.75),
         widthBar = 0.2,
         col = col[1],alpha= 0.5,
         alpha1= 0.8
)
## ====     20.5.9 EFFECT OF SEX ON P0#### 
plot(-1000,xlim=c(0, 2), ylim=c(-1,1),
     ylab="IndCov(sex)",
     xlab="",
     xaxt="n")
abline(h=0)
tmp <- myResults$sims.list$indBetas
col <- c("red","blue","green4")
plotBars(x=tmp[,2],
         at = 1,
         quantile = c(0.0275, 0.975),
         quantile1 = c(0.25, 0.75),
         widthBar = 0.2,
         col = col[1],alpha= 0.5
)
