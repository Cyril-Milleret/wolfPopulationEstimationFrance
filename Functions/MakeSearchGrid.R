#' @title Function to create detectors
#' #'
#' @description
#' \code{MakeSearchGrid} creates  detectors within a defined area (data) and resolution between detectors. 
#' Division precises if subspatial division should be performed (Following PAB model).  
#' 
#' @param data A \code{sp} object. can be a \code{SpatialPolygons} or \code{SpatialPoints} or \code{SpatialLines} or \code{raster}. Spatial Lines objects takes time to compute
#' @param resolution Numeric variable denoting the size of grid cells in units of (\code{polygon}).
#' @param div Numeric variable denoting the number of equally-sized subplots grid cell is to be divided into.
#' @param center Logical variable denoting whether nodes of the resulting search grid are to be centered within their respective subplots.
#' @param plot Logical for whether \code{True} or not (\code{FALSE}) check plots are to be generated.
#' 
#' @author Richard Bischof, \email{richard.bischof@@nmbu.no}
#' @backref R/MakeSearchGrid.R
#' @keywords simul
#'
#' @examples
#' # Make a nested search grid:
#' MakeSearchGrid(...)

MakeSearchGrid <- function(subdetector.r,
                             detResolution=5000,
                             plot = TRUE){
  
  ## if the data is already a raster   
  subdetector <- subdetector.r
  
  
  ## Obtain the resolution of the main detector 
  fact <- detResolution/res(subdetector)[1] 

  ### ==== AGGREGATE TO MAIN DETECTORS  ====
  if(fact>1){
    maindetector.r <- raster::aggregate(subdetector, fact = fact, fun = sum)
  }else{
    maindetector.r <- subdetector
  }
  
  
  ### ==== CREATE POLYGONS FROM RASTER  ====
  ## Main detectors 
  temp.r <- maindetector.r
   #temp.r[!is.na(temp.r)] <- 1:length(temp.r[sum(!is.na(temp.r))])
   temp.r[] <- 1:length(temp.r[])
   
  maindetector.poly <- sf::st_as_sf(stars::st_as_stars(temp.r), 
                                    as_points = FALSE, merge = TRUE)
  
  
  ## Sub-detectors
  temp.r <- subdetector
  #temp.r[!is.na(temp.r)] <- 1:length(temp.r[sum(!is.na(temp.r))])
  temp.r[] <- 1:length(temp.r[])
  
  subdetector.poly <- sf::st_as_sf(stars::st_as_stars(temp.r), 
                                   as_points = FALSE, merge = TRUE)
  
  # plot(habitat.r)
  # plot(subdetector.poly[1:50,]$geometry,add=T,col="red")
  # 
  ### ==== OBTAIN SPATIALPOINTS FROM DETECTORS ====
  ## Main detectors 
  main.detector.xy <- data.frame(crds(maindetector.r,na.rm=F))#[!is.na(maindetector.r[])[,1],]#[!is.na(maindetector.r[])[,1],] #as.data.frame(xyFromCell(maindetector.r, 1:ncell(maindetector.r)))
  main.detector.sf <-  st_as_sf(main.detector.xy, coords = c("x", "y"))
  st_crs(main.detector.sf) <- st_crs(subdetector.r)
  main.detector.sf$main.cell.x <- main.detector.xy[,"x"]
  main.detector.sf$main.cell.y <- main.detector.xy[,"y"]
  main.detector.sf$main.cell.id <- 1:nrow(main.detector.sf)
  
  ## Sub-detectors 
  #detector.xy <- as.data.frame(xyFromCell(subdetector, 1:ncell(subdetector)))
  detector.xy <- data.frame(crds(subdetector,na.rm=F))#[!is.na(maindetector.r[])[,1],] #as.data.frame(xyFromCell(maindetector.r, 1:ncell(maindetector.r)))
  
  sub.detector.sf <-  st_as_sf(detector.xy, coords = c("x", "y"))
  st_crs(sub.detector.sf) <- st_crs(subdetector.r)
 
  
  # col(sub.detector.sf) <- c("x","y")
  sub.detector.sf$Id <- 1:length(sub.detector.sf)
  
  ### ==== SELECT ACTIVE DETECTORS (SEARCHED) ====
  sub.detector.sf$main.cell.id <- as.numeric(st_intersects(sub.detector.sf, maindetector.poly))#over(sub.detector.sf, maindetector.poly)
  
 # merge.df <- st_join(sub.detector.sf, main.detector.sf, by="main.cell.id")
  #colnames(merge.df)[which(colnames(merge.df) %in%"main.cell.id.x" )] <- "main.cell.id"
  
  #merge.df <- merge.df[order(merge.df$Id), ]
  #sub.detector.sf <- merge.df
  
  
sub.detector.sf[sub.detector.sf$main.cell.id%in%4835,]
  
  ## Subset main detectors 
  values.main.r <- terra::extract(maindetector.r, main.detector.sf)
  mainDet <- main.detector.sf$main.cell.id[!is.na(values.main.r[,2])]
  # values.main.r[is.na(values.main.r[,2])] <- 0
   main.detector.sf <- main.detector.sf[!is.na(values.main.r[,2]), ]
  maindetector.poly <- maindetector.poly[!is.na(values.main.r[,2]),]
  ## Subset sub detectors 
  # values.sub.r <- terra::extract(subdetector, sub.detector.sf)
  #values.sub.r[is.na(values.sub.r[,1])] <- 0
  sub.detector.sf <- sub.detector.sf[sub.detector.sf$main.cell.id %in%mainDet, ]
  subdetector.poly <- subdetector.poly[sub.detector.sf$Id,]
  
  
  
  ### ==== MAKE A PLOTTING FUNCTION ====
  if(plot){
     plot(subdetector.r)
    plot(st_geometry(maindetector.poly), add=TRUE, lwd=3)
    plot(st_geometry(subdetector.poly), add=TRUE)
    plot(st_geometry(sub.detector.sf), pch=19, cex=0.6, col=as.numeric(sub.detector.sf$main.cell.id), add=TRUE)
    points(main.cell.y ~ main.cell.x, data=sub.detector.sf, pch=19, cex=1, col=as.numeric(sub.detector.sf$main.cell.id))
  }
  
  ### ==== OUTPUT ====
  out <- list( detector.sf = sub.detector.sf,
               main.detector.sf = main.detector.sf,
               sub.grid.poly = subdetector.poly,
               grid.poly = maindetector.poly,
               maindetector.r = maindetector.r)
}


