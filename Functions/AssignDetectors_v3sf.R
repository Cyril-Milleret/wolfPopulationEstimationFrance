#' @title Detectors Assignment Function
#'
#' @description
#' \code{AssignDetectorsPD} returns a \code{SpatialPointsDataFrame} object with DNA samples
#'  locations and associated detectors'IDs.
#'
#' @param myData A \code{SpatialPointsDataFrame} object containing the filtered data
#' @param myDetector A \code{list} of SpatialPointsDataFrame objects containing the 
#' detectors locations
#' @param mysubDetectors A \code{list} of SpatialPointsDataFrame objects containing
#'  the subdetectors locations for the binomial models
#'
#' @return A \code{SpatialPointsDataFrame} object with the x and y locations of the
#'  different samples with associated detector's Id 

AssignDetectors_v3sf <- function( myData,                              
                              myDetectors,
                              mysubDetectors = NULL,
                              radius = 20000){

   require(RANN)
   require(dplyr)
   myData$record.id<-1:dim(myData)[1]
   
   ## RETRIEVE NUMBER OF YEARS 
   if(is.null(myData$Year)){myData$Year <- 1}
   years <- min(unique(myData$Year)):max(unique(myData$Year))
   
   ## SET UP THE DETECTORS LIST
   # if(!is.list(myDetectors)){myDetectors <- list(myDetectors)}
   if(!sum(class(myDetectors) %in% "list")==1){myDetectors <- list(myDetectors)}
   
   ## CHECK THAT THE NUMBER OF YEARS MATCH THE LENGTH OF THE DETECTORS LIST
   if(length(myDetectors) != length(years)){
      myDetectors <- lapply(1:length(years), function(x) myDetectors[[1]])
      }
   
   ## IF PAB PROCESS
   if(!is.null(mysubDetectors)){
      ## SET UP THE SUB-DETECTORS LIST
      # if(!is.list(mysubDetectors)){mysubDetectors <- list(mysubDetectors)}
     if(!sum(class(mysubDetectors) %in% "list")==1){mysubDetectors <- list(mysubDetectors)}
     
      if(length(mysubDetectors) != length(years)){
         mysubDetectors <- lapply(1:length(years), function(x) mysubDetectors[[1]])
      }
      ## INITIALIZE DETECTORS AND SUB-DETECTORS ID COLUMNS
      myData$sub.detector <- NA 
      myData$Detector <- NA
      for(t in 1:length(years)){
         
         #--[RB] using sf package:
         myDetectors.sf<- mysubDetectors[[t]]
         myData.sf<- myData[myData$Year == years[t],]

         # get coordinate matrices
         detector_coords <- do.call(rbind, st_geometry(myDetectors.sf))
         detection_coords <- do.call(rbind, st_geometry(myData.sf))
         if(!is.null(detection_coords)){
        
         # fast nearest neighbour search
         closest <- nn2(detector_coords,data.frame(detection_coords[,1],detection_coords[,2]),  k = 1, searchtype = "radius", radius = radius)
         closest$nn.idx[closest$nn.idx==0]<-NA
         myData$sub.detector[myData$Year == years[t]] <-closest$nn.idx#myDetectors.sf[closest$nn.idx,]$Id
         
         main.cell.id <- mysubDetectors[[t]]$main.cell.id[closest$nn.idx]     
         myData$Detector[myData$Year == years[t]] <- unlist(
            lapply(1:length(main.cell.id), function(x){
               out<-NA
               if(!is.na(main.cell.id[x]))out<-which(myDetectors[[t]]$main.cell.id %in% main.cell.id[x])
               return(out)
               }))
           }#if
         }#t
   }else{
      myData$Detector <- NA
      for(t in 1:length(years)){
         #--[RB] using sf package:
         myDetectors.sf<- myDetectors[[t]]
         myData.sf<- myData[myData$Year == years[t],]
         
         # get coordinate matrices
         detector_coords <- do.call(rbind, st_geometry(myDetectors.sf))
         detection_coords <- do.call(rbind, st_geometry(myData.sf))
         
         if(!is.null(detection_coords)){
         # fast nearest neighbour search
         closest <- nn2(detector_coords,data.frame(detection_coords[,1],detection_coords[,2]),  k = 1, searchtype = "radius", radius = radius)

         closest$nn.idx[closest$nn.idx==0]<-NA
         myData$Detector[myData$Year == years[t]] <-closest$nn.idx#myDetectors.sf[closest$nn.idx,]$main.cell.id

         }#if 
      }
   }
   
   myData <- myData[!is.na(myData$Detector),]
   myData$Id <- as.character(myData$Id)
   
   ## Output return
   if(!is.null(mysubDetectors)){                                  
      n.trials <- lapply(mysubDetectors,function(x)table(x$main.cell.id))
      return(list(myData.sp = myData, n.trials = n.trials))
   } else {
      return(myData)
   }
}
