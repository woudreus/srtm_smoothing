library(raster)

#top left coordinate with accurate decimal places
TLxcoordinate <- 277461.51
TLYcoordinate <- 664838.56
#windowsize 
windowsize <- 2999
#movingwindow should overlap 
movingwindowsize <- 2500
# use accurate decimal places
pixelsize <- 30.8048

#input rasters
raster1 <- raster("C:\\Users\\playfairs\\Desktop\\noisevar\\cl_SRTM.tif")
raster2 <- raster("C:\\Users\\playfairs\\Desktop\\noisevar\\nwExt_res_cl_log_lee_int_palsarHH.tif") 
raster3 <- raster("C:\\Users\\playfairs\\Desktop\\noisevar\\cl_Hansen.tif") 
raster4 <- raster("C:\\Users\\playfairs\\Desktop\\noisevar\\GLC30m_bin.tif") 
raster5 <- raster("C:\\Users\\playfairs\\Desktop\\noisevar\\nwExt_res_cl_treeheight.tif") 



maxextenty <- (TLYcoordinate - (raster1@nrows) * pixelsize) + pixelsize
maxextentx <- (TLxcoordinate + (raster1@ncols) * pixelsize) + pixelsize

for( i in 6: ceiling((raster1@ncols/windowsize )) ){

windowextentX <- i * (movingwindowsize * pixelsize)

for( j in 0: ceiling((raster1@nrows/windowsize )) ){

  #load DLL 
  dyn.load("C:\\Users\\playfairs\\Desktop\\CUDA\\treeOffset\\x64\\Debug\\treeOffset.dll")
  
windowextentY <- j * (movingwindowsize * pixelsize)

print(paste("this is row i", i, "column j", j))



#Input window extent  
#this method ensures that the input extent window is always a fixed ize
currentExtent    <- NULL
currentExtent[1] <-   TLxcoordinate + windowextentX
currentExtent[2] <-  (TLxcoordinate + ( (windowsize+1) *pixelsize) ) + windowextentX
currentExtent[3] <-   TLYcoordinate - windowextentY
currentExtent[4] <-  (TLYcoordinate - ( (windowsize+1) *pixelsize) ) - windowextentY

#fix raster extents
if(i == ceiling((raster1@ncols/windowsize )) ){
  currentExtent[1] <- maxextentx -  ( (windowsize+2) *pixelsize) 
  currentExtent[2] <- maxextentx 
}

if(j == ceiling((raster1@nrows/windowsize)) ){
  currentExtent[3] <- maxextenty + ( (windowsize+1) *pixelsize)
  currentExtent[4] <- maxextenty  
}

currentExtent <- as(extent(currentExtent), "SpatialPolygons")

#output window extent
#this method ensures that the output pixel size is equal to the source pixelsize 
outputExtent <- NULL
outputExtent[1] <- TLxcoordinate + windowextentX
outputExtent[2] <- (TLxcoordinate + ( (windowsize) *pixelsize) ) + windowextentX
outputExtent[3] <- TLYcoordinate - windowextentY
outputExtent[4] <- (TLYcoordinate - ( (windowsize) *pixelsize) )  - windowextentY 


#fix raster extents
if(i == ceiling((raster1@ncols/windowsize)) ){
  outputExtent[1] <- maxextentx -  ( (windowsize+2) *pixelsize) 
  outputExtent[2] <- maxextentx 
}

if(j == ceiling((raster1@nrows/windowsize)) ){
  outputExtent[3] <- maxextenty + ( (windowsize+1) *pixelsize)
  outputExtent[4] <- maxextenty  
}

outputExtent <- as(extent(outputExtent), "SpatialPolygons")

#clip window extent
#this method removes edge effects
clipExtent <- NULL
clipExtent[1] <- (TLxcoordinate + ( 50 * pixelsize)) + windowextentX
clipExtent[2] <- (TLxcoordinate + ( (windowsize - 50) * pixelsize) ) + windowextentX
clipExtent[3] <- (TLYcoordinate - ( 50 * pixelsize)) - windowextentY
clipExtent[4] <- (TLYcoordinate - ( (windowsize - 50) * pixelsize) ) - windowextentY 


#fix raster extents
if(i == ceiling((raster1@ncols/windowsize)) ){
  clipExtent[1] <- maxextentx -  ( (windowsize+2) *pixelsize) 
  clipExtent[2] <- maxextentx 
}

if(j == ceiling((raster1@nrows/windowsize)) ){
  clipExtent[3] <- maxextenty + ( (windowsize+1) *pixelsize)
  clipExtent[4] <- maxextenty  
}


clipExtent <- as(extent(clipExtent), "SpatialPolygons")

SRTMcrop <- crop(raster1, currentExtent, snap = 'near')  
PALSHHcrop <- crop(raster2, currentExtent, snap = 'near')  
HANScrop <- crop(raster3, currentExtent, snap = 'near')  
THeightcrop <- crop(raster5, currentExtent, snap = 'near')  
NONFcrop <- crop(raster4, currentExtent, snap = 'near')  

SRTMcrs <- crs(raster1)

#convert inputs to 1D arrays
DTM.1D <- array(as.matrix(SRTMcrop))
PALSHH.1D <- array(as.matrix(PALSHHcrop))
HANS.1D   <- array(as.matrix(HANScrop))
NONF.1D <- array(as.matrix(NONFcrop))
TRH.1D <- array(as.matrix(THeightcrop)) 

#remove NA's
DTM.1D[is.na(DTM.1D)] <- 0 
PALSHH.1D[is.na(PALSHH.1D)]<- 0
HANS.1D[is.na(HANS.1D)]<- 0
NONF.1D[is.na(NONF.1D)]<- 0
TRH.1D[is.na(TRH.1D)]<- 0

#number of elements near input data
n <- length(DTM.1D)
ncolInput <- ncol(SRTMcrop)
nrowInput <- nrow(SRTMcrop)
#load CUDA DLL

#RUN CUDA DLL
resultTreeOffset <-  .C("treeOffset", as.double(DTM.1D), as.double(PALSHH.1D), as.double(HANS.1D), as.double(NONF.1D),
                      as.double(TRH.1D), DEM = double(length= n), as.integer(ncolInput),  as.integer(n))

#summary(resultNoiseVar$aggregatedNoiseVar)
DEM <- raster(matrix(resultTreeOffset$DEM, nrow = nrowInput, ncol = ncolInput))

crs(DEM) <- SRTMcrs

DEM <-  setExtent(DEM, alignExtent(SRTMcrop, DEM, snap = "near"), snap = FALSE, keepres = FALSE)


DEM <- crop(DEM, clipExtent)  


uuid <-  floor(runif(1, 100,999))

writeRaster(DEM, paste0("C:\\Users\\playfairs\\Desktop\\noisevar\\treeoffset\\",uuid,"_", i,"_",j,"DEM.tif"),
            format = "GTiff", overwrite = TRUE, dataType= 'INT2S')

rm(resultTreeOffset)
gc()
#Release CUDA DLL
dyn.unload("C:\\Users\\playfairs\\Desktop\\CUDA\\treeOffset\\x64\\Debug\\treeOffset.dll")
  }
}


.rs.restartR()
