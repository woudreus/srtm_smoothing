library(raster)
#input rasters

#top left coordinate with accurate decimal places
TLxcoordinate <- 282082.24
TLYcoordinate <- 660217.85
#windowsize 
windowsize <- 2999
#movingwindow should overlap 
movingwindowsize <- 2500
# use accurate decimal places
pixelsize <- 30.8048

raster1 <- raster('V:\\GPS - Rapporten\\Playfair\\noisevar\\FOTO\\input_25112017\\srtmutm_corrected.tif') 

maxextenty <- (TLYcoordinate - (raster1@nrows) * pixelsize) + pixelsize
maxextentx <- (TLxcoordinate + (raster1@ncols) * pixelsize) + pixelsize

for( i in 0: ceiling((raster1@ncols/windowsize )) ){

  windowextentX <- i * (movingwindowsize * pixelsize)

  for( j in 0: ceiling((raster1@nrows/windowsize )) ){

       windowextentY <- j * (movingwindowsize * pixelsize)
    
    print(paste("this is row i", i, "column j", j))
    
    #load CUDA DLL
    dyn.load("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\noiseVariance\\x64\\Debug\\noiseVariance.dll")
    
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
      clipExtent[1] <- maxextentx -  ( (windowsize+2) *pixelsize) + ( 50 * pixelsize)
      clipExtent[2] <- maxextentx - ( 50 * pixelsize)
    }
    
    if(j == ceiling((raster1@nrows/windowsize)) ){
      clipExtent[3] <- maxextenty + ( (windowsize+2) *pixelsize)  - ( 50 * pixelsize)
      clipExtent[4] <- maxextenty  + ( 50 * pixelsize)
    }
    
    clipExtent <- as(extent(clipExtent), "SpatialPolygons")

SRTMcrop <- crop(raster1, currentExtent, snap = 'near')  

SRTMcrs <- crs(raster1)

#convert inputs to 1D arrays
DEM.1D <- array(as.matrix(SRTMcrop))

#remove NA's
DEM.1D[is.na(DEM.1D)] <- 0 

#number of elements in input data
n <- length(DEM.1D)
ncolInput <- ncol(SRTMcrop)
nrowInput <- nrow(SRTMcrop)
lengthout <- (as.integer(ncolInput)/5)^2

#RUN CUDA DLL
resultNoiseVar <-  .C("NoiseVar", as.double(DEM.1D), aggregatedNoiseVar = double(length= n),
                      as.integer(n), as.integer(ncolInput), as.integer(lengthout) )

#summary(resultNoiseVar$aggregatedNoiseVar)
resampledNoiseVariance <- raster(matrix(resultNoiseVar$aggregatedNoiseVar, nrow = nrowInput, ncol = ncolInput))

crs(resampledNoiseVariance) <- SRTMcrs

resampledNoiseVariance <-  setExtent(resampledNoiseVariance, alignExtent(SRTMcrop,resampledNoiseVariance, snap = "near"), 
                           snap = FALSE, keepres =FALSE)

resampledNoiseVariance <- crop(resampledNoiseVariance, clipExtent)

uuid <-  floor(runif(1, 100,999))

writeRaster(resampledNoiseVariance, paste0("V:\\GPS - Rapporten\\Playfair\\noisevar\\noisemap\\",uuid,"_", i,"_",j,"_w.tif"),
            format = "GTiff", overwrite = TRUE, dataType= 'INT2S')

dyn.unload("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\noiseVariance\\x64\\Debug\\noiseVariance.dll")
 }
}

#Release DLL
rm(resultNoiseVar)
gc()
.rs.restartR()

