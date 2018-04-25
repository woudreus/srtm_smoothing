library(raster)
#input rasters

#top left coordinate with accurate decimal places
TLxcoordinate <- 280541.997
TLYcoordinate <- 661758.084

#windowsize 
windowsize <-3455

#movingwindow should overlap 
movingwindowsize <- 2500

# use accurate decimal places
pixelsize <- 30.8048

#input random rasters
raster1 <- raster("F:\\GPS - Rapporten\\Playfair\\noisevar\\DEM_Gallant.tif") 
raster2 <- raster("F:\\GPS - Rapporten\\Playfair\\noisevar\\noisemap_v2.tif") 

maxextenty <- (TLYcoordinate - (raster2@nrows) * pixelsize) + pixelsize
maxextentx <- (TLxcoordinate + (raster2@ncols) * pixelsize) + pixelsize

#for( i in 0: ceiling((raster1@ncols/windowsize )) ){
i=3
windowextentX <- i * (movingwindowsize * pixelsize)

#for( j in 0: ceiling((raster1@nrows/windowsize )) ){
j=3
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
  currentExtent[1] <- (maxextentx - (5*pixelsize) ) -  ( (windowsize+1) *pixelsize) 
  currentExtent[2] <- (maxextentx - (5*pixelsize) )
}

if(j == ceiling((raster1@nrows/windowsize)) ){
  currentExtent[3] <- (maxextenty+(5*pixelsize) ) + ( (windowsize+1) *pixelsize)
  currentExtent[4] <- (maxextenty+(5*pixelsize) )  
}

currentExtent <- as(extent(currentExtent), "SpatialPolygons")

#output window extent
#this method ensures that the output pixel size is equal to the source pixelsize 
outputExtent <- NULL
outputExtent[1] <- TLxcoordinate + windowextentX
outputExtent[2] <- (TLxcoordinate + ( (windowsize+1) *pixelsize) ) + windowextentX
outputExtent[3] <- TLYcoordinate - windowextentY
outputExtent[4] <- (TLYcoordinate - ( (windowsize+1) *pixelsize) ) - windowextentY 


#fix raster extents
if(i == ceiling((raster1@ncols/windowsize)) ){
  outputExtent[1] <- (maxextentx - (5*pixelsize) ) -  ( (windowsize+1) *pixelsize) 
  outputExtent[2] <- (maxextentx - (5*pixelsize) )
}

if(j == ceiling((raster1@nrows/windowsize)) ){
  outputExtent[3] <- (maxextenty + (5*pixelsize) ) + ( (windowsize+1) *pixelsize)
  outputExtent[4] <- (maxextenty + (5*pixelsize) ) 
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
  clipExtent[1] <- (maxextentx - (55*pixelsize) ) - (( (windowsize) *pixelsize)  - ( 60 * pixelsize))
  clipExtent[2] <- (maxextentx - (55*pixelsize) ) 
}

if(j == ceiling((raster1@nrows/windowsize)) ){
  clipExtent[3] <- (maxextenty + (55*pixelsize) )  + (( (windowsize) *pixelsize) - ( 60 * pixelsize))
  clipExtent[4] <- (maxextenty + (55*pixelsize) ) 
}

clipExtent <- as(extent(clipExtent), "SpatialPolygons")

#crop rasters####    
zcrop  <-  crop(raster1, currentExtent, snap = 'near')  
v0crop  <-  crop(raster2, currentExtent, snap = 'near')  

vgcrop  <- raster( matrix(0, ncol = 3456, nrow = 3456)) 
ncrop  <-  raster( matrix(1, ncol = 3456, nrow = 3456))

SRTMcrs <- crs(raster1)

iniNcol <- 3456
nArray <- (3456)^2

#limit V0 to 5############
v0.1D <- array(as.matrix(v0crop))

#v0.1D[v0.1D] <- ( ( (v0.1D[v0.1D >= 5]  - min(v0.1D)) / max(v0.1D) ) * 5)+5
v0.1D[v0.1D >=4.5] <- 4.5
#v0.1D <- 13^2*(0.5)^v0.1D
#v0.1D <- array(matrix(2,nrow = 3456, ncol = 3456))
hist(v0.1D, density = 5, breaks = 300, col = "blue" )
summary(v0.1D)

# v0.1D <- ( (v0.1D  

#v0.1D[v0.1D > 5 ] <-  5
#dyn.load("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\cuGauss\\x64\\Debug\\cuGauss.dll")
#resultV0Gauss <- .C("cuGauss", as.double(v0.1D), varianceOut = double(length = nArray), as.integer(iniNcol), as.integer(nArray) )
#v0.1D <- resultV0Gauss$varianceOut
#dyn.unload("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\cuGauss\\x64\\Debug\\cuGauss.dll")

v0crop <- raster(matrix(v0.1D, nrow = iniNcol))

##################

#convert inputs to 1D arrays
z.1D <-   array(as.matrix(zcrop))
w.1D <-   1/v0.1D

wsq.1D <- w.1D^2
vg.1D <-   array(as.matrix(vgcrop))
n.1D <-   array(as.matrix(ncrop))

#remove NA
z.1D[is.na(z.1D)] <- 0 
z.1D[z.1D < 0 ] <- 0 

w.1D[is.na(w.1D)] <- 0
w.1D[w.1D < 0 ] <- 0

wsq.1D[is.na(wsq.1D)] <- 0
wsq.1D[wsq.1D < 0 ] <- 0


writeRaster(zcrop, paste0(tempdir(), "\\Zvalue_res_0.tif"), format = "GTiff", overwrite =TRUE)
writeRaster(v0crop, paste0(tempdir(), "\\Vvalue_res_0.tif"), format = "GTiff", overwrite =TRUE)
#writeRaster(raster(matrix(z.1D, nrow = 3456)), paste0(tempdir(), "\\Zvalue_res_0.tif"), format = "GTiff", overwrite =TRUE)
#writeRaster(raster(matrix(v0.1D, nrow = 3456)), paste0(tempdir(), "\\Vvalue_res_0.tif"), format = "GTiff", overwrite =TRUE)
count <-0



repeat{
  
  count  <- count+1
  if (count>1){ break}
  
  
  #number of elements in input data
  n <- length(z.1D)
  ncolInput <- sqrt(length(z.1D))
  nrowInput <- sqrt(length(z.1D))
  nout <- n/9
  ncolOutput <- sqrt(length(z.1D))/3
  nrowOutput <- sqrt(length(z.1D))/3
  #load CUDA DLL
  dyn.load("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\RCUDA2\\x64\\Debug\\RCUDA2.dll")
  
  
  #RUN CUDA DLL
  result <-  .C("RCUDA2", as.double(z.1D), as.double(w.1D),  as.double(wsq.1D), as.double(vg.1D), as.double(n.1D),
                          as.integer(n), as.integer(ncolInput), as.integer(nout), 
                          outz = double(length = nout), outw = double(length = nout), outwsq = double(length = nout),
                          outn = double(length = nout), outvg = double(length = nout),
                          outneff = double(length = nout), outvm = double(length = nout), outvstat = double(length = nout) )
  
  #Release DLL
  dyn.unload("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\RCUDA2\\x64\\Debug\\RCUDA2.dll")
  
  summary(result$outz)
  
  z.1D   <- array(result$outz)
  w.1D   <- array(result$outw)
  wsq.1D <- array(result$outwsq)
  
  n.1D   <- array(result$outn)
  vg.1D  <- array(result$outvg)
  neff.1D <- array(result$outneff)
  vstat.1D <- array(result$outvstat)
  vm.1D <- array(result$outvm)
  
  summary(z.1D)
  neff.1D[neff.1D <= 1.01] <- 1.01
  chisq <- as.array(apply(neff.1D, 1, function(x){ qchisq(0.05, (x - 1), lower.tail = FALSE)}))
  
  #load dll ####
  dyn.load("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\xcritcalc\\x64\\Debug\\xcritcalc.dll")
  
  
  resultv <- .C("xcrit", as.double(vm.1D), as.double(vg.1D), as.double(vstat.1D), as.double(chisq), outv = double(length = nout),
                as.integer(ncolOutput), as.integer(nout), now = TRUE )
  #release dll
  dyn.unload("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\xcritcalc\\x64\\Debug\\xcritcalc.dll")  
  
  chisq[is.na(chisq)] <- 0
  vm.1D[is.na(vm.1D)] <- 0
  vm.1D[is.infinite(vm.1D)] <- 0
  vstat.1D[is.infinite(vstat.1D)] <- 0
  vstat.1D[is.na(vstat.1D)] <- 0
  vg.1D[is.na(vg.1D)] <- 0
  
  v.1D <- array(resultv$outv)
  
#  writeRaster(raster(matrix(chisq,nrow = nrowOutput)), paste0(tempdir(), "\\chisq_",count,".tif"), format = "GTiff", dataType = "INT2S",overwrite =TRUE)
  
#  writeRaster(raster(matrix(vm.1D,nrow = nrowOutput)), paste0(tempdir(), "\\VM_",count,".tif"), format = "GTiff",dataType = "INT2S", overwrite =TRUE)
  
#  writeRaster(raster(matrix(vstat.1D,nrow = nrowOutput)), paste0(tempdir(), "\\VStat_",count,".tif"), format = "GTiff",dataType = "INT2S", overwrite =TRUE)
  
#  writeRaster(raster(matrix(vg.1D,nrow = nrowOutput)), paste0(tempdir(), "\\Vg_",count,".tif"), format = "GTiff", dataType = "INT2S",overwrite =TRUE)
  
#  writeRaster(raster(matrix(neff.1D,nrow = nrowOutput)), paste0(tempdir(), "\\neff_",count,".tif"), format = "GTiff", dataType = "INT2S",overwrite =TRUE)
  
#  writeRaster(raster(matrix(n.1D,nrow = nrowOutput)), paste0(tempdir(), "\\n_res_",count,".tif"), format = "GTiff",dataType = "INT2S", overwrite =TRUE)
  
#  writeRaster(raster(matrix(z.1D, nrow = nrowOutput)), paste0(tempdir(), "\\z_",count,".tif"), format = "GTiff",dataType = "INT2S", overwrite =TRUE)
  
# writeRaster(raster(matrix(w.1D,nrow = nrowOutput)), paste0(tempdir(), "\\w_",count,".tif"), format = "GTiff", dataType = "INT2S",overwrite =TRUE)
  
#  writeRaster(raster(matrix(wsq.1D,nrow = nrowOutput)), paste0(tempdir(), "\\wsq_",count,".tif"), format = "GTiff",dataType = "INT2S", overwrite =TRUE)
  
 # writeRaster(raster(matrix(z.1D,nrow = nrowOutput)), paste0(tempdir(), "\\v_",count,".tif"), format = "GTiff", dataType = "INT2S",overwrite =TRUE)
  
 # writeRaster(raster(matrix(v.1D,nrow = nrowOutput)), paste0(tempdir(), "\\v_",count,".tif"), format = "GTiff", dataType = "INT2S",overwrite =TRUE)
  #remove NA
  #treat V

 # v.1D[v.1D >= 5] <- (( (v.1D[v.1D >= 5]  - 5) / (max(v.1D) - min(v.1D) )) * 5) + max(v.1D[v.1D < 5 ])
#  v.1D[v.1D >= 100] <- 100
# v.1D[v.1D >= 5] <-4^2*(0.5)^v.1D[v.1D >= 5]
#   v.1D <- 0.5*v.1D^0.25
  hist(v.1D, density = 5, breaks = 300, col = "blue" )
  summary(v.1D)
  
#  dyn.load("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\cuGauss\\x64\\Debug\\cuGauss.dll")
#  resultVGaussRepeat <- .C("cuGauss", as.double(v.1D), vOut = double(length = nout), as.integer(ncolOutput), as.integer(nout) )
#  v.1D <- resultVGaussRepeat$vOut
#  dyn.unload("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\cuGauss\\x64\\Debug\\cuGauss.dll")
  
  ########
  
  z.1D[is.na(z.1D)] <- 0 
  z.1D[z.1D < 0 ] <- 0 
  
  w.1D[is.na(w.1D)] <- 0
  w.1D[w.1D < 0 ] <- 0
  
  wsq.1D[is.na(wsq.1D)] <- 0
  wsq.1D[wsq.1D < 0 ] <- 0
  
  vg.1D[is.na(vg.1D)] <-0
  
  zbar <- raster(matrix(z.1D, nrow = nrowOutput, ncol = ncolOutput))
  vbar <- raster(matrix(v.1D, nrow = nrowOutput, ncol = ncolOutput))
  
  crs(zbar) <- SRTMcrs
  crs(vbar) <- SRTMcrs
  
  zbar <-  setExtent(zbar, alignExtent(zcrop,zbar, snap = "near"), 
                    snap = FALSE, keepres = FALSE)
  
  vbar <-  setExtent(vbar, alignExtent(v0crop,vbar, snap = "near"), 
                    snap = FALSE, keepres =FALSE)
  
  #zbar <- crop(zbar, clipExtent)
  #vbar <- crop(vbar, clipExtent)
  
  writeRaster(zbar, paste0(tempdir(), "\\Zvalue_res_",count,".tif"), format = "GTiff", overwrite =TRUE)
  writeRaster(vbar, paste0(tempdir(), "\\Vvalue_res_",count,".tif"), format = "GTiff", overwrite =TRUE)
  
  
  
}
count <- count - 1
zbar.1D <- array(as.matrix(zbar))
vbar.1D <- array(as.matrix(vbar))

repeat{
  
  
  count <- count - 1
  
  rasterC1 <- raster(paste0(tempdir(), "\\Zvalue_res_",count,".tif"))
  rasterC2 <- raster(paste0(tempdir(), "\\Vvalue_res_",count,".tif"))  
 
  z.1D <- array(as.matrix(rasterC1)) 
  v.1D <- array(as.matrix(rasterC2))
  
  z.1D[is.na(z.1D)] <- 0 
  z.1D[z.1D < 0 ] <- 0 
  
  zbar.1D[is.na(zbar.1D)] <- 0
  zbar.1D[ zbar.1D < 0 ] <- 0
  
  vbar.1D[is.na(vbar.1D)] <- 5
  vbar.1D[vbar.1D <= 0.1636] <- 0.1636
  
  v.1D[is.na(v.1D)] <- 5
  v.1D[v.1D <= 0.1636 ] <- 0.1636
  
  ncolF <- ncol(rasterC1)
  nrowF <- nrow(rasterC1)
  nF <- length(rasterC1)
  
 dyn.load("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\upscaling\\x64\\Debug\\upscaling.dll")
 upscaleResultZbar <- .C("bicubUpscaling", as.double(zbar.1D), zUpscaled = double(length = nF), as.integer(ncolF), as.integer(nF), now = TRUE )
 dyn.unload("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\upscaling\\x64\\Debug\\upscaling.dll")
  
  
 dyn.load("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\upscaling\\x64\\Debug\\upscaling.dll")
 upscaleResultVbar  <- .C("bicubUpscaling", as.double(vbar.1D), vUpscaled = double(length = nF), as.integer(ncolF), as.integer(nF), now = TRUE )
 dyn.unload("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\upscaling\\x64\\Debug\\upscaling.dll")
  
  
  zbar.1D <- array(upscaleResultZbar$zUpscaled)
  vbar.1D <- array(upscaleResultVbar$vUpscaled)
#   zbar.1D <- array(disaggregate(raster(matrix(zbar.1D, ncol = ncolF/3, nrow = nrowF/3)), fact = 3, method ='bilinear'))
#   vbar.1D <- array(disaggregate(raster(matrix(vbar.1D, ncol = ncolF/3, nrow = nrowF/3)), fact = 3, method ='bilinear'))
   zbar.1D[is.na(zbar.1D)] <- 0
   zbar.1D[ zbar.1D < 0 ] <- 0
   
   vbar.1D[is.na(vbar.1D)] <- 5
   vbar.1D[vbar.1D <= 0.1636] <- 0.1636
   
  dyn.load( "C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\gallantResampling\\x64\\Debug\\gallantResampling.dll")
  
  resultResampling <-  .C("resamplingFun", as.double(zbar.1D), as.double(vbar.1D), as.double(z.1D), as.double(v.1D), resampledZs = double(length= nF),
                          resampledVs = double(length= nF), as.integer(nF), as.integer(ncolF), now = TRUE )
  
  rm(zbar.1D)
  rm(vbar.1D)
  zbar.1D <- array(as.matrix(resultResampling$resampledZs))
  vbar.1D <- array(as.matrix(resultResampling$resampledVs))
  
  rm(resultResampling)
  gc()
  #Release DLL
  dyn.unload("C:\\Users\\playfairs\\Desktop\\LOTS_O_FILES\\CUDA\\gallantResampling\\x64\\Debug\\gallantResampling.dll")
  
  if (count== 0){break}
}

DEM <- raster(matrix(zbar.1D, nrow = 3456, ncol = 3456))
#DEM <- (DEM - DEM@data@min)/(raster1@data@max - raster1@data@min+50) * raster1@data@max

crs(DEM) <- SRTMcrs

DEM <-  setExtent(DEM, alignExtent(zcrop, DEM, snap = "in"), snap = FALSE, keepres =FALSE)

DEM <- crop(DEM, clipExtent)  

uuid <-  floor(runif(1, 100,999))

writeRaster(DEM, paste0("F:\\GPS - Rapporten\\Playfair\\noisevar\\gallantSmoothing\\",uuid,"_", i,"_",j,"DEMGS.tif"),
            format = "GTiff", overwrite = TRUE)

 }
}


#rm(list =ls())
#gc()
#.rs.restartR()



