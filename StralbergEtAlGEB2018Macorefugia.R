#This code has five major parts, following Stralberg et al. (2018):
#1. Generate climate-change refugia indices for individual species based on species distribution model predictions
#2. Generate a multi-species climate-change refugia index based on species-specific refugia layers 
#3. Calculate topographic metrics
#4. Run quantile regression models to attribute sources of variation in a refugia index
#5. Conduct a limiting factor analysis based on climate variables for a refugia index

#First install and load the required packages:

library(raster) #raster manipulation
library(rgdal) #spatial data library
library(yaImpute) # for ann function
library(gdata) # for nobs function
library(quantreg) #for qr function
library(dplyr) #for sample_frac function
library(rasterVis) #for raster-based levelplot function; requires lattice, latticeExtra and gridExtra

#Next, define paths, set the working directory, and obtain a list of species with SDMs:

#ref = output directory
#species = SDM directory
#ecor = study area extent raster
#sample = sample SDM raster
#climfut = future climate data brick (CMIP5 from McKenney et al. available on AdaptWest website)
#tfacet = land facet layer available on AdaptWest website
#elev = 10-km resolution elevation raster
#id = 10-km resolution unique ID raster

setwd(species)
spec <- list.dirs(species)
sample <- sample*0 + 1

#1a.
#A key component of the index is a distance-decay function based on a fat-tailed dispersal kernel, 
#which accommodates rare long-distance dispersal events, and has been invoked to explain the rapid 
#post-glacial recolonization of trees across northern North America (Clark et al. 1998). The following 
#two functions are used to calculate (1) an individual probability value [ranging from 0 to 1] from a 
#negative exponential distribution for a given distance (x), alpha, and shape (c); and (2) the mean 
#dispersal distance for a given alpha and c. For a fat-tailed distribution, c = 0.5. 

fattail <- function(x, alpha, c) {
  left <- (c/(2*alpha*gamma(1/c)))
  right <- exp(-1*((abs(x/alpha))^c))
  result <- left*right
  return(right)
}

ftmean <- function(alpha, c) {
  result <-(alpha*gamma(2/c))/gamma(1/c)
  return(result)
}

#1b.
#The following code loops through species, calculates study area prevalence (mean probability of occurrence), 
#writes those values to a data frame, and then generates binary presence/absence rasters for current and future 
#time periods using prevalence as a probability threshold:

prevtab <- data.frame("spec" = gsub(paste(species,"/",sep=""),"",spec), "prev" = 0)
for (i in 2:length(spec)){
  velstack<-emptystack
  p <- list.files(spec[i],pattern="current.tif$")
  pres <- raster(paste(spec[i],"\\",p,sep=""))
  presc <- crop(pres,ecor)
  prescm <- mask(presc,ecor)
  prev <- cellStats(prescm,'mean')
  prevtab[i,2] <- prev
  present<-as.data.frame(rasterToPoints(prescm))
  names(present)[3] <- "prob"
  presprev <- prescm
  values(presprev) <- ifelse(values(presprev) > prev, 1, 0) 
  writeRaster(presprev,filename=paste(spec[i],"presprev",sep=""),format="GTiff", overwrite=TRUE)
  f <- list.files(spec[i],pattern="future.tif$")
  for (j in 1:4) {
    fut <- raster(paste(spec[i],"/",f[j],sep=""))
    futc <- crop(fut,ecor)
    futcm <- mask(futc,ecor)
    future<-as.data.frame(rasterToPoints(futcm))
    names(future)[3] <- "prob"	
    futprev <- futcm
    values(futprev) <- ifelse(values(futprev) > prev, 1, 0) 
    writeRaster(futprev,filename=paste(gsub(paste(species,"/",sep=""),"",f[j]),"futprev",sep=""),format="GTiff", overwrite=TRUE)
  }
}

#1c.
#This core bit of code calculates biotic velocity (Carroll et al. 2015) for each species based on the 
#nearest-analog velocity algorithm defined by Hamann et al. (2015), and then applies the distance-decay 
#function to obtain an index ranging from 0 to 1. For a fat-tailed distribution, c = 0.5, and alpha = 8333.33 
#results in a mean migration rate of 500 m/year (50 km/century). Velocity values are averaged over four GCMs, 
#with NA values first converted to zeros. 

emptystack<-sample
emptystack<-addLayer(emptystack,emptystack)
emptystack<-dropLayer(emptystack,c(1,2))
for (i in 2:length(spec)){
  refstack<-emptystack
  futprevstack<-emptystack
  p <- list.files(spec[i],pattern="presprev.tif$")
  presprev <- raster(paste(spec[i],"//",p,sep=""))
  present<-as.data.frame(rasterToPoints(presprev))
  names(present)[3] <- "prev"
  f <- list.files(spec[i],pattern="futprev.tif")
  for (j in 1:4) {
    futprev <- raster(paste(spec[i],"\\",f[j],sep=""))
    future <- as.data.frame(rasterToPoints(futprev))
    names(future)[3] <- "prev"	
    p.xy<-cbind(seq(1,length(present$x),1),present$x,present$y,present$prev)
    f.xy<-cbind(seq(1,length(future$x),1),future$x,future$y,future$prev)
    p.xy2<-p.xy[p.xy[,4]>0.1,1:3,drop=FALSE]
    f.xy2<-f.xy[f.xy[,4]>0.1,1:3,drop=FALSE]
    if(nrow(f.xy)>0){
      d.ann <- as.data.frame(ann(
        as.matrix(p.xy2[,-1,drop=FALSE]),
        as.matrix(f.xy2[,-1,drop=FALSE]),
        k=1, verbose=F)$knnIndexDist)
      d1b <- as.data.frame(cbind(f.xy2, round(sqrt(d.ann[,2]))))
      names(d1b) <- c("ID","X","Y","bvel")
    } else {
      print(spec[i])
    }
    f.xy <- as.data.frame(f.xy)
    names(f.xy) <- c("ID","X","Y","Pres")
    d1b<-merge(f.xy,d1b,by=c("ID","X","Y"),all.x=TRUE)
    d1b$fat <- fattail(d1b$bvel, 8333.3335, 0.5) 
    sppref<-rasterFromXYZ(d1b[,c(2,3,6)])
    sppref[is.na(sppref[])] <- 0
    refstack<-addLayer(refstack,sppref)
    futprevstack<-addLayer(futprevstack,futprev)
  }
  futprevmean <- overlay(futprevstack, fun='mean')
  writeRaster(futprevmean,filename=paste(gsub(paste(species,"/",sep=""),"",spec[i]),"futprev",sep=""),format="GTiff",overwrite=TRUE)
  refmean <- overlay(refstack, fun='mean')
  writeRaster(refmean,filename=paste(ref,gsub(paste(species,"/",sep=""),"",spec[i]),"refugia",sep=""),format="GTiff",overwrite=TRUE)
}

#2.
#The next step is to calculate the multi-species refugia index. The following code reads in refugia index layers 
#for individual species, weights them according to the projected future proportional change in distribution 
#(thereby downweighting increasing species), adds them to a raster stack, and then calculates a weighted mean as the  
#refugia index. A table of change proportions is also generated. Zeros are replaced with NA values outside of the  
#combined current and future distribution of a species.

emptystack<-sample
emptystack<-addLayer(emptystack,emptystack)
emptystack<-dropLayer(emptystack,c(1,2))
velstack<-emptystack
refstack <-emptystack
chgtab <- data.frame(index=1:length(spec),chg=0)
for (i in 2:length(spec)){
  s <- list.dirs(species,full.names=FALSE)
  d <- grep(pattern=substr(spec[i],1,12),s,value=TRUE)
  p <- list.files(d,pattern="presprev.tif$")
  f <- list.files(d,pattern="futprev.tif$")
  current<-raster(paste(d,"/",p,sep=""))
  future<-raster(paste(d,"/",f,sep=""))
  curfut <- stack(current)
  curfut <- addLayer(curfut, future)
  combo <- calc(curfut,fun=max,na.rm=TRUE) 
  combo[combo==0] <- NA
  cp <- cellStats(future,sum)/cellStats(current,sum)
  chgtab$chg[i] <- cp
  chgtab$spec[i] <- d
  cp <- ifelse(cp<0.5,0.5,cp) 
  sppref <- raster(paste(ref,"/",gsub("_presprev","_refugia",p),sep=""))
  sppref <- sppref*combo
  sppref[sppref==0] <- NA
  refstack<-addLayer(refstack,sppref/cp)
}
refstackmean<-calc(refstack,fun=mean,na.rm=TRUE) 
writeRaster(refstackmean,filename=paste(ref,"_multispeciesrefugia",sep=""), format="GTiff", ,overwrite=TRUE)

#3a.
#The following code can be used to calculate local topographic metrics (landform proportions) from an land facet layer 
#(here at 100-m resolution). First, set the R temp directory to a location large enough to hold very large rasters, by 
#creating a file named ".Renviron" containing the specified temp directory:
#write("TMP = '<your-desired-tempdir>'", file=file.path(Sys.getenv('R_USER'), '.Renviron')).

#Next, create a landform layer based on the first digit of the land facet code 
#(see https://adaptwest.databasin.org/pages/adaptwest-landfacets) and calculate the proportion of each of 9 landform 
#types within each cell of a 10-km ID grid.

lf <- floor(tfacet/1000)
id <- raster::crop(id,lf)
lf <- raster::extend(lf,id)
lf <- raster::crop(lf,id)
id10 <- disaggregate(id,fact=10)
id100 <- disaggregate(id10,fact=10)
lfc <- raster::crop(lf,id100)
idcrop <- raster::crop(id100,lfc)
idpts <- rasterToPoints(id)
idpts <- as.data.frame(idpts)
names(idpts)[3] <- "ID"
lfx <- crosstab(lfc,idcrop)
lfxdf <- as.data.frame(lfx)
names(lfxdf)[1] <- "landform"
names(lfxdf)[2] <- "ID"
names(lfxdf)[3] <- "freq"
lfxw <- reshape(lfxdf, timevar = "landform", idvar = c("ID"), direction = "wide")
lfxw2 <- merge(idpts,lfxw)
lf1 <- rasterize(as.matrix(lfxw2[,2:3]), elev, lfxw2[,4])
lf2 <- rasterize(as.matrix(lfxw2[,2:3]), elev, lfxw2[,5])
lf3 <- rasterize(as.matrix(lfxw2[,2:3]), elev, lfxw2[,6])
lf4 <- rasterize(as.matrix(lfxw2[,2:3]), elev, lfxw2[,7])
lf5 <- rasterize(as.matrix(lfxw2[,2:3]), elev, lfxw2[,8])
lf6 <- rasterize(as.matrix(lfxw2[,2:3]), elev, lfxw2[,9])
lf7 <- rasterize(as.matrix(lfxw2[,2:3]), elev, lfxw2[,10])
lf8 <- rasterize(as.matrix(lfxw2[,2:3]), elev, lfxw2[,11])
lf9 <- rasterize(as.matrix(lfxw2[,2:3]), elev, lfxw2[,12])

#3b.
#The following code can be used to calculate landscape topographic metrics from an elevation layer 
#(here at 10-km resolution). 
tpi<-terrain(elev,opt='tpi',neighbors=8)
rough<-terrain(elev,opt='roughness',neighbors=8)

#3c.
#This code can be used to calculate regional topographic position indices from an elevation layer 
#(here at 10-km resolution). 
elevmean9 <- focal(elev, w=matrix(1,nrow=9,ncol=9),na.rm=TRUE,fun=mean) 
elevmean9 <- mask(elevmean9,elev)
tpi9 <- elev/(elevmean9+0.5)
elevmean25 <- focal(elev, w=matrix(1,nrow=9,ncol=9),na.rm=TRUE,fun=mean) 
elevmean25 <- mask(elevmean25,elev)
tpi25 <- elev/(elevmean25+0.5)
elevmean81 <- focal(elev, w=matrix(1,nrow=9,ncol=9),na.rm=TRUE,fun=mean) 
elevmean81 <- mask(elevmean81,elev)
tpi81 <- elev/(elevmean81+0.5)

#3d.
#And here is code to calculate the north-south corridor index from a 10-km DEM
elevdd <- projectRaster(elev,crs=geo.prj)
elev5xdd <- aggregate(elevdd, fact=5, fun=mean)
elev10xdd <- aggregate(elevdd, fact=10, fun=mean)
elev5x <- aggregate(elev,fact=5,fun=mean)
elev10x <- aggregate(elev,fact=10,fun=mean)
flow5xdd <- terrain(elev5xdd, opt='flowdir', neighbors=8)
flow10xdd <- terrain(elev10xdd, opt='flowdir', neighbors=8)
flow5xxdd <- log2(flow5xdd)+1
flow10xxdd <- log2(flow10xdd)+1
ns5xxdd <- -1*abs(flow5xxdd-4) +4
ns10xxdd <- -1*abs(flow10xxdd-4) +4
slope5xdd <- terrain(elev5xdd, opt='slope', neighbors=8)
slope5x <- terrain(elev5x, opt='slope',neighbors=8)
slope10xdd <- terrain(elev10xdd, opt='slope', neighbors=8)
slope10x <- terrain(elev10x, opt='slope',neighbors=8)
slope <- terrain(elev, opt='slope',neighbors=8)
aspectdd <- terrain(elevdd, opt='aspect',neighbors=8)
aspect<-projectRaster(aspectdd,elev)

ns5xx <- projectRaster(ns5xxdd,elev)
slope5xx <- raster::resample(slope5x,elev)
nscorr <- ns5xx*slope5xx
nscorr2 <- nscorr*ns5xx


#4a. 
#Combine continental/topographic and future climate layers in a raster stack and crops them to match the refugia
#layer extent.

xy <- rasterToPoints(dem)
m <- as.matrix(xy[,1:2])
lat <- rasterize(m, dem, field=m[,2])
topo <- raster::stack(elev,lat,coastd,facetd,tpi9,tpi25,tpi81,nscorr2,tpi,rough,lf1,lf2,lf3,lf4)
names(topo) <- c("elev","lat","coastd","facetd","tpi9","tpi25","tpi81","nscorr2","tpi","rough","lf1","lf2","lf3","lf4")
topo1 <- crop(topo,refugia)
ecor <- crop(ecor,refugia)
ecor <- mask(ecor,refugia)
combo <- raster::stack(refugia,ecor,topo1)
names(combo)[[1]] <- "refugia"
names(combo)[[2]] <- "ecor"
climfut <- crop(climfut,combo[[1]])
comboclim <- stack(combo, climfut)

#4b.
#Extract environmental data to a regular sample of points covering approximately 10% of the study region.
set.seed(999)
ecosamp <- sampleRegular(ecor, size=44000, cells=TRUE, xy=TRUE, na.rm=TRUE)
ecosamp <- na.omit(ecosamp)
ecosamp <- as.data.frame(ecosamp)
combosamp <- extract(comboclim, ecosamp[,1:2])
combosamp <- as.data.frame(combosamp)
combosamp <- na.omit(combosamp)
combosamp2 <- cbind(combosamp[,1],combosamp[,3:140])
names(combosamp2)[1] <- "refugia"
comboscale <- scale(combosamp2[,3:139])
comboscale2 <- cbind(combosamp2[,1:2],comboscale)

#4c.
#The following code can be used to build and analyze quantile regression models (Cade and Noon 2003) for nested and separate models
#containing continental ("cont"), topographic ("top") at regional ("reg"), landscape ("land"), and local "loc" scales,
#and climatic ("clim") variables. First the nested models:
base <- rq(refugia ~ 1, data=comboscale2, c(0.50,0.75,0.90,0.99,0.999))
cont <- rq(refugia ~ lat + coastd , data=comboscale2, c(0.50,0.75,0.90,0.99,0.999))
contopreg <- rq(refugia ~ lat + coastd + nscorr2 + tpi81 + tpi25 + tpi9, data=comboscale2, c(0.50,0.75,0.90,0.99,0.999))
contopregland <- rq(refugia ~ lat + coastd + nscorr2 + tpi81  + tpi25 + tpi9+ tpi  + rough + elev , data=comboscale2, c(0.50,0.75,0.90,0.99,0.999))
contop <- rq(refugia ~ lat + coastd + nscorr2 + tpi81 + tpi25 + tpi9 + tpi + rough +  elev + lf1 + lf2 + lf3 + lf4 + facetd, data=comboscale2, c(0.50,0.75,0.90,0.99,0.999))
contopclim <- rq(refugia ~ lat + coastd + nscorr2 + tpi81  + tpi25 + tpi9 + tpi  + rough +  elev+ lf1 + lf2 + lf3 + lf4 + facetd + cmi_sum + bio_01 + bio_02 + bio_10 + bio_11 + bio_12 + bio_18, data=comboscale2, c(0.50,0.75,0.90,0.99,0.999))
#And the remaining separate models:
topreg <- rq(refugia ~ nscorr2 + tpi81 + tpi25 + tpi9, data=comboscale2, c(0.50,0.75,0.90,0.99,0.999))
topland <- rq(refugia ~ tpi + rough +  elev, data=comboscale2, c(0.50,0.75,0.90,0.99,0.999))
toploc <- rq(refugia ~ lf1 + lf2 + lf3 + lf4 + facetd, data=comboscale2, c(0.50,0.75,0.90,0.99,0.999))
clim <- rq(refugia ~ cmi_sum + bio_01 + bio_02 + bio_10 + bio_11 + bio_12 + bio_18, data=comboscale2, c(0.50,0.75,0.90,0.99,0.999))
#Pseudo-R2 values:
contrq <- (1-(cont$rho/base$rho))
contopregrq <- (1-(contopreg$rho/cont$rho))
contopreglandrq <- (1-(contopregland$rho/contopreg$rho))
contoprq <- (1-(contop$rho/contopregland$rho))
contopclimrq <- (1-(contopclim$rho/contop$rho))
topregrq <- 1-(topreg$rho/base$rho)
toplandrq <- 1-(topland$rho/base$rho)
toplocrq <- 1-(toploc$rho/base$rho)
climrq <- 1-(clim$rho/base$rho)

#5a.
#The final portions of code can be used to conduct a limiting factor analysis, sensu Greenberg et al.(2015)
#but using quantile regression analysis (Cade and Noon 2003). First, reselect the relevant climate layers for analysis:
combo <- raster::stack(refugia,ecor,topo1)
climfut <- crop(climfut,combo[[1]])
comboclim <- stack(combo, climfut)
clim7 <- dropLayer(climfut,c(2:9,13:17,19:31,33:125))
clim7 <- crop(clim7, refugia)
clim7 <- mask(clim7, refugia)
names(clim7) <- c("TempAnn","TempWarm","TempCold","PrecipAnn","PrecipWarm","CMI")
combosamp3 <- combosamp2[,c(1,16,25:27,33,47)]
names(combosamp3) <- c("refugia","TempAnn","TempWarm","TempCold","PrecipAnn","PrecipWarm","CMI")

#5b.
#Next, generate scatterplots for each set of univariate relationships, develop and predict univariate regression models
#for various quantile values, and plot regression lines for each quantile. These plots can be used to visualize 
#variable relationships for each quantile value.
emptystack<-climfut[[1]]
emptystack<-addLayer(emptystack,emptystack)
emptystack<-dropLayer(emptystack,c(1,2))
lf999 <- lf99 <- lf90 <- lf50 <- lf75 <- emptystack
pdf("birdclimplots_6var.pdf")
par(mfrow = c(3, 3))
par(cex = 0.6)
par(mar = c(2, 0, 2, 0), oma = c(4, 4, 0, 0.5))
par(tcl = -0.25)
par(mgp = c(3, 1, 0))
for (i in 2:7) {
  plot(combosamp3[,i], combosamp3$refugia, main=names(combosamp3)[i])
  f <- as.formula(paste("refugia ~ ",names(combosamp3)[i],sep=""))
  lf.q99 <- rq(f, data=combosamp3, 0.99)
  lf.ref99 <- raster::predict(clim7, lf.q99)
  lf99 <- addLayer(lf99, lf.ref99)
  abline(reg=lf.q99, col="blue")
  lf.q999 <- rq(f, data=combosamp3, 0.999)
  lf.ref999 <- raster::predict(clim7, lf.q999)
  lf999 <- addLayer(lf999, lf.ref999)
  abline(reg=lf.q999, col="purple")
  lf.q75 <- rq(f, data=combosamp3, 0.75)
  lf.ref75 <- raster::predict(clim7, lf.q75)
  lf75 <- addLayer(lf75, lf.ref75)
  abline(reg=lf.q75, col="orange")
  lf.q90 <- rq(f, data=combosamp3, 0.90)
  lf.ref90 <- raster::predict(clim7, lf.q90)
  lf90 <- addLayer(lf90, lf.ref90)
  abline(reg=lf.q90, col="green")
  lf.q50 <- rq(f, data=combosamp3, 0.50)
  lf.ref50 <- raster::predict(clim7, lf.q50)
  lf50 <- addLayer(lf50, lf.ref50)
  abline(reg=lf.q50, col="red")		
}
dev.off()

#5c. 
#Finally, calculate limiting factors for each pixel based on bootstrap resampling of quantile regression models,
#and generate raster stacks containing limiting factors (variables yielding minimum predictions) for each quantile value
blf999 <- blf99 <- blf90 <- blf75 <- blf50 <- emptystack
for (j in 1:100) {
  boot <- sample_frac(combosamp3,size=1,replace=TRUE)
  lf99 <- emptystack
  lf75 <- emptystack
  lf90 <- emptystack
  lf50 <- emptystack
  lf999 <- emptystack
  for (i in 2:7) {
    f <- as.formula(paste("refugia ~ ",names(boot)[i],sep=""))
    minMax = range(boot[,i])
    xVals = seq(minMax[1], minMax[2], len = 100) 
    nd <- data.frame(lin=xVals,quad=xVals*xVals)
    names(nd) <- c(names(boot)[i])
    lf.q999 <- rq(f, data=boot, 0.999)
    lf.ref999 <- raster::predict(clim7, lf.q999)
    lf999 <- addLayer(lf999, lf.ref999)
    lf.q <- rq(f, data=boot, 0.99)
    lf.ref <- raster::predict(clim7, lf.q)
    lf99 <- addLayer(lf99, lf.ref)
    lf.q75 <- rq(f, data=boot, 0.75)
    lf.ref75 <- raster::predict(clim7, lf.q75)
    lf75 <- addLayer(lf75, lf.ref75)
    lf.q90 <- rq(f, data=boot, 0.90)
    lf.ref90 <- raster::predict(clim7, lf.q90)
    lf90 <- addLayer(lf90, lf.ref90)
    lf.q50 <- rq(f, data=boot, 0.50)
    lf.ref50 <- raster::predict(clim7, lf.q50)
    lf50 <- addLayer(lf50, lf.ref50)
  }
  names(lf50) <- names(lf75) <- names(lf90) <- names(lf99) <- names(lf999) <- names(boot[,2:7])
  lf50.b <- raster::which.min(lf50)
  lf90.b <- raster::which.min(lf90)
  lf99.b <- raster::which.min(lf99)
  lf75.b <- raster::which.min(lf75)
  lf999.b <- raster::which.min(lf999)
  blf50 <- addLayer(blf50, lf50.b)
  blf90 <- addLayer(blf90, lf90.b)	
  blf99 <- addLayer(blf99, lf99.b)	
  blf75 <- addLayer(blf75, lf75.b)
  blf999 <- addLayer(blf999, lf999.b)
}

#5d. 
#Use the Mode and Mode_count functions, defined below (from Lucas Fortini on stackoverflow.com), to calcuate the most
#frequent limiting factor variable (mode), and the frequency of that variable (mode_count).

Mode <- function(x) { 
  ux <- unique(x)
  ux=ux[!is.na(ux)]
  ux[which.max(tabulate(match(x, ux)))]
}

Mode_count <- function(x) {
  ux <- unique(x)    
  ux=ux[!is.na(ux)]
  mode=ux[which.max(tabulate(match(x, ux)))]
  sum(x==mode, na.rm=T)
}

blf50.mode <- ratify(calc(blf50, fun=Mode))
rat <- levels(blf50.mode)[[1]]
rat$names <- c('Temp Mean Annual', 'Temp Warm Month', 'Temp Cold Month', 'Precip Annual', 'Precip Warm Month', 'Moisture Index')
rat$codes <- c('BIO1','BIO10', 'BIO11','BIO12','BIO18','CMI')
levels(blf50.mode) <- rat
blf90.mode <- ratify(calc(blf90, fun=Mode))
levels(blf90.mode) <- rat
blf99.mode <- ratify(calc(blf99, fun=Mode))
levels(blf99.mode) <- rat
blf75.mode <- ratify(calc(blf75, fun=Mode))
levels(blf75.mode) <- rat
mode <- stack(blf50.mode, blf75.mode, blf90.mode, blf99.mode)
mode <- setZ(mode, c(50,75,90,99))
names(mode) <- c("percentile 50", "percentile 75", "percentile 90", "percentile 99")

blf50.freq <- calc(blf50, fun=Mode_count)
blf90.freq <- calc(blf90, fun=Mode_count)
blf99.freq <- calc(blf99, fun=Mode_count)
blf75.freq <- calc(blf75, fun=Mode_count)
freq <- stack(blf50.freq, blf75.freq, blf90.freq, blf99.freq)
freq <- mask(freq,refugia)
freq <- setZ(freq, c(50,75,90,99))
names(freq) <- c("percentile 50", "percentile 75", "percentile 90", "percentile 99")

myTheme=rasterTheme(region=brewer.pal('Set1', n=6))
freqTheme=rasterTheme(region=brewer.pal('YlGnBu',n=5))

pdf("elf_6var_mode_freq.pdf", width=11, height=5)
plot1 <- levelplot(mode, att='names',scales=list(draw=FALSE), layout=c(5,1), par.settings=myTheme)
par(mar=c(0,0,0,10))
plot2 <- levelplot(freq, scales=list(draw=FALSE), layout=c(5,1), par.settings=freqTheme)
grid.arrange(plot1,plot2, nrow=2)
dev.off()

pdf("elf_6var_mode.pdf", width=8, height=8)
plot1 <- levelplot(mode, att='names',scales=list(draw=FALSE), layout=c(2,3), par.settings=myTheme, between=list(0,0))
plot1
dev.off()



# References
# Cade, B. S., and B. R. Noon. 2003. A gentle introduction to quantile regression for ecologists. Frontiers in Ecology and the Environment 1:412-420.
# Carroll, C., J. J. Lawler, D. R. Roberts, and A. Hamann. 2015. Biotic and climatic velocity identify contrasting areas of vulnerability to climate change. PLoS ONE 10:e0140486.
# Clark, J. S., C. Fastie, G. Hurtt, S. T. Jackson, C. Johnson, G. A. King, M. Lewis, J. Lynch, S. Pacala, C. Prentice, E. W. Schupp, I. I. I. T. Webb, and P. Wyckoff. 1998. Reid's Paradox of Rapid Plant MigrationDispersal theory and interpretation of paleoecological records. BioScience 48:13-24.
# Greenberg, J. A., M. J. Santos, S. Z. Dobrowski, V. C. Vanderbilt, and S. L. Ustin. 2015. Quantifying environmental limiting factors on tree cover using geospatial data. PLoS ONE 10:e0114648.
# Hamann, A., D. Roberts, Q. Barber, C. Carroll, and S. Nielsen. 2015. Velocity of climate change algorithms for guiding conservation and management. Global Change Biology 21:997-1004.
# Stralberg, D., C. Carroll, J. H. Pedlar, C. B. Wilsey, D. W. McKenney, and S. E. Nielsen. 2018. Macrorefugia for North American trees and songbirds: Climatic limiting factors and multi-scale topographic influences. Global Ecology and Biogeography. https://doi.org/10.1111/geb.12731 