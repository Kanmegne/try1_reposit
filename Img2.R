

########################################################################################################
#################               DIRECTORIES             ################################################

setwd("D:/PhD/Data/imageProcessing/buyo_gueyo/result1")
dir1 = "F:/shapefile/Sentinel2/S2B_MSIL1C_20180222T104029_N0206_R008_T29NRG_20180223T184036.SAFE/GRANULE/L1C_T29NRG_A005039_20180222T105258/IMG_DATA"
dir2 = "D:/PhD/Data/imageProcessing/Inputs/Landcover"

####################################################################################################
##############################        LIBRARIES               ##########################################
library(raster) 
library(rgdal)
library(rgeos)
library(caTools)
library(sf)
########################################################################################################
#####################           SHAPEFILES             ################################################

civ            = shapefile("D:/PhD/Data/imageProcessing/Inputs/Shapefiles/CIV_adm0.shp")
roi            = shapefile("D:/PhD/Data/imageProcessing/buyo_gueyo/Gueyo.shp")
polygon1       = shapefile("D:/PhD/Data/imageProcessing/buyo_gueyo/polygon_Gueyo.shp")
polygon2       = shapefile("D:/PhD/Data/imageProcessing/buyo_gueyo/polygon_buyo.shp")
polygon3       = shapefile("D:/PhD/Data/imageProcessing/buyo_gueyo/others_shape.shp")
polyg_test     = polygon1 + polygon3
trainData      = polygon1 + polygon2 + polygon3

writeOGR(trainData, ".", "trainData", driver = "ESRI Shapefile", overwrite_layer = T)
trainData      = shapefile("D:/PhD/Data/imageProcessing/buyo_gueyo/result/trainData.shp")                                           #Reload shapefile
trainData      = spTransform(trainData, crs(proj4string(b)))


#########################################################################################################
########################          MASK DATA         #####################################################

landcover = list.files(dir2, pattern = ".tif$", full.names = T)

bare_class        = projectRaster((crop(raster(unlist(landcover[1])),roi)),crs = crs(proj4string(trainData)))
forest_class      = projectRaster((crop(raster(unlist(landcover[5])),roi)),crs = crs(proj4string(trainData)))
urban_class       = projectRaster((crop(raster(unlist(landcover[11])),roi)),crs = crs(proj4string(trainData)))
water_perm_class  = projectRaster((crop(raster(unlist(landcover[12])),roi)),crs = crs(proj4string(trainData)))
water_seas_class  = projectRaster((crop(raster(unlist(landcover[13])),roi)),crs = crs(proj4string(trainData)))

#Creating a mask cote d'ivoire
mask_bare      = rasterToPolygons(bare_class >=5, dissolve = T)
mask_forest    = rasterToPolygons((forest_class==2)+(forest_class==4), dissolve = T)
mask_urban     = rasterToPolygons(urban_class >=10, dissolve = T)
mask_Gueyo     = mask_bare + mask_forest +mask_urban 

mask_vect  = gBuffer(mask_Gueyo, byid = T, width = 0)        # to prevent unmatching geometry
cadre  = as(extent(mask_Gueyo), "SpatialPolygons")
cadre  = gBuffer(cadre, byid = T, width = 0)        # to prevent unmatching geometry

crs(cadre) = crs(mask_Gueyo)
cadre  = mask(cadre, mask_Gueyo)
mask       = erase(cadre, mask_vect)

#######################################################################################################
#########################        IMAGEs FOR GUEYO        ##############################################
# A- Sentinel-1
s1 = projectRaster(brick("D:/PhD/Data/imageProcessing/buyo_gueyo/sentinel1/subset_S1A_IW_GRDH_1SDV_20180211T184233_20180211T184258_020563_0232F8_F069_Orb_Cal_Spk_TC.tif"),crs = crs(proj4string(b)))
s1 = crop(s1, mask)

VH = resample(mask(s1[[1]], mask_gueyo, inverse=T), spectral, method="bilinear", progress='text') 
VV = resample(mask(s1[[2]], mask_gueyo, inverse=T), spectral, method="bilinear", progress='text')  

# Texture parameters generation with GLCM
library(glcm)
text1 = glcm(VH,n_grey = 32, window = c(5,5), shift = c(1,1), statistics = c("correlation","entropy", "variance","contrast"))
text2 = glcm(VV,n_grey = 32, window = c(5,5), shift = c(1,1), statistics = c("correlation","entropy", "variance","contrast"))

spatial = stack(VH,VV,text1$glcm_correlation, text1$glcm_variance, text1$glcm_contrast, text1$glcm_entropy
                , text2$glcm_correlation, text2$glcm_variance, text2$glcm_contrast, text2$glcm_entropy)
names(spatial) =c("VH","VV","VH_correlation", "VH_entropy", "VH_variance", "VH_contrast", "VV_correlation", "VV_entropy", "VV_variance","VV_contrast")
writeRaster(spatial, "spatial_gueyo.tif", overwrite=T)

# B- Sentinel-2 
setwd(dir1)
B      = crop(raster('T29NRG_20180222T104029_B02.jp2'), mask_gueyo)
G      = crop(raster('T29NRG_20180222T104029_B03.jp2'), mask_gueyo)
R      = crop(raster('T29NRG_20180222T104029_B04.jp2'), mask_gueyo)
RE     = resample((crop(raster('T29NRG_20180222T104029_B05.jp2'),mask_gueyo)),B, method="bilinear", progress="text")
NIR    = crop(raster('T29NRG_20180222T104029_B08.jp2'), mask_gueyo)

# vegetation indices
NDVI       = mask(((NIR-R)/(NIR+R)), mask_gueyo, inverse=T)
GLI        = mask(((2*G-R-B)/(2*G+R+B)), mask_gueyo, inverse =T)
EVI        = mask((2.5*((NIR-R)/(NIR+6*R*B+1))), mask_gueyo, inverse=T)
SAVI       = mask((((NIR/R)/(NIR+R+0.5))*(1+0.5)), mask_gueyo, inverse=T)
MSAVI      = mask(((2*NIR+1-sqrt((2*NIR+1)^2-8*(NIR-R)))/2), mask_gueyo, inverse=T)
TCARI      = mask((3*(RE-R)-0.2*(RE-G)*(RE/R)), mask_gueyo, inverse=T)
VARI       = mask(((G-R)/(G+R-B)), mask_gueyo, inverse=T)

spectral        = stack(B,G,R,RE,NIR, NDVI,GLI,EVI,SAVI,MSAVI,TCARI,VARI)
names(spectral) = c("B","G","R","RE","NIR", "NDVI","GLI","EVI","SAVI","MSAVI","TCARI","VARI")
writeRaster(spectral, "Img_gueyo.tif", overwrite=T, progress='text')

# C- Sentinel 1 + Sentinel 2
Sentinels        = stack(spectral, spatial)
names(Sentinels)  = c("B","G","R","RE","NIR", "NDVI","GLI","EVI","SAVI","MSAVI","TCARI","VARI","VH","VV","VH_correlation", "VH_entropy", "VH_variance", "VH_entropy", "VV_correlation", "VV_entropy", "VV_variance","VV_entropy")
writeRaster(Sentinels, "combined_s1_s2.tif")


##########################################################################################################
#######################         IMAGE TO DATA           ############################################
# Reload image data
# sentinel 1
spatial = brick("spatial_gueyo.tif")
names(spatial) =c("VH","VV","VH_correlation", "VH_entropy", "VH_variance", "VH_contrast", "VV_correlation", "VV_entropy", "VV_variance","VV_contrast")

#sentinel 2
spectral         = projectRaster(brick("spect_Img_gueyo.tif"), crs=crs(spatial))
names(spectral)  = c("B","G","R","RE","NIR", "NDVI","GLI","EVI","SAVI","MSAVI","TCARI","VARI")

# sentinel 1 and 2 combined
Sentinels = brick("combined_s1_s2.tif")
names(Sentinels)  = c("B","G","R","RE","NIR", "NDVI","GLI","EVI","SAVI","MSAVI","TCARI","VARI","VH","VV","VH_correlation", "VH_entropy", "VH_variance", "VH_contrast", "VV_correlation", "VV_entropy", "VV_variance","VV_contrast")

#plot image
plotRGB(spectral, r=8,g=4,b=3, stretch="lin")                          

#1.spectral features extraction
#for (i in 1:length(unique(trainData[[responseCol]]))){
df_s1         = data.frame(matrix(vector(), nrow=0, ncol = length(names(spatial))+1))
df_s2         = data.frame(matrix(vector(), nrow=0, ncol = length(names(spectral))+1))
df_s3         = data.frame(matrix(vector(), nrow=0, ncol = length(names(Sentinels))+1))
responseCol   = "AFS"
for (i in 1:length(unique(trainData[[responseCol]]))){
category    = unique(trainData[[responseCol]])[i]
categoryMap = trainData[trainData[[responseCol]]== category,]
dataset1     = sapply(extract(spatial, categoryMap), function(x){cbind(x, class = rep(category, nrow(x)))})
dataset2     = sapply(extract(spectral, categoryMap), function(x){cbind(x, class = rep(category, nrow(x)))})
dataset3     = sapply(extract(Sentinels, categoryMap), function(x){cbind(x, class = rep(category, nrow(x)))})
df1          = do.call("rbind",dataset1)
df2          = do.call("rbind",dataset2)
df3          = do.call("rbind",dataset3)
df_s1       = rbind(df_s1, df1)
df_s2       = rbind(df_s2, df2)
df_s3       = rbind(df_s3, df3)
}
write.csv(df_s1, "Gueyo_data_s1.csv")
write.csv(df_s2, "gueyo_data_s2.csv")
write.csv(df_s3, "gueyo_data_s3.csv")

###############################################################################################################
#                                         Matrix Data Gueyo
# load data image data matrix
s1_data         = na.omit(read.csv("gueyo_data_s1.csv", stringsAsFactors = T)[,-1])
s2_data         = na.omit(read.csv("gueyo_data_s2.csv", stringsAsFactors = T)[,-1])
s3_data         = na.omit(read.csv("gueyo_data_s3.csv", stringsAsFactors = T)[,-1])

names(s1_data)  = c("VH","VV","VH_correlation", "VH_entropy", "VH_variance", "VH_contrast", "VV_correlation", "VV_entropy", "VV_variance","VV_contrast","class")
names(s2_data)  = c("B","G","R","RE","NIR", "NDVI","GLI","EVI","SAVI","MSAVI","TCARI","VARI","class")
names(s3_data)  = c("B","G","R","RE","NIR", "NDVI","GLI","EVI","SAVI","MSAVI","TCARI","VARI","VH","VV","VH_correlation", "VH_entropy", "VH_variance", "VH_contrast", "VV_correlation", "VV_entropy", "VV_variance","VV_contrast","class")

#subset and processing over 1000 samples to reduce the computing load
nsamples = 1000
set.seed(515)
s1_data = s1_data[sample(1:nrow(s1_data),nsamples),]
s2_data = s2_data[sample(1:nrow(s2_data),nsamples),]
s3_data = s3_data[sample(1:nrow(s3_data),nsamples),]


#II. sample partitionning (change the code according to the dataset you are working with "s1_..s2_..)
# split training data 2/3
library(caTools)
split         = sample.split(s3_data, SplitRatio = 2/3)
train_set     = subset(s3_data, split==TRUE)
test_set      = subset(s3_data, split==FALSE)
train_labels  = train_set[,length(s3_data)]
test_labels   = test_set[,length(s3_data)]

############################################################################################################
#####################              MODEL TRAINING              ##############################################
library(randomForest)
library(caret)
# parameter tuning
set.seed(515)
#sentinel 1
ctrl           = trainControl(method = "repeatedcv", number = 10, repeats = 10)
grid_rf_model  = expand.grid(.mtry=c(2,3,6,10))
tuned         = train(class~., data = train_set, method="rf", metric= "Kappa", trControl = ctrl, tuneGrid= grid_rf_model)
tuned

# Model Evaluation 
rf_pred        = predict(tuned, test_set)
result1        = confusionMatrix(rf_pred, test_labels)
result1

# Image classification
beginCluster()
map1 = clusterR(spatial, raster::predict, args = list(model=tuned))
endCluster()
writeRaster(map1, "map1.tif", overwrite=T)


# sentinel 2
grid_rf_model2  = expand.grid(.mtry=c(2,4,8,12))
tuned2          = train(class~., data = train_set, method="rf", metric= "Kappa", trControl = ctrl, tuneGrid= grid_rf_model2)
tuned2

# Model Evaluation 
rf_pred2      = predict(tuned2, test_set)
result2       = confusionMatrix(rf_pred2, test_labels)
result2

# Image classification
beginCluster()
map2 = clusterR(spectral, raster::predict, args = list(model=tuned2))
endCluster()
writeRaster(map2, "map2.tif", overwrite=T)


# sentinel 3
grid_rf_model3  = expand.grid(.mtry=c(3,5,10,22))
tuned3          = train(class~., data = train_set, method="rf", metric= "Kappa", trControl = ctrl, tuneGrid= grid_rf_model3)
tuned3

# Model Evaluation 
rf_pred3      = predict(tuned3, test_set)
result3       = confusionMatrix(rf_pred3, test_labels)
result3

# Image classification
beginCluster()
map3 = clusterR(Sentinels, raster::predict, args = list(model=tuned3))
endCluster()
writeRaster(map3, "map3.tif", overwrite=T)

###############################################
#### Plot classified map
plot(map2.smth, col=c("green","blue","red", "yellow"), legend=F, axes=F, main="classification s2")
legend("right", legend = c("cocoa","farm", "palm", "rubber"), fill=c("green","blue","red", "yellow"), border = F, bty = "n")

# smooth
map2.smth = focal(map2, w=matrix(1,3,3), fun=modal)
writeRaster(map2.smth, "smooth_gueyo.tif", overwrite=T)

# zonal statistics give the total number of pixel per class
zones = zonal(map2.smth, map2.smth, fun= 'sum', digits=3, na.rm=T)


# reprojection classified map
ras = projectRaster(map2.smth,crs = crs(B))
plot(ras, col=c("green","blue","red", "yellow"), legend=F, axes=F, main="agroforestry systems")
legend("right", legend = c("cocoa","farm", "palm", "rubber"), fill=c("green","blue","red", "yellow"), border = F, bty = "n")
writeRaster(ras, "zstat.tif", overwrite=T)

