#load in libraries
library(splines)
library(foreach)
library(nlme)
library(mgcv)
library(Metrics)
library(MASS)
library(sm) 
library(mgcv.helper)
library(tibble)
library(dplyr)
library(fmsb)
library(randomForest)
library(KRLS)
library(matrixStats)
#remove variables from prior run
rm(list = ls())
#define model concurvity threshold:
concurvity_threshold = 0.4

##=========================pick best model based on BurnArea_GAM_v0*.R
#read in data

#0-1100m:
Climate_Fire_Data = read.table("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter1.csv",sep=",")
#define variables:
WYs = Climate_Fire_Data[1:37,1]
BAs = Climate_Fire_Data[1:37,2]
BAs_Z1 = (BAs*0.000247105) /(10^6);
SWEIs_Z1 = Climate_Fire_Data[1:37,3]
WinPRCP_Z1 = Climate_Fire_Data[1:37,4]
WinTMP_Z1 = Climate_Fire_Data[1:37,5]
WinterPDSI_Z1 = Climate_Fire_Data[1:37,6]
Spring_PRCP_Z1 = Climate_Fire_Data[1:37,7] 
Spring_TMP_Z1 = Climate_Fire_Data[1:37,8] 
Spring_PDSI_Z1 = Climate_Fire_Data[1:37,9] 
Spring_SWEI_Z1 = Climate_Fire_Data[1:37,10] 
Spring_VPD_Z1 = Climate_Fire_Data[1:37,11] 
Winter_VPD_Z1 = Climate_Fire_Data[1:37,12] 
MODIS_BA_Z1 = Climate_Fire_Data[1:37,13] 
MODIS_BA_Z1 = (MODIS_BA_Z1*0.000247105) /(10^6);
Spring_ET_Z1 = Climate_Fire_Data[1:37,14] 
Spring_PET_Z1 = Climate_Fire_Data[1:37,15] 
Spring_PETminusET_Z1 = Climate_Fire_Data[1:37,16] 
Winter_ET_Z1 = Climate_Fire_Data[1:37,17] 
Winter_PET_Z1 = Climate_Fire_Data[1:37,18] 
Winter_PETminusET_Z1 = Climate_Fire_Data[1:37,19] 
Winter_DroughtArea_Z1 = Climate_Fire_Data[1:37,33] 
Spring_DroughtArea_Z1 = Climate_Fire_Data[1:37,34] 

Z1_ID = rep(1,37)

#relate MODIS and MTBS BA linearly to allow MODIS to gap fill MTBS 2020 BA:
BA_MTBS <- BAs_Z1[18:36]
BA_MODIS <- MODIS_BA_Z1[18:36]
BA_DF = data.frame(BA_MTBS=BA_MTBS,BA_MODIS=BA_MODIS)
BA_lm <- lm(data=BA_DF,BA_MTBS~BA_MODIS)
data_2020 = data.frame(BA_MODIS = MODIS_BA_Z1[37])
BA_2020_MTBS=predict(BA_lm,newdata=data_2020)
BAs_Z1[37] = BA_2020_MTBS
if (BAs_Z1[37] < 0){
  BAs_Z1[37] = 0
}

#1100-2200:
Climate_Fire_Data = read.table("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter2.csv",sep=",")
#define variables:
#define variables:
WYs = Climate_Fire_Data[1:37,1]
BAs = Climate_Fire_Data[1:37,2]
BAs_Z2 = (BAs*0.000247105) /(10^6);
SWEIs_Z2 = Climate_Fire_Data[1:37,3]
WinPRCP_Z2 = Climate_Fire_Data[1:37,4]
WinTMP_Z2 = Climate_Fire_Data[1:37,5]
WinterPDSI_Z2 = Climate_Fire_Data[1:37,6]
Spring_PRCP_Z2 = Climate_Fire_Data[1:37,7] 
Spring_TMP_Z2 = Climate_Fire_Data[1:37,8] 
Spring_PDSI_Z2 = Climate_Fire_Data[1:37,9] 
Spring_SWEI_Z2 = Climate_Fire_Data[1:37,10] 
Spring_VPD_Z2 = Climate_Fire_Data[1:37,11] 
Winter_VPD_Z2 = Climate_Fire_Data[1:37,12] 
MODIS_BA_Z2 = Climate_Fire_Data[1:37,13] 
Spring_ET_Z2 = Climate_Fire_Data[1:37,14] 
Spring_PET_Z2 = Climate_Fire_Data[1:37,15] 
Spring_PETminusET_Z2 = Climate_Fire_Data[1:37,16] 
Winter_ET_Z2 = Climate_Fire_Data[1:37,17] 
Winter_PET_Z2 = Climate_Fire_Data[1:37,18] 
Winter_PETminusET_Z2 = Climate_Fire_Data[1:37,19] 
Winter_DroughtArea_Z2 = Climate_Fire_Data[1:37,33] 
Spring_DroughtArea_Z2 = Climate_Fire_Data[1:37,34]  
Z2_ID = rep(2,37)

#relate MODIS and MTBS BA linearly to allow MODIS to gap fill MTBS 2020 BA:
BA_MTBS <- BAs_Z2[18:36]
BA_MODIS <- MODIS_BA_Z2[18:36]
BA_DF = data.frame(BA_MTBS=BA_MTBS,BA_MODIS=BA_MODIS)
BA_lm <- lm(data=BA_DF,BA_MTBS~BA_MODIS)
data_2020 = data.frame(BA_MODIS = MODIS_BA_Z2[37])
BA_2020_MTBS=predict(BA_lm,newdata=data_2020)
BAs_Z2[37] = BA_2020_MTBS
if (BAs_Z2[37] < 0){
  BAs_Z2[37] = 0
}

#2200-3300:
Climate_Fire_Data = read.table("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter3.csv",sep=",")
#define variables:
#define variables:
WYs = Climate_Fire_Data[1:37,1]
BAs = Climate_Fire_Data[1:37,2]
BAs_Z3 = (BAs*0.000247105) /(10^6);
SWEIs_Z3 = Climate_Fire_Data[1:37,3]
WinPRCP_Z3 = Climate_Fire_Data[1:37,4]
WinTMP_Z3 = Climate_Fire_Data[1:37,5]
WinterPDSI_Z3 = Climate_Fire_Data[1:37,6]
Spring_PRCP_Z3 = Climate_Fire_Data[1:37,7] 
Spring_TMP_Z3 = Climate_Fire_Data[1:37,8] 
Spring_PDSI_Z3 = Climate_Fire_Data[1:37,9] 
Spring_SWEI_Z3 = Climate_Fire_Data[1:37,10] 
Spring_VPD_Z3 = Climate_Fire_Data[1:37,11] 
Winter_VPD_Z3 = Climate_Fire_Data[1:37,12] 
MODIS_BA_Z3 = Climate_Fire_Data[1:37,13] 
Spring_ET_Z3 = Climate_Fire_Data[1:37,14] 
Spring_PET_Z3 = Climate_Fire_Data[1:37,15] 
Spring_PETminusET_Z3 = Climate_Fire_Data[1:37,16] 
Winter_ET_Z3 = Climate_Fire_Data[1:37,17] 
Winter_PET_Z3 = Climate_Fire_Data[1:37,18] 
Winter_PETminusET_Z3 = Climate_Fire_Data[1:37,19] 
Winter_DroughtArea_Z3 = Climate_Fire_Data[1:37,33] 
Spring_DroughtArea_Z3 = Climate_Fire_Data[1:37,34]  
Z3_ID = rep(3,37)

#relate MODIS and MTBS BA linearly to allow MODIS to gap fill MTBS 2020 BA:
BA_MTBS <- BAs_Z3[18:36]
BA_MODIS <- MODIS_BA_Z3[18:36]
BA_DF = data.frame(BA_MTBS=BA_MTBS,BA_MODIS=BA_MODIS)
BA_lm <- lm(data=BA_DF,BA_MTBS~BA_MODIS)
data_2020 = data.frame(BA_MODIS = MODIS_BA_Z3[37])
BA_2020_MTBS=predict(BA_lm,newdata=data_2020)
BAs_Z3[37] = BA_2020_MTBS
if (BAs_Z3[37] < 0){
  BAs_Z3[37] = 0
}

idx<-c(1:37)
WYs = rep(WYs[idx],3)
BAs = cbind(t(BAs_Z1[idx]),t(BAs_Z2[idx]),t(BAs_Z3[idx]))
SWEIs = cbind(t(SWEIs_Z1[idx]),t(SWEIs_Z2[idx]),t(SWEIs_Z3[idx]))
WinPRCP = cbind(t(WinPRCP_Z1[idx]),t(WinPRCP_Z2[idx]),t(WinPRCP_Z3[idx]))
WinTMP = cbind(t(WinTMP_Z1[idx]),t(WinTMP_Z2[idx]),t(WinTMP_Z3[idx]))
WinterPDSI = cbind(t(WinterPDSI_Z1[idx]),t(WinterPDSI_Z2[idx]),t(WinterPDSI_Z3[idx]))
Spring_PRCP = cbind(t(Spring_PRCP_Z1[idx]),t(Spring_PRCP_Z2[idx]),t(Spring_PRCP_Z3[idx]))
Spring_TMP = cbind(t(Spring_TMP_Z1[idx]),t(Spring_TMP_Z2[idx]),t(Spring_TMP_Z3[idx]))
Spring_PDSI = cbind(t(Spring_PDSI_Z1[idx]),t(Spring_PDSI_Z2[idx]),t(Spring_PDSI_Z3[idx]))
Spring_SWEI = cbind(t(Spring_SWEI_Z1[idx]),t(Spring_SWEI_Z2[idx]),t(Spring_SWEI_Z3[idx]))
Spring_VPD = cbind(t(Spring_VPD_Z1[idx]),t(Spring_VPD_Z2[idx]),t(Spring_VPD_Z3[idx]))
Winter_VPD = cbind(t(Winter_VPD_Z1[idx]),t(Winter_VPD_Z2[idx]),t(Winter_VPD_Z3[idx]))
Spring_ET = cbind(t(Spring_ET_Z1[idx]),t(Spring_ET_Z2[idx]),t(Spring_ET_Z3[idx]))
Spring_PET = cbind(t(Spring_PET_Z1[idx]),t(Spring_PET_Z2[idx]),t(Spring_PET_Z3[idx]))
Spring_PETminusET = cbind(t(Spring_PETminusET_Z1[idx]),t(Spring_PETminusET_Z2[idx]),t(Spring_PETminusET_Z3[idx]))
Winter_ET = cbind(t(Winter_ET_Z1[idx]),t(Winter_ET_Z2[idx]),t(Winter_ET_Z3[idx]))
Winter_PET = cbind(t(Winter_PET_Z1[idx]),t(Winter_PET_Z2[idx]),t(Winter_PET_Z3[idx]))
Winter_PETminusET = cbind(t(Winter_PETminusET_Z1[idx]),t(Winter_PETminusET_Z2[idx]),t(Winter_PETminusET_Z3[idx]))
Winter_DroughtArea = cbind(t(Winter_DroughtArea_Z1[idx]),t(Winter_DroughtArea_Z2[idx]),t(Winter_DroughtArea_Z3[idx]))
Spring_DroughtArea = cbind(t(Spring_DroughtArea_Z1[idx]),t(Spring_DroughtArea_Z2[idx]),t(Spring_DroughtArea_Z3[idx]))
Winter_Spring_Temp = (Spring_TMP+WinTMP)/2
Winter_Spring_Precip = (Spring_PRCP+WinPRCP)/2
Winter_Spring_VPD = (Winter_VPD+Spring_VPD)/2
Winter_Spring_PET = (Winter_PET+Spring_PET)/2
Winter_Spring_ET = (Winter_ET+Spring_ET)/2

Zs = cbind(t(Z1_ID[idx]),t(Z2_ID[idx]),t(Z3_ID[idx]))

Climate_Fire_DF = data.frame(Spring_PDSI=t(Spring_PDSI),WinPRCP = t(WinPRCP), WinTMP = t(WinTMP),Spring_PRCP = t(Spring_PRCP),Spring_TMP=t(Spring_TMP),Spring_VPD=t(Spring_VPD),Winter_VPD=t(Winter_VPD),Spring_ET=t(Spring_ET),Spring_PET=t(Spring_PET),Winter_ET=t(Winter_ET),Winter_PET=t(Winter_PET),Spring_DroughtArea=t(Spring_DroughtArea),Winter_Spring_Temp=t(Winter_Spring_Temp),Winter_Spring_Precip=t(Winter_Spring_Precip),Winter_Spring_VPD=t(Winter_Spring_VPD),Winter_Spring_PET=t(Winter_Spring_PET),Winter_Spring_ET=t(Winter_Spring_ET))

#loop through all combinations of covariates:
Y = as.vector(BAs)
nvars <- 17

#for  combos of 2 variables:
n_predictors <- 2
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

store_AIC_2vars = c()
store_R_2vars = c()
store_i_2vars=c()
for (i in 1:ncombos){
  col1 <- combo_IDs[1,i]
  col2 <- combo_IDs[2,i]
  
  predictor1 <-Climate_Fire_DF[,col1]
  predictor2 <-Climate_Fire_DF[,col2]
  #PCA combo model:
  Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2)
  df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
  pcs = df.pca$x
  pc1 = pcs[,1]
  pc2 = pcs[,2]
  current_DF_pcs <- data.frame(Y=Y,predictor1=pc1,predictor2=pc2,predictor3=as.vector(Zs))
  mod <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+predictor3,family="gaussian")
  
  X = concurvity(mod,full=FALSE)
  X=X$estimate
  idx <- which(X==1)
  X[idx] = 0
  Max_Concuvrity <- max(X)
  if (Max_Concuvrity <= concurvity_threshold){
    AIC = mod$aic
    R = cor(predict(mod),Y)
    
    store_R_2vars = c(store_R_2vars,R)
    store_AIC_2vars = c(store_AIC_2vars,AIC)
    store_i_2vars = c(store_i_2vars,i)
  }
}

idx_bestmod_2vars <- which(store_AIC_2vars==min(store_AIC_2vars))
store_R_2vars[idx_bestmod_2vars]
idx_bestmod_2vars = store_i_2vars[idx_bestmod_2vars]

#for  combos of 3 variables:
n_predictors <- 3
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

store_AIC_3vars = c()
store_R_3vars = c()
store_i_3vars = c()
for (i in 1:ncombos){
  col1 <- combo_IDs[1,i]
  col2 <- combo_IDs[2,i]
  col3 <- combo_IDs[3,i]
  
  predictor1 <-Climate_Fire_DF[,col1]
  predictor2 <-Climate_Fire_DF[,col2]
  predictor3 <-Climate_Fire_DF[,col3]
  #PCA combo model:
  Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3)
  df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
  pcs = df.pca$x
  pc1 = pcs[,1]
  pc2 = pcs[,2]
  pc3 = pcs[,3]
  current_DF_pcs <- data.frame(Y=Y,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=as.vector(Zs))
  mod <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+predictor4,family="gaussian")
  
  X = concurvity(mod,full=FALSE)
  X=X$estimate
  idx <- which(X==1)
  X[idx] = 0
  Max_Concuvrity <- max(X)
  if (Max_Concuvrity <= concurvity_threshold){
    AIC = mod$aic
    R = cor(predict(mod),Y)
    
    store_R_3vars = c(store_R_3vars,R)
    store_AIC_3vars = c(store_AIC_3vars,AIC)
    store_i_3vars = c(store_i_3vars,i)
  }
}

idx_bestmod_3vars <- which(store_AIC_3vars==min(store_AIC_3vars))
store_R_3vars[idx_bestmod_3vars]
idx_bestmod_3vars = store_i_3vars[idx_bestmod_3vars]


#for  combos of 4 variables:
n_predictors <- 4
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

store_AIC_4vars = c()
store_R_4vars = c()
store_i_4vars=c()
for (i in 1:ncombos){
  col1 <- combo_IDs[1,i]
  col2 <- combo_IDs[2,i]
  col3 <- combo_IDs[3,i]
  col4 <- combo_IDs[4,i]
  
  predictor1 <-Climate_Fire_DF[,col1]
  predictor2 <-Climate_Fire_DF[,col2]
  predictor3 <-Climate_Fire_DF[,col3]
  predictor4 <-Climate_Fire_DF[,col4]
  #PCA combo model:
  Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4)
  df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
  pcs = df.pca$x
  pc1 = pcs[,1]
  pc2 = pcs[,2]
  pc3 = pcs[,3]
  pc4 = pcs[,4]
  current_DF_pcs <- data.frame(Y=Y,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=as.vector(Zs))
  mod <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4)+predictor5,family="gaussian")
  
  X = concurvity(mod,full=FALSE)
  X=X$estimate
  idx <- which(X==1)
  X[idx] = 0
  Max_Concuvrity <- max(X)
  if (Max_Concuvrity <= concurvity_threshold){
    AIC = mod$aic
    R = cor(predict(mod),Y)
    
    store_R_4vars = c(store_R_4vars,R)
    store_AIC_4vars = c(store_AIC_4vars,AIC)
    store_i_4vars = c(store_i_4vars,i)
  }
}

idx_bestmod_4vars <- which(store_AIC_4vars==min(store_AIC_4vars))
store_R_4vars[idx_bestmod_4vars]
idx_bestmod_4vars = store_i_3vars[idx_bestmod_4vars]


#for  combos of 5 variables:
n_predictors <- 5
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

store_AIC_5vars = c()
store_R_5vars = c()
store_i_5vars=c()
for (i in 1:ncombos){
  col1 <- combo_IDs[1,i]
  col2 <- combo_IDs[2,i]
  col3 <- combo_IDs[3,i]
  col4 <- combo_IDs[4,i]
  col5 <- combo_IDs[5,i]
  
  predictor1 <-Climate_Fire_DF[,col1]
  predictor2 <-Climate_Fire_DF[,col2]
  predictor3 <-Climate_Fire_DF[,col3]
  predictor4 <-Climate_Fire_DF[,col4]
  predictor5 <-Climate_Fire_DF[,col5]
  #PCA combo model:
  Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5)
  df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
  pcs = df.pca$x
  pc1 = pcs[,1]
  pc2 = pcs[,2]
  pc3 = pcs[,3]
  pc4 = pcs[,4]
  pc5 = pcs[,5]
  current_DF_pcs <- data.frame(Y=Y,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=as.vector(Zs))
  mod <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4)+s(predictor5)+predictor6,family="gaussian")
  
  X = concurvity(mod,full=FALSE)
  X=X$estimate
  idx <- which(X==1)
  X[idx] = 0
  Max_Concuvrity <- max(X)
  if (Max_Concuvrity <= concurvity_threshold){
    AIC = mod$aic
    R = cor(predict(mod),Y)
    
    store_R_5vars = c(store_R_5vars,R)
    store_AIC_5vars = c(store_AIC_5vars,AIC)
    store_i_5vars = c(store_i_5vars,i)
  }
}

idx_bestmod_5vars <- which(store_AIC_5vars==min(store_AIC_5vars))
store_R_5vars[idx_bestmod_5vars]
idx_bestmod_5vars = store_i_5vars[idx_bestmod_5vars]

#for  combos of 6 variables:
n_predictors <- 6
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

store_AIC_6vars = c()
store_R_6vars = c()
store_i_6vars=c()
for (i in 1:ncombos){
  col1 <- combo_IDs[1,i]
  col2 <- combo_IDs[2,i]
  col3 <- combo_IDs[3,i]
  col4 <- combo_IDs[4,i]
  col5 <- combo_IDs[5,i]
  col6 <- combo_IDs[6,i]
  
  predictor1 <-Climate_Fire_DF[,col1]
  predictor2 <-Climate_Fire_DF[,col2]
  predictor3 <-Climate_Fire_DF[,col3]
  predictor4 <-Climate_Fire_DF[,col4]
  predictor5 <-Climate_Fire_DF[,col5]
  predictor6 <-Climate_Fire_DF[,col6]

  #PCA combo model:
  Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6)
  df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
  pcs = df.pca$x
  pc1 = pcs[,1]
  pc2 = pcs[,2]
  pc3 = pcs[,3]
  pc4 = pcs[,4]
  pc5 = pcs[,5]
  pc6 = pcs[,6]
  current_DF_pcs <- data.frame(Y=Y,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=as.vector(Zs))
  mod <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4)+s(predictor5)+s(predictor6)+predictor7,family="gaussian")
  
  X = concurvity(mod,full=FALSE)
  X=X$estimate
  idx <- which(X==1)
  X[idx] = 0
  Max_Concuvrity <- max(X)
  if (Max_Concuvrity <= concurvity_threshold){
    AIC = mod$aic
    R = cor(predict(mod),Y)
    
    store_R_6vars = c(store_R_6vars,R)
    store_AIC_6vars = c(store_AIC_6vars,AIC)
    store_i_6vars = c(store_i_6vars,i)
  }
}

idx_bestmod_6vars <- which(store_AIC_6vars==min(store_AIC_6vars))
store_R_6vars[idx_bestmod_6vars]
idx_bestmod_6vars = store_i_6vars[idx_bestmod_6vars]

#for  combos of 7 variables:
n_predictors <- 7
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

store_AIC_7vars = c()
store_R_7vars = c()
store_i_7vars = c()
for (i in 1:ncombos){
  col1 <- combo_IDs[1,i]
  col2 <- combo_IDs[2,i]
  col3 <- combo_IDs[3,i]
  col4 <- combo_IDs[4,i]
  col5 <- combo_IDs[5,i]
  col6 <- combo_IDs[6,i]
  col7 <- combo_IDs[7,i]
  
  predictor1 <-Climate_Fire_DF[,col1]
  predictor2 <-Climate_Fire_DF[,col2]
  predictor3 <-Climate_Fire_DF[,col3]
  predictor4 <-Climate_Fire_DF[,col4]
  predictor5 <-Climate_Fire_DF[,col5]
  predictor6 <-Climate_Fire_DF[,col6]
  predictor7 <-Climate_Fire_DF[,col7]
  #PCA combo model:
  Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7)
  df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
  pcs = df.pca$x
  pc1 = pcs[,1]
  pc2 = pcs[,2]
  pc3 = pcs[,3]
  pc4 = pcs[,4]
  pc5 = pcs[,5]
  pc6 = pcs[,6]
  pc7 = pcs[,7]
  current_DF_pcs <- data.frame(Y=Y,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=as.vector(Zs))
  mod <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4)+s(predictor5)+s(predictor6)+s(predictor7)+predictor8,family="gaussian")

  X = concurvity(mod,full=FALSE)
  X=X$estimate
  idx <- which(X==1)
  X[idx] = 0
  Max_Concuvrity <- max(X)
  if (Max_Concuvrity <= concurvity_threshold){
    AIC = mod$aic
    R = cor(predict(mod),Y)
    
    store_R_7vars = c(store_R_7vars,R)
    store_AIC_7vars = c(store_AIC_7vars,AIC)
    store_i_7vars = c(store_i_7vars,i)
  } 
}

idx_bestmod_7vars <- which(store_AIC_7vars==min(store_AIC_7vars))
store_R_7vars[idx_bestmod_7vars]
idx_bestmod_7vars = store_i_7vars[idx_bestmod_7vars]


#for  combos of 8 variables:
n_predictors <- 8
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

store_AIC_8vars = c()
store_R_8vars = c()
store_i_8vars = c()
for (i in 1:ncombos){
  col1 <- combo_IDs[1,i]
  col2 <- combo_IDs[2,i]
  col3 <- combo_IDs[3,i]
  col4 <- combo_IDs[4,i]
  col5 <- combo_IDs[5,i]
  col6 <- combo_IDs[6,i]
  col7 <- combo_IDs[7,i]
  col8 <- combo_IDs[8,i]
  
  predictor1 <-Climate_Fire_DF[,col1]
  predictor2 <-Climate_Fire_DF[,col2]
  predictor3 <-Climate_Fire_DF[,col3]
  predictor4 <-Climate_Fire_DF[,col4]
  predictor5 <-Climate_Fire_DF[,col5]
  predictor6 <-Climate_Fire_DF[,col6]
  predictor7 <-Climate_Fire_DF[,col7]
  predictor8 <-Climate_Fire_DF[,col8]
  #PCA combo model:
  Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8,predictor9=as.vector(Zs))
  df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
  pcs = df.pca$x
  pc1 = pcs[,1]
  pc2 = pcs[,2]
  pc3 = pcs[,3]
  pc4 = pcs[,4]
  pc5 = pcs[,5]
  pc6 = pcs[,6]
  pc7 = pcs[,7]
  pc8 = pcs[,8]
  current_DF_pcs <- data.frame(Y=Y,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8,predictor9=as.vector(Zs))
  mod <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4)+s(predictor5)+s(predictor6)+s(predictor7)+s(predictor8)+predictor9,family="gaussian")
  
  X = concurvity(mod,full=FALSE)
  X=X$estimate
  idx <- which(X==1)
  X[idx] = 0
  Max_Concuvrity <- max(X)
  if (Max_Concuvrity <= concurvity_threshold){
    AIC = mod$aic
    R = cor(predict(mod),Y)
    
    store_R_8vars = c(store_R_8vars,R)
    store_AIC_8vars = c(store_AIC_8vars,AIC)
    store_i_8vars = c(store_i_8vars,i)
  }
}

idx_bestmod_8vars <- which(store_AIC_8vars==min(store_AIC_8vars))
store_R_8vars[idx_bestmod_8vars]
idx_bestmod_8vars = store_i_8vars[idx_bestmod_8vars]


#for  combos of 9 variables:
n_predictors <- 9
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

store_AIC_9vars = c()
store_R_9vars = c()
store_i_9vars = c()
for (i in 1:ncombos){
  col1 <- combo_IDs[1,i]
  col2 <- combo_IDs[2,i]
  col3 <- combo_IDs[3,i]
  col4 <- combo_IDs[4,i]
  col5 <- combo_IDs[5,i]
  col6 <- combo_IDs[6,i]
  col7 <- combo_IDs[7,i]
  col8 <- combo_IDs[8,i]
  col9 <- combo_IDs[9,i]
  
  predictor1 <-Climate_Fire_DF[,col1]
  predictor2 <-Climate_Fire_DF[,col2]
  predictor3 <-Climate_Fire_DF[,col3]
  predictor4 <-Climate_Fire_DF[,col4]
  predictor5 <-Climate_Fire_DF[,col5]
  predictor6 <-Climate_Fire_DF[,col6]
  predictor7 <-Climate_Fire_DF[,col7]
  predictor8 <-Climate_Fire_DF[,col8]
  predictor9 <-Climate_Fire_DF[,col9]
  
  #PCA combo model:
  Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8,predictor9=predictor9)
  df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
  pcs = df.pca$x
  pc1 = pcs[,1]
  pc2 = pcs[,2]
  pc3 = pcs[,3]
  pc4 = pcs[,4]
  pc5 = pcs[,5]
  pc6 = pcs[,6]
  pc7 = pcs[,7]
  pc8 = pcs[,8]
  pc9 = pcs[,9]
  current_DF_pcs <- data.frame(Y=Y,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8,predictor9=pc9,predictor10=as.vector(Zs))
  mod <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4)+s(predictor5)+s(predictor6)+s(predictor7)+s(predictor8)+s(predictor9)+predictor10,family="gaussian")
  
  X = concurvity(mod,full=FALSE)
  X=X$estimate
  idx <- which(X==1)
  X[idx] = 0
  Max_Concuvrity <- max(X)
  if (Max_Concuvrity <= concurvity_threshold){
    AIC = mod$aic
    R = cor(predict(mod),Y)
    
    store_R_9vars = c(store_R_9vars,R)
    store_AIC_9vars = c(store_AIC_9vars,AIC)
    store_i_9vars = c(store_i_9vars,i)
  }
}

idx_bestmod_9vars <- which(store_AIC_9vars==min(store_AIC_9vars))
store_R_9vars[idx_bestmod_9vars]
idx_bestmod_9vars = store_i_9vars[idx_bestmod_9vars]


#for  combos of 10 variables:
n_predictors <- 10
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

store_AIC_10vars = c()
store_R_10vars = c()
store_i_10vars = c()
for (i in 1:ncombos){
  col1 <- combo_IDs[1,i]
  col2 <- combo_IDs[2,i]
  col3 <- combo_IDs[3,i]
  col4 <- combo_IDs[4,i]
  col5 <- combo_IDs[5,i]
  col6 <- combo_IDs[6,i]
  col7 <- combo_IDs[7,i]
  col8 <- combo_IDs[8,i]
  col9 <- combo_IDs[9,i]
  col10 <- combo_IDs[10,i]
  
  predictor1 <-Climate_Fire_DF[,col1]
  predictor2 <-Climate_Fire_DF[,col2]
  predictor3 <-Climate_Fire_DF[,col3]
  predictor4 <-Climate_Fire_DF[,col4]
  predictor5 <-Climate_Fire_DF[,col5]
  predictor6 <-Climate_Fire_DF[,col6]
  predictor7 <-Climate_Fire_DF[,col7]
  predictor8 <-Climate_Fire_DF[,col8]
  predictor9 <-Climate_Fire_DF[,col9]
  predictor10 <-Climate_Fire_DF[,col10]
  
  #PCA combo model:
  Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8,predictor9=predictor9,predictor10=predictor10)
  df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
  pcs = df.pca$x
  pc1 = pcs[,1]
  pc2 = pcs[,2]
  pc3 = pcs[,3]
  pc4 = pcs[,4]
  pc5 = pcs[,5]
  pc6 = pcs[,6]
  pc7 = pcs[,7]
  pc8 = pcs[,8]
  pc9 = pcs[,9]
  pc10 = pcs[,10]
  current_DF_pcs <- data.frame(Y=Y,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8,predictor9=pc9,predictor10=pc10,predictor11=as.vector(Zs))
  mod <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4)+s(predictor5)+s(predictor6)+s(predictor7)+s(predictor8)+s(predictor9)+s(predictor10)+predictor11,family="gaussian")
  
  X = concurvity(mod,full=FALSE)
  X=X$estimate
  idx <- which(X==1)
  X[idx] = 0
  Max_Concuvrity <- max(X)
  if (Max_Concuvrity <= concurvity_threshold){
    AIC = mod$aic
    R = cor(predict(mod),Y)
    
    store_R_10vars = c(store_R_10vars,R)
    store_AIC_10vars = c(store_AIC_10vars,AIC)
    store_i_10vars = c(store_i_10vars,i)
  }
}

idx_bestmod_10vars <- which(store_AIC_10vars==min(store_AIC_9vars))
store_R_10vars[idx_bestmod_10vars]
idx_bestmod_10vars = store_i_10vars[idx_bestmod_10vars]

#for  combos of 11 variables:
n_predictors <- 11
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

store_AIC_11vars = c()
store_R_11vars = c()
store_i_11vars = c()
for (i in 1:ncombos){
  col1 <- combo_IDs[1,i]
  col2 <- combo_IDs[2,i]
  col3 <- combo_IDs[3,i]
  col4 <- combo_IDs[4,i]
  col5 <- combo_IDs[5,i]
  col6 <- combo_IDs[6,i]
  col7 <- combo_IDs[7,i]
  col8 <- combo_IDs[8,i]
  col9 <- combo_IDs[9,i]
  col10 <- combo_IDs[10,i]
  col11 <- combo_IDs[11,i]
  
  predictor1 <-Climate_Fire_DF[,col1]
  predictor2 <-Climate_Fire_DF[,col2]
  predictor3 <-Climate_Fire_DF[,col3]
  predictor4 <-Climate_Fire_DF[,col4]
  predictor5 <-Climate_Fire_DF[,col5]
  predictor6 <-Climate_Fire_DF[,col6]
  predictor7 <-Climate_Fire_DF[,col7]
  predictor8 <-Climate_Fire_DF[,col8]
  predictor9 <-Climate_Fire_DF[,col9]
  predictor10 <-Climate_Fire_DF[,col10]
  predictor11 <-Climate_Fire_DF[,col11]
  
  #PCA combo model:
  Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8,predictor9=predictor9,predictor10=predictor10,predictor11=predictor11)
  df.pca <- prcomp(Covariates_DF, center = TRUE,scale. = TRUE)
  pcs = df.pca$x
  pc1 = pcs[,1]
  pc2 = pcs[,2]
  pc3 = pcs[,3]
  pc4 = pcs[,4]
  pc5 = pcs[,5]
  pc6 = pcs[,6]
  pc7 = pcs[,7]
  pc8 = pcs[,8]
  pc9 = pcs[,9]
  pc10 = pcs[,10]
  pc11 = pcs[,11]
  current_DF_pcs <- data.frame(Y=Y,predictor1=pc1,predictor2=pc2,predictor3=pc3,predictor4=pc4,predictor5=pc5,predictor6=pc6,predictor7=pc7,predictor8=pc8,predictor9=pc9,predictor10=pc10,predictor11=pc11,predictor12=as.vector(Zs))
  mod <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4)+s(predictor5)+s(predictor6)+s(predictor7)+s(predictor8)+s(predictor9)+s(predictor10)+s(predictor11)+predictor12,family="gaussian")
  
  X = concurvity(mod,full=FALSE)
  X=X$estimate
  idx <- which(X==1)
  X[idx] = 0
  Max_Concuvrity <- max(X)
  if (Max_Concuvrity <= concurvity_threshold){
    AIC = mod$aic
    R = cor(predict(mod),Y)
    
    store_R_11vars = c(store_R_11vars,R)
    store_AIC_11vars = c(store_AIC_11vars,AIC)
    store_i_11vars = c(store_i_11vars,i)
  }
}

idx_bestmod_11vars <- which(store_AIC_11vars==min(store_AIC_9vars))
store_R_11vars[idx_bestmod_11vars]
idx_bestmod_11vars = store_i_11vars[idx_bestmod_11vars]

#12 vars - model has more coefficients than data for drop 1

#13 vars - model has more coefficients than data for fit



#define the best models by AIC:

vars=c("PDSI","WinPRCP","WinTMP","Spring-prcp","Spring TMP","Spring VPD","Winter VPD","Spring ET","Spring PET","Winter ET","Winter PET","Spring area","winter+spring temp","winter+spring precip","winter+spring VPD","winter+spring PET")

# 7 variable model is best:

#store best models:
AIC_7vars=store_AIC_7vars
i_7vars=store_i_7vars
store_i=c()
for (b in 1:100){
  idx_bestmod_7vars <- which(AIC_7vars==min(AIC_7vars))
  i_bestmod_7vars = i_7vars[idx_bestmod_7vars]
  store_i =cbind(store_i,i_bestmod_7vars)
  
  AIC_7vars=AIC_7vars[-idx_bestmod_7vars]
  i_7vars=i_7vars[-idx_bestmod_7vars]
}

best_is=store_i
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_PDSImod_PCA_bestmods.csv")
write.table(best_is,file=outfilename,sep=",",row.names = FALSE)


##record results across all variable options:

#2vars:
AIC_2vars=store_AIC_2vars
i_2vars=store_i_2vars
store_i=c()
for (b in 1:100){
  idx_bestmod_2vars <- which(AIC_2vars==min(AIC_2vars))
  i_bestmod_2vars = i_2vars[idx_bestmod_2vars]
  store_i =cbind(store_i,i_bestmod_2vars)
  
  AIC_2vars=AIC_2vars[-idx_bestmod_2vars]
  i_2vars=i_2vars[-idx_bestmod_2vars]
}

best_is=store_i
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_PDSImod_PCA_bestmods_2vars.csv")
write.table(best_is,file=outfilename,sep=",",row.names = FALSE)


#3vars:
AIC_3vars=store_AIC_3vars
i_3vars=store_i_3vars
store_i=c()
for (b in 1:100){
  idx_bestmod_3vars <- which(AIC_3vars==min(AIC_3vars))
  i_bestmod_3vars = i_3vars[idx_bestmod_3vars]
  store_i =cbind(store_i,i_bestmod_3vars)
  
  AIC_3vars=AIC_3vars[-idx_bestmod_3vars]
  i_3vars=i_3vars[-idx_bestmod_3vars]
}

best_is=store_i
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_PDSImod_PCA_bestmods_3vars.csv")
write.table(best_is,file=outfilename,sep=",",row.names = FALSE)

#4vars:
AIC_4vars=store_AIC_4vars
i_4vars=store_i_4vars
store_i=c()
for (b in 1:100){
  idx_bestmod_4vars <- which(AIC_4vars==min(AIC_4vars))
  i_bestmod_4vars = i_4vars[idx_bestmod_4vars]
  store_i =cbind(store_i,i_bestmod_4vars)
  
  AIC_4vars=AIC_4vars[-idx_bestmod_4vars]
  i_4vars=i_4vars[-idx_bestmod_4vars]
}

best_is=store_i
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_PDSImod_PCA_bestmods_4vars.csv")
write.table(best_is,file=outfilename,sep=",",row.names = FALSE)

#5vars:
AIC_5vars=store_AIC_5vars
i_5vars=store_i_5vars
store_i=c()
for (b in 1:100){
  idx_bestmod_5vars <- which(AIC_5vars==min(AIC_5vars))
  i_bestmod_5vars = i_5vars[idx_bestmod_5vars]
  store_i =cbind(store_i,i_bestmod_5vars)
  
  AIC_5vars=AIC_5vars[-idx_bestmod_5vars]
  i_5vars=i_5vars[-idx_bestmod_5vars]
}

best_is=store_i
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_PDSImod_PCA_bestmods_5vars.csv")
write.table(best_is,file=outfilename,sep=",",row.names = FALSE)

#6vars:
AIC_6vars=store_AIC_6vars
i_6vars=store_i_6vars
store_i=c()
for (b in 1:100){
  idx_bestmod_6vars <- which(AIC_6vars==min(AIC_6vars))
  i_bestmod_6vars = i_6vars[idx_bestmod_6vars]
  store_i =cbind(store_i,i_bestmod_6vars)
  
  AIC_6vars=AIC_6vars[-idx_bestmod_6vars]
  i_6vars=i_6vars[-idx_bestmod_6vars]
}

best_is=store_i
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_PDSImod_PCA_bestmods_6vars.csv")
write.table(best_is,file=outfilename,sep=",",row.names = FALSE)

#7vars:
AIC_7vars=store_AIC_7vars
i_7vars=store_i_7vars
store_i=c()
for (b in 1:100){
  idx_bestmod_7vars <- which(AIC_7vars==min(AIC_7vars))
  i_bestmod_7vars = i_7vars[idx_bestmod_7vars]
  store_i =cbind(store_i,i_bestmod_7vars)
  
  AIC_7vars=AIC_7vars[-idx_bestmod_7vars]
  i_7vars=i_7vars[-idx_bestmod_7vars]
}

best_is=store_i
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_PDSImod_PCA_bestmods_7vars.csv")
write.table(best_is,file=outfilename,sep=",",row.names = FALSE)

#8vars:
AIC_8vars=store_AIC_8vars
i_8vars=store_i_8vars
store_i=c()
for (b in 1:100){
  idx_bestmod_8vars <- which(AIC_8vars==min(AIC_8vars))
  i_bestmod_8vars = i_8vars[idx_bestmod_8vars]
  store_i =cbind(store_i,i_bestmod_8vars)
  
  AIC_8vars=AIC_8vars[-idx_bestmod_8vars]
  i_8vars=i_8vars[-idx_bestmod_8vars]
}

best_is=store_i
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_PDSImod_PCA_bestmods_8vars.csv")
write.table(best_is,file=outfilename,sep=",",row.names = FALSE)

#9vars:
AIC_9vars=store_AIC_9vars
i_9vars=store_i_9vars
store_i=c()
for (b in 1:100){
  idx_bestmod_9vars <- which(AIC_9vars==min(AIC_9vars))
  i_bestmod_9vars = i_9vars[idx_bestmod_9vars]
  store_i =cbind(store_i,i_bestmod_9vars)
  
  AIC_9vars=AIC_9vars[-idx_bestmod_9vars]
  i_9vars=i_9vars[-idx_bestmod_9vars]
}

best_is=store_i
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_PDSImod_PCA_bestmods_9vars.csv")
write.table(best_is,file=outfilename,sep=",",row.names = FALSE)


#10vars:
AIC_10vars=store_AIC_10vars
i_10vars=store_i_10vars
store_i=c()
for (b in 1:100){
  idx_bestmod_10vars <- which(AIC_10vars==min(AIC_10vars))
  i_bestmod_10vars = i_10vars[idx_bestmod_10vars]
  store_i =cbind(store_i,i_bestmod_10vars)
  
  AIC_10vars=AIC_10vars[-idx_bestmod_10vars]
  i_10vars=i_10vars[-idx_bestmod_10vars]
}

best_is=store_i
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_PDSImod_PCA_bestmods_10vars.csv")
write.table(best_is,file=outfilename,sep=",",row.names = FALSE)

#11vars:
AIC_11vars=store_AIC_11vars
i_11vars=store_i_11vars
store_i=c()
for (b in 1:100){
  idx_bestmod_11vars <- which(AIC_11vars==min(AIC_11vars))
  i_bestmod_11vars = i_11vars[idx_bestmod_11vars]
  store_i =cbind(store_i,i_bestmod_11vars)
  
  AIC_11vars=AIC_11vars[-idx_bestmod_11vars]
  i_11vars=i_11vars[-idx_bestmod_11vars]
}

best_is=store_i
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_PDSImod_PCA_bestmods_11vars.csv")
write.table(best_is,file=outfilename,sep=",",row.names = FALSE)


