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
Spring_ET_Z1 = Climate_Fire_Data[1:37,14] 
Spring_PET_Z1 = Climate_Fire_Data[1:37,15] 
Spring_PETminusET_Z1 = Climate_Fire_Data[1:37,16] 
Winter_ET_Z1 = Climate_Fire_Data[1:37,17] 
Winter_PET_Z1 = Climate_Fire_Data[1:37,18] 
Winter_PETminusET_Z1 = Climate_Fire_Data[1:37,19] 
Winter_DroughtArea_Z1 = Climate_Fire_Data[1:37,31] 
Spring_DroughtArea_Z1 = Climate_Fire_Data[1:37,32] 

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
Winter_DroughtArea_Z2 = Climate_Fire_Data[1:37,31] 
Spring_DroughtArea_Z2 = Climate_Fire_Data[1:37,32]  
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
Winter_DroughtArea_Z3 = Climate_Fire_Data[1:37,31] 
Spring_DroughtArea_Z3 = Climate_Fire_Data[1:37,32]  
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

Climate_Fire_DF = data.frame(Spring_SWEI=t(Spring_SWEI),WinPRCP = t(WinPRCP), WinTMP = t(WinTMP),Spring_PRCP = t(Spring_PRCP),Spring_TMP=t(Spring_TMP),Spring_VPD=t(Spring_VPD),Winter_VPD=t(Winter_VPD),Spring_ET=t(Spring_ET),Spring_PET=t(Spring_PET),Winter_ET=t(Winter_ET),Winter_PET=t(Winter_PET),Spring_DroughtArea=t(Spring_DroughtArea),Winter_Spring_Temp=t(Winter_Spring_Temp),Winter_Spring_Precip=t(Winter_Spring_Precip),Winter_Spring_VPD=t(Winter_Spring_VPD),Winter_Spring_PET=t(Winter_Spring_PET),Winter_Spring_ET=t(Winter_Spring_ET))


#loop through all combinations of covariates:


#2 variable model:
best_mods_idx = read.csv("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_bestmods_2vars.csv",sep=",",header=TRUE)
best_mods_idx=as.numeric(best_mods_idx)
nmods = length(best_mods_idx)
#for  combos of 2 variables:
nvars<- 17
n_predictors <- 2
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

SWEI_iter=0 #count which of the best mods include SWEI as predictor
good_mod_iter=0
Good_mods=c()
store_Taylor = c()
store_AIC=c()
n=length(BAs)
index = 1:n
for (j in 1:nmods){
  i <- best_mods_idx[j]
  Y = as.vector(BAs)
  
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
  Fit_AIC = mod$aic
  store_AIC = c(store_AIC,Fit_AIC)
  #GAM drop 1 year
  yest_gam=1
  standar_error=1
  Y_iter=0
  for(Y in 1984:2020){
    Y_iter = Y_iter+1
    #drop 1 year
    IDX_drop <-which(WYs == Y)
    index1=index[-IDX_drop]
    dropped_DF=Climate_Fire_DF[index1,]
    
    #PCA combo model:
    Zs_dropped=as.vector(Zs)
    Zs_dropped = Zs_dropped[index1]
    
    Y_dropped = as.vector(BAs)
    Y_dropped = Y_dropped[index1]
    
    dropped_DF_pcs <- data.frame(Y=Y_dropped,predictor1=pc1[index1],predictor2=pc2[index1],predictor3=as.vector(Zs_dropped))
    
    #model
    fit_drop_gam = gam(mod$formula,data=dropped_DF_pcs,family="gaussian")
    #now estimate at the point that was dropped
    newdata = data.frame(predictor1=pc1[IDX_drop],predictor2=pc2[IDX_drop],predictor3=as.vector(Zs[IDX_drop]))
    yest_gam[IDX_drop]=predict(fit_drop_gam,newdata=newdata)
    #get the confidence interval:
    yhat <- predict(fit_drop_gam,newdata=newdata,se.fit = TRUE)
    standar_error[Y_iter] <- sum(yhat$se.fit)
  }
  IDX_neg <- which(yest_gam<0)
  if (length(IDX_neg) > 0){
  yest_gam[IDX_neg] <- min(yest_gam[-IDX_neg])
  }
  #sum BAs across all Z bins:
  iter=0
  store_BA_obs=0
  store_BA_est=0
  for (Y in 1984:2020){
    iter <- iter+1
    current_IDX <-which(WYs == Y)
    store_BA_est[iter] = sum(yest_gam[current_IDX])
    store_BA_obs[iter] = sum(BAs[current_IDX])
  }
  #record drop 1-year taylor score:
  Drop1_Ann_R <- cor(store_BA_est,store_BA_obs)
  sigma_f = sd(store_BA_est)
  sigma_r = sd(store_BA_obs)
  Taylor_Score = (4*(1+Drop1_Ann_R)^4)/( ( (sigma_f/sigma_r)+(sigma_r/sigma_f))^2 *(1+1)^4)
  store_Taylor = c(store_Taylor,Taylor_Score)
  
}
output<-rbind(store_AIC,store_Taylor)
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_SWEImod_PCA_2vars_AICfit_Drop1Taylor.csv")
write.table(output,file=outfilename,sep=",",row.names = FALSE)


#3 variable model:
best_mods_idx = read.csv("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_bestmods_3vars.csv",sep=",",header=TRUE)
best_mods_idx=as.numeric(best_mods_idx)
nmods = length(best_mods_idx)
#for  combos of 3 variables:

n_predictors <- 3
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

SWEI_iter=0 #count which of the best mods include SWEI as predictor
good_mod_iter=0
Good_mods=c()
store_Taylor = c()
store_AIC=c()
n=length(BAs)
index = 1:n
for (j in 1:nmods){
  i <- best_mods_idx[j]
  Y = as.vector(BAs)
  
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
  Fit_AIC = mod$aic
  store_AIC = c(store_AIC,Fit_AIC)
  #GAM drop 1 year
  yest_gam=1
  standar_error=1
  Y_iter=0
  for(Y in 1984:2020){
    Y_iter = Y_iter+1
    #drop 1 year
    IDX_drop <-which(WYs == Y)
    index1=index[-IDX_drop]
    dropped_DF=Climate_Fire_DF[index1,]
    
    #PCA combo model:
    Zs_dropped=as.vector(Zs)
    Zs_dropped = Zs_dropped[index1]
    
    Y_dropped = as.vector(BAs)
    Y_dropped = Y_dropped[index1]
    
    dropped_DF_pcs <- data.frame(Y=Y_dropped,predictor1=pc1[index1],predictor2=pc2[index1],predictor3=pc3[index1],predictor4=as.vector(Zs_dropped))
    
    #model
    fit_drop_gam = gam(mod$formula,data=dropped_DF_pcs,family="gaussian")
    #now estimate at the point that was dropped
    newdata = data.frame(predictor1=pc1[IDX_drop],predictor2=pc2[IDX_drop],predictor3=pc3[IDX_drop],predictor4=as.vector(Zs[IDX_drop]))
    yest_gam[IDX_drop]=predict(fit_drop_gam,newdata=newdata)
    #get the confidence interval:
    yhat <- predict(fit_drop_gam,newdata=newdata,se.fit = TRUE)
    standar_error[Y_iter] <- sum(yhat$se.fit)
  }
  IDX_neg <- which(yest_gam<0)
  if (length(IDX_neg) > 0){
    yest_gam[IDX_neg] <- min(yest_gam[-IDX_neg])
  }  
  #sum BAs across all Z bins:
  iter=0
  store_BA_obs=0
  store_BA_est=0
  for (Y in 1984:2020){
    iter <- iter+1
    current_IDX <-which(WYs == Y)
    store_BA_est[iter] = sum(yest_gam[current_IDX])
    store_BA_obs[iter] = sum(BAs[current_IDX])
  }
  #record drop 1-year taylor score:
  Drop1_Ann_R <- cor(store_BA_est,store_BA_obs)
  sigma_f = sd(store_BA_est)
  sigma_r = sd(store_BA_obs)
  Taylor_Score = (4*(1+Drop1_Ann_R)^4)/( ( (sigma_f/sigma_r)+(sigma_r/sigma_f))^2 *(1+1)^4)
  store_Taylor = c(store_Taylor,Taylor_Score)
  
}
output<-rbind(store_AIC,store_Taylor)
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_SWEImod_PCA_3vars_AICfit_Drop1Taylor.csv")
write.table(output,file=outfilename,sep=",",row.names = FALSE)


#4 variable model:
best_mods_idx = read.csv("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_bestmods_4vars.csv",sep=",",header=TRUE)
best_mods_idx=as.numeric(best_mods_idx)
nmods = length(best_mods_idx)
#for  combos of 4 variables:

n_predictors <- 4
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

SWEI_iter=0 #count which of the best mods include SWEI as predictor
good_mod_iter=0
Good_mods=c()
store_Taylor = c()
store_AIC=c()
n=length(BAs)
index = 1:n
for (j in 1:nmods){
  i <- best_mods_idx[j]
  Y = as.vector(BAs)
  
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
  Fit_AIC = mod$aic
  store_AIC = c(store_AIC,Fit_AIC)
  #GAM drop 1 year
  yest_gam=1
  standar_error=1
  Y_iter=0
  for(Y in 1984:2020){
    Y_iter = Y_iter+1
    #drop 1 year
    IDX_drop <-which(WYs == Y)
    index1=index[-IDX_drop]
    dropped_DF=Climate_Fire_DF[index1,]
    
    #PCA combo model:
    Zs_dropped=as.vector(Zs)
    Zs_dropped = Zs_dropped[index1]
    
    Y_dropped = as.vector(BAs)
    Y_dropped = Y_dropped[index1]
    
    dropped_DF_pcs <- data.frame(Y=Y_dropped,predictor1=pc1[index1],predictor2=pc2[index1],predictor3=pc3[index1],predictor4=pc4[index1],predictor5=as.vector(Zs_dropped))
    
    #model
    fit_drop_gam = gam(mod$formula,data=dropped_DF_pcs,family="gaussian")
    #now estimate at the point that was dropped
    newdata = data.frame(predictor1=pc1[IDX_drop],predictor2=pc2[IDX_drop],predictor3=pc3[IDX_drop],predictor4=pc4[IDX_drop],predictor5=as.vector(Zs[IDX_drop]))
    yest_gam[IDX_drop]=predict(fit_drop_gam,newdata=newdata)
    #get the confidence interval:
    yhat <- predict(fit_drop_gam,newdata=newdata,se.fit = TRUE)
    standar_error[Y_iter] <- sum(yhat$se.fit)
  }
  IDX_neg <- which(yest_gam<0)
  if (length(IDX_neg) > 0){
    yest_gam[IDX_neg] <- min(yest_gam[-IDX_neg])
  }  

  #sum BAs across all Z bins:
  iter=0
  store_BA_obs=0
  store_BA_est=0
  for (Y in 1984:2020){
    iter <- iter+1
    current_IDX <-which(WYs == Y)
    store_BA_est[iter] = sum(yest_gam[current_IDX])
    store_BA_obs[iter] = sum(BAs[current_IDX])
  }
  #record drop 1-year taylor score:
  Drop1_Ann_R <- cor(store_BA_est,store_BA_obs)
  sigma_f = sd(store_BA_est)
  sigma_r = sd(store_BA_obs)
  Taylor_Score = (4*(1+Drop1_Ann_R)^4)/( ( (sigma_f/sigma_r)+(sigma_r/sigma_f))^2 *(1+1)^4)
  store_Taylor = c(store_Taylor,Taylor_Score)
  
}
output<-rbind(store_AIC,store_Taylor)
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_SWEImod_PCA_4vars_AICfit_Drop1Taylor.csv")
write.table(output,file=outfilename,sep=",",row.names = FALSE)

#5 variable model:
best_mods_idx = read.csv("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_bestmods_5vars.csv",sep=",",header=TRUE)
best_mods_idx=as.numeric(best_mods_idx)
nmods = length(best_mods_idx)
#for  combos of 5 variables:

n_predictors <- 5
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

SWEI_iter=0 #count which of the best mods include SWEI as predictor
good_mod_iter=0
Good_mods=c()
store_Taylor = c()
store_AIC=c()
n=length(BAs)
index = 1:n
for (j in 1:nmods){
  i <- best_mods_idx[j]
  Y = as.vector(BAs)
  
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
  Fit_AIC = mod$aic
  store_AIC = c(store_AIC,Fit_AIC)
  #GAM drop 1 year
  yest_gam=1
  standar_error=1
  Y_iter=0
  for(Y in 1984:2020){
    Y_iter = Y_iter+1
    #drop 1 year
    IDX_drop <-which(WYs == Y)
    index1=index[-IDX_drop]
    dropped_DF=Climate_Fire_DF[index1,]
    
    #PCA combo model:
    Zs_dropped=as.vector(Zs)
    Zs_dropped = Zs_dropped[index1]
    
    Y_dropped = as.vector(BAs)
    Y_dropped = Y_dropped[index1]
    
    dropped_DF_pcs <- data.frame(Y=Y_dropped,predictor1=pc1[index1],predictor2=pc2[index1],predictor3=pc3[index1],predictor4=pc4[index1],predictor5=pc5[index1],predictor6=as.vector(Zs_dropped))
    
    #model
    fit_drop_gam = gam(mod$formula,data=dropped_DF_pcs,family="gaussian")
    #now estimate at the point that was dropped
    newdata = data.frame(predictor1=pc1[IDX_drop],predictor2=pc2[IDX_drop],predictor3=pc3[IDX_drop],predictor4=pc4[IDX_drop],predictor5=pc5[IDX_drop],predictor6=as.vector(Zs[IDX_drop]))
    yest_gam[IDX_drop]=predict(fit_drop_gam,newdata=newdata)
    #get the confidence interval:
    yhat <- predict(fit_drop_gam,newdata=newdata,se.fit = TRUE)
    standar_error[Y_iter] <- sum(yhat$se.fit)
  }
  IDX_neg <- which(yest_gam<0)
  if (length(IDX_neg) > 0){
    yest_gam[IDX_neg] <- min(yest_gam[-IDX_neg])
  }  
  #sum BAs across all Z bins:
  iter=0
  store_BA_obs=0
  store_BA_est=0
  for (Y in 1984:2020){
    iter <- iter+1
    current_IDX <-which(WYs == Y)
    store_BA_est[iter] = sum(yest_gam[current_IDX])
    store_BA_obs[iter] = sum(BAs[current_IDX])
  }
  #record drop 1-year taylor score:
  Drop1_Ann_R <- cor(store_BA_est,store_BA_obs)
  sigma_f = sd(store_BA_est)
  sigma_r = sd(store_BA_obs)
  Taylor_Score = (4*(1+Drop1_Ann_R)^4)/( ( (sigma_f/sigma_r)+(sigma_r/sigma_f))^2 *(1+1)^4)
  store_Taylor = c(store_Taylor,Taylor_Score)
  
}
output<-rbind(store_AIC,store_Taylor)
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_SWEImod_PCA_5vars_AICfit_Drop1Taylor.csv")
write.table(output,file=outfilename,sep=",",row.names = FALSE)

#6 variable model:
best_mods_idx = read.csv("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_bestmods_6vars.csv",sep=",",header=TRUE)
best_mods_idx=as.numeric(best_mods_idx)
nmods = length(best_mods_idx)
#for  combos of 6 variables:

n_predictors <- 6
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

SWEI_iter=0 #count which of the best mods include SWEI as predictor
good_mod_iter=0
Good_mods=c()
store_Taylor = c()
store_AIC=c()
n=length(BAs)
index = 1:n
for (j in 1:nmods){
  i <- best_mods_idx[j]
  Y = as.vector(BAs)
  
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
  Fit_AIC = mod$aic
  store_AIC = c(store_AIC,Fit_AIC)
  #GAM drop 1 year
  yest_gam=1
  standar_error=1
  Y_iter=0
  for(Y in 1984:2020){
    Y_iter = Y_iter+1
    #drop 1 year
    IDX_drop <-which(WYs == Y)
    index1=index[-IDX_drop]
    dropped_DF=Climate_Fire_DF[index1,]
    
    #PCA combo model:
    Zs_dropped=as.vector(Zs)
    Zs_dropped = Zs_dropped[index1]
    
    Y_dropped = as.vector(BAs)
    Y_dropped = Y_dropped[index1]
    
    dropped_DF_pcs <- data.frame(Y=Y_dropped,predictor1=pc1[index1],predictor2=pc2[index1],predictor3=pc3[index1],predictor4=pc4[index1],predictor5=pc5[index1],predictor6=pc6[index1],predictor7=as.vector(Zs_dropped))
    
    #model
    fit_drop_gam = gam(mod$formula,data=dropped_DF_pcs,family="gaussian")
    #now estimate at the point that was dropped
    newdata = data.frame(predictor1=pc1[IDX_drop],predictor2=pc2[IDX_drop],predictor3=pc3[IDX_drop],predictor4=pc4[IDX_drop],predictor5=pc5[IDX_drop],predictor6=pc6[IDX_drop],predictor7=as.vector(Zs[IDX_drop]))
    yest_gam[IDX_drop]=predict(fit_drop_gam,newdata=newdata)
    #get the confidence interval:
    yhat <- predict(fit_drop_gam,newdata=newdata,se.fit = TRUE)
    standar_error[Y_iter] <- sum(yhat$se.fit)
  }
  IDX_neg <- which(yest_gam<0)
  if (length(IDX_neg) > 0){
    yest_gam[IDX_neg] <- min(yest_gam[-IDX_neg])
  }  
  #sum BAs across all Z bins:
  iter=0
  store_BA_obs=0
  store_BA_est=0
  for (Y in 1984:2020){
    iter <- iter+1
    current_IDX <-which(WYs == Y)
    store_BA_est[iter] = sum(yest_gam[current_IDX])
    store_BA_obs[iter] = sum(BAs[current_IDX])
  }
  #record drop 1-year taylor score:
  Drop1_Ann_R <- cor(store_BA_est,store_BA_obs)
  sigma_f = sd(store_BA_est)
  sigma_r = sd(store_BA_obs)
  Taylor_Score = (4*(1+Drop1_Ann_R)^4)/( ( (sigma_f/sigma_r)+(sigma_r/sigma_f))^2 *(1+1)^4)
  store_Taylor = c(store_Taylor,Taylor_Score)
  
}
output<-rbind(store_AIC,store_Taylor)
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_SWEImod_PCA_6vars_AICfit_Drop1Taylor.csv")
write.table(output,file=outfilename,sep=",",row.names = FALSE)

#7 variable model:
best_mods_idx = read.csv("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_bestmods_7vars.csv",sep=",",header=TRUE)
best_mods_idx=as.numeric(best_mods_idx)
nmods = length(best_mods_idx)
#for  combos of 7 variables:

n_predictors <- 7
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

SWEI_iter=0 #count which of the best mods include SWEI as predictor
good_mod_iter=0
Good_mods=c()
store_Taylor = c()
store_AIC=c()
n=length(BAs)
index = 1:n
for (j in 1:nmods){
  i <- best_mods_idx[j]
  Y = as.vector(BAs)
  
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
  Fit_AIC = mod$aic
  store_AIC = c(store_AIC,Fit_AIC)
  #GAM drop 1 year
  yest_gam=1
  standar_error=1
  Y_iter=0
  for(Y in 1984:2020){
    Y_iter = Y_iter+1
    #drop 1 year
    IDX_drop <-which(WYs == Y)
    index1=index[-IDX_drop]
    dropped_DF=Climate_Fire_DF[index1,]
    
    #PCA combo model:
    Zs_dropped=as.vector(Zs)
    Zs_dropped = Zs_dropped[index1]
    
    Y_dropped = as.vector(BAs)
    Y_dropped = Y_dropped[index1]
    
    dropped_DF_pcs <- data.frame(Y=Y_dropped,predictor1=pc1[index1],predictor2=pc2[index1],predictor3=pc3[index1],predictor4=pc4[index1],predictor5=pc5[index1],predictor6=pc6[index1],predictor7=pc7[index1],predictor8=as.vector(Zs_dropped))
    
    #model
    fit_drop_gam = gam(mod$formula,data=dropped_DF_pcs,family="gaussian")
    #now estimate at the point that was dropped
    newdata = data.frame(predictor1=pc1[IDX_drop],predictor2=pc2[IDX_drop],predictor3=pc3[IDX_drop],predictor4=pc4[IDX_drop],predictor5=pc5[IDX_drop],predictor6=pc6[IDX_drop],predictor7=pc7[IDX_drop],predictor8=as.vector(Zs[IDX_drop]))
    yest_gam[IDX_drop]=predict(fit_drop_gam,newdata=newdata)
    #get the confidence interval:
    yhat <- predict(fit_drop_gam,newdata=newdata,se.fit = TRUE)
    standar_error[Y_iter] <- sum(yhat$se.fit)
  }
  IDX_neg <- which(yest_gam<0)
  if (length(IDX_neg) > 0){
    yest_gam[IDX_neg] <- min(yest_gam[-IDX_neg])
  }  
  #sum BAs across all Z bins:
  iter=0
  store_BA_obs=0
  store_BA_est=0
  for (Y in 1984:2020){
    iter <- iter+1
    current_IDX <-which(WYs == Y)
    store_BA_est[iter] = sum(yest_gam[current_IDX])
    store_BA_obs[iter] = sum(BAs[current_IDX])
  }
  #record drop 1-year taylor score:
  Drop1_Ann_R <- cor(store_BA_est,store_BA_obs)
  sigma_f = sd(store_BA_est)
  sigma_r = sd(store_BA_obs)
  Taylor_Score = (4*(1+Drop1_Ann_R)^4)/( ( (sigma_f/sigma_r)+(sigma_r/sigma_f))^2 *(1+1)^4)
  store_Taylor = c(store_Taylor,Taylor_Score)
  
}
output<-rbind(store_AIC,store_Taylor)
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_SWEImod_PCA_7vars_AICfit_Drop1Taylor.csv")
write.table(output,file=outfilename,sep=",",row.names = FALSE)


#8 variable model:
best_mods_idx = read.csv("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_bestmods_8vars.csv",sep=",",header=TRUE)
best_mods_idx=as.numeric(best_mods_idx)
nmods = length(best_mods_idx)
#for  combos of 8 variables:

n_predictors <- 8
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

SWEI_iter=0 #count which of the best mods include SWEI as predictor
good_mod_iter=0
Good_mods=c()
store_Taylor = c()
store_AIC=c()
n=length(BAs)
index = 1:n
for (j in 1:nmods){
  i <- best_mods_idx[j]
  Y = as.vector(BAs)
  
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
  Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8)
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
  Fit_AIC = mod$aic
  store_AIC = c(store_AIC,Fit_AIC)
  #GAM drop 1 year
  yest_gam=1
  standar_error=1
  Y_iter=0
  for(Y in 1984:2020){
    Y_iter = Y_iter+1
    #drop 1 year
    IDX_drop <-which(WYs == Y)
    index1=index[-IDX_drop]
    dropped_DF=Climate_Fire_DF[index1,]
    
    #PCA combo model:
    Zs_dropped=as.vector(Zs)
    Zs_dropped = Zs_dropped[index1]
    
    Y_dropped = as.vector(BAs)
    Y_dropped = Y_dropped[index1]
    
    dropped_DF_pcs <- data.frame(Y=Y_dropped,predictor1=pc1[index1],predictor2=pc2[index1],predictor3=pc3[index1],predictor4=pc4[index1],predictor5=pc5[index1],predictor6=pc6[index1],predictor7=pc7[index1],predictor8=pc8[index1],predictor9=as.vector(Zs_dropped))
    
    #model
    fit_drop_gam = gam(mod$formula,data=dropped_DF_pcs,family="gaussian")
    #now estimate at the point that was dropped
    newdata = data.frame(predictor1=pc1[IDX_drop],predictor2=pc2[IDX_drop],predictor3=pc3[IDX_drop],predictor4=pc4[IDX_drop],predictor5=pc5[IDX_drop],predictor6=pc6[IDX_drop],predictor7=pc7[IDX_drop],predictor8=pc8[IDX_drop],predictor9=as.vector(Zs[IDX_drop]))
    yest_gam[IDX_drop]=predict(fit_drop_gam,newdata=newdata)
    #get the confidence interval:
    yhat <- predict(fit_drop_gam,newdata=newdata,se.fit = TRUE)
    standar_error[Y_iter] <- sum(yhat$se.fit)
  }
  IDX_neg <- which(yest_gam<0)
  if (length(IDX_neg) > 0){
    yest_gam[IDX_neg] <- min(yest_gam[-IDX_neg])
  }  
  #sum BAs across all Z bins:
  iter=0
  store_BA_obs=0
  store_BA_est=0
  for (Y in 1984:2020){
    iter <- iter+1
    current_IDX <-which(WYs == Y)
    store_BA_est[iter] = sum(yest_gam[current_IDX])
    store_BA_obs[iter] = sum(BAs[current_IDX])
  }
  #record drop 1-year taylor score:
  Drop1_Ann_R <- cor(store_BA_est,store_BA_obs)
  sigma_f = sd(store_BA_est)
  sigma_r = sd(store_BA_obs)
  Taylor_Score = (4*(1+Drop1_Ann_R)^4)/( ( (sigma_f/sigma_r)+(sigma_r/sigma_f))^2 *(1+1)^4)
  store_Taylor = c(store_Taylor,Taylor_Score)
  
}
output<-rbind(store_AIC,store_Taylor)
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_SWEImod_PCA_8vars_AICfit_Drop1Taylor.csv")
write.table(output,file=outfilename,sep=",",row.names = FALSE)

#9 variable model:
best_mods_idx = read.csv("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_bestmods_9vars.csv",sep=",",header=TRUE)
best_mods_idx=as.numeric(best_mods_idx)
nmods = length(best_mods_idx)
#for  combos of 9 variables:

n_predictors <- 9
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

SWEI_iter=0 #count which of the best mods include SWEI as predictor
good_mod_iter=0
Good_mods=c()
store_Taylor = c()
store_AIC=c()
n=length(BAs)
index = 1:n
for (j in 1:nmods){
  i <- best_mods_idx[j]
  Y = as.vector(BAs)
  
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
  Fit_AIC = mod$aic
  store_AIC = c(store_AIC,Fit_AIC)
  #GAM drop 1 year
  yest_gam=1
  standar_error=1
  Y_iter=0
  for(Y in 1984:2020){
    Y_iter = Y_iter+1
    #drop 1 year
    IDX_drop <-which(WYs == Y)
    index1=index[-IDX_drop]
    dropped_DF=Climate_Fire_DF[index1,]
    
    #PCA combo model:
    Zs_dropped=as.vector(Zs)
    Zs_dropped = Zs_dropped[index1]
    
    Y_dropped = as.vector(BAs)
    Y_dropped = Y_dropped[index1]
    
    dropped_DF_pcs <- data.frame(Y=Y_dropped,predictor1=pc1[index1],predictor2=pc2[index1],predictor3=pc3[index1],predictor4=pc4[index1],predictor5=pc5[index1],predictor6=pc6[index1],predictor7=pc7[index1],predictor8=pc8[index1],predictor9=pc9[index1],predictor10=as.vector(Zs_dropped))
    
    #model
    fit_drop_gam = gam(mod$formula,data=dropped_DF_pcs,family="gaussian")
    #now estimate at the point that was dropped
    newdata = data.frame(predictor1=pc1[IDX_drop],predictor2=pc2[IDX_drop],predictor3=pc3[IDX_drop],predictor4=pc4[IDX_drop],predictor5=pc5[IDX_drop],predictor6=pc6[IDX_drop],predictor7=pc7[IDX_drop],predictor8=pc8[IDX_drop],predictor9=pc9[IDX_drop],predictor10=as.vector(Zs[IDX_drop]))
    yest_gam[IDX_drop]=predict(fit_drop_gam,newdata=newdata)
    #get the confidence interval:
    yhat <- predict(fit_drop_gam,newdata=newdata,se.fit = TRUE)
    standar_error[Y_iter] <- sum(yhat$se.fit)
  }
  IDX_neg <- which(yest_gam<0)
  if (length(IDX_neg) > 0){
    yest_gam[IDX_neg] <- min(yest_gam[-IDX_neg])
  }  
  #sum BAs across all Z bins:
  iter=0
  store_BA_obs=0
  store_BA_est=0
  for (Y in 1984:2020){
    iter <- iter+1
    current_IDX <-which(WYs == Y)
    store_BA_est[iter] = sum(yest_gam[current_IDX])
    store_BA_obs[iter] = sum(BAs[current_IDX])
  }
  #record drop 1-year taylor score:
  Drop1_Ann_R <- cor(store_BA_est,store_BA_obs)
  sigma_f = sd(store_BA_est)
  sigma_r = sd(store_BA_obs)
  Taylor_Score = (4*(1+Drop1_Ann_R)^4)/( ( (sigma_f/sigma_r)+(sigma_r/sigma_f))^2 *(1+1)^4)
  store_Taylor = c(store_Taylor,Taylor_Score)
  
}
output<-rbind(store_AIC,store_Taylor)
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_SWEImod_PCA_9vars_AICfit_Drop1Taylor.csv")
write.table(output,file=outfilename,sep=",",row.names = FALSE)


#10 variable model:
best_mods_idx = read.csv("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_bestmods_10vars.csv",sep=",",header=TRUE)
best_mods_idx=as.numeric(best_mods_idx)
nmods = length(best_mods_idx)
#for  combos of 10 variables:

n_predictors <- 10
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

SWEI_iter=0 #count which of the best mods include SWEI as predictor
good_mod_iter=0
Good_mods=c()
store_Taylor = c()
store_AIC=c()
n=length(BAs)
index = 1:n
for (j in 1:nmods){
  i <- best_mods_idx[j]
  Y = as.vector(BAs)
  
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
  Fit_AIC = mod$aic
  store_AIC = c(store_AIC,Fit_AIC)
  #GAM drop 1 year
  yest_gam=1
  standar_error=1
  Y_iter=0
  for(Y in 1984:2020){
    Y_iter = Y_iter+1
    #drop 1 year
    IDX_drop <-which(WYs == Y)
    index1=index[-IDX_drop]
    dropped_DF=Climate_Fire_DF[index1,]
    
    #PCA combo model:
    Zs_dropped=as.vector(Zs)
    Zs_dropped = Zs_dropped[index1]
    
    Y_dropped = as.vector(BAs)
    Y_dropped = Y_dropped[index1]
    
    dropped_DF_pcs <- data.frame(Y=Y_dropped,predictor1=pc1[index1],predictor2=pc2[index1],predictor3=pc3[index1],predictor4=pc4[index1],predictor5=pc5[index1],predictor6=pc6[index1],predictor7=pc7[index1],predictor8=pc8[index1],predictor9=pc9[index1],predictor10=pc10[index1],predictor11=as.vector(Zs_dropped))
    
    #model
    fit_drop_gam = gam(mod$formula,data=dropped_DF_pcs,family="gaussian")
    #now estimate at the point that was dropped
    newdata = data.frame(predictor1=pc1[IDX_drop],predictor2=pc2[IDX_drop],predictor3=pc3[IDX_drop],predictor4=pc4[IDX_drop],predictor5=pc5[IDX_drop],predictor6=pc6[IDX_drop],predictor7=pc7[IDX_drop],predictor8=pc8[IDX_drop],predictor9=pc9[IDX_drop],predictor10=pc10[IDX_drop],predictor11=as.vector(Zs[IDX_drop]))
    yest_gam[IDX_drop]=predict(fit_drop_gam,newdata=newdata)
    #get the confidence interval:
    yhat <- predict(fit_drop_gam,newdata=newdata,se.fit = TRUE)
    standar_error[Y_iter] <- sum(yhat$se.fit)
  }
  IDX_neg <- which(yest_gam<0)
  if (length(IDX_neg) > 0){
    yest_gam[IDX_neg] <- min(yest_gam[-IDX_neg])
  }  
  #sum BAs across all Z bins:
  iter=0
  store_BA_obs=0
  store_BA_est=0
  for (Y in 1984:2020){
    iter <- iter+1
    current_IDX <-which(WYs == Y)
    store_BA_est[iter] = sum(yest_gam[current_IDX])
    store_BA_obs[iter] = sum(BAs[current_IDX])
  }
  #record drop 1-year taylor score:
  Drop1_Ann_R <- cor(store_BA_est,store_BA_obs)
  sigma_f = sd(store_BA_est)
  sigma_r = sd(store_BA_obs)
  Taylor_Score = (4*(1+Drop1_Ann_R)^4)/( ( (sigma_f/sigma_r)+(sigma_r/sigma_f))^2 *(1+1)^4)
  store_Taylor = c(store_Taylor,Taylor_Score)
  
}
output<-rbind(store_AIC,store_Taylor)
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_SWEImod_PCA_10vars_AICfit_Drop1Taylor.csv")
write.table(output,file=outfilename,sep=",",row.names = FALSE)

#11 variable model:
best_mods_idx = read.csv("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_bestmods_11vars.csv",sep=",",header=TRUE)
best_mods_idx=as.numeric(best_mods_idx)
nmods = length(best_mods_idx)
#for  combos of 11 variables:

n_predictors <- 11
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)

SWEI_iter=0 #count which of the best mods include SWEI as predictor
good_mod_iter=0
Good_mods=c()
store_Taylor = c()
store_AIC=c()
n=length(BAs)
index = 1:n
for (j in 1:nmods){
  i <- best_mods_idx[j]
  Y = as.vector(BAs)
  
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
  Fit_AIC = mod$aic
  store_AIC = c(store_AIC,Fit_AIC)
  #GAM drop 1 year
  yest_gam=1
  standar_error=1
  Y_iter=0
  for(Y in 1984:2020){
    Y_iter = Y_iter+1
    #drop 1 year
    IDX_drop <-which(WYs == Y)
    index1=index[-IDX_drop]
    dropped_DF=Climate_Fire_DF[index1,]
    
    #PCA combo model:
    Zs_dropped=as.vector(Zs)
    Zs_dropped = Zs_dropped[index1]
    
    Y_dropped = as.vector(BAs)
    Y_dropped = Y_dropped[index1]
    
    dropped_DF_pcs <- data.frame(Y=Y_dropped,predictor1=pc1[index1],predictor2=pc2[index1],predictor3=pc3[index1],predictor4=pc4[index1],predictor5=pc5[index1],predictor6=pc6[index1],predictor7=pc7[index1],predictor8=pc8[index1],predictor9=pc9[index1],predictor10=pc10[index1],predictor11=pc11[index1],predictor12=as.vector(Zs_dropped))
    S=dim(dropped_DF_pcs)
    if(S[1]!=111-3){
      S[1]
      stop('wrong dim')
    }
    #model
    fit_drop_gam = gam(mod$formula,data=dropped_DF_pcs,family="gaussian")
    #now estimate at the point that was dropped
    newdata = data.frame(predictor1=pc1[IDX_drop],predictor2=pc2[IDX_drop],predictor3=pc3[IDX_drop],predictor4=pc4[IDX_drop],predictor5=pc5[IDX_drop],predictor6=pc6[IDX_drop],predictor7=pc7[IDX_drop],predictor8=pc8[IDX_drop],predictor9=pc9[IDX_drop],predictor10=pc10[IDX_drop],predictor11=pc11[IDX_drop],predictor12=as.vector(Zs[IDX_drop]))
    yest_gam[IDX_drop]=predict(fit_drop_gam,newdata=newdata)
    #get the confidence interval:
    yhat <- predict(fit_drop_gam,newdata=newdata,se.fit = TRUE)
    standar_error[Y_iter] <- sum(yhat$se.fit)
  }
  IDX_neg <- which(yest_gam<0)
  if (length(IDX_neg) > 0){
    yest_gam[IDX_neg] <- min(yest_gam[-IDX_neg])
  }  
  #sum BAs across all Z bins:
  iter=0
  store_BA_obs=0
  store_BA_est=0
  for (Y in 1984:2020){
    iter <- iter+1
    current_IDX <-which(WYs == Y)
    store_BA_est[iter] = sum(yest_gam[current_IDX])
    store_BA_obs[iter] = sum(BAs[current_IDX])
  }
  #record drop 1-year taylor score:
  Drop1_Ann_R <- cor(store_BA_est,store_BA_obs)
  sigma_f = sd(store_BA_est)
  sigma_r = sd(store_BA_obs)
  Taylor_Score = (4*(1+Drop1_Ann_R)^4)/( ( (sigma_f/sigma_r)+(sigma_r/sigma_f))^2 *(1+1)^4)
  store_Taylor = c(store_Taylor,Taylor_Score)
  
}
output<-rbind(store_AIC,store_Taylor)
outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_SWEImod_PCA_11vars_AICfit_Drop1Taylor.csv")
write.table(output,file=outfilename,sep=",",row.names = FALSE)
