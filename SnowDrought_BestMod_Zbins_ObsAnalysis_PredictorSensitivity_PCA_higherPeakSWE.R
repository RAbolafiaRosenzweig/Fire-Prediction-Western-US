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
Climate_Fire_Data = read.table("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter1_higherPeakSWE.csv",sep=",")
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
Climate_Fire_Data = read.table("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter2_higherPeakSWE.csv",sep=",")
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
Climate_Fire_Data = read.table("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter3_higherPeakSWE.csv",sep=",")
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
Y = as.vector(BAs)
nvars <- 17 


n_predictors <- 8
ncombos = factorial(nvars)/(factorial(n_predictors)*factorial(nvars-n_predictors))
x=1:nvars
combo_IDs<-combn(x,n_predictors)


#calculate performance removing 1 climate variable at a time:

#load in the list of best models:
best_mod_ids = read.csv("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_bestmods_higherPeakSWE.csv",sep=",",header=TRUE)
best_mod_ids=as.numeric(best_mod_ids)
nmods = length(best_mod_ids)

for (j in 1:nmods){
  i = best_mod_ids[j]
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
  mod_original <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4)+s(predictor5)+s(predictor6)+s(predictor7)+s(predictor8)+predictor9,family="gaussian")
  yhat=predict(mod_original)
  
  iter=0
  yhat_ann=1
  Y_ann=1
  for (year in 1984:2020){
    iter <- iter+1
    current_IDX <-which(WYs == year)
    yhat_ann[iter] = sum(yhat[current_IDX])
    Y_ann[iter] = sum(Y[current_IDX])
  }
  R_fit=cor(yhat_ann,Y_ann)
  PBias_fit = mean(yhat - Y)/mean(Y)
  RMSE_fit = sqrt( sum( (yhat - Y)^2 )  /length(Y))
  
  iter=0
  store_Stats_rm=c()
  for (p in 1:n_predictors){
    if(p==1){
      #Remove predictor 1:
      Covariates_DF = data.frame(predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8)
    }
    if (p==2){
      #Remove predictor 2:
      Covariates_DF = data.frame(predictor1=predictor1,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8)
    }
    if (p==3){
      #Remove predictor 3:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8)
    }
    if (p==4){
      #Remove predictor 4:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8)
    }
    if (p==5){
      #Remove predictor 5:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor6=predictor6,predictor7=predictor7,predictor8=predictor8)
    }
    if (p==6){
      #Remove predictor 6:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor7=predictor7,predictor8=predictor8)
    }
    if (p==7){
      #Remove predictor 7:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor8=predictor8)
    }
    if (p==8){
      #Remove predictor 8:
      Covariates_DF = data.frame(predictor1=predictor1,predictor2=predictor2,predictor3=predictor3,predictor4=predictor4,predictor5=predictor5,predictor6=predictor6,predictor7=predictor7)
    }
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
    mod_rm <- gam(data=current_DF_pcs,Y~s(predictor1)+s(predictor2)+s(predictor3)+s(predictor4)+s(predictor5)+s(predictor6)+s(predictor7)+predictor8,family="gaussian")
    yhat=predict(mod_rm)
    
    for (year in 1984:2020){
      iter <- iter+1
      current_IDX <-which(WYs == year)
      yhat_ann[iter] = sum(yhat[current_IDX])
      Y_ann[iter] = sum(Y[current_IDX])
    }
    R_rm=cor(yhat_ann,Y_ann)
    PBias_rm = mean(yhat - Y)/mean(Y)
    RMSE_rm = sqrt( sum( (yhat - Y)^2 )  /length(Y))
    
    store_Stats_rm=rbind(store_Stats_rm,c(R_rm,PBias_rm,RMSE_rm) )
  }
  #store output: [predictor ID, R2-rm, R-fit(same for all predictors)]
  predictor_IDs=combo_IDs[,i]
  #R2:
  R_rm = store_Stats_rm[,1]
  R2_rm = R_rm^2
  R2_fit = R_fit^2
  R2_fit = rep(R2_fit,times=8)
  #PBIAS:
  PBIAS_rm = store_Stats_rm[,2]
  PBias_fit = rep(PBias_fit,times=8)
  #RMSE:
  RMSE_rm = store_Stats_rm[,3]
  RMSE_fit = rep(RMSE_fit,times=8)
  
  output = data.frame(predictor_IDs=predictor_IDs,R2_rm=R2_rm,R2_fit=R2_fit,PBIAS_rm=PBIAS_rm,PBias_fit=PBias_fit,RMSE_rm=RMSE_rm,RMSE_fit=RMSE_fit)
  outfilename = sprintf("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_drop_predictor_sensitivity_mod%d_higherPeakSWE.csv",i)
  write.table(output,file=outfilename,sep=",",row.names = FALSE)
}

