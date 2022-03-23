clc;clear all;close all;

%load data:

%% Z = 0-1100
Climate_Fire_Data = csvread("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter1.csv");
WYs = Climate_Fire_Data(:,1);
BAs = Climate_Fire_Data(:,2);
BAs_Z1 = (BAs*0.000247105) /(10^6);
SWEIs_Z1 = Climate_Fire_Data(:,3);
WinPRCP_Z1 = Climate_Fire_Data(:,4);
WinTMP_Z1 = Climate_Fire_Data(:,5);
WinterPDSI_Z1 = Climate_Fire_Data(:,6);
Spring_PRCP_Z1 = Climate_Fire_Data(:,7); 
Spring_TMP_Z1 = Climate_Fire_Data(:,8); 
Spring_PDSI_Z1 = Climate_Fire_Data(:,9); 
Spring_SWEI_Z1 = Climate_Fire_Data(:,10); 
Spring_VPD_Z1 = Climate_Fire_Data(:,11); 
Winter_VPD_Z1 = Climate_Fire_Data(:,12); 
MODIS_BA_Z1 = Climate_Fire_Data(:,13); 
MODIS_BA_Z1 = (MODIS_BA_Z1*0.000247105) /(10^6);
Spring_ET_Z1 = Climate_Fire_Data(:,14); 
Spring_PET_Z1 = Climate_Fire_Data(:,15); 
Spring_PETminusET_Z1 = Climate_Fire_Data(:,16);
Winter_ET_Z1 = Climate_Fire_Data(:,17); 
Winter_PET_Z1 = Climate_Fire_Data(:,18); 
Winter_DroughtArea_Z1 = Climate_Fire_Data(:,28);
Spring_DroughtArea_Z1 = Climate_Fire_Data(:,29);

%relate MODIS and MTBS BA linearly to allow MODIS to gap fill MTBS 2020 BA:
BA_MTBS = BAs_Z1(18:36);
BA_MODIS = MODIS_BA_Z1(18:36);
p = polyfit(BA_MODIS,BA_MTBS,1);
BA_lm = @(BA_MODIS) BA_MODIS*p(1) + p(2);
BA_2020_MTBS = BA_lm(MODIS_BA_Z1(37));
BAs_Z1(37) = BA_2020_MTBS;
if (BAs_Z1(37) < 0)
  BAs_Z1(37) = 0;
end

%% Z = 1100-2200
Climate_Fire_Data = csvread("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter2.csv");
WYs = Climate_Fire_Data(:,1);
BAs = Climate_Fire_Data(:,2);
BAs_Z2 = (BAs*0.000247105) /(10^6);
SWEIs_Z2 = Climate_Fire_Data(:,3);
WinPRCP_Z2 = Climate_Fire_Data(:,4);
WinTMP_Z2 = Climate_Fire_Data(:,5);
WinterPDSI_Z2 = Climate_Fire_Data(:,6);
Spring_PRCP_Z2 = Climate_Fire_Data(:,7); 
Spring_TMP_Z2 = Climate_Fire_Data(:,8); 
Spring_PDSI_Z2 = Climate_Fire_Data(:,9); 
Spring_SWEI_Z2 = Climate_Fire_Data(:,10); 
Spring_VPD_Z2 = Climate_Fire_Data(:,11); 
Winter_VPD_Z2 = Climate_Fire_Data(:,12); 
MODIS_BA_Z2 = Climate_Fire_Data(:,13); 
MODIS_BA_Z2 = (MODIS_BA_Z2*0.000247105) /(10^6);
Spring_ET_Z2 = Climate_Fire_Data(:,14); 
Spring_PET_Z2 = Climate_Fire_Data(:,15); 
Spring_PETminusET_Z2 = Climate_Fire_Data(:,16);
Winter_ET_Z2 = Climate_Fire_Data(:,17); 
Winter_PET_Z2 = Climate_Fire_Data(:,18); 
Winter_DroughtArea_Z2 = Climate_Fire_Data(:,28);
Spring_DroughtArea_Z2 = Climate_Fire_Data(:,29);

%relate MODIS and MTBS BA linearly to allow MODIS to gap fill MTBS 2020 BA:
BA_MTBS = BAs_Z2(18:36);
BA_MODIS = MODIS_BA_Z2(18:36);
p = polyfit(BA_MODIS,BA_MTBS,1);
BA_lm = @(BA_MODIS) BA_MODIS*p(1) + p(2);
BA_2020_MTBS = BA_lm(MODIS_BA_Z2(37));
BAs_Z2(37) = BA_2020_MTBS;
if (BAs_Z2(37) < 0)
  BAs_Z2(37) = 0;
end

%% Z = 2200-3300
Climate_Fire_Data = csvread("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter3.csv");
WYs = Climate_Fire_Data(:,1);
BAs = Climate_Fire_Data(:,2);
BAs_Z3 = (BAs*0.000247105) /(10^6);
SWEIs_Z3 = Climate_Fire_Data(:,3);
WinPRCP_Z3 = Climate_Fire_Data(:,4);
WinTMP_Z3 = Climate_Fire_Data(:,5);
WinterPDSI_Z3 = Climate_Fire_Data(:,6);
Spring_PRCP_Z3 = Climate_Fire_Data(:,7); 
Spring_TMP_Z3 = Climate_Fire_Data(:,8); 
Spring_PDSI_Z3 = Climate_Fire_Data(:,9); 
Spring_SWEI_Z3 = Climate_Fire_Data(:,10); 
Spring_VPD_Z3 = Climate_Fire_Data(:,11); 
Winter_VPD_Z3 = Climate_Fire_Data(:,12); 
MODIS_BA_Z3 = Climate_Fire_Data(:,13); 
MODIS_BA_Z3 = (MODIS_BA_Z3*0.000247105) /(10^6);
Spring_ET_Z3 = Climate_Fire_Data(:,14); 
Spring_PET_Z3 = Climate_Fire_Data(:,15); 
Spring_PETminusET_Z3 = Climate_Fire_Data(:,16);
Winter_ET_Z3 = Climate_Fire_Data(:,17); 
Winter_PET_Z3 = Climate_Fire_Data(:,18); 
Winter_DroughtArea_Z3 = Climate_Fire_Data(:,28);
Spring_DroughtArea_Z3 = Climate_Fire_Data(:,29);

%relate MODIS and MTBS BA linearly to allow MODIS to gap fill MTBS 2020 BA:
BA_MTBS = BAs_Z3(18:36);
BA_MODIS = MODIS_BA_Z3(18:36);
p = polyfit(BA_MODIS,BA_MTBS,1);
BA_lm = @(BA_MODIS) BA_MODIS*p(1) + p(2);
BA_2020_MTBS = BA_lm(MODIS_BA_Z3(37));
BAs_Z3(37) = BA_2020_MTBS;
if (BAs_Z3(37) < 0)
  BAs_Z3(37) = 0;
end

%% aggregate data for best predictors:
Spring_DroughtArea = [Spring_DroughtArea_Z1;Spring_DroughtArea_Z2;Spring_DroughtArea_Z3];
WinTMP = [WinTMP_Z1;WinTMP_Z2;WinTMP_Z3];
Winter_ET = [Winter_ET_Z1;Winter_ET_Z2;Winter_ET_Z3];
WinterPDSI = [WinterPDSI_Z1;WinterPDSI_Z2;WinterPDSI_Z3];
SWEIs = [SWEIs_Z1;SWEIs_Z2;SWEIs_Z3];
Spring_ET = [Spring_ET_Z1;Spring_ET_Z2;Spring_ET_Z3];
Spring_PDSI = [Spring_PDSI_Z1;Spring_PDSI_Z2;Spring_PDSI_Z3];
Spring_PET=[Spring_PET_Z1;Spring_PET_Z2;Spring_PET_Z3];
BAs = [BAs_Z1;BAs_Z2;BAs_Z3];
    
vars={'BAs','Spring_DroughtArea','WinTMP','Winter_ET','WinterPDSI','SWEIs','Spring_ET','Spring_PDSI','Spring_PET'};
All_Data = [BAs,Spring_DroughtArea,WinTMP,Winter_ET,WinterPDSI,SWEIs,Spring_ET,Spring_PDSI,Spring_PET];

%% calculate their correlations with eachother and with BA:
nvars = 9;
%rows = correlation with each var
%1 col for each variable
store_R = nan(nvars-1,nvars);
for i=1:nvars
    current_var = All_Data(:,i);
    cols = 1:nvars;
    cols(i) = [];
    for j=1:nvars-1
        comparing_var = All_Data(:,cols(j));
        R = corr(current_var,comparing_var);
        store_R(j,i) = R;
    end
end

%report max cor for each var:
Max_R = max(store_R);

%% report results:
for i=1:nvars
    IDX = 1:nvars;
    IDX(i) =[];
    current_R  = abs(store_R(:,i));
    idx = find(current_R == Max_R(i));
    idx = IDX(idx);
    sprintf('%s correlates highest with %s (R=%.2f)',vars{i},vars{idx},Max_R(i))
end


