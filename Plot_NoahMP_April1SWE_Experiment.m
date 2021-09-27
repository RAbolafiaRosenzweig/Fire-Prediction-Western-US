clc;clear all;close all;

%% goals of this script:
%this script will plot summer SM that corresponds with each of the 3 NoahMP
%simulations:

%% Define BA grid:
latlim = [32 51];
lonlim = [-125 -104];
BA_lat = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200708.nc','XLAT_M');
BA_latvec = BA_lat(:);
BA_lon = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200708.nc','XLONG_M');
BA_lonvec = BA_lon(:);
%remove NaN:
idx_NAN_BAgrid = find(isnan(BA_latvec));
BA_lonvec(idx_NAN_BAgrid)=[];
BA_latvec(idx_NAN_BAgrid)=[];
%trim to bounding box
idx_BB = find(BA_latvec>=latlim(1) & BA_latvec<=latlim(2) & BA_lonvec>=lonlim(1) & BA_lonvec<=lonlim(2));

%% define screening info:
% % DEWS_IDX = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/DEWS_REgions_Indices.mat');
% % IN_pnw = DEWS_IDX.DATA.IN_pnw;
% % IN_ca_nv = DEWS_IDX.DATA.IN_ca_nv;
% % IN_west_mount = DEWS_IDX.DATA.IN_west_mount;
% % IN_missouri = DEWS_IDX.DATA.IN_missouri;
% % 
% % idx_pnw = find(IN_pnw==1);
% % idx_ca_nv = find(IN_ca_nv==1);
% % idx_west_mount = find(IN_west_mount==1);
% % idx_missouri = find(IN_missouri==1);
% % 
% % IDX_DEWS  = sort(unique([idx_pnw;idx_ca_nv;idx_west_mount;idx_missouri]));
% % 
% % %% calculate SWEI drought index values from peak SWE (Hunning and AghaKouchak, 2020)
% % UA_peakSWE_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_PeakSWE/';
% % IDX_NN_UAswe = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_UAswe.mat');
% % IDX_NN_UAswe = IDX_NN_UAswe.IDX_NN_UAswe;
% % idx_nan = find(isnan(IDX_NN_UAswe));
% % IDX_NN_UAswe(idx_nan)= [];
% % store_peak_swe=[];
% % for year = 1984:2020
% %     infilename = sprintf('UA_PeakSWE_%04d.mat',year);
% %     UA_PeakSWE = load([UA_peakSWE_dir,infilename]);
% %     UA_PeakSWE = UA_PeakSWE.Data.current_peak_SWE;
% %     idx = find(UA_PeakSWE == 0);
% %     UA_PeakSWE(idx) = rand(length(idx),1)*10^-10;
% %     store_peak_swe = [store_peak_swe,UA_PeakSWE(:)];
% % end
% % store_peak_swe = store_peak_swe(IDX_NN_UAswe,:);
% % store_peak_swe_BB = store_peak_swe(idx_BB,:);
% % 
% % %elevation data:
% % Elevation_coarse=load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Elevation_coarse.mat');
% % Elevation_coarse = Elevation_coarse.Elevation_coarse;
% % IDX_NN_Z = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_Elevation.mat');
% % IDX_NN_Z = IDX_NN_Z.IDX_NN_Elevation;
% % idx_nan = find(isnan(IDX_NN_Z));
% % IDX_NN_Z(idx_nan)=[];
% % Elevation_vec = Elevation_coarse(IDX_NN_Z);
% % Elevation_vec_BB = Elevation_vec(idx_BB);
% % Z_thresh_min = 0;
% % Z_thresh_max = 3300;
% % 
% % %% load in MODIS landcover type:
% % %from Get_MODIS_LC_class.m:
% % MODIS_LC = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/MODIS_Landcover/Gridded_MCD12Q1_Landcover.mat');
% % MODIS_LC = MODIS_LC.MCD12Q1.LC_class;
% % MODIS_LC = MODIS_LC(:);
% % MODIS_LC(idx_nan) = [];
% % MODIS_LC_BB = MODIS_LC(idx_BB);
% % %define index of veg types to keep:
% % idx_ENF = find(MODIS_LC_BB == 1);
% % idx_EBF = find(MODIS_LC_BB == 2);
% % idx_DNF = find(MODIS_LC_BB == 3);
% % idx_DBF = find(MODIS_LC_BB == 4);
% % idx_MF = find(MODIS_LC_BB == 5);
% % idx_CS = find(MODIS_LC_BB == 6);
% % idx_OS = find(MODIS_LC_BB == 7);
% % idx_WSav = find(MODIS_LC_BB == 8);
% % idx_Sav = find(MODIS_LC_BB == 9);
% % idx_grass = find(MODIS_LC_BB == 10);
% % disp('to justify this - make a pie chart of BA over time period in the snowy domain by veg type')
% % Total_LC_index_BB =[idx_ENF;idx_EBF;idx_DNF;idx_DBF;idx_MF;idx_WSav;idx_Sav;idx_grass];
% % Total_LC_index_BB = sort(Total_LC_index_BB);
% % idx_LC = Total_LC_index_BB;
% % 
% % %screen data based on elevation:
% % idx_Zthresh = find(Elevation_vec_BB>=Z_thresh_min & Elevation_vec_BB<=Z_thresh_max);
% % %combined screening:
% % peakSWE_thresh = 100; %mm
% % idx_PeakSWEthresh = find(max(store_peak_swe_BB(:,1:end)') >= peakSWE_thresh);
% % IDX_screened = intersect(idx_Zthresh,idx_PeakSWEthresh);
% % IDX_screened = intersect(IDX_screened,idx_LC);
% % IDX_screened = intersect(IDX_screened,IDX_DEWS);
% % save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_screened.mat','IDX_screened');

IDX_screened = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_screened.mat','IDX_screened');
IDX_screened = IDX_screened.IDX_screened;

%% Bring Noah-MP data to BA grid:
% % NoahMP_LAT = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NoahMPData/RawForcing/2010060103.LDASIN_DOMAIN1.nc','XLAT');
% % NoahMP_LAT = double(NoahMP_LAT);
% % NoahMP_LON = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NoahMPData/RawForcing/2010060103.LDASIN_DOMAIN1.nc','XLONG');
% % NoahMP_LON = double(NoahMP_LON);
% % 
% % BA_lat = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200708.nc','XLAT_M');
% % BA_lon = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200708.nc','XLONG_M');
% % 
% % NoahMP_GRID = [NoahMP_LAT(:)';NoahMP_LON(:)'];
% % BA_GRID = [BA_lat(:)';BA_lon(:)'];
% % IDX_NN_NoahMP=nearestneighbour(BA_GRID,NoahMP_GRID);
% % 
% % save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_NoahMP.mat','IDX_NN_NoahMP');

IDX_NN_NoahMP = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_NoahMP.mat');
IDX_NN_NoahMP=IDX_NN_NoahMP.IDX_NN_NoahMP;
idx_nan = find(isnan(IDX_NN_NoahMP));
IDX_NN_NoahMP(idx_nan)= [];
IDX_NN_NoahMP_BB = IDX_NN_NoahMP(idx_BB);

%% load in NoahMP basecase outputs:
% % Basecase_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NoahMPData/BasecaseOutputs/';
% % Basecase_SWE = ncread([Basecase_dir,'daily_2D_WY2012_Basecase_WY2012.nc'],'SWE_nanmean');
% % Basecase_SM = ncread([Basecase_dir,'daily_2D_WY2012_Basecase_WY2012.nc'],'SOIL_M');
% % Basecase_LH = ncread([Basecase_dir,'daily_2D_WY2012_Basecase_WY2012.nc'],'LH');
% % Basecase_QH = ncread([Basecase_dir,'daily_2D_WY2012_Basecase_WY2012.nc'],'HFX');
% % Basecase_EDIR = ncread([Basecase_dir,'daily_2D_WY2012_Basecase_WY2012.nc'],'EDIR');
% % Basecase_ETRAN = ncread([Basecase_dir,'daily_2D_WY2012_Basecase_WY2012.nc'],'ETRAN');
% % 
% % %% get spatially-averaged daily mean of each variables for the basecase:
% % ndays = 366;
% % dates = datenum([2011 10 1]):datenum([2012 9 30]);
% % 
% % store_basecase_data = [];
% % for d = 1:ndays
% %     %SWE:
% %     current_SWE = squeeze(Basecase_SWE(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_SWE < -1*10^30);
% %     current_SWE(idx) = NaN;
% %     current_SWE = current_SWE(IDX_NN_NoahMP_BB);
% %     current_SWE = current_SWE(IDX_screened);
% %     mean_SWE = nanmean(current_SWE);
% %     
% %     %SM:
% %     current_SM = squeeze(Basecase_SM(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_SM < -1*10^30);
% %     current_SM(idx) = NaN;
% %     current_SM = current_SM(IDX_NN_NoahMP_BB);
% %     current_SM = current_SM(IDX_screened);
% %     mean_SM = nanmean(current_SM);
% %     
% %     %LH:
% %     current_LH = squeeze(Basecase_LH(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_LH < -1*10^30);
% %     current_LH(idx) = NaN;
% %     current_LH = current_LH(IDX_NN_NoahMP_BB);
% %     current_LH = current_LH(IDX_screened);
% %     mean_LH = nanmean(current_LH);
% %     
% %     %QH:
% %     current_QH = squeeze(Basecase_QH(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_QH < -1*10^30);
% %     current_QH(idx) = NaN;
% %     current_QH = current_QH(IDX_NN_NoahMP_BB);
% %     current_QH = current_QH(IDX_screened);
% %     mean_QH = nanmean(current_QH);
% %     
% %     %ETRAN:
% %     current_ETRAN = squeeze(Basecase_ETRAN(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_ETRAN < -1*10^30);
% %     current_ETRAN(idx) = NaN;
% %     current_ETRAN = current_ETRAN(IDX_NN_NoahMP_BB);
% %     current_ETRAN = current_ETRAN(IDX_screened);
% %     mean_ETRAN = nanmean(current_ETRAN);
% %     
% %     %EDIR:
% %     current_EDIR = squeeze(Basecase_EDIR(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_EDIR < -1*10^30);
% %     current_EDIR(idx) = NaN;
% %     current_EDIR = current_EDIR(IDX_NN_NoahMP_BB);
% %     current_EDIR = current_EDIR(IDX_screened);
% %     mean_EDIR = nanmean(current_EDIR);
% %     
% %     store_basecase_data = [store_basecase_data;dates(d),mean_SWE,mean_SM,mean_LH,mean_ETRAN,mean_EDIR,mean_QH];
% % end
% % save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NoahMPData/BasecaseOutputs/store_basecase_data.mat','store_basecase_data');

store_basecase_data = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NoahMPData/BasecaseOutputs/store_basecase_data.mat');
store_basecase_data = store_basecase_data.store_basecase_data;

%% load in NoahMP DecrSWE outputs:
% % DecrSWE_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NoahMPData/DecrSWEOutputs/';
% % DecrSWE_SWE = ncread([DecrSWE_dir,'daily_2D_WY2012_Decr_April1SWE_WY2012.nc'],'SWE_nanmean');
% % DecrSWE_SM = ncread([DecrSWE_dir,'daily_2D_WY2012_Decr_April1SWE_WY2012.nc'],'SOIL_M');
% % DecrSWE_LH = ncread([DecrSWE_dir,'daily_2D_WY2012_Decr_April1SWE_WY2012.nc'],'LH');
% % DecrSWE_QH = ncread([DecrSWE_dir,'daily_2D_WY2012_Decr_April1SWE_WY2012.nc'],'HFX');
% % DecrSWE_EDIR = ncread([DecrSWE_dir,'daily_2D_WY2012_Decr_April1SWE_WY2012.nc'],'EDIR');
% % DecrSWE_ETRAN = ncread([DecrSWE_dir,'daily_2D_WY2012_Decr_April1SWE_WY2012.nc'],'ETRAN');
% % 
% % %% get spatially-averaged daily mean of each variables for the DecrSWE:
% % ndays = 366;
% % dates = datenum([2011 10 1]):datenum([2012 9 30]);
% % 
% % store_DecrSWE_data = [];
% % for d = 1:ndays
% %     %SWE:
% %     current_SWE = squeeze(DecrSWE_SWE(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_SWE < -1*10^30);
% %     current_SWE(idx) = NaN;
% %     current_SWE = current_SWE(IDX_NN_NoahMP_BB);
% %     current_SWE = current_SWE(IDX_screened);
% %     mean_SWE = nanmean(current_SWE);
% %     
% %     %SM:
% %     current_SM = squeeze(DecrSWE_SM(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_SM < -1*10^30);
% %     current_SM(idx) = NaN;
% %     current_SM = current_SM(IDX_NN_NoahMP_BB);
% %     current_SM = current_SM(IDX_screened);
% %     mean_SM = nanmean(current_SM);
% %     
% %     %LH:
% %     current_LH = squeeze(DecrSWE_LH(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_LH < -1*10^30);
% %     current_LH(idx) = NaN;
% %     current_LH = current_LH(IDX_NN_NoahMP_BB);
% %     current_LH = current_LH(IDX_screened);
% %     mean_LH = nanmean(current_LH);
% %     
% %     %QH:
% %     current_QH = squeeze(DecrSWE_QH(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_QH < -1*10^30);
% %     current_QH(idx) = NaN;
% %     current_QH = current_QH(IDX_NN_NoahMP_BB);
% %     current_QH = current_QH(IDX_screened);
% %     mean_QH = nanmean(current_QH);
% %     
% %     %ETRAN:
% %     current_ETRAN = squeeze(DecrSWE_ETRAN(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_ETRAN < -1*10^30);
% %     current_ETRAN(idx) = NaN;
% %     current_ETRAN = current_ETRAN(IDX_NN_NoahMP_BB);
% %     current_ETRAN = current_ETRAN(IDX_screened);
% %     mean_ETRAN = nanmean(current_ETRAN);
% %     
% %     %EDIR:
% %     current_EDIR = squeeze(DecrSWE_EDIR(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_EDIR < -1*10^30);
% %     current_EDIR(idx) = NaN;
% %     current_EDIR = current_EDIR(IDX_NN_NoahMP_BB);
% %     current_EDIR = current_EDIR(IDX_screened);
% %     mean_EDIR = nanmean(current_EDIR);
% %     
% %     store_DecrSWE_data = [store_DecrSWE_data;dates(d),mean_SWE,mean_SM,mean_LH,mean_ETRAN,mean_EDIR,mean_QH];
% % end
% % 
% % save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NoahMPData/DecrSWEOutputs/store_DecrSWE_data.mat','store_DecrSWE_data');

store_DecrSWE_data = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NoahMPData/DecrSWEOutputs/store_DecrSWE_data.mat');
store_DecrSWE_data = store_DecrSWE_data.store_DecrSWE_data;

%% load in NoahMP IncrSWE outputs:
% % IncrSWE_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NoahMPData/IncrSWEOutputs/';
% % IncrSWE_SWE = ncread([IncrSWE_dir,'daily_2D_WY2012_Incr_April1SWE_WY2012.nc'],'SWE_nanmean');
% % IncrSWE_SM = ncread([IncrSWE_dir,'daily_2D_WY2012_Incr_April1SWE_WY2012.nc'],'SOIL_M');
% % IncrSWE_LH = ncread([IncrSWE_dir,'daily_2D_WY2012_Incr_April1SWE_WY2012.nc'],'LH');
% % IncrSWE_QH = ncread([IncrSWE_dir,'daily_2D_WY2012_Incr_April1SWE_WY2012.nc'],'HFX');
% % IncrSWE_EDIR = ncread([IncrSWE_dir,'daily_2D_WY2012_Incr_April1SWE_WY2012.nc'],'EDIR');
% % IncrSWE_ETRAN = ncread([IncrSWE_dir,'daily_2D_WY2012_Incr_April1SWE_WY2012.nc'],'ETRAN');
% % 
% % %% get spatially-averaged daily mean of each variables for the IncrSWE:
% % ndays = 366;
% % dates = datenum([2011 10 1]):datenum([2012 9 30]);
% % 
% % store_IncrSWE_data = [];
% % for d = 1:ndays
% %     %SWE:
% %     current_SWE = squeeze(IncrSWE_SWE(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_SWE < -1*10^30);
% %     current_SWE(idx) = NaN;
% %     current_SWE = current_SWE(IDX_NN_NoahMP_BB);
% %     current_SWE = current_SWE(IDX_screened);
% %     mean_SWE = nanmean(current_SWE);
% %     
% %     %SM:
% %     current_SM = squeeze(IncrSWE_SM(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_SM < -1*10^30);
% %     current_SM(idx) = NaN;
% %     current_SM = current_SM(IDX_NN_NoahMP_BB);
% %     current_SM = current_SM(IDX_screened);
% %     mean_SM = nanmean(current_SM);
% %     
% %     %LH:
% %     current_LH = squeeze(IncrSWE_LH(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_LH < -1*10^30);
% %     current_LH(idx) = NaN;
% %     current_LH = current_LH(IDX_NN_NoahMP_BB);
% %     current_LH = current_LH(IDX_screened);
% %     mean_LH = nanmean(current_LH);
% %     
% %     %QH:
% %     current_QH = squeeze(IncrSWE_QH(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_QH < -1*10^30);
% %     current_QH(idx) = NaN;
% %     current_QH = current_QH(IDX_NN_NoahMP_BB);
% %     current_QH = current_QH(IDX_screened);
% %     mean_QH = nanmean(current_QH);
% %     
% %     %ETRAN:
% %     current_ETRAN = squeeze(IncrSWE_ETRAN(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_ETRAN < -1*10^30);
% %     current_ETRAN(idx) = NaN;
% %     current_ETRAN = current_ETRAN(IDX_NN_NoahMP_BB);
% %     current_ETRAN = current_ETRAN(IDX_screened);
% %     mean_ETRAN = nanmean(current_ETRAN);
% %     
% %     %EDIR:
% %     current_EDIR = squeeze(IncrSWE_EDIR(:,:,d));
% %     %replace fill value with NaN:
% %     idx = find(current_EDIR < -1*10^30);
% %     current_EDIR(idx) = NaN;
% %     current_EDIR = current_EDIR(IDX_NN_NoahMP_BB);
% %     current_EDIR = current_EDIR(IDX_screened);
% %     mean_EDIR = nanmean(current_EDIR);
% %     
% %     store_IncrSWE_data = [store_IncrSWE_data;dates(d),mean_SWE,mean_SM,mean_LH,mean_ETRAN,mean_EDIR,mean_QH];
% % end
% % save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NoahMPData/IncrSWEOutputs/store_IncrSWE_data.mat','store_IncrSWE_data');

store_IncrSWE_data = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NoahMPData/IncrSWEOutputs/store_IncrSWE_data.mat');
store_IncrSWE_data = store_IncrSWE_data.store_IncrSWE_data;

%% plot data:
%basecase data:
dates = store_basecase_data(:,1);
basecase_SWE = store_basecase_data(:,2);
basecase_SM = store_basecase_data(:,3);
basecase_LH = store_basecase_data(:,4);
basecase_ETRAN = store_basecase_data(:,5);
basecase_EDIR = store_basecase_data(:,6);
basecase_QH = store_basecase_data(:,7);

basecase_ET = basecase_ETRAN+basecase_EDIR;
basecase_BR = basecase_QH./basecase_LH;

%decreased SWE data:
DecrSWE_SWE = store_DecrSWE_data(:,2);
DecrSWE_SM = store_DecrSWE_data(:,3);
DecrSWE_LH = store_DecrSWE_data(:,4);
DecrSWE_ETRAN = store_DecrSWE_data(:,5);
DecrSWE_EDIR = store_DecrSWE_data(:,6);
DecrSWE_QH = store_DecrSWE_data(:,7);

DecrSWE_ET = DecrSWE_ETRAN+DecrSWE_EDIR;
DecrSWE_BR = DecrSWE_QH./DecrSWE_LH;

%increased SWE data:
IncrSWE_SWE = store_IncrSWE_data(:,2);
IncrSWE_SM = store_IncrSWE_data(:,3);
IncrSWE_LH = store_IncrSWE_data(:,4);
IncrSWE_ETRAN = store_IncrSWE_data(:,5);
IncrSWE_EDIR = store_IncrSWE_data(:,6);
IncrSWE_QH = store_IncrSWE_data(:,7);

IncrSWE_ET = IncrSWE_ETRAN+IncrSWE_EDIR;
IncrSWE_BR = IncrSWE_QH./IncrSWE_LH;

f=figure;
f.Position = [-1840         -20        1766         804];
pp
subplot(3,1,1)
hold on
p1=plot(dates,DecrSWE_SWE-basecase_SWE,'-r','linewidth',2.5);
p2=plot(dates,IncrSWE_SWE-basecase_SWE,'-b','linewidth',2.5);
xticks([datenum([2012 4 1]),datenum([2012 5 1]),datenum([2012 6 1]),datenum([2012 7 1]),datenum([2012 8 1]),datenum([2012 9 1])]);
datetick('x','yyyy/mm','keepticks','keeplimits')
plot(linspace(0,datenum([2012 10 1]),1000),linspace(0,0,1000),'--k','linewidth',2.5)
set(gca,'fontsize',20)
ylabel({'Departure from'; 'Reference SWE (mm)'},'fontsize',20)
grid on
legend([p1 p2],{'-30% April 1 SWE','+30% April 1 SWE'},'fontsize',20)
xlim([ datenum([2012 6 1]) datenum([2012 9 30]) ] )

subplot(3,1,2)
hold on
p1=plot(dates,DecrSWE_SM-basecase_SM,'-r','linewidth',2.5);
p2=plot(dates,IncrSWE_SM-basecase_SM,'-b','linewidth',2.5);
xticks([datenum([2012 4 1]),datenum([2012 5 1]),datenum([2012 6 1]),datenum([2012 7 1]),datenum([2012 8 1]),datenum([2012 9 1])]);
datetick('x','yyyy/mm','keepticks','keeplimits')
plot(linspace(0,datenum([2012 10 1]),1000),linspace(0,0,1000),'--k','linewidth',2.5)
set(gca,'fontsize',20)
ylabel({'Departure from'; 'Reference SM'; '(mm^{3}/mm^{3})'},'fontsize',20)
grid on
xlim([ datenum([2012 6 1]) datenum([2012 9 30]) ] )
subplot(3,1,3)
hold on
p1=plot(dates,DecrSWE_BR - basecase_BR,'-r','linewidth',2.5);
p2=plot(dates,IncrSWE_BR - basecase_BR,'-b','linewidth',2.5);
xticks([datenum([2012 4 1]),datenum([2012 5 1]),datenum([2012 6 1]),datenum([2012 7 1]),datenum([2012 8 1]),datenum([2012 9 1])]);
datetick('x','yyyy/mm','keepticks','keeplimits')
plot(linspace(0,datenum([2012 10 1]),1000),linspace(0,0,1000),'--k','linewidth',2.5)
set(gca,'fontsize',20)
ylabel({'Departure from'; 'Reference B.R.'},'fontsize',20)
grid on
xlim([ datenum([2012 6 1]) datenum([2012 9 30]) ] )

saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/NoahMP_Experiment_DeltaVars_Timeseries.png')

f=figure;

subplot(2,1,1)
hold on
p1=plot(dates,DecrSWE_SM-basecase_SM,'-r','linewidth',2.5);
p2=plot(dates,IncrSWE_SM-basecase_SM,'-b','linewidth',2.5);
xticks([datenum([2012 4 1]),datenum([2012 5 1]),datenum([2012 6 1]),datenum([2012 7 1]),datenum([2012 8 1]),datenum([2012 9 1])]);
datetick('x','yyyy/mm','keepticks','keeplimits')
plot(linspace(0,datenum([2012 10 1]),1000),linspace(0,0,1000),'--k','linewidth',2.5)
set(gca,'fontsize',20)
ylabel({'Departure from'; 'Reference SM'; '(mm^{3}/mm^{3})'},'fontsize',20)
leg=legend([p1 p2],{'-30% April 1 SWE simulation','+30% April 1 SWE simulation'},'fontsize',24);
grid on
box on
xlim([ datenum([2012 6 1]) datenum([2012 9 30]) ] )

subplot(2,1,2)
hold on
p1=plot(dates,DecrSWE_BR - basecase_BR,'-r','linewidth',2.5);
p2=plot(dates,IncrSWE_BR - basecase_BR,'-b','linewidth',2.5);
xticks([datenum([2012 4 1]),datenum([2012 5 1]),datenum([2012 6 1]),datenum([2012 7 1]),datenum([2012 8 1]),datenum([2012 9 1])]);
datetick('x','yyyy/mm','keepticks','keeplimits')
plot(linspace(0,datenum([2012 10 1]),1000),linspace(0,0,1000),'--k','linewidth',2.5)
set(gca,'fontsize',20)
ylabel({'Departure from'; 'Reference B.R.'},'fontsize',20)
grid on
box on
xlim([ datenum([2012 6 1]) datenum([2012 9 30]) ] )

f.Position = [-1840         305        1766         479];
leg.Position = [0.7024    0.5859    0.2019    0.1284];
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/NoahMP_Experiment_DeltaVars_Timeseries_SM_and_BR.eps','epsc')

%% calculate difference in melt out:
idx_decrease_DOD=find(DecrSWE_SWE<1);
idx=find(idx_decrease_DOD >150);
idx_decrease_DOD = idx_decrease_DOD(idx(1))

idx_increase_DOD=find(IncrSWE_SWE<1);
idx=find(idx_increase_DOD >150);
idx_increase_DOD = idx_increase_DOD(idx(1))
%% literature to relate LH results to the VPD conclusion:
%source: https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1002/2016JD025855
%discusses how limited ET from limited SM can increase VPD
%increasing VPD increases ET

%source: file:///Users/abolafia/Downloads/[15588432%20-%20Journal%20of%20Applied%20Meteorology%20and%20Climatology]%20Climatology,%20Variability,%20and%20Trends%20in%20the%20U.S.%20Vapor%20Pressure%20Deficit,%20an%20Important%20Fire-Related%20Meteorological%20Quantity.pdf
%limited SM reduces ET which increases VPD

%source: https://agupubs.onlinelibrary.wiley.com/doi/pdfdirect/10.1029/2019MS001790:
%increasing VPD increases ET

%source: https://journals.ametsoc.org/view/journals/apme/54/6/jamc-d-14-0321.1.xml
%"As the soil dries out, incoming solar radiation needs to be increasingly balanced by sensible and longwave radiative heat loss, and less by evapotranspiration. This requires an increase in surface temperature and less moisture flux from the surface to the atmosphere, both effects that increase VPD."
