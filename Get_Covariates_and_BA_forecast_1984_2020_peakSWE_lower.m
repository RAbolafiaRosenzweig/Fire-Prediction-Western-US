clc;clear all;close all;

%% This script considers the relationship between burn area and snow drought
%the study domain is the western CONUS for areas of high elevation, high
%LAI, and annual SWE cycles


%% define bounding box for area of analysis: western CONUS:
latlim = [32 51];
lonlim = [-125 -104];

%% Define thresholds for study domain
Z_thresh_mins = 0:1100:2200;
Z_thresh_maxs = 1100:1100:3300;

%arbitrary
peakLAI_thresh = 0.0;

%arbitrary
peakSWE_thresh = 50; %mm

%% Define BA grid:
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

%% trim by  DEWS regions for domain separation:
% % DEWS_areas = shaperead('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Raw/DEWS_Regions/DEWS/DEWS_AllRegions202103.shp');
% % DEWS_PNW = DEWS_areas(6);
% % DEWS_CA_NV = DEWS_areas(2);
% % DEWS_West_mount = DEWS_areas(3);
% % DEWS_missouri = DEWS_areas(5);
% %
% % BA_coords_BB = [BA_latvec(idx_BB),BA_lonvec(idx_BB)];
% % IN_pnw =inpolygon(BA_coords_BB(:,2),BA_coords_BB(:,1),DEWS_PNW.X,DEWS_PNW.Y);
% % IN_ca_nv =inpolygon(BA_coords_BB(:,2),BA_coords_BB(:,1),DEWS_CA_NV.X,DEWS_CA_NV.Y);
% % IN_west_mount =inpolygon(BA_coords_BB(:,2),BA_coords_BB(:,1),DEWS_West_mount.X,DEWS_West_mount.Y);
% % IN_missouri =inpolygon(BA_coords_BB(:,2),BA_coords_BB(:,1),DEWS_missouri.X,DEWS_missouri.Y);
% % DATA.IN_pnw = IN_pnw;
% % DATA.IN_ca_nv = IN_ca_nv;
% % DATA.IN_west_mount = IN_west_mount;
% % DATA.IN_missouri = IN_missouri;
% % outfilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/DEWS_REgions_Indices.mat');
% % save(outfilename,'DATA');

DEWS_IDX = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/DEWS_REgions_Indices.mat');
IN_pnw = DEWS_IDX.DATA.IN_pnw;
IN_ca_nv = DEWS_IDX.DATA.IN_ca_nv;
IN_west_mount = DEWS_IDX.DATA.IN_west_mount;
IN_missouri = DEWS_IDX.DATA.IN_missouri;

idx_pnw = find(IN_pnw==1);
idx_ca_nv = find(IN_ca_nv==1);
idx_west_mount = find(IN_west_mount==1);
idx_missouri = find(IN_missouri==1);

IDX_DEWS  = sort(unique([idx_pnw;idx_ca_nv;idx_west_mount;idx_missouri]));
%% define peak SWE directory (Get_UA_Annual_Peak_SWE.m):
UA_peakSWE_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_PeakSWE/';
%get nearestneighbor index to match with BA grid (NN_UA_and_PRISM_to_BAgrid.m)
IDX_NN_UAswe = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_UAswe.mat');
IDX_NN_UAswe = IDX_NN_UAswe.IDX_NN_UAswe;
idx_nan = find(isnan(IDX_NN_UAswe));
IDX_NN_UAswe(idx_nan)= [];

%% calculate SWEI drought index values from peak SWE (Hunning and AghaKouchak, 2020)
store_peak_swe=[];
for year = 1984:2020
    infilename = sprintf('UA_PeakSWE_%04d.mat',year);
    UA_PeakSWE = load([UA_peakSWE_dir,infilename]);
    UA_PeakSWE = UA_PeakSWE.Data.current_peak_SWE;
    idx = find(UA_PeakSWE == 0);
    UA_PeakSWE(idx) = rand(length(idx),1)*10^-10;
    store_peak_swe = [store_peak_swe,UA_PeakSWE(:)];
end
store_peak_swe = store_peak_swe(IDX_NN_UAswe,:);
store_peak_swe_BB = store_peak_swe(idx_BB,:);

% % %rank the peak SWE values across the years for each pixel:
% % [~,ii] = sort(store_peak_swe_BB,2);
% % [~,i] = sort(ii,2);
% % N=size(i,2);
% % p = (i - 0.44)./(N+0.12);
% %
% % %calculate nonparametric SWEI by transforming the empirical probability, p,
% % %to the standard normal distribution:
% % store_SWEI = NaN(length(p),37);
% % for i = 1:length(p)
% %     current_pixel_p = p(i,:);
% %     store_SWEI(i,:) = norminv(current_pixel_p);
% % end
% % %export:
% % save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_PeakSWE_SWEI/UA_SWEI_Hindcast.mat','store_SWEI', '-v7.3');
store_SWEI = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_PeakSWE_SWEI/UA_SWEI_Hindcast.mat','store_SWEI');
store_SWEI = store_SWEI.store_SWEI;

%% calculate SWEI drought index values from spring SWE (Hunning and AghaKouchak, 2020)
%from Get_UA_Total_Spring_SWE.m
% % store_spring_swe=[];
% % for year = 1984:2020
% %     infilename = sprintf('UA_SpringSWE_%04d.mat',year);
% %     UA_PeakSWE = load([UA_peakSWE_dir,infilename]);
% %     UA_PeakSWE = UA_PeakSWE.Data.current_spring_SWE;
% %     idx = find(UA_PeakSWE == 0);
% %     UA_PeakSWE(idx) = rand(length(idx),1)*10^-10;
% %     store_spring_swe = [store_spring_swe,UA_PeakSWE(:)];
% % end
% % store_spring_swe = store_spring_swe(IDX_NN_UAswe,:);
% % store_spring_swe = store_spring_swe(idx_BB,:);

% % % % %rank the peak SWE values across the years for each pixel:
% % [~,ii] = sort(store_spring_swe,2);
% % [~,i] = sort(ii,2);
% % N=size(i,2);
% % p = (i - 0.44)./(N+0.12);
% %
% % %calculate nonparametric SWEI by transforming the empirical probability, p,
% % %to the standard normal distribution:
% % store_SWEI = NaN(length(p),37);
% % for i = 1:length(p)
% %     current_pixel_p = p(i,:);
% %     store_SWEI(i,:) = norminv(current_pixel_p);
% % end
%export:
% % save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_SpringSWE_SWEI/UA_Spring_SWEI_Hindcast.mat','store_SWEI', '-v7.3');
store_SWEI_spring = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_SpringSWE_SWEI/UA_Spring_SWEI_Hindcast.mat','store_SWEI');
store_SWEI_spring = store_SWEI_spring.store_SWEI;

%% Load in elevation grid: (from Map_elevation_to_BA_product.m)
Elevation_coarse=load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Elevation_coarse.mat');
Elevation_coarse = Elevation_coarse.Elevation_coarse;
IDX_NN_Z = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_Elevation.mat');
IDX_NN_Z = IDX_NN_Z.IDX_NN_Elevation;
idx_nan = find(isnan(IDX_NN_Z));
IDX_NN_Z(idx_nan)=[];
Elevation_vec = Elevation_coarse(IDX_NN_Z);
Elevation_vec_BB = Elevation_vec(idx_BB);

%% define temp and precip directory (Get_Annual_PRISM_data.m)
PRISM_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/PRISM_Annual/';
%get nearestneighbor index to match with BA grid (NN_UA_and_PRISM_to_BAgrid.m)
IDX_NN_PRISM = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_PRISM.mat');
IDX_NN_PRISM = IDX_NN_PRISM.IDX_NN_PRISM;
idx_nan = find(isnan(IDX_NN_PRISM));
IDX_NN_PRISM(idx_nan)=[];
IDX_NN_PRISM_BB = IDX_NN_PRISM(idx_BB);

%get nearestneighbor index to match with BA grid (NN_NLDAS_to_BAgrid.m)
IDX_NN_NLDAS = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_NLDAS.mat');
IDX_NN_NLDAS = IDX_NN_NLDAS.IDX_NN_NLDAS;
idx_nan = find(isnan(IDX_NN_NLDAS));
IDX_NN_NLDAS(idx_nan)= [];
IDX_NN_NLDAS_BB = IDX_NN_NLDAS(idx_BB);

%% load in MODIS landcover type:
%from Get_MODIS_LC_class.m:
MODIS_LC = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/MODIS_Landcover/Gridded_MCD12Q1_Landcover.mat');
MODIS_LC = MODIS_LC.MCD12Q1.LC_class;
MODIS_LC = MODIS_LC(:);
MODIS_LC(idx_nan) = [];
MODIS_LC_BB = MODIS_LC(idx_BB);
%define index of veg types to keep:
idx_ENF = find(MODIS_LC_BB == 1);
idx_EBF = find(MODIS_LC_BB == 2);
idx_DNF = find(MODIS_LC_BB == 3);
idx_DBF = find(MODIS_LC_BB == 4);
idx_MF = find(MODIS_LC_BB == 5);
idx_CS = find(MODIS_LC_BB == 6);
idx_OS = find(MODIS_LC_BB == 7);
idx_WSav = find(MODIS_LC_BB == 8);
idx_Sav = find(MODIS_LC_BB == 9);
idx_grass = find(MODIS_LC_BB == 10);
disp('to justify this - make a pie chart of BA over time period in the snowy domain by veg type')
Total_LC_index_BB =[idx_ENF;idx_EBF;idx_DNF;idx_DBF;idx_MF;idx_WSav;idx_Sav;idx_grass];
Total_LC_index_BB = sort(Total_LC_index_BB);

%% define merged BA product directory (/Users/abolafia/Drought_Fire_Snow/tools/Merged_BA)
MTBS_BA_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Monthly_gridded_MTBS/';
MODIS_BA_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Monthly_gridded_MCD64/';

%% Compare summer (June-September) burn area with SWEI:
%screen based on peak SWE:
%     idx_PeakSWEthresh = find(max(store_peak_swe_BB(:,1:end)') >= peakSWE_thresh);
idx_PeakSWEthresh = find(max(store_peak_swe_BB(:,1:end)') >= peakSWE_thresh);
%screen based on landcover type:
idx_LC = Total_LC_index_BB;

for Ziter = 1:length(Z_thresh_mins)
    Ziter
    %z-threshold based on Westerling et al. (2006)
    % % % Z_thresh_min = 1680; %m
    % % % Z_thresh_max = 2590; %m
    Z_thresh_min = Z_thresh_mins(Ziter); %m
    Z_thresh_max = Z_thresh_maxs(Ziter); %m
    
    %define modis grid cell size:
    MOD_gridcell_area = 463.3127^2; %m2
    
    %screen data based on elevation:
    if Ziter < length(Z_thresh_mins)
        idx_Zthresh = find(Elevation_vec_BB>=Z_thresh_min & Elevation_vec_BB<Z_thresh_max);
    else
        idx_Zthresh = find(Elevation_vec_BB>=Z_thresh_min & Elevation_vec_BB<=Z_thresh_max);
    end
    %combined screening:
    IDX_screened = intersect(idx_Zthresh,idx_PeakSWEthresh);
    IDX_screened = intersect(IDX_screened,idx_LC);
    IDX_screened = intersect(IDX_screened,IDX_DEWS);
    
    %initialize output (WY SWEI and total summer BA):
    store_covariates = [];
    for WY = 1984:2020
        WY
        %get burn area:
        store_BurnFraction_MTBS = [];
        store_BurnFraction_MODIS = [];
        if WY < 2020
            for month = 6:9
                BA_filename = sprintf('Gridded_MTBS_%04d%02d.mat',WY,month);
                MTBS_Data = load([MTBS_BA_dir,BA_filename]);
                BurnFraction = MTBS_Data.Gridded_MTBS.BurnFraction;
                BurnFraction_vec = BurnFraction(:);
                BurnFraction_vec(idx_NAN_BAgrid) =[];
                store_BurnFraction_MTBS = [store_BurnFraction_MTBS,BurnFraction_vec(idx_BB)];
            end
        end
        
        if WY>=2001 %get MODIS BA:
            for month = 6:9
                BA_filename = sprintf('Gridded_MCD64_%04d%02d.mat',WY,month);
                MCD64_Data = load([MODIS_BA_dir,BA_filename]);
                BurnFraction = MCD64_Data.MCD64_burndate;
                idx = find(BurnFraction>0);
                BurnFraction(idx)=1;
                idx = find(BurnFraction<0);
                BurnFraction(idx) = 0;
                BurnFraction_vec = BurnFraction(:);
                BurnFraction_vec(idx_NAN_BAgrid) =[];
                store_BurnFraction_MODIS = [store_BurnFraction_MODIS,BurnFraction_vec(idx_BB)];
            end
        end
        
        %get total MTBS BA:
        if WY < 2020
            Total_fraction = nanmax(store_BurnFraction_MTBS')';
            Total_fraction_screened = Total_fraction(IDX_screened);
            Total_BA = nansum(Total_fraction_screened.*MOD_gridcell_area);
        else
            Total_BA = NaN;
        end
        
        %get total MODIS BA:
        if WY >= 2001
            Total_fraction_MODIS = nanmax(store_BurnFraction_MODIS')';
            Total_fraction_screened_MODIS = Total_fraction_MODIS(IDX_screened);
            Total_BA_MODIS = nansum(Total_fraction_screened_MODIS.*MOD_gridcell_area);
        else
            Total_BA_MODIS = NaN;
        end
        
        %Get VPD: (Get_Annual_PRISM_VPD.m)
        VPD_filename = sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/PRISM_Annual/Annual_PRISM_VPD_WY%04d.mat',WY);
        VPD_Data = load(VPD_filename);
        
        Spring_VPD = VPD_Data.Data.VPD_spring;
        Winter_VPD = VPD_Data.Data.VPD_winter;
        Summer_VPD = VPD_Data.Data.VPD_summer;
        
        Spring_VPD = Spring_VPD(:);
        Winter_VPD = Winter_VPD(:);
        Summer_VPD = Summer_VPD(:);
        
        %perform NN matching:
        Spring_VPD = Spring_VPD(IDX_NN_PRISM_BB);
        Winter_VPD = Winter_VPD(IDX_NN_PRISM_BB);
        Summer_VPD = Summer_VPD(IDX_NN_PRISM_BB);
        
        %screen:
        Spring_VPD = Spring_VPD(IDX_screened);
        Winter_VPD = Winter_VPD(IDX_screened);
        Summer_VPD = Summer_VPD(IDX_screened);
        
        %get SWEI - peak SWE:
        current_SWEI = squeeze(store_SWEI(:,WY-1983));
        current_SWEI_screened = current_SWEI(IDX_screened);
        SWEI_WY_mean = nanmean(current_SWEI_screened);
        
        %get SWEI - spring SWE:
        current_SWEI_spring = squeeze(store_SWEI_spring(:,WY-1983));
        current_SWEI_spring_screened = current_SWEI_spring(IDX_screened);
        springSWEI_WY_mean = nanmean(current_SWEI_spring_screened);
        
        %Get precipitation and temperature data (outputs from Get_Annual_PRISM_data.m):
        %these data are used to categorize snow droughts as "dry" or "hot"
        PRISM_filename = sprintf('Annual_PRISM_met_WY%04d.mat',WY);
        PRSIM_data = load([PRISM_dir,PRISM_filename]);
        
        Winter_Precip = PRSIM_data.Data.total_prcp_winter;
        Winter_Precip = Winter_Precip(:);
        Winter_Precip = Winter_Precip(IDX_NN_PRISM_BB);
        
        Spring_Precip = PRSIM_data.Data.total_prcp_spring;
        Spring_Precip = Spring_Precip(:);
        Spring_Precip = Spring_Precip(IDX_NN_PRISM_BB);
        
        Summer_Precip = PRSIM_data.Data.total_prcp_summer;
        Summer_Precip = Summer_Precip(:);
        Summer_Precip = Summer_Precip(IDX_NN_PRISM_BB);
        
        Winter_Temp = PRSIM_data.Data.mean_temp_winter;
        Winter_Temp = Winter_Temp(:);
        Winter_Temp = Winter_Temp(IDX_NN_PRISM_BB);
        
        Spring_Temp = PRSIM_data.Data.mean_temp_spring;
        Spring_Temp = Spring_Temp(:);
        Spring_Temp = Spring_Temp(IDX_NN_PRISM_BB);
        
        Summer_Temp = PRSIM_data.Data.mean_temp_summer;
        Summer_Temp = Summer_Temp(:);
        Summer_Temp = Summer_Temp(IDX_NN_PRISM_BB);
        
        Winter_Precip_screened = Winter_Precip(IDX_screened);
        Spring_Precip_screened = Spring_Precip(IDX_screened);
        Summer_Precip_screened = Summer_Precip(IDX_screened);
        
        Winter_temp_screened = Winter_Temp(IDX_screened);
        Spring_Temp_screened = Spring_Temp(IDX_screened);
        Summer_Temp_screened = Summer_Temp(IDX_screened);
        %Get PDSI:
        PDSI_filename = sprintf('Annual_PRISM_PDSI_WY%04d.mat',WY);
        PDSI_data = load([PRISM_dir,PDSI_filename]);
        
        Winter_PDSI = PDSI_data.Data.pdsi_winter;
        Winter_PDSI = Winter_PDSI(IDX_NN_PRISM_BB);
        
        Spring_PDSI = PDSI_data.Data.pdsi_spring;
        Spring_PDSI = Spring_PDSI(IDX_NN_PRISM_BB);
        
        Summer_PDSI = PDSI_data.Data.pdsi_summer;
        Summer_PDSI = Summer_PDSI(IDX_NN_PRISM_BB);
        
        Winter_PDSI_screened = Winter_PDSI(IDX_screened);
        Spring_PDSI_screened = Spring_PDSI(IDX_screened);
        Summer_PDSI_screened = Summer_PDSI(IDX_screened);
        
        %calculate total winter and spring area effected by drought:
        idx_pdsi_drought_winter = find(Winter_PDSI_screened < -0.5); %based on Alley 1983 & ul drought monitor (https://droughtmonitor.unl.edu/About/AbouttheData/DroughtClassification.aspx)
        idx_pdsi_drought_spring = find(Spring_PDSI_screened < -0.5); %based on ul drought monitor
        idx_spring_SWEI_drought = find(current_SWEI_spring_screened < -0.5);
        idx_peak_SWEI_drought = find(current_SWEI_screened < -0.5);
        idx_pdsi_summer_drought = find(Summer_PDSI_screened < -0.5);
        
        SWEI_winter_drought_area = length(idx_peak_SWEI_drought);
        SWEI_spring_drought_area = length(idx_spring_SWEI_drought);
        PDSI_winter_drought_area = length(idx_pdsi_drought_winter);
        PDSI_spring_drought_area = length(idx_pdsi_drought_spring);
        PDSI_summer_drought_area = length(idx_pdsi_summer_drought);
        
        winter_drought_area_total = length(unique([idx_peak_SWEI_drought;idx_pdsi_drought_winter']));
        spring_drought_area_total = length(unique([idx_spring_SWEI_drought;idx_pdsi_drought_spring']));
        summer_drought_area = PDSI_summer_drought_area;
        
        %Get PET and ET NLDAS data:
        ET_filename = sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NLDAS_ET_PET/Annual/Annual_ET_PET_WY%04d.mat',WY);
        ET_Data = load(ET_filename);
        
        Spring_ET = ET_Data.Data.ET_spring;
        Spring_PET = ET_Data.Data.PET_spring;
        Spring_PETminusET = ET_Data.Data.PETminusET_spring;
        
        Winter_ET = ET_Data.Data.ET_winter;
        Winter_PET = ET_Data.Data.PET_winter;
        Winter_PETminusET = ET_Data.Data.PETminusET_winter;
        
        summer_ET = ET_Data.Data.ET_summer;
        summer_PET = ET_Data.Data.PET_summer;
        summer_PETminusET = ET_Data.Data.PETminusET_summer;
        
        Spring_ET = Spring_ET(IDX_NN_NLDAS_BB);
        Spring_PET = Spring_PET(IDX_NN_NLDAS_BB);
        Spring_PETminusET = Spring_PETminusET(IDX_NN_NLDAS_BB);
        Winter_ET = Winter_ET(IDX_NN_NLDAS_BB);
        Winter_PET = Winter_PET(IDX_NN_NLDAS_BB);
        Winter_PETminusET = Winter_PETminusET(IDX_NN_NLDAS_BB);
        summer_ET = summer_ET(IDX_NN_NLDAS_BB);
        summer_PET = summer_PET(IDX_NN_NLDAS_BB);
        summer_PETminusET = summer_PETminusET(IDX_NN_NLDAS_BB);
        
        Spring_ET=Spring_ET(IDX_screened);
        Spring_PET=Spring_PET(IDX_screened);
        Spring_PETminusET=Spring_PETminusET(IDX_screened);
        Winter_ET=Winter_ET(IDX_screened);
        Winter_PET=Winter_PET(IDX_screened);
        Winter_PETminusET=Winter_PETminusET(IDX_screened);
        summer_ET=summer_ET(IDX_screened);
        summer_PET=summer_PET(IDX_screened);
        summer_PETminusET=summer_PETminusET(IDX_screened);
        
        Spring_ET = nanmean(Spring_ET);
        Spring_PET = nanmean(Spring_PET);
        Spring_PETminusET = nanmean(Spring_PETminusET);
        Winter_ET = nanmean(Winter_ET);
        Winter_PET = nanmean(Winter_PET);
        Winter_PETminusET = nanmean(Winter_PETminusET);
        summer_ET = nanmean(summer_ET);
        summer_PET = nanmean(summer_PET);
        summer_PETminusET = nanmean(summer_PETminusET);
        
        %get Summer SM (Get_Summer_NLDAS_SM.m):
        infilename_sm = sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NLDAS_Summer_SM/NLDAS_summer_SM_WY%04d.mat',WY);
        Summer_SM = load(infilename_sm);
        Summer_SM = Summer_SM.Data.summer_SM;
        Summer_SM = Summer_SM(IDX_NN_NLDAS_BB);
        Summer_SM=Summer_SM(IDX_screened);
        Summer_SM = nanmean(Summer_SM);
        
        
        
        %store outputs:
        store_covariates = [store_covariates;WY,SWEI_WY_mean,Total_BA,nanmean(Winter_Precip_screened),nanmean(Winter_temp_screened),nanmean(Winter_PDSI_screened),nanmean(Spring_Precip_screened),nanmean(Spring_Temp_screened),nanmean(Spring_PDSI_screened),springSWEI_WY_mean,nanmean(Spring_VPD),nanmean(Winter_VPD),Total_BA_MODIS,Spring_ET,Spring_PET,Spring_PETminusET,Winter_ET,Winter_PET,Winter_PETminusET,nanmean(Summer_VPD),nanmean(Summer_PDSI_screened),nanmean(Summer_Precip_screened),nanmean(Summer_Temp_screened),summer_ET,summer_PET,summer_PETminusET,Summer_SM,winter_drought_area_total,spring_drought_area_total,summer_drought_area,SWEI_winter_drought_area,SWEI_spring_drought_area,PDSI_winter_drought_area,PDSI_spring_drought_area];
    end
    %this data has the annual time series across all zins:
    %     outfilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/BA_covariates_Obs_forecast_All_Zbins.mat',Ziter);
    %this is data for each Zbin
    outfilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/BA_covariates_Obs_forecast_Ziter%d_higherPeakSWE.mat',Ziter);
    save(outfilename,'store_covariates');
    store_covariates = load(outfilename);
    store_covariates = store_covariates.store_covariates;
     %export variables:
    WYs = store_covariates(:,1);
    SWEIs = store_covariates(:,2);
    BAs = store_covariates(:,3);
    Winter_PRCP = store_covariates(:,4);
    Winter_temp = store_covariates(:,5);
    Winter_PDSI = store_covariates(:,6);
    Spring_PRCP = store_covariates(:,7);
    Spring_TMP = store_covariates(:,8);
    Spring_PDSI = store_covariates(:,9);
    Spring_SWEI = store_covariates(:,10);
    Spring_VPD = store_covariates(:,11);
    Winter_VPD= store_covariates(:,12);
    MODIS_BA = store_covariates(:,13);
    Spring_ET = store_covariates(:,14);
    Spring_PET = store_covariates(:,15);
    Spring_PETminusET = store_covariates(:,16);
    Winter_ET = store_covariates(:,17);
    Winter_PET = store_covariates(:,18);
    Winter_PETminusET = store_covariates(:,19);
    %get summer variables as well for further analysis:
    Summer_VPD = store_covariates(:,20);
    Summer_PDSI = store_covariates(:,21);
    Summer_PRCP = store_covariates(:,22);
    Summer_Temp = store_covariates(:,23);
    Summer_ET = store_covariates(:,24);
    Summer_PET= store_covariates(:,25);
    Summer_PETminusET = store_covariates(:,26);
    Summer_SM = store_covariates(:,27);
    Winter_Drought_area = store_covariates(:,28);
    Spring_Drought_area = store_covariates(:,29);
    Summer_Drought_area = store_covariates(:,30);
    Winter_SWEI_Drought_area = store_covariates(:,31);
    Spring_SWEI_Drought_area = store_covariates(:,32);
    Winter_PDSI_Drought_area = store_covariates(:,33);
    Spring_PDSI_Drought_area = store_covariates(:,34);
    
    Annual_Fire_Climate_Variables = [WYs,BAs,SWEIs,Winter_PRCP,Winter_temp,Winter_PDSI,Spring_PRCP,Spring_TMP,Spring_PDSI,Spring_SWEI,Spring_VPD,Winter_VPD,MODIS_BA,Spring_ET,Spring_PET,Spring_PETminusET,Winter_ET,Winter_PET,Winter_PETminusET,Summer_VPD,Summer_PDSI,Summer_PRCP,Summer_Temp,Summer_ET,Summer_PET,Summer_PETminusET,Summer_SM,Winter_Drought_area,Spring_Drought_area,Summer_Drought_area,Winter_SWEI_Drought_area,Spring_SWEI_Drought_area,Winter_PDSI_Drought_area,Spring_PDSI_Drought_area];
    outfilename_csv=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter%d_lowerPeakSWE.csv',Ziter);
    dlmwrite(outfilename_csv,Annual_Fire_Climate_Variables,'delimiter',',','precision','%.5f');
end
