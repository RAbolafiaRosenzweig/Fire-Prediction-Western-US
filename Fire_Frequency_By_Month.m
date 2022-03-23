clc;clear all;close all;

% % %% define bounding box for area of analysis: western CONUS:
% % latlim = [32 51];
% % lonlim = [-125 -104];
% % 
% % %% Define thresholds for study domain
% % Z_thresh_min = 0;
% % Z_thresh_max = 3300;
% % % % %by zbins:
% % % % Z_thresh_max1 = 1100;
% % % % Z_thresh_max2 = 2200;
% % % % Z_thresh_max3 = 3300;
% % 
% % %arbitrary
% % peakLAI_thresh = 0.0;
% % 
% % %arbitrary
% % peakSWE_thresh = 100; %mm
% % 
% % %% Define BA grid:
% % BA_lat = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200708.nc','XLAT_M');
% % BA_latvec = BA_lat(:);
% % BA_lon = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200708.nc','XLONG_M');
% % BA_lonvec = BA_lon(:);
% % %remove NaN:
% % idx_NAN_BAgrid = find(isnan(BA_latvec));
% % BA_lonvec(idx_NAN_BAgrid)=[];
% % BA_latvec(idx_NAN_BAgrid)=[];
% % %trim to bounding box
% % idx_BB = find(BA_latvec>=latlim(1) & BA_latvec<=latlim(2) & BA_lonvec>=lonlim(1) & BA_lonvec<=lonlim(2));
% % 
% % %% trim by  DEWS regions for domain separation:
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
% % %% define peak SWE directory (Get_UA_Annual_Peak_SWE.m):
% % UA_peakSWE_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_PeakSWE/';
% % %get nearestneighbor index to match with BA grid (NN_UA_and_PRISM_to_BAgrid.m)
% % IDX_NN_UAswe = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_UAswe.mat');
% % IDX_NN_UAswe = IDX_NN_UAswe.IDX_NN_UAswe;
% % idx_nan = find(isnan(IDX_NN_UAswe));
% % IDX_NN_UAswe(idx_nan)= [];
% % 
% % %% calculate SWEI drought index values from peak SWE (Hunning and AghaKouchak, 2020)
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
% % % % %rank the peak SWE values across the years for each pixel:
% % % % [~,ii] = sort(store_peak_swe_BB,2);
% % % % [~,i] = sort(ii,2);
% % % % N=size(i,2);
% % % % p = (i - 0.44)./(N+0.12);
% % % %
% % % % %calculate nonparametric SWEI by transforming the empirical probability, p,
% % % % %to the standard normal distribution:
% % % % store_SWEI = NaN(length(p),37);
% % % % for i = 1:length(p)
% % % %     current_pixel_p = p(i,:);
% % % %     store_SWEI(i,:) = norminv(current_pixel_p);
% % % % end
% % % % %export:
% % % % save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_PeakSWE_SWEI/UA_SWEI_Hindcast.mat','store_SWEI', '-v7.3');
% % store_SWEI = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_PeakSWE_SWEI/UA_SWEI_Hindcast.mat','store_SWEI');
% % store_SWEI = store_SWEI.store_SWEI;
% % 
% % %% calculate SWEI drought index values from spring SWE (Hunning and AghaKouchak, 2020)
% % %from Get_UA_Total_Spring_SWE.m
% % % % store_spring_swe=[];
% % % % for year = 1984:2020
% % % %     infilename = sprintf('UA_SpringSWE_%04d.mat',year);
% % % %     UA_PeakSWE = load([UA_peakSWE_dir,infilename]);
% % % %     UA_PeakSWE = UA_PeakSWE.Data.current_spring_SWE;
% % % %     idx = find(UA_PeakSWE == 0);
% % % %     UA_PeakSWE(idx) = rand(length(idx),1)*10^-10;
% % % %     store_spring_swe = [store_spring_swe,UA_PeakSWE(:)];
% % % % end
% % % % store_spring_swe = store_spring_swe(IDX_NN_UAswe,:);
% % % % store_spring_swe = store_spring_swe(idx_BB,:);
% % 
% % % % % % %rank the peak SWE values across the years for each pixel:
% % % % [~,ii] = sort(store_spring_swe,2);
% % % % [~,i] = sort(ii,2);
% % % % N=size(i,2);
% % % % p = (i - 0.44)./(N+0.12);
% % % %
% % % % %calculate nonparametric SWEI by transforming the empirical probability, p,
% % % % %to the standard normal distribution:
% % % % store_SWEI = NaN(length(p),37);
% % % % for i = 1:length(p)
% % % %     current_pixel_p = p(i,:);
% % % %     store_SWEI(i,:) = norminv(current_pixel_p);
% % % % end
% % %export:
% % % % save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_SpringSWE_SWEI/UA_Spring_SWEI_Hindcast.mat','store_SWEI', '-v7.3');
% % store_SWEI_spring = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_SpringSWE_SWEI/UA_Spring_SWEI_Hindcast.mat','store_SWEI');
% % store_SWEI_spring = store_SWEI_spring.store_SWEI;
% % 
% % %% Load in elevation grid: (from Map_elevation_to_BA_product.m)
% % Elevation_coarse=load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Elevation_coarse.mat');
% % Elevation_coarse = Elevation_coarse.Elevation_coarse;
% % IDX_NN_Z = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_Elevation.mat');
% % IDX_NN_Z = IDX_NN_Z.IDX_NN_Elevation;
% % idx_nan = find(isnan(IDX_NN_Z));
% % IDX_NN_Z(idx_nan)=[];
% % Elevation_vec = Elevation_coarse(IDX_NN_Z);
% % Elevation_vec_BB = Elevation_vec(idx_BB);
% % 
% % %% define temp and precip directory (Get_Annual_PRISM_data.m)
% % PRISM_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/PRISM_Annual/';
% % %get nearestneighbor index to match with BA grid (NN_UA_and_PRISM_to_BAgrid.m)
% % IDX_NN_PRISM = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_PRISM.mat');
% % IDX_NN_PRISM = IDX_NN_PRISM.IDX_NN_PRISM;
% % idx_nan = find(isnan(IDX_NN_PRISM));
% % IDX_NN_PRISM(idx_nan)=[];
% % IDX_NN_PRISM_BB = IDX_NN_PRISM(idx_BB);
% % 
% % %get nearestneighbor index to match with BA grid (NN_NLDAS_to_BAgrid.m)
% % IDX_NN_NLDAS = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_NLDAS.mat');
% % IDX_NN_NLDAS = IDX_NN_NLDAS.IDX_NN_NLDAS;
% % idx_nan = find(isnan(IDX_NN_NLDAS));
% % IDX_NN_NLDAS(idx_nan)= [];
% % IDX_NN_NLDAS_BB = IDX_NN_NLDAS(idx_BB);
% % 
% % %% load in MODIS landcover type:
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
% % 
% % %define forest LC:
% % Forest_LC_index_BB =sort([idx_ENF;idx_EBF;idx_DNF;idx_DBF;idx_MF]);
% % %define grass LC:
% % Grass_LC_index_BB =[idx_grass];
% % %define CS LC:
% % Shrub_LC_index_BB =[idx_CS;idx_OS];
% % %define Sav LC:
% % Sav_LC_index_BB =sort([idx_WSav;idx_Sav]);
% % 
% % %% define merged BA product directory (/Users/abolafia/Drought_Fire_Snow/tools/Merged_BA)
% % MTBS_BA_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Monthly_gridded_MTBS/';
% % MODIS_BA_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Monthly_gridded_MCD64/';
% % 
% % %% Compare summer (June-September) burn area with SWEI:
% % %screen based on peak SWE:
% % idx_PeakSWEthresh = find(max(store_peak_swe_BB(:,1:end)') >= peakSWE_thresh);
% % %screen based on landcover type:
% % idx_LC = Total_LC_index_BB;
% % %define modis grid cell size:
% % MOD_gridcell_area = 463.3127^2; %m2
% % 
% % %screen data based on elevation:
% % idx_Zthresh = find(Elevation_vec_BB>=Z_thresh_min & Elevation_vec_BB<=Z_thresh_max);
% % 
% % %combined screening:
% % IDX_screened = intersect(idx_Zthresh,idx_PeakSWEthresh);
% % IDX_screened = intersect(IDX_screened,idx_LC);
% % IDX_screened = intersect(IDX_screened,IDX_DEWS);
% % 
% % %consider total w US area:
% % ncoords = length(IDX_screened);
% % Z_thresh_min=-10000;
% % Z_thresh_max=10000;
% % idx_Zthresh = find(Elevation_vec_BB>=Z_thresh_min & Elevation_vec_BB<=Z_thresh_max);
% % IDX_screened_total_CONUS = intersect(IDX_DEWS,idx_Zthresh);
% % 
% % sprintf('domain is %.2f percent of total w. US',length(IDX_screened)/length(IDX_screened_total_CONUS)*100)
% % 
% % %consider screened area without elevation threshold:
% % IDX_screened_noZ= intersect(idx_LC,idx_PeakSWEthresh);
% % IDX_screened_noZ= intersect(IDX_screened_noZ,IDX_DEWS);
% % 
% % %consider zbins:
% % idx_Zthresh1 = find(Elevation_vec_BB>=0 & Elevation_vec_BB<1100);
% % idx_Zthresh2 = find(Elevation_vec_BB>=1100 & Elevation_vec_BB<2200);
% % idx_Zthresh3 = find(Elevation_vec_BB>=2200 & Elevation_vec_BB<=3300);
% % idx_Zthresh4 = find(Elevation_vec_BB>3300);
% % 
% % IDX_screened1 = intersect(idx_Zthresh1,idx_PeakSWEthresh);
% % IDX_screened1 = intersect(IDX_screened1,idx_LC);
% % IDX_screened1 = intersect(IDX_screened1,IDX_DEWS);
% % 
% % IDX_screened2 = intersect(idx_Zthresh2,idx_PeakSWEthresh);
% % IDX_screened2 = intersect(IDX_screened2,idx_LC);
% % IDX_screened2 = intersect(IDX_screened2,IDX_DEWS);
% % 
% % IDX_screened3 = intersect(idx_Zthresh3,idx_PeakSWEthresh);
% % IDX_screened3 = intersect(IDX_screened3,idx_LC);
% % IDX_screened3 = intersect(IDX_screened3,IDX_DEWS);
% % 
% % IDX_screened4 = intersect(idx_Zthresh4,idx_PeakSWEthresh);
% % IDX_screened4 = intersect(IDX_screened4,idx_LC);
% % IDX_screened4 = intersect(IDX_screened4,IDX_DEWS);
% % 
% % %bin by LC class:
% % %forest:
% % IDX_screened_forest = intersect(idx_Zthresh,idx_PeakSWEthresh);
% % IDX_screened_forest = intersect(IDX_screened_forest,Forest_LC_index_BB);
% % IDX_screened_forest = intersect(IDX_screened_forest,IDX_DEWS);
% % %grass:
% % IDX_screened_grass = intersect(idx_Zthresh,idx_PeakSWEthresh);
% % IDX_screened_grass = intersect(IDX_screened_grass,Grass_LC_index_BB);
% % IDX_screened_grass = intersect(IDX_screened_grass,IDX_DEWS);
% % %CS:
% % IDX_screened_Shrub = intersect(idx_Zthresh,idx_PeakSWEthresh);
% % IDX_screened_Shrub = intersect(IDX_screened_Shrub,Shrub_LC_index_BB);
% % IDX_screened_Shrub = intersect(IDX_screened_Shrub,IDX_DEWS);
% % %Sav:
% % IDX_screened_Sav = intersect(idx_Zthresh,idx_PeakSWEthresh);
% % IDX_screened_Sav = intersect(IDX_screened_Sav,Sav_LC_index_BB);
% % IDX_screened_Sav = intersect(IDX_screened_Sav,IDX_DEWS);
% % % % 
% % %initialize output (WY SWEI and total summer BA):
% % store_BurnArea = [];
% % store_BurnArea_MODIS=[];
% % WY_months=[10 11 12,1:9];
% % for WY = 1984:2020
% %     WY
% %     for month = 1:12
% %         WY_month = WY_months(month);
% %         if WY < 2020
% %             if WY_month>=10
% %                 BA_filename = sprintf('Gridded_MTBS_%04d%02d.mat',WY-1,WY_month);
% %             else
% %                 BA_filename = sprintf('Gridded_MTBS_%04d%02d.mat',WY,WY_month);
% %             end
% %             if exist([MTBS_BA_dir,BA_filename],'file')>0
% %                 MTBS_Data = load([MTBS_BA_dir,BA_filename]);
% %                 BurnFraction = MTBS_Data.Gridded_MTBS.BurnFraction;
% %                 BurnFraction_vec = BurnFraction(:);
% %                 BurnFraction_vec(idx_NAN_BAgrid) =[];
% %                 BurnFraction_vec = BurnFraction_vec(idx_BB);
% %                 BurnFraction_vec_screened = BurnFraction_vec(IDX_screened);
% %                 Total_BA = nansum(BurnFraction_vec_screened.*MOD_gridcell_area);
% %                 %record results for total w CONUS as well:
% %                 BurnFraction_vec_wCONUS = BurnFraction_vec(IDX_screened_total_CONUS);
% %                 Total_BA_wCONUS = nansum(BurnFraction_vec_wCONUS.*MOD_gridcell_area);
% %                 %record results for Z bins as well:
% %                 %zbin1
% %                 BurnFraction_vec_Z1 = BurnFraction_vec(IDX_screened1);
% %                 Total_BA_Z1 = nansum(BurnFraction_vec_Z1.*MOD_gridcell_area);
% %                 %zbin2
% %                 BurnFraction_vec_Z2 = BurnFraction_vec(IDX_screened2);
% %                 Total_BA_Z2 = nansum(BurnFraction_vec_Z2.*MOD_gridcell_area);
% %                 %zbin3
% %                 BurnFraction_vec_Z3 = BurnFraction_vec(IDX_screened3);
% %                 Total_BA_Z3 = nansum(BurnFraction_vec_Z3.*MOD_gridcell_area);
% %                 %zbin > zbin3
% %                 BurnFraction_vec_Z4 = BurnFraction_vec(IDX_screened4);
% %                 Total_BA_Z4 = nansum(BurnFraction_vec_Z4.*MOD_gridcell_area);
% %                 
% %                 %recrod results by LC:
% %                 %forest
% %                 BurnFraction_vec_forest = BurnFraction_vec(IDX_screened_forest);
% %                 Total_BA_forest = nansum(BurnFraction_vec_forest.*MOD_gridcell_area);
% %                 %grass
% %                 BurnFraction_vec_grass = BurnFraction_vec(IDX_screened_grass);
% %                 Total_BA_grass = nansum(BurnFraction_vec_grass.*MOD_gridcell_area);
% %                 %CS
% %                 BurnFraction_vec_Shrub = BurnFraction_vec(IDX_screened_Shrub);
% %                 Total_BA_Shrub = nansum(BurnFraction_vec_Shrub.*MOD_gridcell_area);
% %                 %Sav LC
% %                 BurnFraction_vec_Sav = BurnFraction_vec(IDX_screened_Sav);
% %                 Total_BA_Sav  = nansum(BurnFraction_vec_Sav.*MOD_gridcell_area);
% %                 
% %                 %store:
% %                 store_BurnArea = [store_BurnArea;WY,WY_month,Total_BA,Total_BA_wCONUS,Total_BA_Z1,Total_BA_Z2,Total_BA_Z3,Total_BA_Z4,Total_BA_forest,Total_BA_grass,Total_BA_Shrub,Total_BA_Sav];
% %             else
% %                 store_BurnArea = [store_BurnArea;WY,WY_month,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
% %             end
% %         end
% %         if WY>=2002
% %             if WY_month>=10
% %                 BA_filename = sprintf('Gridded_MCD64_%04d%02d.mat',WY-1,WY_month);
% %             else
% %                 BA_filename = sprintf('Gridded_MCD64_%04d%02d.mat',WY,WY_month);
% %             end
% %             MCD64_Data = load([MODIS_BA_dir,BA_filename]);
% %             BurnFraction = MCD64_Data.MCD64_burndate;
% %             idx = find(BurnFraction>0);
% %             BurnFraction(idx)=1;
% %             idx = find(BurnFraction<0);
% %             BurnFraction(idx) = 0;
% %             BurnFraction_vec = BurnFraction(:);
% %             BurnFraction_vec(idx_NAN_BAgrid) =[];
% %             BurnFraction_vec = BurnFraction_vec(idx_BB);
% %             BurnFraction_vec_screened = BurnFraction_vec(IDX_screened);
% %             Total_BA = nansum(BurnFraction_vec_screened.*MOD_gridcell_area);
% %             %record results for total w CONUS as well:
% %             BurnFraction_vec_wCONUS = BurnFraction_vec(IDX_screened_total_CONUS);
% %             Total_BA_wCONUS = nansum(BurnFraction_vec_wCONUS.*MOD_gridcell_area);
% %             %record results for Z bins as well:
% %             %zbin1
% %             BurnFraction_vec_Z1 = BurnFraction_vec(IDX_screened1);
% %             Total_BA_Z1 = nansum(BurnFraction_vec_Z1.*MOD_gridcell_area);
% %             %zbin2
% %             BurnFraction_vec_Z2 = BurnFraction_vec(IDX_screened2);
% %             Total_BA_Z2 = nansum(BurnFraction_vec_Z2.*MOD_gridcell_area);
% %             %zbin3
% %             BurnFraction_vec_Z3 = BurnFraction_vec(IDX_screened3);
% %             Total_BA_Z3 = nansum(BurnFraction_vec_Z3.*MOD_gridcell_area);
% %             %zbin > zbin3
% %             BurnFraction_vec_Z4 = BurnFraction_vec(IDX_screened4);
% %             Total_BA_Z4 = nansum(BurnFraction_vec_Z4.*MOD_gridcell_area);
% %             
% %             %recrod results by LC:
% %                 %forest
% %                 BurnFraction_vec_forest = BurnFraction_vec(IDX_screened_forest);
% %                 Total_BA_forest = nansum(BurnFraction_vec_forest.*MOD_gridcell_area);
% %                 %grass
% %                 BurnFraction_vec_grass = BurnFraction_vec(IDX_screened_grass);
% %                 Total_BA_grass = nansum(BurnFraction_vec_grass.*MOD_gridcell_area);
% %                 %CS
% %                 BurnFraction_vec_Shrub = BurnFraction_vec(IDX_screened_Shrub);
% %                 Total_BA_Shrub = nansum(BurnFraction_vec_Shrub.*MOD_gridcell_area);
% %                 %Sav LC
% %                 BurnFraction_vec_Sav = BurnFraction_vec(IDX_screened_Sav);
% %                 Total_BA_Sav  = nansum(BurnFraction_vec_Sav.*MOD_gridcell_area);
% %                 
% %                 %store:
% %                 store_BurnArea_MODIS = [store_BurnArea_MODIS;WY,WY_month,Total_BA,Total_BA_wCONUS,Total_BA_Z1,Total_BA_Z2,Total_BA_Z3,Total_BA_Z4,Total_BA_forest,Total_BA_grass,Total_BA_Shrub,Total_BA_Sav];
% %         end
% %     end
% % end
% % %% relate MODIS to MTBS BA - for screened area:
% % %no data for Jan 2019:
% % idx_nan = find(isnan(store_BurnArea));
% % store_BurnArea(idx_nan) = 0;
% % 
% % MODIS_BA = store_BurnArea_MODIS(1:end-12,3);
% % MTBS_BA = store_BurnArea(217:end,3);
% % p = polyfit(MODIS_BA,MTBS_BA,1);
% % MODIS_to_MTBS_func = @(MODIS_BAs)  MODIS_BAs*p(1) + p(2);
% % MODIS_2020 =store_BurnArea_MODIS(end-11:end,3);
% % MTBS_2020 = MODIS_to_MTBS_func(MODIS_2020);
% % %% relate MODIS to MTBS BA - for total W CONUS:
% % MODIS_BA = store_BurnArea_MODIS(1:end-12,4);
% % MTBS_BA = store_BurnArea(217:end,4);
% % p = polyfit(MODIS_BA,MTBS_BA,1);
% % MODIS_to_MTBS_func = @(MODIS_BAs)  MODIS_BAs*p(1) + p(2);
% % MODIS_2020 =store_BurnArea_MODIS(end-11:end,4);
% % MTBS_2020_wCONUS = MODIS_to_MTBS_func(MODIS_2020);
% % %% relate MODIS to MTBS BA - for Z1:
% % MODIS_BA = store_BurnArea_MODIS(1:end-12,5);
% % MTBS_BA = store_BurnArea(217:end,5);
% % p = polyfit(MODIS_BA,MTBS_BA,1);
% % MODIS_to_MTBS_func = @(MODIS_BAs)  MODIS_BAs*p(1) + p(2);
% % MODIS_2020 =store_BurnArea_MODIS(end-11:end,5);
% % MTBS_2020_Z1 = MODIS_to_MTBS_func(MODIS_2020);
% % %% relate MODIS to MTBS BA - for Z2:
% % MODIS_BA = store_BurnArea_MODIS(1:end-12,6);
% % MTBS_BA = store_BurnArea(217:end,6);
% % p = polyfit(MODIS_BA,MTBS_BA,1);
% % MODIS_to_MTBS_func = @(MODIS_BAs)  MODIS_BAs*p(1) + p(2);
% % MODIS_2020 =store_BurnArea_MODIS(end-11:end,6);
% % MTBS_2020_Z2 = MODIS_to_MTBS_func(MODIS_2020);
% % %% relate MODIS to MTBS BA - for Z3:
% % MODIS_BA = store_BurnArea_MODIS(1:end-12,7);
% % MTBS_BA = store_BurnArea(217:end,7);
% % p = polyfit(MODIS_BA,MTBS_BA,1);
% % MODIS_to_MTBS_func = @(MODIS_BAs)  MODIS_BAs*p(1) + p(2);
% % MODIS_2020 =store_BurnArea_MODIS(end-11:end,7);
% % MTBS_2020_Z3 = MODIS_to_MTBS_func(MODIS_2020);
% % %% relate MODIS to MTBS BA - for Z4:
% % MODIS_BA = store_BurnArea_MODIS(1:end-12,8);
% % MTBS_BA = store_BurnArea(217:end,8);
% % p = polyfit(MODIS_BA,MTBS_BA,1);
% % MODIS_to_MTBS_func = @(MODIS_BAs)  MODIS_BAs*p(1) + p(2);
% % MODIS_2020 =store_BurnArea_MODIS(end-11:end,8);
% % MTBS_2020_Z4 = MODIS_to_MTBS_func(MODIS_2020);
% % 
% % %% relate MODIS to MTBS BA - for forest:
% % MODIS_BA = store_BurnArea_MODIS(1:end-12,9);
% % MTBS_BA = store_BurnArea(217:end,9);
% % p = polyfit(MODIS_BA,MTBS_BA,1);
% % MODIS_to_MTBS_func = @(MODIS_BAs)  MODIS_BAs*p(1) + p(2);
% % MODIS_2020 =store_BurnArea_MODIS(end-11:end,9);
% % MTBS_2020_forest = MODIS_to_MTBS_func(MODIS_2020);
% % 
% % %% relate MODIS to MTBS BA - for grass:
% % MODIS_BA = store_BurnArea_MODIS(1:end-12,10);
% % MTBS_BA = store_BurnArea(217:end,10);
% % p = polyfit(MODIS_BA,MTBS_BA,1);
% % MODIS_to_MTBS_func = @(MODIS_BAs)  MODIS_BAs*p(1) + p(2);
% % MODIS_2020 =store_BurnArea_MODIS(end-11:end,10);
% % MTBS_2020_grass = MODIS_to_MTBS_func(MODIS_2020);
% % 
% % %% relate MODIS to MTBS BA - for shrublands:
% % MODIS_BA = store_BurnArea_MODIS(1:end-12,11);
% % MTBS_BA = store_BurnArea(217:end,11);
% % p = polyfit(MODIS_BA,MTBS_BA,1);
% % MODIS_to_MTBS_func = @(MODIS_BAs)  MODIS_BAs*p(1) + p(2);
% % MODIS_2020 =store_BurnArea_MODIS(end-11:end,11);
% % MTBS_2020_Shrub = MODIS_to_MTBS_func(MODIS_2020);
% % 
% % %% relate MODIS to MTBS BA - for Sav LC:
% % MODIS_BA = store_BurnArea_MODIS(1:end-12,12);
% % MTBS_BA = store_BurnArea(217:end,12);
% % p = polyfit(MODIS_BA,MTBS_BA,1);
% % MODIS_to_MTBS_func = @(MODIS_BAs)  MODIS_BAs*p(1) + p(2);
% % MODIS_2020 =store_BurnArea_MODIS(end-11:end,12);
% % MTBS_2020_Sav_LC = MODIS_to_MTBS_func(MODIS_2020);
% % 
% % store_BurnArea_total_MTBS = [store_BurnArea;linspace(2020,2020,12)',WY_months',MTBS_2020,MTBS_2020_wCONUS,MTBS_2020_Z1,MTBS_2020_Z2,MTBS_2020_Z3,MTBS_2020_Z4,MTBS_2020_forest,MTBS_2020_grass,MTBS_2020_Shrub,MTBS_2020_Sav_LC];
% % save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/MTBS_Snowy_BA_Timeseries.mat','store_BurnArea_total_MTBS');

BA_data = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/MTBS_Snowy_BA_Timeseries.mat');
BA_data = BA_data.store_BurnArea_total_MTBS;
years = BA_data(:,1);
months = BA_data(:,2);
BA = BA_data(:,3);
BA = (BA*0.000247105) /(10^6); %m2 --> millions of acres

%% aggregate to monthly climo and plot:
store_monthly_BA=[];
for m=1:12
    idx = find(months==m);
    BA_month = BA(idx);
    store_monthly_BA=[store_monthly_BA,BA_month];
end

%% plot BA climo data:
f=figure;
monthly_months = 1:12;
hold on
b=boxplot(store_monthly_BA,monthly_months);
set(b,'linew',3)
% xlabel('Month','fontsize',34)
ylabel('Burned Area (millions of acres)','fontsize',34)
set(gca,'fontsize',34)
xlim([0.5 12.5])
ylim([0 2])
xticks(monthly_months)
f.Position=[-1919         -27        1600         793];
xticklabels({'Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sept.','Oct.','Nov.','Dec.'})
grid on
box on
%fill red fire season:
x1_fireseason = [linspace(5.5,5.5,10)];
x2_fireseason = [linspace(9.5,9.5,10)];
y_fill = [linspace(0,2,10)];

plot(x1_fireseason,y_fill,'r')
plot(x2_fireseason,y_fill,'r')

fill_y = [y_fill,fliplr(y_fill)];
fill_x = [x1_fireseason,fliplr(x2_fireseason)];
fill_shape=fill(fill_x,fill_y,'r');
fill_shape.FaceAlpha = 0.2;
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/Snowy_BurnArea_ByMonth_Boxplot.eps','epsc')

%% calculate increasing trend of BA for each month:
f=figure;
summer_BA = sum(store_monthly_BA(:,6:9)');
titles = {'Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sept.','Oct.','Nov.','Dec.'};
years = 1984:2020;
store_trends=[];
store_R=[];
for m=1:12
    current_month_BA = store_monthly_BA(:,m);
    subplot(3,4,m)
    hold on
    plot(years,current_month_BA,'-r','linewidth',2);
    %get trend line (excluding 2020):
    p = polyfit(years(1:end-1),current_month_BA(1:end-1),1);
    y=@(x) p(1)*x + p(2);
    BA_trend = y(years(1:end-1));
    plot(years(1:end-1),BA_trend,'--r','linewidth',1.5)
    TITLE = titles{m};
    title(TITLE,'fontsize',22)
    set(gca,'fontsize',22)
    grid(gca,'minor')
    grid on
    xlim([1984 2020])
    ylim([0 max(current_month_BA)+0.1*max(current_month_BA)])
    %store trend as percent increase/decrease per year for respective months
    store_trends = [store_trends;p(1)/mean(current_month_BA(1:end-1))*100];
    %store this months correlation with summer BA:
    R = corr(summer_BA',current_month_BA);
    store_R = [store_R;R];
end
f.Position = [-1909         -73        1764         923];
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/Burn_area_timeseries_by_month.eps','epsc')
pp
%% Plot correlation between screened area and total w. CONUS:
%trim data to summer months:
idx = find(BA_data(:,2)>=6 & BA_data(:,2)<=9);
BA_data_summer = BA_data(idx,:);

BA_screened = BA_data_summer(:,3);
BA_wCONUS = BA_data_summer(:,4);
BA_screened=(BA_screened*0.000247105) /(10^6); %m2 --> millions of acres
BA_wCONUS=(BA_wCONUS*0.000247105) /(10^6); %m2 --> millions of acres

%aggregate to summer:
[u,~,j] = unique(BA_data_summer(:,1));
Summer_BA_screened = accumarray(j,BA_screened,[],@nansum);
Summer_BA_wCONUS = accumarray(j,BA_wCONUS,[],@nansum);
R_screened_total = corr(Summer_BA_screened,Summer_BA_wCONUS);
R_screened_total

fig=figure;
hold on
plot(Summer_BA_screened,Summer_BA_wCONUS,'ro','markersize',12,'markerfacecolor','r')
p=polyfit(Summer_BA_screened,Summer_BA_wCONUS,1);
f=@(x) p(1)*x + p(2);
best_fit = f(Summer_BA_screened);
plot(Summer_BA_screened,best_fit,'-k','linewidth',3)
xlabel({'Summer burned area for screened domain'; '(millions of acres)'},'fontsize',24)
ylabel({'Summer burned area for total Western US'; '(millions of acres)'},'fontsize',24)
set(gca,'fontsize',24)
fig.Position = [-1358         -11         936         633];
grid on
box on
xlim([0 5])
ylim([0 7])
saveas(fig,'/Users/abolafia/Drought_Fire_Snow/Plots/Screened_TotalwUS_Scatter.eps','epsc')

%% plot pie chart of BA by Z bins:
BA_data_summer_Z1 = nansum(BA_data_summer(:,5));
BA_data_summer_Z2 = nansum(BA_data_summer(:,6));
BA_data_summer_Z3 = nansum(BA_data_summer(:,7));
BA_data_summer_Z4 = nansum(BA_data_summer(:,8));

fig=figure;
X=[BA_data_summer_Z1;BA_data_summer_Z2;BA_data_summer_Z3;BA_data_summer_Z4];
p=pie(X);
title('Burned Area by Elevation','fontsize',24)
set(gca,'fontsize',24)
set(p(2:2:end),'FontSize',24);

labels = {'0-1100m','1100-2200m','2200-3300m','Z > 3300m'};
lgd = legend(labels);
lgd.Location = 'eastoutside';
fig.Position = [-1459         -76         985         604];
set(gca,'fontsize',25)
saveas(fig,'/Users/abolafia/Drought_Fire_Snow/Plots/Pie_Chart_BA_By_Elevation.eps','epsc')

%% plot pie chart of BA by LC bins:
%this plots the total burned area including all base case domain screening,
%for each LC class:
BA_data_summer_forest = nansum(BA_data_summer(:,9));
BA_data_summer_grass = nansum(BA_data_summer(:,10));
BA_data_summer_Shrub = nansum(BA_data_summer(:,11));
BA_data_summer_Sav = nansum(BA_data_summer(:,12));

fig=figure;
X=[BA_data_summer_forest;BA_data_summer_grass;BA_data_summer_Sav;BA_data_summer_Shrub];
p=pie(X);
title('Burned Area by Landcover','fontsize',24)
set(gca,'fontsize',24)
set(p(2:2:end),'FontSize',24);
p(1).FaceColor = [0, 0.5, 0];
p(3).FaceColor = [0.4940, 0.1840, 0.5560];
p(5).FaceColor = [0.75, 0.75, 0];
p(7).FaceColor = [0.8500, 0.3250, 0.0980];

labels = {'Forest','Grass','Savanna','Shrubland'};
lgd = legend(labels);
lgd.Location = 'westoutside';
fig.Position = [-1459         -76         985         604];
set(gca,'fontsize',25)
saveas(fig,'/Users/abolafia/Drought_Fire_Snow/Plots/Pie_Chart_BA_By_LCtype.eps','epsc')

%% plot pie chart of domain by LC:
LC_types_screened = MODIS_LC_BB(IDX_screened);
idx_ENF = find(LC_types_screened == 1);
idx_EBF = find(LC_types_screened == 2);
idx_DNF = find(LC_types_screened == 3);
idx_DBF = find(LC_types_screened == 4);
idx_MF = find(LC_types_screened == 5);
idx_WSav = find(LC_types_screened == 8);
idx_Sav = find(LC_types_screened == 9);
idx_grass = find(LC_types_screened == 10);

L_ENF = length(idx_ENF);
L_EBF = length(idx_EBF);
L_DNF = length(idx_DNF);
L_DBF = length(idx_DBF);
L_MF = length(idx_MF);
L_forest = L_ENF+L_EBF+L_DNF+L_DBF+L_MF;
L_Savanna = length([idx_WSav;idx_Sav]);
L_grass = length(idx_grass);
X = [L_forest,L_grass,L_Savanna];

fig=figure;
p=pie(X);
title('Domain Area by Landcover','fontsize',24)
set(gca,'fontsize',24)
set(p(2:2:end),'FontSize',24);
p(1).FaceColor = [0, 0.5, 0];
p(3).FaceColor = [0.4940, 0.1840, 0.5560];
p(5).FaceColor = [0.75, 0.75, 0];
labels = {'Forest','Grassland','Savanna'};
lgd = legend(labels);
lgd.Location = 'westoutside';
fig.Position = [-1459         -76         985         604];
set(gca,'fontsize',25)
saveas(fig,'/Users/abolafia/Drought_Fire_Snow/Plots/Pie_Chart_DomainArea_By_LCtype.eps','epsc')
