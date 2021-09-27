clc;clear all;close all;

%% NN PRISM data to BA grid:

% % PRISM_lat = linspace(24.08333333,49.91666667,621);
% % PRISM_lon = linspace(-125,-66.49,1405);
% % [PRISM_LAT,PRISM_LON] = meshgrid(PRISM_lat,PRISM_lon);
% % 
% % BA_lat = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200708.nc','XLAT_M');
% % BA_lon = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200708.nc','XLONG_M');
% % 
% % PRISM_GRID = [PRISM_LAT(:)';PRISM_LON(:)'];
% % BA_GRID = [BA_lat(:)';BA_lon(:)'];
% % IDX_NN_PRISM=nearestneighbour(BA_GRID,PRISM_GRID);
% % 
% % save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_PRISM.mat','IDX_NN_PRISM');

%% NN UA SWE data to BA grid
UA_lat = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Raw/UA_SWE/UA_SWE_Depth_WY2001.nc','lat');
UA_lon = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Raw/UA_SWE/UA_SWE_Depth_WY2001.nc','lon');
[UA_LAT,UA_LON] = meshgrid(UA_lat,UA_lon);
UA_LAT = double(UA_LAT);
UA_LON = double(UA_LON);

BA_lat = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200708.nc','XLAT_M');
BA_lon = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200708.nc','XLONG_M');

UA_GRID = [UA_LAT(:)';UA_LON(:)'];
BA_GRID = [BA_lat(:)';BA_lon(:)'];

IDX_NN_UAswe=nearestneighbour(BA_GRID,UA_GRID);

save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_UAswe.mat','IDX_NN_UAswe');