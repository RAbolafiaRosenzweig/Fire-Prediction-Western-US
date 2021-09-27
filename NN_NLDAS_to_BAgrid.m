clc;clear all;close all;

%% NN NLDAS data to BA grid:

NLDAS_lat = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NLDAS/Daily_Data/NLDAS_19980515.nc','lat');
NLDAS_lat = double(NLDAS_lat);
NLDAS_lon = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NLDAS/Daily_Data/NLDAS_19980515.nc','lon');
NLDAS_lon = double(NLDAS_lon);

[NLDAS_LAT,NLDAS_LON] = meshgrid(NLDAS_lat,NLDAS_lon);

BA_lat = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200708.nc','XLAT_M');
BA_lon = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Merged_netcdf/Merged_BurnArea_200708.nc','XLONG_M');

NLDAS_GRID = [NLDAS_LAT(:)';NLDAS_LON(:)'];
BA_GRID = [BA_lat(:)';BA_lon(:)'];
IDX_NN_NLDAS=nearestneighbour(BA_GRID,NLDAS_GRID);

save('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/IDX_NN_NLDAS.mat','IDX_NN_NLDAS');

