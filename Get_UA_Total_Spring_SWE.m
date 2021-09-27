clc;clear all;close all;

%% in this script output a matrix for each year recording the total Aprial-June SWE for each grid cell
%read in LAT and LON before the loop - only done once:
UA_lat = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Raw/UA_SWE/UA_SWE_Depth_WY2001.nc','lat');
UA_lon = ncread('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Raw/UA_SWE/UA_SWE_Depth_WY2001.nc','lon');
[UA_LAT,UA_LON] = meshgrid(UA_lat,UA_lon);
in_UA_directory = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Raw/UA_SWE/';
for year=1983:2020
    year
    
    infilename = sprintf('UA_SWE_Depth_WY%04d.nc',year);
    current_SWE =ncread([in_UA_directory,infilename],'SWE');
    annual_dates = datenum([year-1 10 1]):datenum([year 9 30]);
    idx_spring = find(annual_dates>=datenum([year 3 1]) & annual_dates<datenum([year 6 1]));
    current_SWE_spring = current_SWE(:,:,idx_spring);
    total_spring_SWE = nansum(current_SWE,3);
    %output peak SWE:
    Data.current_spring_SWE=total_spring_SWE;
    Data.LAT=UA_LAT;
    Data.LON=UA_LON;
    outputfilename = sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_PeakSWE/UA_SpringSWE_%04d',year);
    save(outputfilename,'Data');
end
