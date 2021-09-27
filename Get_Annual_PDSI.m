clc;clear all;close all;

%% Script goal:
%this script take monthly time series of PDSI data and converts it to
%average seasonal PDSI:

%% load in monthly pdsi data:
%from Get_Monthly_PDSI_from_PRISM.m
PDSI = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/PRISM_Monthly/Monthly_PRISM_PDSI.mat');
PDSI = PDSI.PDSI;
%% define dates of interest:
PDSI_dates = datenum([1981 10 1]):datenum([2020 9 30]);
PDSI_dates = datevec(PDSI_dates);
idx=find(PDSI_dates(:,3)==1);
PDSI_datevec = PDSI_dates(idx,:);
PDSI_dates = datenum(PDSI_datevec);
ndates = length(PDSI_dates);

%read in LAT and LON before the loop - only done once:
PRISM_lat = linspace(24.08333333,49.91666667,621);
PRISM_lon = linspace(-125,-66.49,1405);
[PRISM_LAT,PRISM_LON] = meshgrid(PRISM_lat,PRISM_lon);

%% get data:
for WY = 1982:2020 
    WY
    
    %initialize outputs:
    mean_annual_pdsi = nan(1405,621);
    mean_pdsi_spring = nan(1405,621);
    mean_pdsi_winter = nan(1405,621);
    
    %define date index for this year:
    idx_this_year = find(PDSI_dates>=datenum([WY-1 10 1]) & PDSI_dates<= datenum([WY 9 30]));
    PDSI_this_year = PDSI(:,idx_this_year);
    PDSI_dates_this_year = PDSI_dates(idx_this_year);
    
    annual_PDSI=[];
    pdsi_spring=[];
    pdsi_winter=[];
    for i = 1 : 1405*621
        %get annual data:
        annual_PDSI(i) = nanmean(PDSI_this_year(i,:));
        
        time_period_winter = 2:5;
        pdsi_winter(i) = nanmean(PDSI_this_year(i,time_period_winter));
        
        time_period_spring = 6:8;
        pdsi_spring(i) = nanmean(PDSI_this_year(i,time_period_spring));
        idx=find(isnan(pdsi_spring));
        
        time_period_summer = 9:12;
        pdsi_summer(i) = nanmean(PDSI_this_year(i,time_period_summer));
    end
    
    %export:
    Data.annual_PDSI=annual_PDSI;
    Data.pdsi_summer = pdsi_summer;
    Data.pdsi_spring = pdsi_spring;
    Data.pdsi_winter = pdsi_winter;
    Data.LAT=PRISM_LAT;
    Data.LON=PRISM_LON;
    outputfilename = sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/PRISM_Annual/Annual_PRISM_PDSI_WY%04d.mat',WY);
    save(outputfilename,'Data');
end