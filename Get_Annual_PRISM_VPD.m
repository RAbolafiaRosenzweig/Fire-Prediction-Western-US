clc;clear all;close all;

%% in this script output a matrix for each year recording the:
%cumulative precip up to peak swe, cumulative annual precip,
%mean temp up to peak swe, annual mean temp:

%% define directories:
vpd_dir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/PRISM_Daily/VPD/';
VPD_fieldname = 'VPD';

%read in LAT and LON before the loop - only done once:
PRISM_lat = linspace(24.08333333,49.91666667,621);
PRISM_lon = linspace(-125,-66.49,1405);
[PRISM_LAT,PRISM_LON] = meshgrid(PRISM_lat,PRISM_lon);

%% define dates of interest:
start_date = datenum([1982 10 1]);
end_date = datenum([2020 9 30]);
datelist_daily = start_date:end_date;
datevec_daily = datevec(datelist_daily);
ndates = length(datelist_daily);
WY_end_idx = find(datevec_daily(:,2)==9 & datevec_daily(:,3)==30);

%% get data:
WY_iter = 1;
store_VPD=[];
for d=1:ndates
    year = datevec_daily(d,1);
    month = datevec_daily(d,2);
    day = datevec_daily(d,3);
    [year month day]

    vpd_filename = sprintf('VPD_PRISM_%04d%02d%02d.nc',year,month,day);
    
    if exist([vpd_dir,vpd_filename],'file')>0
        VPD = ncread([vpd_dir,vpd_filename],VPD_fieldname);
        store_VPD = cat(3,store_VPD,VPD);
    else
        store_VPD = cat(3,store_VPD,nan(1405,621));
    end
    
            
    idx = find(WY_end_idx == d);
    if idx == WY_iter
        %initialize outputs:
        VPD_spring = nan(1405,621);
        VPD_winter = nan(1405,621);
        
        daysinyear = size(store_VPD);
        daysinyear = daysinyear(3);
        
        if daysinyear == 365
            tim_period_summer = 244:daysinyear;%June - Sept.
            tim_period_spring = 152:243;%March - May
            tim_period_winter = 32:151;%Nov - Feb
        else
            tim_period_summer = 245:daysinyear;%June - Sept.
            tim_period_spring = 153:244;%March - May
            tim_period_winter = 32:152;%Nov - Feb
        end
        
        for i = 1 : 1405*621
            [r,c]=ind2sub(size(store_VPD),i);
            %get spring data:
            vpd_mean = nanmean(store_VPD(r,c,tim_period_spring));
            VPD_spring(r,c) = vpd_mean;
            
            %get Winter data:
            vpd_mean = nanmean(store_VPD(r,c,tim_period_winter));
            VPD_winter(r,c) = vpd_mean;
            
            %get summer data:
            vpd_mean = nanmean(store_VPD(r,c,tim_period_summer));
            VPD_summer(r,c) = vpd_mean;
        end
        
        %export:
        Data.VPD_spring = VPD_spring;
        Data.VPD_winter = VPD_winter;
        Data.VPD_summer = VPD_summer;
        Data.LAT=PRISM_lat;
        Data.LON=PRISM_lon;
        outputfilename = sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/PRISM_Annual/Annual_PRISM_VPD_WY%04d.mat',year);
        save(outputfilename,'Data');
        
        %reset:
        WY_iter = WY_iter+1;
        store_VPD=[];
        store_temp=[];
    end
    
end
