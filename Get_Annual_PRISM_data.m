clc;clear all;close all;

%% in this script output a matrix for each year recording the:

%% define directories:
precip_dir = '/glade/scratch/abolafia/Drought_Fire_Snow/DataSets/PRISM/precip/daily_nc/';
precip_fieldname = 'RAINRATE_Pri';

temp_dir = '/glade/scratch/abolafia/Drought_Fire_Snow/DataSets/PRISM/surfTemp/daily_nc/';
temp_fieldname = 'tmean';

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
store_prcp=[];
store_temp=[];
for d=1:ndates
    year = datevec_daily(d,1);
    month = datevec_daily(d,2);
    day = datevec_daily(d,3);
    [year month day]
    
    prcp_filename = sprintf('Rainrate_PRISM_%04d%02d%02d.nc',year,month,day);
    temp_filename = sprintf('SurfTemp_PRISM_%04d%02d%02d.nc',year,month,day);
    
    if exist([precip_dir,prcp_filename],'file')>0
        prcp = ncread([precip_dir,prcp_filename],precip_fieldname);
        prcp = prcp.*(3600*24); %mm/day
        store_prcp = cat(3,store_prcp,prcp);
    else
        disp('missing prcp file')
        store_prcp = cat(3,store_prcp,nan(1405,621));
    end
    
    if exist([temp_dir,temp_filename],'file')>0
        temp = ncread([temp_dir,temp_filename],temp_fieldname);
        store_temp = cat(3,store_temp,temp);
    else
        disp('missing temp file')
        store_temp = cat(3,store_temp,nan(1405,621));
    end
    
    
    idx = find(WY_end_idx == d);
    if idx == WY_iter
        daysinyear = size(store_temp);
        daysinyear = daysinyear(3);
        
        %initialize outputs:
        total_annual_prcp = nansum(store_prcp,3);
        mean_annual_temp = nanmean(store_temp,3);
        
        total_prcp_summer = nan(1405,621);
        mean_temp_summer = nan(1405,621);
        
        total_prcp_spring = nan(1405,621);
        mean_temp_spring = nan(1405,621);
        
        total_prcp_winter = nan(1405,621);
        mean_temp_winter = nan(1405,621);
        
        WY_dates = datenum([year-1 10 1]):datenum([year 9 30]);
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
            
            [r,c]=ind2sub(size(store_prcp),i);
            
            %get summer data:
            prcp_sum = nansum(store_prcp(r,c,tim_period_summer));
            temp_mean = nanmean(store_temp(r,c,tim_period_summer));
            total_prcp_summer(r,c) = prcp_sum;
            mean_temp_summer(r,c) = temp_mean;
            
            %get spring data:
            prcp_sum = nansum(store_prcp(r,c,tim_period_spring));
            temp_mean = nanmean(store_temp(r,c,tim_period_spring));
            total_prcp_spring(r,c) = prcp_sum;
            mean_temp_spring(r,c) = temp_mean;
            
            %get Winter data:
            prcp_sum = nansum(store_prcp(r,c,tim_period_winter));
            temp_mean = nanmean(store_temp(r,c,tim_period_winter));
            total_prcp_winter(r,c) = prcp_sum;
            mean_temp_winter(r,c) = temp_mean;
        end
        
        %export:
        %         Data.total_annual_prcp=total_annual_prcp;
        %         Data.mean_annual_temp = mean_annual_temp;
        Data.total_prcp_summer = total_prcp_summer;
        Data.mean_temp_summer = mean_temp_summer;
        Data.total_prcp_winter = total_prcp_winter;
        Data.mean_temp_winter = mean_temp_winter;
        Data.total_prcp_spring = total_prcp_spring;
        Data.mean_temp_spring = mean_temp_spring;
        Data.LAT=PRISM_lat;
        Data.LON=PRISM_lon;
        
        outputfilename = sprintf('/glade/scratch/abolafia/Drought_Fire_Snow/AnalysisData/PRISM_Annual/Annual_PRISM_met_WY%04d.mat',year);
        save(outputfilename,'Data');
        
        %reset:
        WY_iter = WY_iter+1;
        store_prcp=[];
        store_temp=[];
    end
end
