clc;clear all;close all;

% % % % %% Define PRISM grid:
% % PRISM_lat = linspace(24.08333333,49.91666667,621);
% % PRISM_lon = linspace(-125,-66.49,1405);
% % [PRISM_LAT,PRISM_LON] = meshgrid(PRISM_lat,PRISM_lon);
% % 
% % %% define input directories:
% % precip_dir = '/glade/scratch/abolafia/Drought_Fire_Snow/DataSets/PRISM/precip/daily_nc/';
% % precip_fieldname = 'RAINRATE_Pri';
% % 
% % temp_dir = '/glade/scratch/abolafia/Drought_Fire_Snow/DataSets/PRISM/surfTemp/daily_nc/';
% % temp_fieldname = 'tmean';
% % 
% % %% define dates of interest:
% % start_date = datenum([1981 10 1]);
% % end_date = datenum([2020 9 30]);
% % datelist_daily = start_date:end_date;
% % datevec_daily = datevec(datelist_daily);
% % ndates = length(datelist_daily);
% % WY_end_idx = find(datevec_daily(:,2)==9 & datevec_daily(:,3)==30);
% % 
% % 
% % %% Get PRISM data:
% % last_month = datevec_daily(1,2);
% % 
% % store_prcp=[];
% % store_temp=[];
% % store_monthly_prcp=[];
% % store_monthly_temp=[];
% % month_iter=0;
% % for d=1:ndates
% %     year = datevec_daily(d,1);
% %     month = datevec_daily(d,2);
% %     day = datevec_daily(d,3);
% %     
% %     prcp_filename = sprintf('Rainrate_PRISM_%04d%02d%02d.nc',year,month,day);
% %     temp_filename = sprintf('SurfTemp_PRISM_%04d%02d%02d.nc',year,month,day);
% %     
% %     upcoming_date = datenum([year month day])+1;
% %     upcoming_month = datevec(upcoming_date); upcoming_month=upcoming_month(2);
% %     if upcoming_month == month
% %         if exist([precip_dir,prcp_filename],'file')>0
% %             prcp = ncread([precip_dir,prcp_filename],precip_fieldname);
% %             prcp = prcp.*(3600*24); %mm/day
% %             store_prcp = cat(3,store_prcp,prcp);
% %         else
% %             store_prcp = cat(3,store_prcp,nan(1405,621));
% %             disp('missing prcp')
% %         end
% %         
% %         if exist([temp_dir,temp_filename],'file')>0
% %             temp = ncread([temp_dir,temp_filename],temp_fieldname);
% %             store_temp = cat(3,store_temp,temp);
% %         else
% %             store_temp = cat(3,store_temp,nan(1405,621));
% %             disp('missing temp')
% %         end
% %     else
% %         disp('new month')
% %         %store monthly total PRCP (mm/month):
% %         monthly_prcp = nansum(store_prcp,3);
% %         store_monthly_prcp = [store_monthly_prcp,monthly_prcp(:)];
% %         store_prcp=[];
% %         idx_nan_precip = find(isnan(monthly_prcp(:)));
% %         
% %         %store average temp (C)
% %         monthly_temp = nanmean(store_temp,3);
% %         store_monthly_temp = [store_monthly_temp,monthly_temp(:)];
% %         store_temp=[];
% %         idx_nan_temp = find(isnan(monthly_temp(:)));
% % % %          assert(length(idx_nan_precip) == length(idx_nan_temp),'NaN in temp');
% %         
% %         [year month day]
% %         month_iter = month_iter+1;
% %     end
% %     %for detecting month changes
% %     last_month = month;
% % end
% % 
% % Data.monthly_prcp=store_monthly_prcp;
% % Data.monhtly_temp = store_monthly_temp;
% % Data.dates = datelist_daily;
% % Data.LAT=PRISM_LAT;
% % Data.LON=PRISM_LON;
% % outputfilename = '/glade/scratch/abolafia/Drought_Fire_Snow/AnalysisData/PRISM_PDSI/Monthly_PRISM_met.mat';
% % save(outputfilename,'Data', '-v7.3');

Data = load('/glade/scratch/abolafia/Drought_Fire_Snow/AnalysisData/PRISM_PDSI/Monthly_PRISM_met.mat');

%% Define input PDSI constant parameters:
% Specify which dimension of the temperature and precipitation data is the
% monthly time step. By default, pdsi treats the first dimension as time.
dim = 2;

% years: The first and last year of the T and P data. A two element vector.
years = [1982 , 2020];

% cafecYears: The first and last year of the calibration period used to
%    compute CAFEC normalization. A two element vector.
cafecYears = [1982 , 2020];

%show progress bar
showprogress = true;

%% calculate PDSI:
lats = Data.Data.LAT;
lats = lats(:);

T = Data.Data.monhtly_temp;
idx_nan_T = find(isnan(T));

P = Data.Data.monthly_prcp;
idx_nan_P = find(isnan(P));

idx_nan = unique([idx_nan_T;idx_nan_P]);
idx_nan = sort(idx_nan);

%place hold NaN values with 0 for now, then replace with NaN after PDSIfunc:
nmonths = length(P(1,:));
P(idx_nan) = rand(length(idx_nan),1).*10^-8;
T(idx_nan) = rand(length(idx_nan),1).*10^-8;
% lats(idx_nan,:)=[];

%define soil water parameters:
% awcs: The available water capacity of the surface layer for each site
%    (in mm). A common default value is 25.4 mm. A numeric array. Must have
%    the same size as lats.
awcs = linspace(25.4,25.4,length(P(:,1)));
awcs = awcs';
% awcu: The available water capacity of the underlying layer for each
%    site (in mm). A common default value is 127 mm. A numeric array. Must
%    have the same size as lats.
awcu = linspace(127,127,length(P(:,1)));
awcu = awcu';

%calculate PDSI
[X, Xm, Z, PE] = pdsi(T, P, years, lats, awcs, awcu, cafecYears, dim, showprogress);
PDSI = X;
PDSI(idx_nan) = NaN;
Modified_PDSI = Xm;
Modified_PDSI(idx_nan) = NaN;
%export data:
outputfilename = '/glade/scratch/abolafia/Drought_Fire_Snow/AnalysisData/PRISM_PDSI/Monthly_PRISM_PDSI.mat';
save(outputfilename,'PDSI', '-v7.3');
disp('finished')