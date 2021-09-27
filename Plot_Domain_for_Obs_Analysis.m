clc;clear all;close all;

%% define bounding box for area of analysis: western CONUS:
latlim = [32 51];
lonlim = [-125 -104];

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

store_SWEI = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_PeakSWE_SWEI/UA_SWEI_Hindcast.mat','store_SWEI');
store_SWEI = store_SWEI.store_SWEI;

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
Z_thresh_min=0;
Z_thresh_max=3300;
idx_Zthresh = find(Elevation_vec_BB>=Z_thresh_min & Elevation_vec_BB<=Z_thresh_max);

%arbitrary
peakSWE_thresh = 100; %mm
idx_PeakSWEthresh = find(max(store_peak_swe_BB(:,1:end)') >= peakSWE_thresh);
%screen based on landcover type:
idx_LC = Total_LC_index_BB;

IDX_screened = intersect(idx_Zthresh,idx_PeakSWEthresh);
IDX_screened = intersect(IDX_screened,idx_LC);
IDX_screened=intersect(IDX_screened,IDX_DEWS);

IDX_screened_SWEthreshOnly = intersect(IDX_DEWS,idx_PeakSWEthresh);

%define domain:
DOMAIN_X = BA_lonvec(idx_BB);
DOMAIN_Y = BA_latvec(idx_BB);

DOMAIN_X = DOMAIN_X(IDX_screened);
DOMAIN_Y = DOMAIN_Y(IDX_screened);

f=figure;
sz=1;
geoscatter(DOMAIN_Y,DOMAIN_X,sz,[0.9 0.9 0.9],'.');
hold on

%show state lines:
states=shaperead('usastatehi', 'UseGeoCoords', true);
for i=1:49
    lat = states(i).Lat;
    lon = states(i).Lon;
    geoplot(lat,lon,'LineWidth',1.5,'color','k','linewidth',2);
end

% % %Show DEWS regions:
% % DEWS_areas = shaperead('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Raw/DEWS_Regions/DEWS/DEWS_AllRegions202103.shp');
% % DEWS_PNW = DEWS_areas(6);
% % DEWS_CA_NV = DEWS_areas(2);
% % DEWS_West_mount = DEWS_areas(3);
% % DEWS_missouri = DEWS_areas(5);
% % geoplot(DEWS_PNW.Y,DEWS_PNW.X,'LineWidth',8,'color','k','linewidth',5);
% % geoplot(DEWS_CA_NV.Y,DEWS_CA_NV.X,'LineWidth',8,'color','k','linewidth',3);
% % geoplot(DEWS_West_mount.Y,DEWS_West_mount.X,'LineWidth',5,'color','r','linewidth',5);
% % geoplot(DEWS_missouri.Y,DEWS_missouri.X,'LineWidth',5,'color','r','linewidth',5);

%backdrop the satellite image
geobasemap satellite
set(gca,'fontsize',25)
f.Position =[-1916        -110         864         915];
% % f.PaperOrientation='portrait';
% % f.PaperPosition = [0   0   12.0000*0.83   12.7083*0.83];
geolimits([32 49.1],[-125 -104])

%% calculate what percent of total western US is included in domain:
ncoords = length(IDX_screened);
Z_thresh_min=-10000;
Z_thresh_max=10000;
idx_Zthresh = find(Elevation_vec_BB>=Z_thresh_min & Elevation_vec_BB<=Z_thresh_max);
IDX_screened_wUS = intersect(IDX_DEWS,idx_Zthresh);
total_coords = length(IDX_screened_wUS);
sprintf('%.2f of pixels are in western US',ncoords/total_coords*100)

%% calculate the percent of burned area covered by region of interest:
store_SWEI_spring = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/UA_SpringSWE_SWEI/UA_Spring_SWEI_Hindcast.mat','store_SWEI');
store_SWEI_spring = store_SWEI_spring.store_SWEI;

MOD_gridcell_area = 463.3127^2; %m2
store_SWEI=[];

store_BA_MTBS=[];
store_BurnFraction_screened_MTBS=[];

store_BA_MODIS=[];
store_BurnFraction_screened_MODIS=[];

store_BurnFraction_Total_MTBS=[];
store_BurnFraction_Total_MODIS=[];
for WY = 1984:2020
    WY
    %get burn area:
    store_BurnFraction = [];
    store_BurnFraction_MODIS=[];
    if WY < 2020
        for month = 6:9
            BA_filename = sprintf('Gridded_MTBS_%04d%02d.mat',WY,month);
            MTBS_Data = load([MTBS_BA_dir,BA_filename]);
            BurnFraction = MTBS_Data.Gridded_MTBS.BurnFraction;
            BurnFraction_vec = BurnFraction(:);
            BurnFraction_vec(idx_NAN_BAgrid) =[];
            store_BurnFraction = [store_BurnFraction,BurnFraction_vec(idx_BB)];
        end
    end
    
    if WY>2001 %get MODIS BA:
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
    %get MTBS fractions:
    if WY < 2020
        Total_fraction = nanmax(store_BurnFraction')';
        Total_fraction_western_US = Total_fraction(IDX_screened_wUS);
        Total_fraction_screened = Total_fraction(IDX_screened);
        Total_BA_western_US = nansum(Total_fraction_western_US.*MOD_gridcell_area);
        Total_BA_screened = nansum(Total_fraction_screened.*MOD_gridcell_area);
        
        Total_fraction_screened_SWEthresh_only = Total_fraction(IDX_screened_SWEthreshOnly);
        Total_BA_screened_SWEthresh_only = nansum(Total_fraction_screened_SWEthresh_only.*MOD_gridcell_area);
        
        %store annual BA:
        store_BA_MTBS = [store_BA_MTBS;Total_BA_western_US,Total_BA_screened,Total_BA_screened_SWEthresh_only];
        store_BurnFraction_screened_MTBS = [store_BurnFraction_screened_MTBS,Total_fraction_screened];
        store_BurnFraction_Total_MTBS = [store_BurnFraction_Total_MTBS,Total_fraction_western_US];
    end
    
    %get MODIS fractions:
    if WY > 2001
        Total_fraction = nanmax(store_BurnFraction_MODIS')';
        Total_fraction_western_US = Total_fraction(IDX_screened_wUS);
        Total_fraction_screened = Total_fraction(IDX_screened);
        Total_BA_western_US = nansum(Total_fraction_western_US.*MOD_gridcell_area);
        Total_BA_screened = nansum(Total_fraction_screened.*MOD_gridcell_area);
        
        Total_fraction_screened_SWEthresh_only = Total_fraction(IDX_screened_SWEthreshOnly);
        Total_BA_screened_SWEthresh_only = nansum(Total_fraction_screened_SWEthresh_only.*MOD_gridcell_area);
        
        %store annual BA:
        store_BA_MODIS = [store_BA_MODIS;Total_BA_western_US,Total_BA_screened,Total_BA_screened_SWEthresh_only];
        store_BurnFraction_screened_MODIS = [store_BurnFraction_screened_MODIS,Total_fraction_screened];
        store_BurnFraction_Total_MODIS = [store_BurnFraction_Total_MODIS,Total_fraction_western_US];
    end
    
    %get SWEI - spring SWE:
    current_SWEI_spring = squeeze(store_SWEI_spring(:,WY-1983));
    current_SWEI_spring_screened = current_SWEI_spring(IDX_screened);
    springSWEI_WY_mean = nanmean(current_SWEI_spring_screened);
    store_SWEI = [store_SWEI;springSWEI_WY_mean];
end
store_BurnFraction_screened = [store_BurnFraction_screened_MTBS,store_BurnFraction_screened_MODIS(:,end)];

%report percent change in domain size based on screening by elevation and vegetation:
sprintf('The domain area is reduced by %.2f percent from veg and Z screening',length(IDX_screened)/length(IDX_screened_SWEthreshOnly)*100)
%report corresponding change in total BA:
BA_screened = sum([store_BA_MTBS(:,2);store_BA_MODIS(end,2)]);
BA_SWE_screened_only = sum([store_BA_MTBS(:,3);store_BA_MODIS(end,3)]);
sprintf('Burned area is reduced by %.2f percent from veg and Z screening',BA_screened/BA_SWE_screened_only*100)

AllTime_Fraction=nanmax(store_BurnFraction_screened,[],2);
idx_burned = find(AllTime_Fraction>0);
l1=length(idx_burned);
cmap=hot(11);
cmap = flipud(cmap);
cmap(1,:) = [0.9 0.9 0.9];
colormap(cmap);
C=AllTime_Fraction(idx_burned);
geoscatter(DOMAIN_Y(idx_burned),DOMAIN_X(idx_burned),sz,C,'.');
c=colorbar;
ylabel(c, 'Burn Fraction')

% % %plot bounding box of idaho future projection region:
% % study_latlim = [44 46];
% % study_lonlim = [-116 -114];
% % gp1=geoplot([study_latlim(1) study_latlim(1)],[study_lonlim(1) study_lonlim(2)],'b-','linewidth',3);
% % gp2=geoplot([study_latlim(2) study_latlim(2)],[study_lonlim(1) study_lonlim(2)],'b-','linewidth',3);
% % gp3=geoplot([study_latlim(1) study_latlim(2)],[study_lonlim(1) study_lonlim(1)],'b-','linewidth',3);
% % gp4=geoplot([study_latlim(1) study_latlim(2)],[study_lonlim(2) study_lonlim(2)],'b-','linewidth',3);
% % 

saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/SnowySpatialDomain_ObsAnalysis.eps','epsc')

%% create scatter plot of BA in domain vs total western US:
%use lm for the screened domain to convert MODIS to MTBS BA:
screened_BA_MTBS = store_BA_MTBS(19:36,2);
screened_BA_MODIS = store_BA_MODIS(1:end-1,2);

p=polyfit(screened_BA_MODIS,screened_BA_MTBS,1);
lm_BA_screened = @(BA_MODIS) p(1)*BA_MODIS + p(2);
BA_2020_screened = lm_BA_screened(store_BA_MODIS(end,2));

BA_screened = [store_BA_MTBS(:,2);BA_2020_screened];
BA_screened = (BA_screened*0.000247105) /(10^6); %m2 --> millions of acres

%use lm for the entire western US domain to convert MODIS to MTBS BA:
wUS_BA_MTBS = store_BA_MTBS(19:36,1);

wUS_BA_MODIS = store_BA_MODIS(1:end-1,1);

p=polyfit(wUS_BA_MODIS,wUS_BA_MTBS,1);
lm_BA_wUS = @(BA_MODIS) p(1)*BA_MODIS + p(2);
BA_2020_wUS = lm_BA_wUS(store_BA_MODIS(end,1));

BA_wUS = [store_BA_MTBS(:,1);BA_2020_wUS];
BA_wUS = (BA_wUS*0.000247105) /(10^6); %m2 --> millions of acres

a=sum(BA_screened)/sum(BA_wUS);
sprintf('Study domain includes %.2f of burned area from 1984-2020',a*100)


f=figure;
hold on
plot(BA_screened,BA_wUS,'ro','markersize',12,'markerfacecolor','r')
p=polyfit(BA_screened,BA_wUS,1);
y=@(x) p(1)*x+p(2);
plot(BA_screened,y(BA_screened),'-k','linewidth',2.5)
xlabel({'Burned Area in Study Domain';'(millions of acres)'},'fontsize',34)
ylabel({'Burned Area in total  western US';'(millions of acres)'},'fontsize',34)
set(gca,'fontsize',34)
xlim([0 max(BA_screened)])
ylim([0 9])
f.Position=[194    67   962   738];
grid on
box on
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/Screened_TotalwUS_Scatter.eps','epsc')
%export time series of both:
WesternUS_BA.screened = BA_screened;
WesternUS_BA.total = BA_wUS;
outfilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/BA_time_series_screened_and_total.mat');
save(outfilename,'WesternUS_BA');
pp

%% plot pie chart of BA by vegetation type:
LC_types_screened = MODIS_LC_BB(IDX_screened);
idx_ENF = find(LC_types_screened == 1);
idx_EBF = find(LC_types_screened == 2);
idx_DNF = find(LC_types_screened == 3);
idx_DBF = find(LC_types_screened == 4);
idx_MF = find(LC_types_screened == 5);
idx_CS = find(LC_types_screened == 6);
idx_WSav = find(LC_types_screened == 8);
idx_Sav = find(LC_types_screened == 9);
idx_grass = find(LC_types_screened == 10);

L_ENF = length(idx_ENF);
L_EBF = length(idx_EBF);
L_DNF = length(idx_DNF);
L_DBF = length(idx_DBF);
L_MF = length(idx_MF);
L_CS = length(idx_CS);
L_Savanna = length([idx_WSav;idx_Sav]);
L_grass = length(idx_grass);
Pie_Data = [L_ENF,L_EBF,L_DNF,L_DBF,L_MF,L_Savanna,L_grass];
labels = {'ENF','EBF','DNF','DBF','MF','Savanna','Grass'};
f=figure;
pie(Pie_Data,labels);

%% include spatial plot that show burn fraction for entire western CONUS as well:
store_BurnFraction_W_US = [store_BurnFraction_Total_MTBS,store_BurnFraction_Total_MODIS(:,end)];

f=figure;
sz=1;
DOMAIN_X = BA_lonvec(idx_BB);
DOMAIN_Y = BA_latvec(idx_BB);
DOMAIN_X = DOMAIN_X(IDX_screened);
DOMAIN_Y = DOMAIN_Y(IDX_screened);
geoscatter(DOMAIN_Y,DOMAIN_X,sz,[0.9 0.9 0.9],'.');
hold on

%show state lines:
states=shaperead('usastatehi', 'UseGeoCoords', true);
for i=1:49
    lat = states(i).Lat;
    lon = states(i).Lon;
    geoplot(lat,lon,'LineWidth',1.5,'color','k','linewidth',2);
end

%backdrop the satellite image
geobasemap satellite
set(gca,'fontsize',25)
f.Position =[-1916        -110         864         915];
% % f.PaperOrientation='portrait';
% % f.PaperPosition = [0   0   12.0000*0.83   12.7083*0.83];
geolimits([32 49.1],[-125 -104])

AllTime_Fraction=nanmax(store_BurnFraction_W_US,[],2);
idx_burned = find(AllTime_Fraction>0);
l2=length(idx_burned);
cmap=hot(11);
cmap = flipud(cmap);
colormap(cmap);
C=AllTime_Fraction(idx_burned);
DOMAIN_X = BA_lonvec(idx_BB);
DOMAIN_Y = BA_latvec(idx_BB);
DOMAIN_X = DOMAIN_X(IDX_screened_wUS);
DOMAIN_Y = DOMAIN_Y(IDX_screened_wUS);

geoscatter(DOMAIN_Y(idx_burned),DOMAIN_X(idx_burned),sz,C,'.');
c=colorbar;
ylabel(c, 'Burn Fraction')

saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/SnowySpatialDomain_ObsAnalysis_Total_W-US_BA.png')

