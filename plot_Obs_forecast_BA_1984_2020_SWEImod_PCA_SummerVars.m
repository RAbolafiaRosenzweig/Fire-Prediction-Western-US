clc;clear all;close all

%% Plot data from the best GAMs (BurnArea_GAM_NLDAS_Forecasting_Ziter_validation.R)
%from BurnArea_GAM_NLDAS_Forecasting_Ziter_bestmod_select2.R we determined
%that GAMs with 10 predictors outperform those with less (in general)
%% load in data:
Good_mods = readtable('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Good_Mods_SWEImod_PCA_SummerVars.csv');
Good_mods = table2array(Good_mods);
nmods = length(Good_mods);
years = 1984:2020;

nbins = 3; %3 Z bins x 4 regions:
WYs = repmat(1984:2020,1,nbins);

%% fit

%initialize annual prediction and standard error data:
store_predictions=[];
store_se =[];
for i=1:nmods
    mod_id = Good_mods(i);
    %get current prediction
    infilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_GAM_Zbins_mod%d_SWEImod_PCA_SummerVars.csv',mod_id);
    current_mod_data = readtable(infilename);
    current_mod_data = table2array(current_mod_data);
    %aggregate to annual_data:
    [u,~,j] = unique(WYs);
    Annual_BA_predict = accumarray(j,current_mod_data(:,3),[],@sum);
    Annual_BA_obs = accumarray(j,current_mod_data(:,1),[],@sum);
    
    %constrain predicted BA to 0
    idx = find(Annual_BA_predict<0);
    Annual_BA_predict(idx)=0;
    store_predictions = [store_predictions,Annual_BA_predict];
    %record the standard error along with predition
    infilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_mod%d_SWEImod_PCA_SummerVars.csv',mod_id);
    se = readtable(infilename);
    se = table2array(se);
    store_se = [store_se,se];
end
Annual_BA_obs = [Annual_BA_obs];

ensemble_mean = median(store_predictions')';
ensemble_mean_bestfit = ensemble_mean;
%ensemble_mean = store_predictions;

se = max(store_se')';
%record ensemble upper and lower bounds based on se from gam predictions:
ensemble_upr = prctile(store_predictions',97.5)';
ensemble_lwr = prctile(store_predictions',2.5)';

%% plot data:
lb=ensemble_mean-ensemble_lwr;
ub=ensemble_upr-ensemble_mean;
idx = find(ensemble_mean - lb < 0);
lb(idx) = ensemble_mean(idx);
f=figure;
hold on
shadedErrorBar(years, ensemble_mean, [ub';lb'], 'k-',1);
p1=plot(years, ensemble_mean,'-k','linewidth',3);
%include trend line:
p=polyfit(years,ensemble_mean,1);
y = @(x) p(1)*x + p(2);
plot(years,y(years),'--k','linewidth',2)
sprintf('bestfit increasing trend: %.2f/year',p(1)/mean(ensemble_mean)*100)

p2 = plot(years,Annual_BA_obs,'-r','linewidth',3);
%include trend line:
p=polyfit(years,Annual_BA_obs,1);
y = @(x) p(1)*x + p(2);
plot(years,y(years),'--r','linewidth',2)
sprintf('obs increasing trend: %.2f/year',p(1)/mean(Annual_BA_obs)*100)
ylim([0 6])
legend([p1 p2],{'Predicted','Observed'},'fontsize',25,'location','northwest')
ylabel('Burned area (millions of acres)','fontsize',25)
set(gca,'fontsize',25)
xlim([min(years) max(years)])
f.Position =[1         313        1440         492];
grid on
xticks([1984:2:2020])
xtickangle(45)
box on
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/BA_predict_1984_2020_Obs_GAM_timeseries_SWEImod_PCA_SummerVars_fit.eps','epsc')

%scatter plot:
f=figure;
hold on
cmap = flipud(hot(37));
colormap(cmap)
sz = 250;
c=1984:2020;
p1=scatter(Annual_BA_obs,ensemble_mean,sz,c,'filled',...
    'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'MarkerEdgeColor','k');
%include 1:1 line
plot(linspace(0,6,10),linspace(0,6,10),'--k','linewidth',4)

%format:
set(gca,'fontsize',25)
xlabel('Observed Burned Area (millions of acres)','fontsize',25)
ylabel('Predicted Burned Area (millions of acres)','fontsize',25)
% title('1984 - 2020 Fire Season Burned Area','fontsize',25)
set(gca,'fontsize',25)
%R and RMSE
R = corr(Annual_BA_obs,ensemble_mean,'rows','complete')
RMSE = sqrt( mean((Annual_BA_obs-ensemble_mean).^2) )
%percent of time +/- anomalies are correctly predicted:
Annual_BA_obs_anom = Annual_BA_obs - mean(Annual_BA_obs);
ensemble_mean_anom = ensemble_mean - mean(ensemble_mean);
idx=find( (Annual_BA_obs_anom>0 & ensemble_mean_anom>0) | (Annual_BA_obs_anom<0 & ensemble_mean_anom<0) | (Annual_BA_obs_anom==0 & ensemble_mean_anom==0));
sprintf('correct anomaly sign predicted %.2f percent of the time',length(idx)/length(ensemble_mean_anom)*100)
xlim([0 5.1])
ylim([0 5.1])
c = colorbar;
c.Ticks = 1984:2:2020;
grid on
f.Position = [139    58   887   746];
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/BA_predict_1984_2020_Obs_GAM_scatter_SWEImod_PCA_SummerVars_fit.eps','epsc')

%%consider if a lm can be used to relate climate conditions in study domain to total western US BA:
BA_data = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/BA_time_series_screened_and_total.mat');
Screened_BA = BA_data.WesternUS_BA.screened;
Total_BA = BA_data.WesternUS_BA.total;
p = polyfit(Screened_BA,Total_BA,1);
screened_to_total = @(x) p(1)*x + p(2);
screnned_converted_to_total = screened_to_total(ensemble_mean);
disp('best fit correlation with total western us ba:')
corr(ensemble_mean,Total_BA)

%% Drop1 results

%initialize annual prediction and standard error data:
store_predictions=[];
store_se =[];
for i=1:nmods
    mod_id = Good_mods(i);
    %get current prediction
    infilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_GAM_Zbins_mod%d_SWEImod_PCA_SummerVars.csv',mod_id);
    current_mod_data = readtable(infilename);
    current_mod_data = table2array(current_mod_data);
    %aggregate to annual_data:
    [u,~,j] = unique(WYs);
    Annual_BA_predict = accumarray(j,current_mod_data(:,2),[],@sum);
    Annual_BA_obs = accumarray(j,current_mod_data(:,1),[],@sum);
    
    %constrain predicted BA to 0
    idx = find(Annual_BA_predict<0);
    Annual_BA_predict(idx)=0;
    store_predictions = [store_predictions,Annual_BA_predict];
    %record the standard error along with predition
    infilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_mod%d_SWEImod_PCA_SummerVars.csv',mod_id);
    se = readtable(infilename);
    se = table2array(se);
    store_se = [store_se,se];
end
Annual_BA_obs = [Annual_BA_obs];

ensemble_mean = median(store_predictions')';
ensemble_mean_drop1 = ensemble_mean;

%ensemble_mean = store_predictions;

se = max(store_se')';
%record ensemble upper and lower bounds based on se from gam predictions:
ensemble_upr = prctile(store_predictions',97.5)';
ensemble_lwr = prctile(store_predictions',2.5)';
%% plot data:
lb=ensemble_mean-ensemble_lwr;
ub=ensemble_upr-ensemble_mean;
idx = find(ensemble_mean - lb < 0);
lb(idx) = ensemble_mean(idx);
f=figure;
hold on
shadedErrorBar(years, ensemble_mean, [ub';lb'], 'k-',1);
p1=plot(years, ensemble_mean,'-k','linewidth',3);
%compute trend line (including 2020):
p=polyfit(years(1:end),ensemble_mean(1:end),1);
y = @(x) p(1)*x + p(2);
sprintf('drop1 prediction increasing trend (including 2020): %.2f/year',p(1)/mean(ensemble_mean(1:end))*100)

%include trend line (excluding 2020):
p=polyfit(years(1:end-1),ensemble_mean(1:end-1),1);
y = @(x) p(1)*x + p(2);
plot(years(1:end-1),y(years(1:end-1)),'--k','linewidth',2)
sprintf('drop1 prediction increasing trend (excluding 2020): %.2f/year',p(1)/mean(ensemble_mean(1:end-1))*100)

p2 = plot(years,Annual_BA_obs,'-r','linewidth',3);

%compute trend line (including 2020):
p=polyfit(years(1:end),Annual_BA_obs(1:end),1);
y = @(x) p(1)*x + p(2);
sprintf('obs increasing trend: %.2f/year (including 2020)',p(1)/mean(Annual_BA_obs)*100)

%include trend line (excluding 2020):
p=polyfit(years(1:end-1),Annual_BA_obs(1:end-1),1);
y = @(x) p(1)*x + p(2);
plot(years(1:end-1),y(years(1:end-1)),'--r','linewidth',2)
sprintf('obs increasing trend: %.2f/year (excluding 2020)',p(1)/mean(Annual_BA_obs(1:end-1))*100)

ylim([0 6])
legend([p1 p2],{'Predicted','Observed'},'fontsize',25,'location','northwest')
ylabel('Burned area (millions of acres)','fontsize',25)
set(gca,'fontsize',25)
xlim([min(years) max(years)])
f.Position =[1         313        1440         492];
grid on
xticks([1984:2:2020])
xtickangle(45)
box on
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/BA_predict_1984_2020_Obs_GAM_timeseries_SWEImod_PCA_SummerVars_drop1.eps','epsc')

%scatter plot:
f=figure;
hold on
cmap = flipud(hot(37));
colormap(cmap)
sz = 250;
c=1984:2020;
p1=scatter(Annual_BA_obs,ensemble_mean,sz,c,'filled',...
    'MarkerFaceAlpha',1,'MarkerEdgeAlpha',1,'MarkerEdgeColor','k');
%include 1:1 line
plot(linspace(0,6,10),linspace(0,6,10),'--k','linewidth',4)

%format:
set(gca,'fontsize',25)
xlabel('Observed Burned Area (millions of acres)','fontsize',25)
ylabel('Predicted Burned Area (millions of acres)','fontsize',25)
% title('1984 - 2020 Fire Season Burned Area','fontsize',25)
set(gca,'fontsize',25)
%R and RMSE
R = corr(Annual_BA_obs,ensemble_mean,'rows','complete')
RMSE = sqrt( mean((Annual_BA_obs-ensemble_mean).^2) )
PBIAS =  mean( ((ensemble_mean-Annual_BA_obs)./mean(Annual_BA_obs)) ) * 100;
%percent of time +/- anomalies are correctly predicted:
Annual_BA_obs_anom = Annual_BA_obs - mean(Annual_BA_obs);
ensemble_mean_anom = ensemble_mean - mean(ensemble_mean);
idx=find( (Annual_BA_obs_anom>0 & ensemble_mean_anom>0) | (Annual_BA_obs_anom<0 & ensemble_mean_anom<0) | (Annual_BA_obs_anom==0 & ensemble_mean_anom==0));
sprintf('correct anomaly sign predicted %.2f percent of the time',length(idx)/length(ensemble_mean_anom)*100)
xlim([0 5.1])
ylim([0 5.1])
c = colorbar;
c.Ticks = 1984:2:2020;
grid on
f.Position = [139    58   887   746];
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/BA_predict_1984_2020_Obs_GAM_scatter_SWEImod_PCA_SummerVars_drop1.eps','epsc')

%%consider if a lm can be used to relate climate conditions in study domain to total western US BA:
%drop 1:
BA_data = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/BA_time_series_screened_and_total.mat');
Screened_BA = BA_data.WesternUS_BA.screened;
Total_BA = BA_data.WesternUS_BA.total;
%%consider if a lm can be used to relate climate conditions in study domain to total western US BA:
%drop 1:
BA_data = load('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/Merged_BA/Supplementary/BA_time_series_screened_and_total.mat');
Screened_BA = BA_data.WesternUS_BA.screened;
Total_BA = BA_data.WesternUS_BA.total;

disp('drop 1 correlation with total western us ba:')
corr(ensemble_mean,Total_BA)
