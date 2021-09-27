clc;clear all;close all


%% load in data:
Good_mods = readtable('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Good_Mods_PDSImod_PCA.csv');
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
    infilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_GAM_Zbins_mod%d_PDSImod_PCA.csv',mod_id);
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
    infilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_mod%d_PDSImod_PCA.csv',mod_id);
    se = readtable(infilename);
    se = table2array(se);
    store_se = [store_se,se];
end
Annual_BA_obs = [Annual_BA_obs];

ensemble_mean = median(store_predictions')';
%ensemble_mean = store_predictions;

se = max(store_se')';
%record ensemble upper and lower bounds based on se from gam predictions:
ensemble_upr = max(store_predictions')';
ensemble_lwr = min(store_predictions')';

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

p2 = plot(years,Annual_BA_obs,'-r','linewidth',3);
%include trend line:
p=polyfit(years,Annual_BA_obs,1);
y = @(x) p(1)*x + p(2);
plot(years,y(years),'--r','linewidth',2)

legend([p1 p2],{'Predicted','Observed'},'fontsize',25,'location','northwest')
ylabel('Burned area (millions of acres)','fontsize',25)
set(gca,'fontsize',25)
xlim([min(years) max(years)])
f.Position =[1         313        1440         492];
grid on
xticks([1984:2:2020])
xtickangle(45)
box on
ylim([0 6])
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/BA_predict_1984_2020_Obs_GAM_timeseries_PDSImod_PCA_fit.eps','epsc')

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
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/BA_predict_1984_2020_Obs_GAM_scatter_PDSImod_PCA_fit.eps','epsc')



%% Drop1 results

%initialize annual prediction and standard error data:
store_predictions=[];
store_se =[];
for i=1:nmods
    mod_id = Good_mods(i);
    %get current prediction
    infilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_GAM_Zbins_mod%d_PDSImod_PCA.csv',mod_id);
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
    infilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_mod%d_PDSImod_PCA.csv',mod_id);
    se = readtable(infilename);
    se = table2array(se);
    store_se = [store_se,se];
end
Annual_BA_obs = [Annual_BA_obs];

ensemble_mean = median(store_predictions')';
%ensemble_mean = store_predictions;

se = max(store_se')';
%record ensemble upper and lower bounds based on se from gam predictions:
ensemble_upr = max(store_predictions')';
ensemble_lwr = min(store_predictions')';
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

p2 = plot(years,Annual_BA_obs,'-r','linewidth',3);
%include trend line:
p=polyfit(years,Annual_BA_obs,1);
y = @(x) p(1)*x + p(2);
plot(years,y(years),'--r','linewidth',2)

legend([p1 p2],{'Predicted','Observed'},'fontsize',25,'location','northwest')
ylabel('Burned area (millions of acres)','fontsize',25)
set(gca,'fontsize',25)
xlim([min(years) max(years)])
f.Position =[1         313        1440         492];
grid on
xticks([1984:2:2020])
xtickangle(45)
box on
ylim([0 6])
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/BA_predict_1984_2020_Obs_GAM_timeseries_PDSImod_PCA_drop1.eps','epsc')

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
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/BA_predict_1984_2020_Obs_GAM_scatter_PDSImod_PCA_drop1.eps','epsc')

