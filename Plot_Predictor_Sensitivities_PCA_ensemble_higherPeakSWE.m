clc;clear all;close all;

%% define variables:
vars={'SWEI','Winter Precip.','Winter Temp.','Spring Precip.','Spring Temp.','Spring VPD','Winter VPD','Spring ET','Spring PET','Winter ET','Winter PET','Snow Drought Area','Winter-Spring Temp.','Winter-Spring Precip.','Winter-Spring VPD','Winter-Spring PET','Winter-Spring ET'};


%% load results from top models:
%[best predictors,dropped predictor R2 fit, R2 fit]
best_mod_ids = readtable('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_bestmods_higherPeakSWE.csv');
best_mod_ids = table2array(best_mod_ids);
nmods = length(best_mod_ids);

store_R2_bestfit=[];
store_pbias_bestfit=[];
store_rmse_bestfit=[];
store_mod_data=[];
for i=1:nmods
    mod_id = best_mod_ids(i);
    infilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_Forecast_SE_Zbins_SWEImod_PCA_drop_predictor_sensitivity_mod%d_higherPeakSWE.csv',mod_id);
    mod_data = readtable(infilename);
    mod_data=table2array(mod_data);
    
    %get outputs:
    %predictor ID:
    predictor_IDs = mod_data(:,1);
    %R2:
    dropped_R2 = mod_data(:,2);
    R2_bestfit = mod_data(:,3);
    %PBIAS:
    dropped_pbias = mod_data(:,4);
    pbias_bestfit = mod_data(:,5);
    %RMSE:
    dropped_rmse = mod_data(:,6);
    rmse_bestfit = mod_data(:,7);
    
    %store data:
    store_mod_data=[store_mod_data;mod_data];
    store_R2_bestfit = [store_R2_bestfit;R2_bestfit(1)];
    store_pbias_bestfit = [store_pbias_bestfit;pbias_bestfit(1)];
    store_rmse_bestfit = [store_rmse_bestfit;rmse_bestfit(1)];
end

%% record number of times each predictor appears in 20 best models:
unique_predictors = unique(store_mod_data(:,1));
unique_predictors = sort(unique_predictors);
store_predictor_count=[];
for i=1:length(unique_predictors)
   idx=find(store_mod_data(:,1) == unique_predictors(i));
   store_predictor_count = [store_predictor_count;unique_predictors(i),length(idx)];
end

%% create bar plot showing number of appearances for each predictor in 20 best models:
[s,i] = sort(store_predictor_count(:,2),'descend');
sorted_predictor_ids = store_predictor_count(i,1);
sorted_counts = s;
barplot_categories=vars(i);
X=categorical(barplot_categories);
X = reordercats(X,barplot_categories);

f=figure;
subplot(3,1,1)
hold on
b=bar(X,sorted_counts);
b.BarWidth=0.5;
b.FaceColor = 'r';
b.EdgeColor = 'k';
b.LineWidth = 1.5;
set(gca,'fontsize',20)
grid on
ylabel({'Appearances in';'best models'},'fontsize',20)

%% create a boxplot plot showing the ratio of summary stats explained after removing each predictor:
%only consider predictors that appear in at least 10 models:
idx_over10 = find(sorted_counts>=0);
predictor_IDS = sorted_predictor_ids(idx_over10);

%% R2:
R2_ratio = store_mod_data(:,2)./store_mod_data(:,3);

store_R2_ratios = nan(100,length(idx_over10));

for i = 1:length(idx_over10)
    IDX = find(store_mod_data(:,1) == predictor_IDS(i));
    current_R2_ratio = R2_ratio(IDX);
    store_R2_ratios(1:length(current_R2_ratio),i) = current_R2_ratio;
end

subplot(3,1,3)
labels=barplot_categories(idx_over10);
hold on
b=boxplot(store_R2_ratios,'Notch','off','Labels',labels,'whisker',1);
set(b,'linew',3)
set(gca,'fontsize',20)
ylabel({'r^{2} ratio'},'fontsize',20)
xtickangle(20)
ylim([0.6 1.1])
xlim([0.5 length(idx_over10)+0.5])
grid on
plot(linspace(0.5,length(idx_over10)+0.5,10),linspace(nanmedian(store_R2_ratios(:,3)),nanmedian(store_R2_ratios(:,3)),10),'--r')
plot(linspace(0.5,length(idx_over10)+0.5,10),linspace(prctile(store_R2_ratios(:,3),25),prctile(store_R2_ratios(:,3),25),10),'--b')
plot(linspace(0.5,length(idx_over10)+0.5,10),linspace(prctile(store_R2_ratios(:,3),75),prctile(store_R2_ratios(:,3),75),10),'--b')

%% rmse:
rmse_ratio = store_mod_data(:,6)./store_mod_data(:,7);

store_rmse_ratios = nan(100,length(idx_over10));

for i = 1:length(idx_over10)
    IDX = find(store_mod_data(:,1) == predictor_IDS(i));
    current_rmse_ratio = rmse_ratio(IDX);
    store_rmse_ratios(1:length(current_rmse_ratio),i) = current_rmse_ratio;
end

subplot(3,1,2)
labels=barplot_categories(idx_over10);
hold on
b=boxplot(store_rmse_ratios,'Notch','off','Labels',labels,'whisker',1);
set(b,'linew',3)
set(gca,'fontsize',20)
ylabel({'RMSE ratio'},'fontsize',20)
xtickangle(20)
ylim([0.85 1.6])
xlim([0.5 length(idx_over10)+0.5])
grid on
plot(linspace(0.5,length(idx_over10)+0.5,10),linspace(nanmedian(store_rmse_ratios(:,3)),nanmedian(store_rmse_ratios(:,3)),10),'--r')
plot(linspace(0.5,length(idx_over10)+0.5,10),linspace(prctile(store_rmse_ratios(:,3),25),prctile(store_rmse_ratios(:,3),25),10),'--b')
plot(linspace(0.5,length(idx_over10)+0.5,10),linspace(prctile(store_rmse_ratios(:,3),75),prctile(store_rmse_ratios(:,3),75),10),'--b')

f.Position=[-1919        -144        1868         949];
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/BA_predict_1984_2020_Obs_SWEImod_PCA_predictor_sensitivity_ensemble_higherPeakSWE.eps','epsc')
