clc;clear all;close all;

%% SWEI mods

store_Taylor_Error=[];
for npredictors = 2:11
    infilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_SWEImod_PCA_%dvars_AICfit_Drop1Taylor.csv',npredictors);
    opts = detectImportOptions(infilename);
    data=readtable(infilename,opts);
    data = table2array(data);
    taylor=data(2,:);
    store_Taylor_Error = [store_Taylor_Error;median(taylor)];
end
SWEI_Drop1_Taylor = store_Taylor_Error;
idx=find(store_Taylor_Error==max(store_Taylor_Error));
sprintf('%d variables are optimal for SWEI models',idx+1)


%% PDSI mods
store_Taylor_Error=[];
for npredictors = 2:11
    infilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_PDSImod_PCA_%dvars_AICfit_Drop1Taylor.csv',npredictors);
    opts = detectImportOptions(infilename);
    data=readtable(infilename,opts);
    data = table2array(data);
    taylor=data(2,:);
    store_Taylor_Error = [store_Taylor_Error;median(taylor)];
end
PDSI_Drop1_Taylor = store_Taylor_Error;

idx=find(store_Taylor_Error==max(store_Taylor_Error));
sprintf('%d variables are optimal for PDSI models',idx+1)

%% SWEI+PDSI mods
store_Taylor_Error=[];
for npredictors = 2:11
    infilename=sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Best_Obs_Forecast_outputs/BAs_1984_2020_Obs_SWEI_PDSI_mod_PCA_%dvars_AICfit_Drop1Taylor.csv',npredictors);
    opts = detectImportOptions(infilename);
    data=readtable(infilename,opts);
    data = table2array(data);
    taylor=data(2,:);
    store_Taylor_Error = [store_Taylor_Error;median(taylor)];
end
SWEI_PDSI_Drop1_Taylor = store_Taylor_Error;

idx=find(store_Taylor_Error==max(store_Taylor_Error));
sprintf('%d variables are optimal for SWEI+PDSI models',idx+1)

%% plot results:

f=figure;
hold on
nvars=2:11;
p1=plot(nvars,SWEI_Drop1_Taylor,'bo','markersize',15,'markerfacecolor','b');
p2=plot(nvars,PDSI_Drop1_Taylor,'ro','markersize',15,'markerfacecolor','r');
% p3=plot(nvars,SWEI_PDSI_Drop1_Taylor,'ko','markersize',15,'markerfacecolor','k');

set(gca,'fontsize',24)
xlabel('Number of Climate Predictors','fontsize',24)
ylabel('Taylor Score (S)','fontsize',24)
leg=legend([p1 p2],{'SWEI ensemble','PDSI ensemble'},'fontsize',24);
grid on
box on
f.Position=[-1809         493         650         288];
xlim([1.5 11.1])
leg.Position=[0.5703    0.2732    0.3238    0.2135];

saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/N_Predictor_Selection_Taylor_Score.eps','epsc')

%for PDSI there is minimal improvement in drop1 score from using more than
%7 climate predictors.
%for SWEI there is degraded skill from using more than 7 predictors