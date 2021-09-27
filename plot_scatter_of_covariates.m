clc;clear all;close all;

%load data:

%% Z = 0-1100
Climate_Fire_Data = csvread("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter1.csv");
WYs = Climate_Fire_Data(:,1);
BAs = Climate_Fire_Data(:,2);
BAs_Z1 = (BAs*0.000247105) /(10^6);
SWEIs_Z1 = Climate_Fire_Data(:,3);
WinPRCP_Z1 = Climate_Fire_Data(:,4);
WinTMP_Z1 = Climate_Fire_Data(:,5);
WinterPDSI_Z1 = Climate_Fire_Data(:,6);
Spring_PRCP_Z1 = Climate_Fire_Data(:,7); 
Spring_TMP_Z1 = Climate_Fire_Data(:,8); 
Spring_PDSI_Z1 = Climate_Fire_Data(:,9); 
Spring_SWEI_Z1 = Climate_Fire_Data(:,10); 
Spring_VPD_Z1 = Climate_Fire_Data(:,11); 
Winter_VPD_Z1 = Climate_Fire_Data(:,12); 
MODIS_BA_Z1 = Climate_Fire_Data(:,13); 
MODIS_BA_Z1 = (MODIS_BA_Z1*0.000247105) /(10^6);
Spring_ET_Z1 = Climate_Fire_Data(:,14); 
Spring_PET_Z1 = Climate_Fire_Data(:,15); 
Spring_PETminusET_Z1 = Climate_Fire_Data(:,16);
Winter_ET_Z1 = Climate_Fire_Data(:,17); 
Winter_PET_Z1 = Climate_Fire_Data(:,18); 
Winter_DroughtArea_Z1 = Climate_Fire_Data(:,31);
Spring_DroughtArea_Z1 = Climate_Fire_Data(:,32);
Spring_DroughtArea_PDSI_Z1 = Climate_Fire_Data(:,29);
%summer vars:
Summer_VPD_Z1 = Climate_Fire_Data(1:37,20);
Summer_PRCP_Z1 = Climate_Fire_Data(1:37,22);
Summer_Temp_Z1 = Climate_Fire_Data(1:37,23); 
Summer_ET_Z1 = Climate_Fire_Data(1:37,24); 
Summer_PET_Z1 = Climate_Fire_Data(1:37,25); 
Summer_SM_Z1 = Climate_Fire_Data(1:37,27); 

%relate MODIS and MTBS BA linearly to allow MODIS to gap fill MTBS 2020 BA:
BA_MTBS = BAs_Z1(18:36);
BA_MODIS = MODIS_BA_Z1(18:36);
p = polyfit(BA_MODIS,BA_MTBS,1);
BA_lm = @(BA_MODIS) BA_MODIS*p(1) + p(2);
BA_2020_MTBS = BA_lm(MODIS_BA_Z1(37));
BAs_Z1(37) = BA_2020_MTBS;
if (BAs_Z1(37) < 0)
  BAs_Z1(37) = 0;
end

%% Z = 1100-2200
Climate_Fire_Data = csvread("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter2.csv");
WYs = Climate_Fire_Data(:,1);
BAs = Climate_Fire_Data(:,2);
BAs_Z2 = (BAs*0.000247105) /(10^6);
SWEIs_Z2 = Climate_Fire_Data(:,3);
WinPRCP_Z2 = Climate_Fire_Data(:,4);
WinTMP_Z2 = Climate_Fire_Data(:,5);
WinterPDSI_Z2 = Climate_Fire_Data(:,6);
Spring_PRCP_Z2 = Climate_Fire_Data(:,7); 
Spring_TMP_Z2 = Climate_Fire_Data(:,8); 
Spring_PDSI_Z2 = Climate_Fire_Data(:,9); 
Spring_SWEI_Z2 = Climate_Fire_Data(:,10); 
Spring_VPD_Z2 = Climate_Fire_Data(:,11); 
Winter_VPD_Z2 = Climate_Fire_Data(:,12); 
MODIS_BA_Z2 = Climate_Fire_Data(:,13); 
MODIS_BA_Z2 = (MODIS_BA_Z2*0.000247105) /(10^6);
Spring_ET_Z2 = Climate_Fire_Data(:,14); 
Spring_PET_Z2 = Climate_Fire_Data(:,15); 
Spring_PETminusET_Z2 = Climate_Fire_Data(:,16);
Winter_ET_Z2 = Climate_Fire_Data(:,17); 
Winter_PET_Z2 = Climate_Fire_Data(:,18); 
Winter_DroughtArea_Z2 = Climate_Fire_Data(:,31);
Spring_DroughtArea_Z2 = Climate_Fire_Data(:,32);
Spring_DroughtArea_PDSI_Z2 = Climate_Fire_Data(:,29);
%summer vars:
Summer_VPD_Z2 = Climate_Fire_Data(1:37,20);
Summer_PRCP_Z2 = Climate_Fire_Data(1:37,22);
Summer_Temp_Z2 = Climate_Fire_Data(1:37,23); 
Summer_ET_Z2 = Climate_Fire_Data(1:37,24); 
Summer_PET_Z2 = Climate_Fire_Data(1:37,25); 
Summer_SM_Z2 = Climate_Fire_Data(1:37,27); 

%relate MODIS and MTBS BA linearly to allow MODIS to gap fill MTBS 2020 BA:
BA_MTBS = BAs_Z2(18:36);
BA_MODIS = MODIS_BA_Z2(18:36);
p = polyfit(BA_MODIS,BA_MTBS,1);
BA_lm = @(BA_MODIS) BA_MODIS*p(1) + p(2);
BA_2020_MTBS = BA_lm(MODIS_BA_Z2(37));
BAs_Z2(37) = BA_2020_MTBS;
if (BAs_Z2(37) < 0)
  BAs_Z2(37) = 0;
end

%% Z = 2200-3300
Climate_Fire_Data = csvread("/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/AnlysisData/Annual_Fire_Climate_Obs_forecast_Ziter3.csv");
WYs = Climate_Fire_Data(:,1);
BAs = Climate_Fire_Data(:,2);
BAs_Z3 = (BAs*0.000247105) /(10^6);
SWEIs_Z3 = Climate_Fire_Data(:,3);
WinPRCP_Z3 = Climate_Fire_Data(:,4);
WinTMP_Z3 = Climate_Fire_Data(:,5);
WinterPDSI_Z3 = Climate_Fire_Data(:,6);
Spring_PRCP_Z3 = Climate_Fire_Data(:,7); 
Spring_TMP_Z3 = Climate_Fire_Data(:,8); 
Spring_PDSI_Z3 = Climate_Fire_Data(:,9); 
Spring_SWEI_Z3 = Climate_Fire_Data(:,10); 
Spring_VPD_Z3 = Climate_Fire_Data(:,11); 
Winter_VPD_Z3 = Climate_Fire_Data(:,12); 
MODIS_BA_Z3 = Climate_Fire_Data(:,13); 
MODIS_BA_Z3 = (MODIS_BA_Z3*0.000247105) /(10^6);
Spring_ET_Z3 = Climate_Fire_Data(:,14); 
Spring_PET_Z3 = Climate_Fire_Data(:,15); 
Spring_PETminusET_Z3 = Climate_Fire_Data(:,16);
Winter_ET_Z3 = Climate_Fire_Data(:,17); 
Winter_PET_Z3 = Climate_Fire_Data(:,18); 
Winter_DroughtArea_Z3 = Climate_Fire_Data(:,31);
Spring_DroughtArea_Z3 = Climate_Fire_Data(:,32);
Spring_DroughtArea_PDSI_Z3 = Climate_Fire_Data(:,29);

%summer vars:
Summer_VPD_Z3 = Climate_Fire_Data(1:37,20);
Summer_PRCP_Z3 = Climate_Fire_Data(1:37,22);
Summer_Temp_Z3 = Climate_Fire_Data(1:37,23); 
Summer_ET_Z3 = Climate_Fire_Data(1:37,24); 
Summer_PET_Z3 = Climate_Fire_Data(1:37,25); 
Summer_SM_Z3 = Climate_Fire_Data(1:37,27); 

%relate MODIS and MTBS BA linearly to allow MODIS to gap fill MTBS 2020 BA:
BA_MTBS = BAs_Z3(18:36);
BA_MODIS = MODIS_BA_Z3(18:36);
p = polyfit(BA_MODIS,BA_MTBS,1);
BA_lm = @(BA_MODIS) BA_MODIS*p(1) + p(2);
BA_2020_MTBS = BA_lm(MODIS_BA_Z3(37));
BAs_Z3(37) = BA_2020_MTBS;
if (BAs_Z3(37) < 0)
  BAs_Z3(37) = 0;
end


%% calculate weighted average of predictors, weighted by burned area:
BAs = [BAs_Z1;BAs_Z2;BAs_Z3];
BA_weights = [sum(BAs_Z1),sum(BAs_Z2),sum(BAs_Z3)]./sum(BAs);
Total_BAs = (BAs_Z1+BAs_Z2+BAs_Z3);

%% aggregate data for predictors and BA:
MOD_gridcell_area = 463.3127^2;

Spring_DroughtArea = (Spring_DroughtArea_Z1*BA_weights(1) + Spring_DroughtArea_Z2*BA_weights(2) + Spring_DroughtArea_Z3*BA_weights(3));
Spring_DroughtArea = ( (Spring_DroughtArea*MOD_gridcell_area)*0.000247105 )/(10^6);%# of pixels --> millions of acres

Spring_DroughtArea_PDSI = (Spring_DroughtArea_PDSI_Z3*BA_weights(1) + Spring_DroughtArea_PDSI_Z3*BA_weights(2) + Spring_DroughtArea_PDSI_Z3*BA_weights(3));
Spring_DroughtArea_PDSI = ( (Spring_DroughtArea_PDSI*MOD_gridcell_area)*0.000247105 )/(10^6);%# of pixels --> millions of acres

SWEIs = (SWEIs_Z1*BA_weights(1) + SWEIs_Z2*BA_weights(2) + SWEIs_Z3*BA_weights(3));
Spring_PDSI = (Spring_PDSI_Z1*BA_weights(1) + Spring_PDSI_Z2*BA_weights(2) + Spring_PDSI_Z3*BA_weights(3));

Winter_VPD = (Winter_VPD_Z1*BA_weights(1) + Winter_VPD_Z2*BA_weights(2) + Winter_VPD_Z3*BA_weights(3));
Spring_VPD = (Spring_VPD_Z1*BA_weights(1) + Spring_VPD_Z2*BA_weights(2) + Spring_VPD_Z3*BA_weights(3));
Winter_Spring_VPD = (Winter_VPD+Spring_VPD)./2;
Spring_PET = (Spring_PET_Z1*BA_weights(1) + Spring_PET_Z2*BA_weights(2) + Spring_PET_Z3*BA_weights(3));
Winter_PET = (Winter_PET_Z1*BA_weights(1) + Winter_PET_Z2*BA_weights(2) + Winter_PET_Z3*BA_weights(3));
Winter_Spring_PET = (Spring_PET+Winter_PET)./2;
WinTMP = (WinTMP_Z1*BA_weights(1) + WinTMP_Z2*BA_weights(2) + WinTMP_Z3*BA_weights(3));
Spring_PRCP = (Spring_PRCP_Z1*BA_weights(1) + Spring_PRCP_Z2*BA_weights(2) + Spring_PRCP_Z3*BA_weights(3));
Spring_ET = (Spring_ET_Z1*BA_weights(1) + Spring_ET_Z2*BA_weights(2) + Spring_ET_Z3*BA_weights(3));
Winter_ET = (Winter_ET_Z1*BA_weights(1) + Winter_ET_Z2*BA_weights(2) + Winter_ET_Z3*BA_weights(3));
Winter_Spring_ET = (Spring_ET+Winter_ET)./2;
Spring_TMP = (Spring_TMP_Z1*BA_weights(1) + Spring_TMP_Z2*BA_weights(2) + Spring_TMP_Z3*BA_weights(3));
Winter_Spring_TMP = (WinTMP+Spring_TMP)./2;
WinPRCP = (WinPRCP_Z1*BA_weights(1) + WinPRCP_Z2*BA_weights(2) + WinPRCP_Z3*BA_weights(3));
Winter_Spring_PRCP = (WinPRCP+Spring_PRCP);

All_Data = [Spring_DroughtArea,Winter_VPD,Winter_PET,Winter_ET,WinTMP,WinPRCP,Spring_DroughtArea_PDSI,...
    SWEIs,Spring_VPD,Spring_PET,Spring_ET,Spring_TMP,Spring_PRCP,Spring_PDSI,...
    Winter_Spring_VPD,Winter_Spring_PET,Winter_Spring_ET,Winter_Spring_TMP,Winter_Spring_PRCP];
vars={'Snow Drought Area','Winter VPD','Winter PET','Winter ET','Winter Temp.','Winter Prcp.','PDSI Drought Area',...
    'SWEI','Spring VPD','Spring PET','Spring ET','Spring Temp.','Spring Prcp.','Spring PDSI',...
    'Winter-Spring VPD','Winter-Spring PET','Winter-Spring ET','Winter-Spring Temp.','Winter-Spring Prcp.'};

%summer vars:
Summer_VPD = (Summer_VPD_Z1*BA_weights(1) + Summer_VPD_Z2*BA_weights(2) + Summer_VPD_Z3*BA_weights(3));
Summer_PRCP = (Summer_PRCP_Z1*BA_weights(1) + Summer_PRCP_Z2*BA_weights(2) + Summer_PRCP_Z3*BA_weights(3));
Summer_Temp = (Summer_Temp_Z1*BA_weights(1) + Summer_Temp_Z2*BA_weights(2) + Summer_Temp_Z3*BA_weights(3));
Summer_ET = (Summer_ET_Z1*BA_weights(1) + Summer_ET_Z2*BA_weights(2) + Summer_ET_Z3*BA_weights(3));
Summer_PET = (Summer_PET_Z1*BA_weights(1) + Summer_PET_Z2*BA_weights(2) + Summer_PET_Z3*BA_weights(3));
Summer_SM = (Summer_SM_Z1*BA_weights(1) + Summer_SM_Z2*BA_weights(2) + Summer_SM_Z3*BA_weights(3));

%% calculate their correlations with eachother and with BA:
nvars = 19;
%rows = correlation with each var
%1 col for each variable
store_stats=[];
f=figure;
f.Position = [-1908        -141        1779         946];
for i=1:nvars
    current_var = All_Data(:,i);
    [R,p] = corr(Total_BAs,current_var);
    store_stats = [store_stats;R,p];
    R = round(R,2);
    if i<15
        subplot(3,7,i)
    else
       subplot(3,7,i+1) 
    end
    hold on
    XLABEL = vars{i};
    plot(current_var,Total_BAs,'ko','markersize',6,'markerfacecolor','k')
    xlabel(XLABEL,'fontsize',18)
    ylabel('Burned Area','fontsize',18)
    set(gca,'fontsize',18)
    grid on
    box on
    
    %color r text by significance (p<0.1)
    if p < 0.1
        text(median(current_var),0.9*max(Total_BAs),['r = ',num2str(R)],'FontSize',18,'Color','b');
    else
        text(median(current_var),0.9*max(Total_BAs),['r = ',num2str(R)],'FontSize',18,'Color','r');
    end
    pfit=polyfit(current_var,Total_BAs,1);
    y=@(x) pfit(1)*x + pfit(2);
    
    if p<0.1
        plot(current_var,y(current_var),'-b','linewidth',1)
    else
       plot(current_var,y(current_var),'-r','linewidth',1)
    end
    
    xlim([ min([current_var]) max([current_var])])
    ylim([min(Total_BAs) max(Total_BAs)])
end
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/BA_antecedent_predictor_scatter.eps','epsc')

%% scatter plots for summer vars:
%summer vars:
Summer_VPD = (Summer_VPD_Z1*BA_weights(1) + Summer_VPD_Z2*BA_weights(2) + Summer_VPD_Z3*BA_weights(3));
Summer_PRCP = (Summer_PRCP_Z1*BA_weights(1) + Summer_PRCP_Z2*BA_weights(2) + Summer_PRCP_Z3*BA_weights(3));
Summer_Temp = (Summer_Temp_Z1*BA_weights(1) + Summer_Temp_Z2*BA_weights(2) + Summer_Temp_Z3*BA_weights(3));
Summer_ET = (Summer_ET_Z1*BA_weights(1) + Summer_ET_Z2*BA_weights(2) + Summer_ET_Z3*BA_weights(3));
Summer_PET = (Summer_PET_Z1*BA_weights(1) + Summer_PET_Z2*BA_weights(2) + Summer_PET_Z3*BA_weights(3));
Summer_SM = (Summer_SM_Z1*BA_weights(1) + Summer_SM_Z2*BA_weights(2) + Summer_SM_Z3*BA_weights(3));

All_Data = [Summer_SM,Summer_VPD,Summer_PET,Summer_ET,Summer_Temp,Summer_PRCP];
vars={'Soil Moisture','VPD','PET','ET','Temp.','Prcp.'};

%% calculate their correlations with eachother and with BA:
nvars = 6;
%rows = correlation with each var
%1 col for each variable
store_stats=[];
f=figure;
for i=1:nvars
    current_var = All_Data(:,i);
    [R,p] = corr(Total_BAs,current_var);
    store_stats = [store_stats;R,p];
    R = round(R,2);
    subplot(2,3,i)
    hold on
    XLABEL = vars{i};
    plot(current_var,Total_BAs,'ko','markersize',6,'markerfacecolor','k')
    xlabel(XLABEL,'fontsize',18)
    ylabel('Burned Area','fontsize',18)
    set(gca,'fontsize',18)
    grid on
    box on
    
    %color r text by significance (p<0.1)
    if p < 0.1
        text(prctile(current_var,25),0.9*max(Total_BAs),['r = ',num2str(R)],'FontSize',18,'Color','b');
    else
        text(prctile(current_var,25),0.9*max(Total_BAs),['r = ',num2str(R)],'FontSize',18,'Color','r');
    end
    pfit=polyfit(current_var,Total_BAs,1);
    y=@(x) pfit(1)*x + pfit(2);
    
    if p<0.1
        plot(current_var,y(current_var),'-b','linewidth',1)
    else
       plot(current_var,y(current_var),'-r','linewidth',1)
    end
    
    xlim([ min([current_var]) max([current_var])])
    ylim([min(Total_BAs) max(Total_BAs)])
end
pp
f.Position = [-1919         172        1399         641];
saveas(f,'/Users/abolafia/Drought_Fire_Snow/Plots/BA_summer_predictor_scatter.eps','epsc')

pp
%% plot annual temp, snow drought and BA data:
SWEIs_prctiles = SWEIs./max(SWEIs);
Spring_DroughtArea_prctiles = Spring_DroughtArea./max(Spring_DroughtArea);
Winter_Spring_TMP_prctiles = Winter_Spring_TMP./max(Winter_Spring_TMP);
Winter_Spring_PRCP_prctiles = Winter_Spring_PRCP./max(Winter_Spring_PRCP);
Total_BAs_prctile = Total_BAs./max(Total_BAs);

fig=figure;
years = 1984:2020;
hold on
p_BA = plot(years, Total_BAs_prctile,'-r','linewidth',3);
p_prcp = plot(years, Winter_Spring_PRCP_prctiles,'-b','linewidth',3);
p_tmp = plot(years, Winter_Spring_TMP_prctiles,'-','linewidth',3,'color',[0.9290, 0.6940, 0.1250]);
p_drought_area = plot(years, Spring_DroughtArea_prctiles,'-','linewidth',3,'color',[0.4940, 0.1840, 0.5560]);

%include best fit BA line:
p=polyfit(years, Total_BAs_prctile,1);
y=@(x) p(1)*x+p(2);
plot(years,y(years),'--r','linewidth',1)

%include best fit PRCP line:
p=polyfit(years, Winter_Spring_PRCP_prctiles,1);
y=@(x) p(1)*x+p(2);
plot(years,y(years),'--b','linewidth',1)

%include best fit temp line:
p=polyfit(years, Winter_Spring_TMP_prctiles,1);
y=@(x) p(1)*x+p(2);
plot(years,y(years),'--','linewidth',1,'color',[0.9290, 0.6940, 0.1250])

%include best fit drought area line:
p=polyfit(years, Spring_DroughtArea_prctiles,1);
y=@(x) p(1)*x+p(2);
plot(years,y(years),'--','linewidth',1,'color',[0.4940, 0.1840, 0.5560])

fig.Position = [-1961         217        1937         588];
grid on
box on
leg=legend([p_BA,p_prcp,p_tmp,p_drought_area],{'Burned Area','Precip,','Temp.','Snow Drought Area'},'fontsize',24);
xlim([1984 2020])
set(gca,'fontsize',24)
ylabel('Percentiles','fontsize',24)