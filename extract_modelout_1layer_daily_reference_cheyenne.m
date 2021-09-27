%%% extract model outputs with one vertical layer (surface variable) daily instantaneous or mena
%%% Cenlin He, 05/01/2020
%Updated by Ronnie 01/05/2021

clc;clear all;close all;
store_FSA=[];
%dy = 17; % number of years since 2000
tpn = 'Basecase_WY2012'; % no, tp, snotel
res = 'res4';

for dy=11:11
    % initialize variables
    year_r = dy + [2000,2000,2000,2001,2001,2001,2001,2001,2001,2001,2001,2001];
    mon_ind = [10,11,12,1,2,3,4,5,6,7,8,9];
    yr_num = unique(year_r);
    
    dir_path = ['/glade/scratch/abolafia/Drought_Fire_Snow/NoahMP/Outputs/',tpn,'/'];
    Time = 1;
    south_north = 601;
    snow_layers = 3;
    west_east = 601;
    hrday = 24;
    %define soil layer parameters:
    soil_layers = 4;
    soil_depth = 0.1+0.3+0.6+1;
    soil_layer1 = 0.1/soil_depth;
    soil_layer2 = 0.3/soil_depth;
    soil_layer3 = 0.6/soil_depth;
    soil_layer4 = 1.0/soil_depth;
    
    if (yr_num(2) == 2004 || yr_num(2) == 2008 || yr_num(2) == 2012 || yr_num(2) == 2016 || yr_num(2) == 2020)
        nday_in_mon = [31,29,31,30,31,30,31,31,30,31,30,31];
    else
        nday_in_mon = [31,28,31,30,31,30,31,31,30,31,30,31];
    end
    ndays = nday_in_mon(mon_ind);
    tdays = sum(ndays);
    thrs = tdays * hrday;
    
    
    FSA = NaN(west_east,south_north,tdays);
    FIRA = NaN(west_east,south_north,tdays);
    GRDFLX = NaN(west_east,south_north,tdays);
    HFX = NaN(west_east,south_north,tdays);
    LH = NaN(west_east,south_north,tdays);
    SNOWH_m = NaN(west_east,south_north,tdays);
    SWE_m = NaN(west_east,south_north,tdays);
    ACSNOW = NaN(west_east,south_north,tdays);
    ACSNOM = NaN(west_east,south_north,tdays);
    ALBEDO = NaN(west_east,south_north,tdays);
    SOIL_M = NaN(west_east,south_north,tdays);
    SOIL_T = NaN(west_east,south_north,tdays);
    ETRAN = NaN(west_east,south_north,tdays);
    EDIR = NaN(west_east,south_north,tdays);
    
    day_ct = 0;
    for imon = 1:length(mon_ind)
        yy_out = num2str(year_r(imon));
        if mon_ind(imon) < 10
            mon_out = ['0',num2str(mon_ind(imon))];
        else
            mon_out = num2str(mon_ind(imon));
        end
        nday_m = ndays(imon);
        
        for iday = 1:nday_m
            day_ct = day_ct + 1
            if iday < 10
                day_out = ['0',num2str(iday)];
            else
                day_out = num2str(iday);
            end
            fnamein=strings(24,1);
            
            FSA_hr = NaN(west_east,south_north,24);
            FIRA_hr = NaN(west_east,south_north,24);
            GRDFLX_hr = NaN(west_east,south_north,24);
            HFX_hr = NaN(west_east,south_north,24);
            LH_hr = NaN(west_east,south_north,24);
            ETRAN_hr = NaN(west_east,south_north,24);
            EDIR_hr = NaN(west_east,south_north,24);
            SNOWH_m_hr = NaN(west_east,south_north,24);
            SWE_m_hr = NaN(west_east,south_north,24);
            ACSNOW_hr = NaN(west_east,south_north,24);
            ACSNOM_hr = NaN(west_east,south_north,24);
            ALBEDO_hr = NaN(west_east,south_north,24);
            SOIL_M_hr = NaN(west_east,south_north,24);
            SOIL_T_hr = NaN(west_east,south_north,24);
            
            for ihr = 0:23
                if ihr < 10
                    hr_out = ['0',num2str(ihr)];
                else
                    hr_out = num2str(ihr);
                end
                % print out date
                date_out = [yy_out,mon_out,day_out,hr_out];
                % read in hrly data file
                fnamein(ihr+1,1) = [date_out,'.LDASOUT_DOMAIN1']
                FSA_hr(:,:,ihr+1) = ncread([dir_path,char(fnamein(ihr+1,1))],'FSA');
                FIRA_hr(:,:,ihr+1) = ncread([dir_path,char(fnamein(ihr+1,1))],'FIRA');
                GRDFLX_hr(:,:,ihr+1) = ncread([dir_path,char(fnamein(ihr+1,1))],'GRDFLX');
                HFX_hr(:,:,ihr+1) = ncread([dir_path,char(fnamein(ihr+1,1))],'HFX');
                LH_hr(:,:,ihr+1) = ncread([dir_path,char(fnamein(ihr+1,1))],'LH');
                ETRAN_hr(:,:,ihr+1) = ncread([dir_path,char(fnamein(ihr+1,1))],'ETRAN');
                EDIR_hr(:,:,ihr+1) = ncread([dir_path,char(fnamein(ihr+1,1))],'EDIR');
                SNOWH_hr(:,:,ihr+1) = ncread([dir_path,char(fnamein(ihr+1,1))],'SNOWH');
                SWE_hr(:,:,ihr+1) = ncread([dir_path,char(fnamein(ihr+1,1))],'SNEQV');
                ACSNOW_hr(:,:,ihr+1) = ncread([dir_path,char(fnamein(ihr+1,1))],'ACSNOW');
                ACSNOM_hr(:,:,ihr+1) = ncread([dir_path,char(fnamein(ihr+1,1))],'ACSNOM');
                ALBEDO_hr(:,:,ihr+1) = ncread([dir_path,char(fnamein(ihr+1,1))],'ALBEDO');
                %get soil data layer by layer:
                SM_all = ncread([dir_path,char(fnamein(ihr+1,1))],'SOIL_M');
                SM1 = squeeze(SM_all(:,1,:));
                SM2 = squeeze(SM_all(:,2,:));
                SM3 = squeeze(SM_all(:,3,:));
                SM4 = squeeze(SM_all(:,4,:));
                SOIL_M_hr(:,:,ihr+1) = SM1.*soil_layer1 + SM2.*soil_layer2 + SM3.*soil_layer3 + SM4.*soil_layer4;
                
                ST_all = ncread([dir_path,char(fnamein(ihr+1,1))],'SOIL_T');
                ST1 = squeeze(ST_all(:,1,:));
                ST2 = squeeze(ST_all(:,2,:));
                ST3 = squeeze(ST_all(:,3,:));
                ST4 = squeeze(ST_all(:,4,:));
                SOIL_T_hr(:,:,ihr+1) = ST1.*soil_layer1 + ST2.*soil_layer2 + ST3.*soil_layer3 + ST4.*soil_layer4;
            end % ihr
            
            % compute daily nanmean or instanteneous
            FSA(:,:,day_ct) = nanmean(FSA_hr(:,:,:),3);
            FIRA(:,:,day_ct) = nanmean(FIRA_hr(:,:,:),3);
            GRDFLX(:,:,day_ct) = nanmean(GRDFLX_hr(:,:,:),3);
            HFX(:,:,day_ct) = nanmean(HFX_hr(:,:,:),3);
            LH(:,:,day_ct) = nanmean(LH_hr(:,:,:),3);
            ETRAN(:,:,day_ct) = nanmean(ETRAN_hr(:,:,:),3);
            EDIR(:,:,day_ct) = nanmean(EDIR_hr(:,:,:),3);
            SWE_m(:,:,day_ct) = nanmean(SWE_hr(:,:,:),3);
            SNOWH_m(:,:,day_ct) = nanmean(SNOWH_hr(:,:,:),3);
            ACSNOW(:,:,day_ct) = ACSNOW_hr(:,:,1);
            ACSNOM(:,:,day_ct) = ACSNOM_hr(:,:,1);
            ALBEDO(:,:,day_ct) = nanmean(ALBEDO_hr(:,:,:),3);
            SOIL_M(:,:,day_ct) = nanmean(SOIL_M_hr(:,:,:),3);
            SOIL_T(:,:,day_ct) = nanmean(SOIL_T_hr(:,:,:),3);
            
        end % iday
    end % imon
    
    day_ct % print total days
    
    
    % save to nc file
    'start saving data....'
    
    diroutname=['/glade/scratch/abolafia/Drought_Fire_Snow/NoahMP/Outputs/Basecase_WY2012/analysis/',tpn];
    if exist(diroutname,'dir')==0
        CMD_mk=['mkdir -p ',diroutname];
        system(CMD_mk);
    end
    fnameout = ['/glade/scratch/abolafia/Drought_Fire_Snow/NoahMP/Outputs/Basecase_WY2012/analysis/',tpn,'/daily_2D_WY',num2str(yr_num(2)),'_',tpn,'.nc'];
    
    %remove file if it already exists
    if exist(fnameout,'file')>0
        CMD_rm=['rm ',fnameout];
        system(CMD_rm);
    end
    nccreate(fnameout,'west_east','Datatype','single','Dimensions', {'west_east',west_east});
    ncwrite(fnameout,'west_east', single(1:west_east));
    nccreate(fnameout,'south_north','Datatype','single','Dimensions', {'south_north',south_north});
    ncwrite(fnameout,'south_north', single(1:south_north));
    nccreate(fnameout,'Time','Datatype','single','Dimensions', {'Time',tdays});
    ncwrite(fnameout,'Time',single(1:tdays));
    
    
    nccreate(fnameout,'FSA','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'FSA', single(FSA));
    ncwriteatt(fnameout,'FSA','units','W/m2');
    
    nccreate(fnameout,'FIRA','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'FIRA', single(FIRA));
    ncwriteatt(fnameout,'FIRA','units','W/m2');
    
    nccreate(fnameout,'GRDFLX','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'GRDFLX', single(GRDFLX));
    ncwriteatt(fnameout,'GRDFLX','units','W/m2');
    
    nccreate(fnameout,'HFX','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'HFX', single(HFX));
    ncwriteatt(fnameout,'HFX','units','W/m2');
    
    nccreate(fnameout,'LH','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'LH', single(LH));
    ncwriteatt(fnameout,'LH','units','W/m2');
    
    nccreate(fnameout,'ETRAN','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'ETRAN', single(ETRAN));
    ncwriteatt(fnameout,'ETRAN','units','kg/m2/s');
    
    nccreate(fnameout,'EDIR','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'EDIR', single(EDIR));
    ncwriteatt(fnameout,'EDIR','units','kg/m2/s');
    
    nccreate(fnameout,'SNOWH_nanmean','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'SNOWH_nanmean', single(SNOWH_m));
    ncwriteatt(fnameout,'SNOWH_nanmean','units','m');
    
    nccreate(fnameout,'SWE_nanmean','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'SWE_nanmean', single(SWE_m));
    ncwriteatt(fnameout,'SWE_nanmean','units','mm');
    
    nccreate(fnameout,'ACSNOW','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'ACSNOW', single(ACSNOW));
    ncwriteatt(fnameout,'ACSNOW','units','mm');
    
    nccreate(fnameout,'ACSNOM','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'ACSNOM', single(ACSNOM));
    ncwriteatt(fnameout,'ACSNOM','units','mm');
    
    nccreate(fnameout,'ALBEDO','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'ALBEDO', single(ALBEDO));
    ncwriteatt(fnameout,'ALBEDO','units','-');
    
    nccreate(fnameout,'SOIL_M','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'SOIL_M', single(SOIL_M));
    ncwriteatt(fnameout,'SOIL_M','units','-');
    
    nccreate(fnameout,'SOIL_T','Datatype','single',...
        'Dimensions', {'west_east',west_east,'south_north',south_north,'Time',tdays});
    ncwrite(fnameout,'SOIL_T', single(SOIL_T));
    ncwriteatt(fnameout,'SOIL_T','units','-');
    
end % dy iyear

'end of program'



