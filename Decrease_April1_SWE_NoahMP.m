clc;clear all;close all

%% Define what percent increase you want to induce on April 1 SWE:
Pct_Change = 30/100;

%% Define original restart file:
original_filename = '/Users/abolafia/Drought_Fire_Snow/Data/NoahMP/RESTART.2012040100_DOMAIN1.nc';
west_east = 601;
south_north = 601;
%% define new restart file:
new_filename = '/Users/abolafia/Drought_Fire_Snow/Data/NoahMP/RESTART.2012040100_DOMAIN1_30pctLowerSWE.nc';
if exist(new_filename,'file') > 0
    CMD_rm = ['rm ',new_filename];
    system(CMD_rm);
end
%% Load in original SWE, snow depth and snice data:
SNEQV_orignal = ncread(original_filename,'SNEQV');
SNICE_original = ncread(original_filename,'SNICE');
SNOWH_original = ncread(original_filename,'SNOWH');

%% remove SWE and snow depth fields from the new restart file:
CMD_ncks = ['/opt/local/bin/ncks -x -v SNEQV,SNICE,SNOWH ',original_filename,' ',new_filename];
system(CMD_ncks);

%% write updated snow variables:

%write new SWE variable:
SNEQV_new = SNEQV_orignal - (SNEQV_orignal.*Pct_Change);
delta_SWE = SNEQV_new - SNEQV_orignal;
nccreate(new_filename,'SNEQV','Datatype','single','Dimensions', {'south_north',south_north,'west_east',west_east,'Time',1});
ncwrite(new_filename,'SNEQV', single(SNEQV_new));

%write SNICE variable (add addtional SWE to SNICE)
SNICE_new = SNICE_original;
SNICE_new(:,3,:) = squeeze(SNICE_original(:,3,:)) + delta_SWE;
fill_value = min(SNICE_new(:));
%prevent SNICE from being negative:
idx_layer_3 = find(SNICE_new(:,3,:) < 0 & SNICE_new(:,3,:)>fill_value);
if length(idx_layer_3) > 0
    for i=1:length(idx_layer_3)
        [r,c] = ind2sub(size(SNEQV_new),idx_layer_3(i));
        site_SNICE3 = SNICE_new(r,3,c);
        neg_value = site_SNICE3;
        SNICE_new(r,3,c) = 0;
        SNICE_new(r,2,c) = SNICE_new(r,2,c)  + neg_value;
        if SNICE_new(r,2,c) < 0
            neg_value = SNICE_new(r,2,c);
            SNICE_new(r,2,c) = 0;
            SNICE_new(r,1,c) = SNICE_new(r,1,c) + neg_value;
            if SNICE_new(r,1,c) < 0
                SNICE_new(r,1,c) = 0;
            end
        end
    end
end

idx_check =find(SNICE_new>fill_value);
SNICE_check = SNICE_new(idx_check);
assert(nanmin(SNICE_check) >= 0 ,'negative SNICE');
nccreate(new_filename,'SNICE','Datatype','single','Dimensions', {'south_north',south_north,'snow_layers',3,'west_east',west_east,'Time',1});
ncwrite(new_filename,'SNICE', single(SNICE_new));

%wrtie new snowH variable:
SNOWH_new = SNOWH_original - (SNOWH_original.*Pct_Change);
nccreate(new_filename,'SNOWH','Datatype','single','Dimensions', {'south_north',south_north,'west_east',west_east,'Time',1});
ncwrite(new_filename,'SNOWH', single(SNOWH_new));