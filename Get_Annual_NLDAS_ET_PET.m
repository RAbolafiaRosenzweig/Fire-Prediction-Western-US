clc;clear all;close all;

indir = '/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NLDAS_ET_PET/Aggregated_Monthly_Grids/';

for year = 1984:2021
    year
    store_ET_winter = [];
    store_PET_winter = [];
    store_PETminusET_winter = [];
    %Nov - Dec previous year
    for month = 11:12
        infilename = sprintf('NLDAS_ET_PET_%04d%02d.mat',year-1,month);
        if exist([indir,infilename],'file') > 0
            Data = load([indir,infilename]);
            Data = Data.NLDAS;
            ET = Data.ET;
            PET = Data.PET;
            PET_minus_ET = PET-ET;
            
            store_ET_winter = cat(3,store_ET_winter,ET);
            store_PET_winter = cat(3,store_PET_winter,PET);
            store_PETminusET_winter = cat(3,store_PETminusET_winter,PET_minus_ET);
        end
    end
    %jan-feb this year
    for month = 1:2
        infilename = sprintf('NLDAS_ET_PET_%04d%02d.mat',year,month);
        if exist([indir,infilename],'file') > 0
            Data = load([indir,infilename]);
            Data = Data.NLDAS;
            ET = Data.ET;
            PET = Data.PET;
            PET_minus_ET = PET-ET;
            
            store_ET_winter = cat(3,store_ET_winter,ET);
            store_PET_winter = cat(3,store_PET_winter,PET);
            store_PETminusET_winter = cat(3,store_PETminusET_winter,PET_minus_ET);
        end
    end
    ET_winter = nanmean(store_ET_winter,3);
    PET_winter = nanmean(store_PET_winter,3);
    PETminusET_winter = nanmean(store_PETminusET_winter,3);
    
    %March - may
    store_ET_spring = [];
    store_PET_spring = [];
    store_PETminusET_spring = [];
    for month =3:5
        infilename = sprintf('NLDAS_ET_PET_%04d%02d.mat',year,month);
        if exist([indir,infilename],'file') > 0
            Data = load([indir,infilename]);
            Data = Data.NLDAS;
            ET = Data.ET;
            PET = Data.PET;
            PET_minus_ET = PET-ET;
            
            store_ET_spring = cat(3,store_ET_spring,ET);
            store_PET_spring = cat(3,store_PET_spring,PET);
            store_PETminusET_spring = cat(3,store_PETminusET_spring,PET_minus_ET);
        end
    end
    ET_spring = nanmean(store_ET_spring,3);
    PET_spring = nanmean(store_PET_spring,3);
    PETminusET_spring = nanmean(store_PETminusET_spring,3);
    
    %June - September
    store_ET_summer = [];
    store_PET_summer = [];
    store_PETminusET_summer = [];
    for month =6:9
        infilename = sprintf('NLDAS_ET_PET_%04d%02d.mat',year,month);
        if exist([indir,infilename],'file') > 0
            Data = load([indir,infilename]);
            Data = Data.NLDAS;
            ET = Data.ET;
            PET = Data.PET;
            PET_minus_ET = PET-ET;
            
            store_ET_summer = cat(3,store_ET_summer,ET);
            store_PET_summer = cat(3,store_PET_summer,PET);
            store_PETminusET_summer = cat(3,store_PETminusET_summer,PET_minus_ET);
        end
    end
    ET_summer = nanmean(store_ET_summer,3);
    PET_summer = nanmean(store_PET_summer,3);
    PETminusET_summer = nanmean(store_PETminusET_summer,3);
    
    Data=[];
    Data.ET_winter=ET_winter;
    Data.PET_winter = PET_winter;
    Data.PETminusET_winter = PETminusET_winter;
    Data.ET_spring = ET_spring;
    Data.PET_spring = PET_spring;
    Data.PETminusET_spring = PETminusET_spring;
    Data.ET_summer = ET_summer;
    Data.PET_summer = PET_summer;
    Data.PETminusET_summer = PETminusET_summer;
    
    outputfilename = sprintf('/Volumes/Pruina_External_Elements/DroughtFireSnow/Data/NLDAS_ET_PET/Annual/Annual_ET_PET_WY%04d.mat',year);
    save(outputfilename,'Data','-v7.3');
end