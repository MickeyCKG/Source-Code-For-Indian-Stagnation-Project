%% This is the Matlab module to read and regrid CMIP6 daily output
% Mi Zhou @ Princeton University, 2023.02

clc
clear
tic

%% Control panel
% model and annual day format【attention when copy scripts!】
model_name='NorESM2_LM';
feb_on=0; % set to 1 if Feb=29days
calendar_360=0; % set to 1 if 360 days in each year

% scenario and time【attention when copy scripts!】
scenario_group=2;
GFDL_server=0; % read files from GFDL server
GFDL_path=[];
convert_output=1;
if scenario_group==1 % historical and ssp scenarios
    scenario_name={'historical','ssp126','ssp370','ssp585'};
    yr_st=[1995,2015,2015,2015];
    yr_ed=[2014,2099,2099,2099];
    ens_name={'_r1i1p1f1_gn_'};
    var_name={'tas','uas','vas','sfcWind','pr','ta','ua','va'};
    var_name_output={'tas','sfcWind','pr','ta','ws'};
elseif scenario_group==2 % CO2 forcing scenarios
    scenario_name={'piControl','1pctCO2'};
    yr_st=[1600,1];
    yr_ed=[1749,150];
    ens_name={'_r1i1p1f1_gn_'};
    var_name={'tas','uas','vas','sfcWind','pr','ta'};
    var_name_output={'tas','sfcWind','pr','ta'};
end

% domain configure
lon_range=[60,98]; 
lat_range=[4 40]; 
regrid_fine_scale=20;
re_lon=[63:1:95];
re_lat=[6:1:38];

% temporal resolution and vertical layer 
read_3d_form=[0,0,0,0,0,2,1,1]; % 0 for 2d variables, 1 for day (850hPa=2,500hPa=4), 2 for Eday (925hPa=2,850hPa=3)
start_layer=[0,0,0,0,0,2,4,4];
end_layer=[0,0,0,0,0,2,4,4];

%% Read in output
% use function/module
disp(['<<< Reading ',model_name,' Group ',num2str(scenario_group)])
[output,nlon,nlat,yr_flag,file_read] = Read_CMIP6_daily(model_name,ens_name,var_name,var_name_output,read_3d_form,start_layer,end_layer,lon_range,lat_range,...
    scenario_name,yr_st,yr_ed,feb_on,calendar_360,GFDL_server,GFDL_path,convert_output);

% check whether each ensemble_member have identical numbers of data
for i=1:1:size(output,1)
    for j=1:1:size(output,2)
        tmp=output{i,j};
        tmp(~isnan(tmp))=1;
        tmp(isnan(tmp))=0;
        tmp=squeeze(sum(sum(sum(tmp,4),3),2));
        if length(unique(tmp))==1 && tmp(1)>0
            disp(['identical numbers of data for ',scenario_name{i},' ',var_name_output{j}])
        elseif std(tmp)/mean(tmp)<5e-4
            disp(['trivial variations in the numbers of data for ',scenario_name{i},' ',var_name_output{j}])
        else
            disp(['std/mean in the numbers of data for ',scenario_name{i},' ',var_name_output{j},' =',num2str(std(tmp)/mean(tmp),'%5.5f')])
        end
    end
end

save(['./Matlab_output/',model_name,'_',num2str(scenario_group),'.mat'],'output','nlat','nlon','model_name','var_name_output','yr_flag','yr_st','yr_ed','scenario_name','read_3d_form','start_layer','end_layer','GFDL_server','file_read','-v7.3')

toc
disp('<<< Done Part I')
disp(' ')

%% Regrid output
% use function/module
[rlat,rlon]=meshgrid(re_lat,re_lon);
[regrid_grid,regrid_fraction]=Regrid_Coef(nlon,nlat,rlon,rlat,regrid_fine_scale,1);

output_regrid=cell(length(scenario_name),length(var_name_output)); % scenarios X variables ('tas','sfcWind','pr','ta')
for i=1:1:size(output,1) % dimension of scenarios
    for j=1:1:size(output,2) % dimension of variables
        thisdata=output{i,j}; % raw data: ensembles X lon X lat X days
        if ~isempty(thisdata)
            output_regrid{i,j}=nan(size(output{i,j},1),length(re_lon),length(re_lat),size(output{i,j},4)); % regridded data: ensembles X r_lon X r_lat X days
            for k=1:1:size(output{i,j},1) % dimension of ensemble_members
                for x=1:1:size(output_regrid{i,j},2) % dimension of regridded lon
                    for y=1:1:size(output_regrid{i,j},3) % dimension of regridded lat
                        tmp=zeros(size(thisdata,4),1); 
                        for z=1:1:size(regrid_grid{x,y},1) % dimension of original grids that go into this regridded grid [x,y]
                            tmp=tmp+squeeze(thisdata(k,regrid_grid{x,y}(z,1),regrid_grid{x,y}(z,2),:))*regrid_fraction{x,y}(z);
                        end
                        output_regrid{i,j}(k,x,y,:)=tmp;
                    end
                end
            end
        end
        toc
        disp(['Done regrid for: scenario= ',scenario_name{i},', var= ',var_name_output{j}])
    end
end

% check available data for output_regird by scenarios, by variables, and by ensembles
data_ava=zeros(size(output_regrid,1),size(output_regrid,2),2,size(output_regrid{1,1},1)); % scenarios X variables X metrics X ensembels
data_sta=zeros(size(output_regrid,1),size(output_regrid,2),2,size(output_regrid{1,1},1)); % scenarios X variables X metrics X ensembels
for i=1:1:size(output_regrid,1) % dimension of scenarios
    for j=1:1:size(output_regrid,2) % dimension of variables
        if ~isempty(output_regrid{i,j})
            for k=1:1:size(output_regrid{1,1},1) % dimension of ensembles
                thisdata=squeeze(output_regrid{i,j}(k,:,:,:));
                data_ava(i,j,1,k)=length(find(~isnan(thisdata)));
                data_ava(i,j,2,k)=length(find(thisdata==0));
                thisdata=reshape(thisdata,1,[]);
                data_sta(i,j,1,k)=nanmean(thisdata);
                data_sta(i,j,2,k)=nanstd(thisdata);
            end
        end
    end
end

save(['./Matlab_output/',model_name,'_regrid_',num2str(scenario_group),'.mat'],'output_regrid','rlat','rlon','model_name','var_name_output','yr_flag','yr_st','yr_ed','scenario_name',...
    'read_3d_form','start_layer','end_layer','GFDL_server','file_read','feb_on','calendar_360','data_ava','data_sta','-v7.3')

toc
disp('<<< Done Part II')
disp(' ')