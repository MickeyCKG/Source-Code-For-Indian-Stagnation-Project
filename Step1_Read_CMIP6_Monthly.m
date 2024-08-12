%% This is the Matlab module to read and regrid CMIP6 daily output
% Mi Zhou @ Princeton University, 2023.02

clc
clear
tic

%% Control panel
% scenario
scenario_name='historical'
% note that ssp370_lowNTCF will automatically read in ssp370
% note that read_global is only for monthly read-in, not for daily read-in

% model and variables
switch(scenario_name)
    case 'piControl'
        model_name={'GFDL_ESM4','UKESM1','NorESM2_LM','NorESM2_MM','CanESM5','MRI_ESM2','MIROC6_ES2L','INM_CM4_8'}; 
        start_yr=[1,1960,1600,1200,5550,1850,1850,1850];
        end_yr=[150,2109,1749,1349,5699,1999,1999,1999];
        num_ensemble=[1,1,1,1,1,1,1,1,1];

        % variables
        var_name={'tas','ta'};
        read_global=[1,0];

        % vertical layer 
        read_3d_form=[0,1]; % 0 for 2-d variables, 1 for 3-d variables
        start_layer=[0,2]; % 1=1000hPa, 2=925hPa, 3=850hPa, 4=700hPa, 5=600hPa, 6=500hPa
        end_layer=[0,2];
        WindSpeed_id=[];

    case '1pctCO2'
        model_name={'GFDL_ESM4','UKESM1','NorESM2_LM','NorESM2_MM','CanESM5','MRI_ESM2','MIROC6_ES2L','INM_CM4_8'}; 
        start_yr=[1,1850,1,1,1850,1850,1850,1850];
        end_yr=[150,1999,150,150,1999,1999,1999,1999];
        num_ensemble=[1,1,1,1,1,1,1,1];

        % variables
        var_name={'tas','ta','hurs','hur','huss','hus','pr','uas','vas','sfcWind','ua','va','ws','psl'};
        read_global=[1,0,0,0,0,0,0,0,0,0,0,0,0,0];

        % vertical layer 
        read_3d_form=[0,1,0,1,0,1,0,0,0,0,1,1,1,0]; % 0 for 2-d variables, 1 for 3-d variables
        start_layer=[0,2,0,2,0,2,0,0,0,0,2,2,2,0]; % 1=1000hPa, 2=925hPa, 3=850hPa, 4=700hPa, 5=600hPa, 6=500hPa
        end_layer=[0,6,0,6,0,6,0,0,0,0,6,6,6,0];
        WindSpeed_id=[10,13];

    case {'historical','ssp370','ssp585'}  
        % note that all models, except CanESM5, have aerosols
        model_name={'GFDL_ESM4','UKESM1','NorESM2_LM','NorESM2_MM','CanESM5','MRI_ESM2_atm','MIROC6_ES2L','INM_CM4_8','MRI_ESM2_aer'}; % read aerosol and meteorology seperately for MRI_ESM2 since they have different resolution!
        if strcmp(scenario_name,'historical')
            start_yr=[1995,1995,1995,1995,1995,1995,1995,1995,1995];
            end_yr=[2014,2014,2014,2014,2014,2014,2014,2014,2014];
        else
            start_yr=[2015,2015,2015,2015,2015,2015,2015,2015,2015];
            end_yr=[2099,2099,2099,2099,2099,2099,2099,2099,2099];
        end
        num_ensemble=[1,3,1,1,3,1,1,1,1];

        % variables
        var_name={'tas','ta','pr','sfcWind','mmrpm2p5','mmrbc','mmroa','mmrsoa','mmrdust','mmrss','mmrso4','mmrnh4','mmrno3','calculated_pm25','ps'};
        read_global=[1,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

        % vertical layer 
        read_3d_form=[0,1,0,0,1,1,1,1,1,1,1,1,1,0,0]; % 0 for 2-d variables, 1 for 3-d variables
        start_layer=[0,2,0,0,1,1,1,1,1,1,1,1,1,0,0]; % (not for aerosol) 1=1000hPa, 2=925hPa, 3=850hPa, 4=700hPa, 5=600hPa, 6=500hPa
        end_layer=[0,2,0,0,1,1,1,1,1,1,1,1,1,0,0];
        WindSpeed_id=[];

    case {'ssp370_lowNTCF'}
        model_name={'GFDL_ESM4','UKESM1','NorESM2_LM','MRI_ESM2_atm','MRI_ESM2_aer'}; % read aerosol and meteorology seperately for MRI_ESM2 since they have different resolution!
        start_yr=[2015,2015,2015,2015,2015];
        end_yr=[2054,2054,2054,2054,2054];
        num_ensemble=[1,3,3,1,1];

        % variables
        var_name={'tas','ta','pr','sfcWind','mmrpm2p5','mmrbc','mmroa','mmrsoa','mmrdust','mmrss','mmrso4','ps'};
        read_global=[1,0,0,0,0,0,0,0,0,0,0,0];

        % vertical layer 
        read_3d_form=[0,1,0,0,1,1,1,1,1,1,1,0]; % 0 for 2-d variables, 1 for 3-d variables
        start_layer=[0,2,0,0,1,1,1,1,1,1,1,0]; % (not for aerosol) 1=1000hPa, 2=925hPa, 3=850hPa, 4=700hPa, 5=600hPa, 6=500hPa
        end_layer=[0,2,0,0,1,1,1,1,1,1,1,0];
        WindSpeed_id=[];

end

% domain configure
lon_range=[50,110]; % read
lat_range=[0 45]; % read
regrid_fine_scale=20;
re_lon=[55:1:105]; % interp
re_lat=[5:1:40]; % interp
[rlat,rlon]=meshgrid(re_lat,re_lon);
num_days=[31,28,31,30,31,30,31,31,30,31,30,31];

%% Read and Regrid

output_regrid=cell(sum(num_ensemble),length(var_name)); % models (total ensembles) X variables, regridded
CMIP_global=cell(1,length(var_name)); % variables, averaged

st_ensemble=1;
for model=1:1:length(model_name)

    this_var_name=var_name;
    % deal with the fact that meteorology and aerosols are read seperately in MRI_ESM2 due to their different model resolutions
    if strcmp(model_name{model},'MRI_ESM2_atm')
        model_name{model}='MRI_ESM2';
        this_var_name(5:length(var_name))={'nan'};
    elseif strcmp(model_name{model},'MRI_ESM2_aer')
        model_name{model}='MRI_ESM2';
        this_var_name(1:4)={'nan'};
    end

    % deal with different ensemble numbers for SSP370 NorESM2_LM
    if strcmp(model_name{model},'NorESM2_LM') && strcmp(scenario_name,'ssp370')
        special_suffix='_monthly_1ensemble';
    else
        special_suffix=[];
    end

    output=cell(num_ensemble(model),length(this_var_name)); % (ensemble) X variables
    ensemble_id=[st_ensemble:1:st_ensemble+num_ensemble(model)-1];
    ensemble2model(ensemble_id)=model;
    st_ensemble=st_ensemble+num_ensemble(model);

    % READ, use function/module
    disp(' ')
    disp(['<<< Reading ',scenario_name,' ',model_name{model}])
    [output,global_var,nlon,nlat,ensemble_list{model},return_nothing]=Read_CMIP6_monthly(output,model_name{model},lon_range,lat_range,scenario_name,this_var_name,start_yr(model),end_yr(model),read_global,...
            read_3d_form,start_layer,end_layer,num_ensemble(model),special_suffix);

    if return_nothing==1
        continue % next model, please!
    end

    % read in additional scenario
    if strcmp(scenario_name,'ssp370_lowNTCF') % read in 'ssp370' as well
        tmp_output=cell(num_ensemble(model),length(this_var_name)); % (ensemble) X variables
        [tmp_output,tmp_global_var,nlon,nlat,ensemble_list{model}]=Read_CMIP6_monthly(tmp_output,model_name{model},lon_range,lat_range,'ssp370',this_var_name,start_yr(model),end_yr(model),read_global,...
               read_3d_form,start_layer,end_layer,num_ensemble(model));

        % merge 'ssp370' (after ssp370_lowNTCF)
        for i=1:1:size(tmp_output,1)
            for j=1:1:size(tmp_output,2)
                if ~isempty(tmp_output{i,j})
                    output{i,j}=cat(length(size(tmp_output{i,j})),output{i,j},tmp_output{i,j});
                else
                    output{i,j}=cat(length(size(output{i,j})),output{i,j},nan(size(output{i,j})));
                end
            end
        end
    end

    disp(' ')
    toc
    disp(['<<< Done Reading ',model_name{model}])

    % retrieve global mean from raw data
    for var=1:1:size(global_var,2) % dimension of variables
        this_flag=read_global(var);

        if this_flag==1
            for en=1:1:num_ensemble(model)
                if ~isempty(global_var{en,var})
                    num_data=size(global_var{en,var},2);
                    CMIP_global{var}(ensemble_id(en),1:num_data)=global_var{en,var};

                    if strcmp(scenario_name,'ssp370_lowNTCF') % merge 'ssp370' (after ssp370_lowNTCF)
                        CMIP_global{var}(ensemble_id(en),num_data+1:2*num_data)=tmp_global_var{en,var};
                    end

                else
                    CMIP_global{var}(ensemble_id(en),:)=nan;
                end
            end
        end
    end

    if ~isempty(WindSpeed_id)
        for en=1:1:num_ensemble(model)
            if isempty(output{en,WindSpeed_id(1)}) % if sfcWind is unavailable, calculate it from uas and vas
                output{en,WindSpeed_id(1)}=sqrt(output{en,WindSpeed_id(1)-2}.^2+output{en,WindSpeed_id(1)-1}.^2);
            end
            output{en,WindSpeed_id(2)}=sqrt(output{en,WindSpeed_id(2)-2}.^2+output{en,WindSpeed_id(2)-1}.^2); % additional variable for wind speed on pressure levels
        end
    end

    % REGRID, use function/module
    [regrid_grid,regrid_fraction]=Regrid_Coef(nlon,nlat,rlon,rlat,regrid_fine_scale,1);

    for en=1:1:num_ensemble(model) % dimension of ensembles
        model_list{ensemble_id(en)}=[model_name{model},'_',ensemble_list{model}{en}];
        for var=1:1:size(output,2) % dimension of variables
            this_flag=read_3d_form(var);

            if this_flag==0 % 2-D variable

                thisdata=output{en,var}; % raw data: lon X lat X months
                if ~isempty(thisdata)
                    output_regrid{ensemble_id(en),var}=nan(length(re_lon),length(re_lat),size(thisdata,3)); % regridded data: r_lon X r_lat X months
                    for x=1:1:size(output_regrid{ensemble_id(en),var},1) % dimension of regridded lon
                        for y=1:1:size(output_regrid{ensemble_id(en),var},2) % dimension of regridded lat
                            tmp=zeros(size(thisdata,3),1); 
                            for z=1:1:size(regrid_grid{x,y},1) % dimension of original grids that go into this regridded grid [x,y]
                                tmp=tmp+squeeze(thisdata(regrid_grid{x,y}(z,1),regrid_grid{x,y}(z,2),:))*regrid_fraction{x,y}(z);
                            end
                            output_regrid{ensemble_id(en),var}(x,y,:)=tmp;
                        end
                    end
                end

            else % 3-D variable

                thisdata=output{en,var}; % raw data: lon X lat X layers X months
                if ~isempty(thisdata)
                    output_regrid{ensemble_id(en),var}=nan(length(re_lon),length(re_lat),size(thisdata,3),size(thisdata,4)); % regridded data: r_lon X r_lat X layers X months
                    for x=1:1:size(output_regrid{ensemble_id(en),var},1) % dimension of regridded lon
                        for y=1:1:size(output_regrid{ensemble_id(en),var},2) % dimension of regridded lat
                            for p=1:1:size(output_regrid{ensemble_id(en),var},3) % dimension of layers
                                tmp=zeros(size(thisdata,4),1); 
                                for z=1:1:size(regrid_grid{x,y},1) % dimension of original grids that go into this regridded grid [x,y]
                                    tmp=tmp+squeeze(thisdata(regrid_grid{x,y}(z,1),regrid_grid{x,y}(z,2),p,:))*regrid_fraction{x,y}(z);
                                end
                                output_regrid{ensemble_id(en),var}(x,y,p,:)=tmp;
                            end
                        end
                    end
                end

            end
    
        end
    end

    disp(' ')
    toc
    disp(['<<< Done Regriding ',model_name{model}])
    disp(' ')

end

% Merge MRI_ESM2_atm with MRI_ESM2_aer
if  strcmp(scenario_name,'ssp370_lowNTCF') || strcmp(scenario_name,'historical') || strcmp(scenario_name,'ssp370') || strcmp(scenario_name,'ssp585') 

    % get id of MRI_ESM2 in current model lists
    mri_ens_id=[];
    mri_mod_id=[];
    this_ensemble=0;
    for model=1:1:length(model_name)
        this_ensemble=this_ensemble+num_ensemble(model);
        if contains(model_name{model},'MRI_ESM2')
            mri_ens_id=[mri_ens_id,this_ensemble];
            mri_mod_id=[mri_mod_id,model];
        end
    end

    % merge them into the first id of MRI_ESM2
    for var=1:1:size(output_regrid,2) % dimension of variables
        if isempty(output_regrid{mri_ens_id(1),var}) && ~isempty(output_regrid{mri_ens_id(2),var})
            output_regrid{mri_ens_id(1),var}=output_regrid{mri_ens_id(2),var};
        end
        if ~isempty(CMIP_global{var}) && size(CMIP_global{var},1)>=mri_ens_id(2)
            CMIP_global{var}(mri_ens_id(2),:)=[];
        end
    end
    output_regrid(mri_ens_id(2),:)=[];
    num_ensemble(mri_mod_id(2))=[];
    model_name(mri_mod_id(2))=[];
    model_name{mri_mod_id(1)}='MRI_ESM2';

end

% save gridded results
% save(['Matlab_monthly_output/',scenario_name,'_monthly_after_step1.mat'],'rlat','rlon','output_regrid','CMIP_global','model_name','var_name','num_ensemble','ensemble_list','-v7.3')
save([scenario_name,'_monthly_after_step1.mat'],'rlat','rlon','output_regrid','CMIP_global','model_name','var_name','num_ensemble','ensemble_list','-v7.3')

toc
disp('<<< Done Saving Rawdata')
disp(' ')
