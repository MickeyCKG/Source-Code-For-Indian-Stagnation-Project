function [output,global_var,nlon,nlat,ensemble_list,return_nothing] = Read_CMIP6_monthly(output,model_name,lon_range,lat_range,scenario_name,var_name,start_yr,end_yr,read_global,...
    read_3d_form,start_layer,end_layer,num_ensemble,special_suffix)

    % Mi Zhou @ Princeton, 2023/02

    % flags
    lat_lon_flag=0;
    initiate_flag=zeros(size(output));
    global_var=cell(size(output));
    return_nothing=1;

    if num_ensemble>1 % multiple ensembles are read
        ensemble_list={};
    else
        ensemble_list={'default'};
    end

    folderpath=['./',model_name,'/',scenario_name,special_suffix,'/'];

    if ~exist(folderpath,'dir')
        disp(['<<< Warning: no files available in ',scenario_name, ' scenario from ',model_name])
        return
    end

    filename=dir(fullfile(folderpath));

    for file=1:1:length(filename) % loop over all qualified files and read data

        thisfile=[folderpath,filename(file).name];
        if contains(filename(file).name,'Amon') || contains(filename(file).name,'AERmon') % this is the desired monthly file

            fullname=filename(file).name;
            marker=strfind(fullname,'_');
            var_in_file=fullname(1:marker-1);
           
            ensemble_id=1;
            if num_ensemble>1 % which ensemble is this file?

                this_ensemble=regexp(fullname,'r.i.p.f.','match');
                this_ensemble=this_ensemble{1};

                new_ensemble=1;
                for en=1:1:length(ensemble_list)
                    if strcmp(this_ensemble,ensemble_list{en})
                        new_ensemble=0;
                        ensemble_id=en;
                    end
                end

                if new_ensemble==1
                    ensemble_list{length(ensemble_list)+1}=this_ensemble;
                    ensemble_id=length(ensemble_list);
                end               
            end

            % which var is included in this file?
            var_id=[];
            for var=1:1:length(var_name)
                if strcmp(var_name{var},var_in_file)                   
                    var_id=var;
                    break
                end
            end

            if isempty(var_id) % this file contains var we do not need
                continue
            end

            % what years are included in this file?
            marker=strfind(thisfile,'12.nc');
            start_yr_file=str2double(thisfile(marker-11:marker-8));
            end_yr_file=str2double(thisfile(marker-4:marker-1));

            % calculate the start/end month of this scenario (relative to year 1)
            start_time=(start_yr-1)*12+1;
            end_time=end_yr*12;
    
            % calculate the start/end month of this file (relative to year 1)
            start_time_file=(start_yr_file-1)*12+1;
            end_time_file=end_yr_file*12;

            % caltulate the intersection months
            inter_time=intersect([start_time:end_time],[start_time_file:end_time_file]);
            num_time=length(inter_time); % number of days shall be read from this file
            if num_time==0
                continue
            end
            start_time_scenario=inter_time(1)-start_time+1; % start_time in relative to the start month of this scenario
            start_time_file=inter_time(1)-start_time_file+1; % start_time in relative to the start month of this file

            disp(' ')
            disp([' Reading file: ',thisfile])
            disp(['ensemble = ',num2str(ensemble_id)])

            return_nothing=0;

            % read in geo information
            if lat_lon_flag==0
                lat=ncread(thisfile,'lat',[1],[inf]);
                lon=ncread(thisfile,'lon',[1],[inf]);
                start_x=abs(lon-lon_range(1));
                start_x=find(start_x==min(start_x),1);
                start_y=abs(lat-lat_range(1));
                start_y=find(start_y==min(start_y),1);
                end_x=abs(lon-lon_range(2));
                end_x=find(end_x==min(end_x),1);
                end_y=abs(lat-lat_range(2));
                end_y=find(end_y==min(end_y),1);
                range_x=end_x-start_x+1;
                range_y=end_y-start_y+1;
                [nlat,nlon]=meshgrid(lat(start_y:end_y),lon(start_x:end_x));
                lat_lon_flag=1;
            end

            % initiate output{scenario_id,var_id}: lon X lat X (layer) X month
            num_layer=end_layer(var_id)-start_layer(var_id)+1;
            if initiate_flag(ensemble_id,var_id)==0 
                if num_layer==1 % 2-D variables
                    output{ensemble_id,var_id}=nan(range_x,range_y,end_time-start_time+1);
                else % 3-D variables
                    output{ensemble_id,var_id}=nan(range_x,range_y,num_layer,end_time-start_time+1);
                end

                if read_global(var_id)==1
                    global_var{ensemble_id,var_id}=nan(1,end_time-start_time+1);
                end           

                initiate_flag(ensemble_id,var_id)=1;
            end

            % read in data
            if read_3d_form(var_id)==0 % 2-D variables
                output{ensemble_id,var_id}(:,:,start_time_scenario:start_time_scenario+num_time-1)=ncread(thisfile,var_name{var_id},[start_x,start_y,start_time_file],[range_x,range_y,num_time]);   

                if read_global(var_id)==1
                    lat_weight=cos(lat/180*pi);
                    tmp_global=ncread(thisfile,var_name{var_id},[1,1,start_time_file],[inf,inf,num_time]); % lon X lat X month
                    tmp_global=squeeze(mean(tmp_global,1)); % lat X month

                    for t=1:1:size(tmp_global,2) % dimension of months
                        global_var{ensemble_id,var_id}(start_time_scenario-1+t)=sum(squeeze(tmp_global(:,t)).*lat_weight)/sum(lat_weight);
                    end
                end

            else % 3-D variables
                if num_layer>1 % read in multiple layers
                    output{ensemble_id,var_id}(:,:,:,start_time_scenario:start_time_scenario+num_time-1)=squeeze(ncread(thisfile,var_name{var_id},[start_x,start_y,start_layer(var_id),start_time_file],[range_x,range_y,num_layer,num_time]));
                else % read in single layer
                    output{ensemble_id,var_id}(:,:,start_time_scenario:start_time_scenario+num_time-1)=squeeze(ncread(thisfile,var_name{var_id},[start_x,start_y,start_layer(var_id),start_time_file],[range_x,range_y,num_layer,num_time]));

                    if read_global(var_id)==1
                        lat_weight=cos(lat/180*pi);
                        tmp_global=squeeze(ncread(thisfile,var_name{var_id},[1,1,start_layer(var_id),start_time_file],[inf,inf,num_layer,num_time])); % lon X lat X month
                        tmp_global=squeeze(nanmean(tmp_global,1)); % lat X month
    
                        for t=1:1:size(tmp_global,2) % dimension of months
                            tmp_data=0;
                            tmp_weight=0;
                            for l=1:1:size(tmp_global,1) % dimension of latitudes
                                if ~isnan(tmp_global(l,t))
                                    tmp_data=tmp_global(l,t)*lat_weight(l);
                                    tmp_weight=tmp_weight+lat_weight(l);
                                end
                            end
                            global_var{ensemble_id,var_id}(start_time_scenario-1+t)=sum(tmp_data)/sum(tmp_weight);
                        end
                    end

                end

            end

        end
    end

    if return_nothing==1
        nlon=[];
        nlat=[];
        ensemble_list=[];
    end

end