function [output,nlon,nlat,yr_flag,file_read] = Read_CMIP6_daily(model_name,ens_name,var_name,var_name_output,read_3d_form,start_layer,end_layer,lon_range,lat_range,...
    scenario_name,yr_st,yr_ed,feb_on,calendar_360,GFDL_server,GFDL_path,convert_output)

    % Mi Zhou @ Princeton, 2023/02

    file_read={};
    % Initialize: scenarios X variables
    % ('tas','uas','vas','sfcWind','pr','ta'), for scenario_group=1, add 'ua' and 'va'
    output=cell(length(scenario_name),length(var_name));   
    % Initialize: scenarios X variables ('tas','sfcWind','pr','ta'), for scenario_group=1, add 'ws'
    output2=cell(length(scenario_name),length(var_name_output));   

    % do not read uas and vas when GFDL_server=1 and convert_output=1
    if GFDL_server==1 && convert_output==1
        var_name{2}='wwwwwwww';
        var_name{3}='wwwwwwww';
    end

    % flags
    output_flag=zeros(size(output));
    lat_lon_flag=0;
    yr_flag=cell(size(output));

    % search for all qualified files
    n_file=0;
    n_file_read=0;
    allfiles={};
    allfiles2={};
    scenario_record=[];
    for i=1:1:length(scenario_name) % dimension of specified scenarios

        if GFDL_server==1 % read files in GFDL server
            folderpath=GFDL_path{i};
        else % read files in Dellamodel_name
            folderpath=['./',model_name,'/',scenario_name{i},'/'];
        end
        if ~exist(folderpath,'dir')
            disp(['<<< Warning: no files available in ',scenario_name{i}, ' scenario from ',model_name])
            continue
        end

        filename=dir(fullfile(folderpath));
        for file=1:1:length(filename)
            thisfile=[folderpath,filename(file).name];
            if contains(filename(file).name,'.nc') % this is a CMIP6 output file
                n_file=n_file+1;
                allfiles{n_file}=thisfile; % path + filename
                allfiles2{n_file}=filename(file).name; % filename only
                scenario_record(n_file)=i;
            end
        end

    end 

    num_format1=n_file;

    disp(' ')
    disp(['find ',num2str(n_file),' files !'])
    disp(' ')
    
    % loop through all qualified files
    for file=1:1:length(allfiles)
        thisfile=allfiles{file}; % path+filename
        thisfile2=allfiles2{file}; % filename only
        find_var=0;
        find_ens=0;

        % which ensemble-member is this file?
        if file<=num_format1 && GFDL_server~=1 % format follows CMIP6 website, files in Della
            for i=1:1:length(ens_name)
                if contains(thisfile,ens_name{i})
                    ens_id=i;
                    find_ens=1;
                    break
                end
            end      

        end

        % which scenario is this file?
        scenario_id=scenario_record(file);

        % which var is included in this file?
        if file<=num_format1 % files from CMIP6 website
            for i=1:1:length(var_name)
                if (contains(thisfile2,['.',var_name{i},'.']) && GFDL_server==1) || (contains(thisfile2,[var_name{i},'_']) && GFDL_server~=1)
                    var_id=i;
                    find_var=1;
                    break
                end
            end
        end

        % what years are included in this file?
        marker=strfind(thisfile,'0101-');
        start_yr_file=str2double(thisfile(marker-4:marker-1));
        end_yr_file=str2double(thisfile(marker+5:marker+8));
        end_day=str2double(thisfile(marker+11:marker+12));

        % make sure Eday or day file are chosen according to read_3d_form
        % only check for files containing 3D variables w/ format following CMIP6 website
        if GFDL_server~=1 && file<=num_format1 && find_var==1 && read_3d_form(var_id)~=0
            if (read_3d_form(var_id)==1 && contains(thisfile2,'Eday')) || (read_3d_form(var_id)==2 && ~contains(thisfile2,'Eday'))
                continue
            end
        end

        % range of years for this scenario
        start_yr=yr_st(scenario_id);
        end_yr=yr_ed(scenario_id);

        inter_yr=intersect([start_yr:end_yr],[start_yr_file:end_yr_file]);
        if find_var==0 || find_ens==0 || isempty(inter_yr)
            continue
        end

        if calendar_360==1 % 360 days in each year
            % calculate the start/end day of this scenario (relative to year 1)
            start_time=(start_yr-1)*360+1;
            end_time=end_yr*360;
    
            % calculate the start/end day of this file (relative to year 1)
            start_time_file=(start_yr_file-1)*360+1;
            end_time_file=end_yr_file*360;
        else % 365 or 366 days in each year
            % calculate the start/end day of this scenario (relative to year 1)
            start_time=day2index(start_yr,1,1,feb_on,1,1,1);
            end_time=day2index(end_yr,12,31,feb_on,1,1,1);
    
            % calculate the start/end day of this file (relative to year 1)
            start_time_file=day2index(start_yr_file,1,1,feb_on,1,1,1);
            end_time_file=day2index(end_yr_file,12,31,feb_on,1,1,1)+(end_day-31); % end_day-31 is to deal with the case where 1230 is the end day.
        end

        % caltulate the intersection days/months
        inter_time=intersect([start_time:end_time],[start_time_file:end_time_file]);
        num_time=length(inter_time); % number of days shall be read from this file
        start_time_scenario=inter_time(1)-start_time+1; % start_time in relative to the start day/month of this scenario
        start_time_file=inter_time(1)-start_time_file+1; % start_time in relative to the start day/month of this file

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

        % initiate output{scenario_id,var_id}: ensemble_members X (layers) X lon X lat X days
        if output_flag(scenario_id,var_id)==0 
            num_layer=end_layer(var_id)-start_layer(var_id)+1;
            if num_layer==1
                output{scenario_id,var_id}=nan(length(ens_name),range_x,range_y,end_time-start_time+1);
            else
                output{scenario_id,var_id}=nan(length(ens_name),range_x,range_y,num_layer,end_time-start_time+1);
            end
            output_flag(scenario_id,var_id)=1;
            yr_flag{scenario_id,var_id}=zeros(end_yr-start_yr+1,1);
        end

        disp(' ')
        toc
        disp(['>>> (1) Reading file: ',thisfile])

        % read in data
        if read_3d_form(var_id)==0 % 2-D variables

            if isempty(start_x) % read in Global
                this_data=ncread(thisfile,var_name{var_id},[1,1,start_time_file],[inf,inf,num_time]);
            else % read in Regional
                this_data=ncread(thisfile,var_name{var_id},[start_x,start_y,start_time_file],[range_x,range_y,num_time]);
            end

        else % 3-D variables, layers in day file (1000/850/700...), layers in Eday file (1000/925/850...)

            num_layer=end_layer(var_id)-start_layer(var_id)+1;
            this_data=squeeze(ncread(thisfile,var_name{var_id},[start_x,start_y,start_layer(var_id),start_time_file],[range_x,range_y,num_layer,num_time]));

        end
        if num_layer==1
            output{scenario_id,var_id}(ens_id,:,:,start_time_scenario:start_time_scenario+num_time-1)=this_data;
        else
            output{scenario_id,var_id}(ens_id,:,:,:,start_time_scenario:start_time_scenario+num_time-1)=this_data;
        end

        min_yr=1;
        max_yr=end_yr-start_yr+1;
        yr_flag{scenario_id,var_id}(max(start_yr_file-start_yr+1,min_yr):min(end_yr_file-start_yr+1,max_yr))=yr_flag{scenario_id,var_id}(max(start_yr_file-start_yr+1,min_yr):min(end_yr_file-start_yr+1,max_yr))+1;
    
        n_file_read=n_file_read+1;
        file_read{n_file_read,1}=allfiles2{file};
    end
   
    % additional process
    if convert_output==1

        % contains ua/va?
        read_uava=0;
        uava_id=[];
        for i=1:1:length(var_name)
            if strcmp(var_name{i},'ua')
                read_uava=1;
                uava_id(1)=i;
            elseif strcmp(var_name{i},'va')
                read_uava=1;
                uava_id(2)=i;
            end
        end

        for sce=1:1:size(output,1) % dimension of scenarios
    
            % variables in ('tas','uas','vas','sfcWind','pr','ta')
    
            % deal with surfWind, we have two methods to get it
            if ~isempty(output{sce,4}) % directly read sfcWind is the priority
                output2{sce,2}=output{sce,4};
            elseif ~isempty(output{sce,2}) && ~isempty(output{sce,3}) % indirectly calculate surfWind
                output2{sce,2}=sqrt(output{sce,2}.^2+output{sce,3}.^2);
            end
    
            % variables out ('tas','sfcWind','pr','ta'): ensemble X lon X lat X time (all 2-D)
            output2{sce,1}=output{sce,1};
            output2{sce,3}=output{sce,5};

            if strcmp(model_name,'MRI_ESM2') % 'ta' is T850

                for i=1:1:length(var_name)
                    if strcmp(var_name{i},'ps')
                        ps_id=i;
                        break
                    end
                end

                if ~isempty(output{sce,6})
                    t850=output{sce,6}; 
                    t2=output2{sce,1}; % 'tas'
                    ps=output{sce,ps_id}/1e2; % surface pressure 'ps', the last variable in the list, convert to hPa
                    tmp=t2+(925-ps)./(850-ps).*(t850-t2);
                    tmp(ps<925)=nan;
                    output2{sce,4}=tmp;
                end
            else % 'ta' is already T925
                output2{sce,4}=output{sce,6};
            end
    
            if read_uava==1 % ua and va is included for scenario_group=1
                if ~isempty(output{sce,uava_id(1)}) && ~isempty(output{sce,uava_id(2)})
                    output2{sce,5}=sqrt(output{sce,uava_id(1)}.^2+output{sce,uava_id(2)}.^2); % calculate upper Wind
                end
            end
        end  
        output=output2;
    end

    if isempty(file_read)
        disp('<<<<< Warning: no file read!')
    end
    
end