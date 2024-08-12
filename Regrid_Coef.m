function [regrid_grid,regrid_fraction] = Regrid_Coef(nlon,nlat,rlon,rlat,regrid_fine_scale,method)

    % Mi Zhou @ Princeton, 2022/12
    % NOTE for this function:
    % 【Attention】This regrid function does not apply for Lambert grid
    % {nlon,nlat} is the original coordinates (2-D)
    % {rlon,rlat} is the regridded coordinates (2-D)
    % {fine_lon,fine_lat} is the intermediate very fine coordinates (1-D)
    % We estimate the fraction of {nlon,nlat} that goes into {rlon,rlat} via {fine_lon,fine_lat}
    % We assume that the grid space of {nlon,nlat} is larger than the grid space of {rlon,rlat}
    % method=1, regrid from flux to flux; method=2, regrid from mass to mass
    
    reso_lon_model=nlon(2,1)-nlon(1,1);
    reso_lat_model=nlat(1,2)-nlat(1,1);
    reso_lon_regrid=rlon(2,1)-rlon(1,1);
    reso_lat_regrid=rlat(1,2)-rlat(1,1);
    % (x1-dx/2+dx/2n):dx/n:(xl+dx/2-dx/2n)
    % An example to divide grid C into four (regrid_fine_scale=2) grid c:
    % Devide {nlon,nlat} into {fine_lon,fine_lat}:
    %-----------------
    %|   |   |   |   | 
    %|---c-------c---|
    %|   |   |   |   |
    %--------C-------|
    %|   |   |   |   |
    %|---c-------c---|
    %|   |   |   |   |
    %------------------
    fine_lon=min(min(nlon))-reso_lon_model/2+reso_lon_model/2/regrid_fine_scale:reso_lon_model/regrid_fine_scale:max(max(nlon))+reso_lon_model/2-reso_lon_model/2/regrid_fine_scale; 
    fine_lat=min(min(nlat))-reso_lat_model/2+reso_lat_model/2/regrid_fine_scale:reso_lat_model/regrid_fine_scale:max(max(nlat))+reso_lat_model/2-reso_lat_model/2/regrid_fine_scale;
   
    % dimension of the regridded
    regrid_grid=cell(size(rlon));
    regrid_fraction=cell(size(rlon));

    for i=1:1:length(fine_lon)
        for j=1:1:length(fine_lat)

            % [i2raw,j2raw] is the grid in {nlon,nlat}(the original grid space) this fine_grid belongs to
            i2raw=ceil(i/regrid_fine_scale);
            j2raw=ceil(j/regrid_fine_scale);
%             yingshe1{i,j}=[i2raw,j2raw];
    
            % [x,y] is the grid in {rlon,rlat}(the new grid space) this fine_grid belongs to
            [x,y]=findpoint(fine_lon(i),fine_lat(j),rlon,rlat);

            % if this fine_grid (covered by {nlon,nlat} but not covered by {rlon,rlat}) is too far away from the grid in {rlon,rlat}
            % we don't put this fine_grid into {rlon,rlat}
            if abs(rlon(x,y)-fine_lon(i))>=(reso_lon_regrid/2+reso_lon_model/regrid_fine_scale/2) || abs(rlat(x,y)-fine_lat(j))>=(reso_lat_regrid/2+reso_lat_model/regrid_fine_scale/2)
%                 yingshe2{i,j}=[];
                continue
            end
%             yingshe2{i,j}=[x,y];

            new_grid=1;
            for k=1:1:size(regrid_grid{x,y},1)
                % fraction of grid [i2raw,j2raw] in the original coordinates has already been pointered to grid [x,y] in the regridded coordinates
                % we add one more 'weight' to [i2raw,j2raw]
                if regrid_grid{x,y}(k,1)-i2raw==0 && regrid_grid{x,y}(k,2)-j2raw==0
                    regrid_fraction{x,y}(k)=regrid_fraction{x,y}(k)+1;
                    new_grid=0;
                    break
                end
            end

            % fraction of grid [i2raw,j2raw] in the original coordinates has not been pointered to grid [x,y] in the regridded coordinates
            % we set the current weight to '1'
            if new_grid==1
                % regrid_grid archives the grids in the original coordinates that overlap with the grids in the regirdded coordinates
                regrid_grid{x,y}(size(regrid_grid{x,y},1)+1,1:2)=[i2raw,j2raw]; 
                % regrid_fraction archives the 'weight' of the grids in the original coordinates that overlap with the grids in the regirdded coordinates
                regrid_fraction{x,y}(size(regrid_fraction{x,y},1)+1,1)=1;
            end

        end
    end

    for i=1:1:size(regrid_fraction,1)
        for j=1:1:size(regrid_fraction,2)
            if method==1 % flux, regridded flux is determined by weighted-mean original fluxes

                % for example, if there are 4 original grids going into the regridded grids, and their 'weight' is [1,1,1,1],
                % we estimate their individual contribution to the regridded grid to be [0.25,0.25,0.25,0.25] by using the code below:
                regrid_fraction{i,j}=regrid_fraction{i,j}/sum(regrid_fraction{i,j});

            elseif method==2 % emission/mass, regridded mass is determined by the sum of fraction of original mass

                % for example, if there are 4 original grids going into the regridded grids, and their 'weight' is [1,1,1,1],
                % we estimate the fraction of the original grids that go into the regridded grid to be [1,1,1,1]/regrid_fine_scale^2,
                % note that each original grid is divided to regrid_fine_scale^2 grids
                regrid_fraction{i,j}=regrid_fraction{i,j}/regrid_fine_scale^2;

            end
        end
    end

end