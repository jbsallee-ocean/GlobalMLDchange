addpath('/net/argos/data/peps/jbslod/Data/Routine-Matlab/Fast_marching/toolbox_fast_marching/')

do_bathy=1;

if do_bathy==1
% % %% % get subsample of world
 latvec = -89.9:0.5:89.9;
 lonvec = 0:0.5:359.75;

 [LON LAT] = meshgrid(lonvec,latvec);
 bat_data = int16(LON'.*NaN);
 bat_data_mean = int16(LON'.*NaN);
 for latstep = 1:length(latvec)
   disp(num2str(latvec(latstep)))
   for lonstep = 1:length(lonvec)
 
       if abs(latvec(latstep)) == 90
           tmp = 360;
           lonmin = -180;
           lonmax = 180;
       else
       %    tmp = 0.17/cosd(latvec(latstep));
            tmp = 0.27/cosd(latvec(latstep));
           lonmin = lonvec(lonstep)-tmp;
           lonmax = lonvec(lonstep)+tmp;
	   if lonmin>180 lonmin=lonmin-360; end
	   if lonmax>180 lonmax=lonmax-360; end
 	   if lonmin<-180 lonmin=lonmin+360; end
	   if lonmax<-180 lonmax=lonmax+360; end
       end
 	lonmin_tmp=min(lonmax,lonmin);
	lonmax_tmp=max(lonmin,lonmax);
	lonmin=lonmin_tmp;
	lonmax=lonmax_tmp;
 
       latmin = max([latvec(latstep)-0.27;-89.95]);
       latmax = min([latvec(latstep)+0.27;89.95]);

       [image_data,vlon,vlat] = extract_gebco(lonmin,lonmax,latmin,latmax);

      %if sum(image_data>-1) > numel(image_data)/5
      if sum(image_data(:)>-1) > numel(image_data(:))/1.25
          bat_data(lonstep,latstep) = 10;
      else
          bat_data(lonstep,latstep) = double(round(median(image_data(image_data(:)<0))));
      end
    end
 end


 %save bat_data_1_etopo.mat bat_data LON LAT
 %save bat_data_0.25_etopo_0.66_land.mat bat_data LON LAT
 save bat_data_0.5_gebco_land.mat bat_data LON LAT% %% % get subsample of world
end

% 
% %%


%%
close all
clear all
%
load bat_data_0.5_gebco_land.mat bat_data
bat_data = bat_data';
bat_dataC = bat_data(:);
clear bat_data

load bat_data_0.5_gebco_land.mat
longrid =LON(1,:);
latgrid =LAT(:,1);

LAT = LAT(:);
LON = LON(:);
bat_data = bat_data';
bat_data = bat_data(:);
ii = find(LON>358);
ii2 = find(LON<2);
LAT = [LAT; LAT(ii); LAT(ii2)];
LON = [LON; LON(ii)-360; LON(ii2)+360];

bat_data = [bat_data; bat_data(ii); bat_data(ii2)];
bat_dataC = [bat_dataC; bat_dataC(ii); bat_dataC(ii2)];

%

F = TriScatteredInterp(double(LON(:)),double(LAT(:)),double(bat_data(:)));
FC = TriScatteredInterp(double(LON(:)),double(LAT(:)),double(bat_dataC(:)));


%
%bat_data = bat_data;
% ==============================================================
%mappingdepth = -1000;
mappingdepth =  0;
FM_fname = 'SO_Gebco_FM_05.nc'
%FM_fname = 'global_0.25_etopo_0.66land.nc';
%FM_fname = 'global_FM_shelf030_0.66land.nc';
%FM_fname = 'test';
% ==============================================================

%max dist in latitudinal degrees

maxdist = 30;

radius = sqrt(((2*maxdist)^2)/2);
close all
m_proj('Azimuthal Equal-area','lat',[-45],'long',[180],'rad',radius,'rect','on');
m_coast
m_grid
scales = get(gca,'xlim');

nsteps_coarse = (2*round(maxdist)+1);
grid_vector_coarse = linspace(scales(1),scales(2),nsteps_coarse);

maxstep = 1/4;
small_scale_reduction = 4;
nsteps = (2*round((maxdist/small_scale_reduction)/maxstep)+1);
grid_vector = linspace(scales(1)/small_scale_reduction,scales(2)/small_scale_reduction,nsteps);

max_dist_out_C = 25*110; % in km !!!!
max_dist_out = 7.5*110;


%latgrid =  61%-80:0.5:90;
%longrid =  180%0:0.5:360;

%max_real_dist = distance_km(0-maxdist,0, [0+maxdist,0]);
[lon1 lat1] = m_xy2ll([0,0],[grid_vector(1) grid_vector(end)]);
max_real_dist = m_idist(lon1(1),lat1(1),lon1(2),lat1(2))/1000;
[lon1 lat1] = m_xy2ll([0,0],[grid_vector_coarse(1) grid_vector_coarse(end)]);
max_real_dist_coarse = m_idist(lon1(1),lat1(1),lon1(2),lat1(2))/1000;
%
withinone = 1.2*1/maxstep;  % currently within 2 degrees ... 8 grid points

mlength = nsteps;
mlength_c = nsteps_coarse;
m2length = median(1:nsteps);
m2length_c = median(1:nsteps_coarse);

maxgridpoints = 3000;% 2400 gridpoinst with maxstep 1/4 == 1,500,000 sqkm == circle 691 km radius
% that should always cover enough area!!
% 3421 == circle radius 825km == 7.5?latitude as original
% do iterative export !!!!!!!!

betaplane = zeros(mlength,mlength);
[betaplane_lat_xy betaplane_lon_xy]= meshgrid(grid_vector,grid_vector);
betaplane_z = zeros(mlength,mlength);
D_dummy = NaN(mlength,mlength) + max_dist_out + 1;
betaplane_f = zeros(mlength,mlength);
betaplane_spd = zeros(mlength,mlength);

Cbetaplane = zeros(mlength_c,mlength_c);
[Cbetaplane_lat_xy Cbetaplane_lon_xy]= meshgrid(grid_vector_coarse,grid_vector_coarse);
Cbetaplane_z = zeros(mlength_c,mlength_c);
CD_dummy = NaN(mlength_c,mlength_c) + max_dist_out + 1;
Cbetaplane_f = zeros(mlength_c,mlength_c);
Cbetaplane_spd = zeros(mlength_c,mlength_c);

start_point = [(mlength+1)/2;(mlength+1)/2];
Cstart_point = [(mlength_c+1)/2;(mlength_c+1)/2];

m2complex =m2length + 1i * m2length;

sttime = now;
%%


%  initialize output netcdf file
% write to a LARGE netcdf file:::::
ncid = netcdf.create(FM_fname,'NC_64BIT_OFFSET');

%Define the dimensions of the variable.
dimid_gridpoints          = netcdf.defDim(ncid,'maingrid',0);
dimid_max_num_gridpoint   = netcdf.defDim(ncid,'subgrid',maxgridpoints-200);
dimid_max_num_gridpoint_C = netcdf.defDim(ncid,'subgrid_coarse',maxgridpoints-500);

% Define a new variable in the file.
latID  = netcdf.defVar(ncid,'lat','float',dimid_gridpoints);
lonID  = netcdf.defVar(ncid,'lon','float',dimid_gridpoints);
dptID  = netcdf.defVar(ncid,'depth','short',dimid_gridpoints);
netcdf.putAtt(ncid,dptID,'_FillValue',int16(9999))
dptID2  = netcdf.defVar(ncid,'depth_C','short',dimid_gridpoints);
netcdf.putAtt(ncid,dptID,'_FillValue',int16(9999))

SUBdisID = netcdf.defVar(ncid,'sub_distance','short',[dimid_max_num_gridpoint dimid_gridpoints]);
netcdf.putAtt(ncid,SUBdisID,'_FillValue',int16(9999))
SUBangID  = netcdf.defVar(ncid,'sub_angle','short',[dimid_max_num_gridpoint dimid_gridpoints]);
netcdf.putAtt(ncid,SUBangID,'_FillValue',int16(9999))
SUBlatID = netcdf.defVar(ncid,'sub_lat','float',[dimid_max_num_gridpoint dimid_gridpoints]);
netcdf.putAtt(ncid,SUBlatID,'_FillValue',NaN('single'))
SUBlonID  = netcdf.defVar(ncid,'sub_lon','float',[dimid_max_num_gridpoint dimid_gridpoints]);
netcdf.putAtt(ncid,SUBlonID,'_FillValue',NaN('single'))

SUB2disID = netcdf.defVar(ncid,'sub_distance_coarse','short',[dimid_max_num_gridpoint_C dimid_gridpoints]);
netcdf.putAtt(ncid,SUB2disID,'_FillValue',int16(9999))
SUB2angID  = netcdf.defVar(ncid,'sub_angle_coarse','short',[dimid_max_num_gridpoint_C dimid_gridpoints]);
netcdf.putAtt(ncid,SUB2angID,'_FillValue',int16(9999))
SUB2latID = netcdf.defVar(ncid,'sub_lat_coarse','float',[dimid_max_num_gridpoint_C dimid_gridpoints]);
netcdf.putAtt(ncid,SUB2latID,'_FillValue',NaN('single'))
SUB2lonID  = netcdf.defVar(ncid,'sub_lon_coarse','float',[dimid_max_num_gridpoint_C dimid_gridpoints]);
netcdf.putAtt(ncid,SUB2lonID,'_FillValue',NaN('single'))

% Leave define mode and enter data mode to write data.
netcdf.endDef(ncid)


nstep = 0;


%%
for ystep = 1:length(latgrid)
    disp(['starting ' num2str(latgrid(ystep)) '?N at ' datestr(now,'HH:MM') ' ---- ' datestr(now-sttime,'HH:MM:SS') ' elapsed'])
    
    
    for xstep = 1:1:length(longrid)
        
        testcase = F(longrid(xstep),latgrid(ystep));
        
        %if testcase < -5
        if testcase < mappingdepth
            
            m_proj('Azimuthal Equal-area','lat',latgrid(ystep),'long',longrid(xstep),'rad',radius,'rect','on');
            
            [betaplane_lon,betaplane_lat] = m_xy2ll(betaplane_lon_xy,betaplane_lat_xy);
            betaplane_lon(betaplane_lon< 0) = betaplane_lon(betaplane_lon <0)+360;
            betaplane_lon(betaplane_lon>=360) = betaplane_lon(betaplane_lon>=360)-360;
            
            [Cbetaplane_lon,Cbetaplane_lat] = m_xy2ll(Cbetaplane_lon_xy,Cbetaplane_lat_xy);
            Cbetaplane_lon(Cbetaplane_lon<0) = Cbetaplane_lon(Cbetaplane_lon<0)+360;
            Cbetaplane_lon(Cbetaplane_lon>=360) = Cbetaplane_lon(Cbetaplane_lon>=360)-360;
            %
            betaplane_z = F(betaplane_lon,betaplane_lat);%-(mappingdepth + 250);
            Cbetaplane_z = FC(Cbetaplane_lon,Cbetaplane_lat);%-(mappingdepth + 250);
            % betaplane_z(betaplane_z<-6500-(mappingdepth + 250))= -6500-(mappingdepth + 250);
            % Cbetaplane_z(Cbetaplane_z<-6500-(mappingdepth + 250))= -6500-(mappingdepth + 250);
            betaplane_z(betaplane_z<-5500)= -5500;
            Cbetaplane_z(Cbetaplane_z<-5500)= -5500;
            
            betaplane_z(betaplane_z<-5 & betaplane_z >-75) = -75;
            Cbetaplane_z(Cbetaplane_z<-5 & Cbetaplane_z >-75) = -75;
            
            Cbetaplane_z(Cbetaplane_z>0) = 0;
            betaplane_z(betaplane_z>0) = 0;
            
            % for bathymetry data out of bounds (N of 80?N:
            betaplane_z(isnan(betaplane_z)) = 0;
            Cbetaplane_z(isnan(Cbetaplane_z)) = 0;
            
            betaplane_z = abs(betaplane_z);
            Cbetaplane_z = abs(Cbetaplane_z);
            %
            % adjust grid and assign
            % if betaplane_z(m2length,m2length) < -20
            %create speedmap using f/H
            
            % general boundary shape:
            betaplane_spd = (1-abs(log(sqrt((betaplane_z(m2length,m2length)))./sqrt((betaplane_z)))));% .* (10/betaplane_z(m2length,m2length)))));
            % betaplane_spd(betaplane_z>betaplane_z(m2length,m2length)) = betaplane_spd(betaplane_z>betaplane_z(m2length,m2length)).^3;
            Cbetaplane_spd = (1-abs(log(sqrt(Cbetaplane_z(m2length_c,m2length_c))./sqrt(Cbetaplane_z))));% .* (10/betaplane_z(m2length,m2length)))));
            % Cbetaplane_spd(Cbetaplane_z>Cbetaplane_z(m2length_c,m2length_c)) = Cbetaplane_spd(Cbetaplane_z>Cbetaplane_z(m2length_c,m2length_c)).^3;
            %
            %  betaplane_spd(betaplane_spd < 0.04) = 0.04;
            % add exception fo shallow water
            %betaplane_spd(betaplane_spd < 0.2 & (abs(betaplane_z - betaplane_z(m2length,m2length))) < 100) = 0.2;
            
            % add equatorial speeds:
            betaplane_spd = betaplane_spd.* exp(-abs(betaplane_lat - betaplane_lat(m2length,m2length))./exp(betaplane_lat(m2length,m2length).^2/56.25)); % 56.25 = 7.5 ^ 2
            betaplane_spd(betaplane_spd < 0.05 ) = 0.05;
            Cbetaplane_spd = Cbetaplane_spd.* exp(-abs(Cbetaplane_lat - Cbetaplane_lat(m2length_c,m2length_c))./exp(Cbetaplane_lat(m2length_c,m2length_c).^2/56.25)); % 56.25 = 7.5 ^ 2
            Cbetaplane_spd(Cbetaplane_spd < 0.05 ) = 0.05;
            
            % remove minor speeds and NAN for arithmetics:
            % betaplane_spd(betaplane_spd<=1e-12) = 1e-4;
            betaplane_spd(isnan(betaplane_spd)) = 1e-3;
            Cbetaplane_spd(isnan(Cbetaplane_spd)) = 1e-3;
            
            %
            %set options and perform fast marching
            % clear options;
            
            % CHANGED TO 3e3 from 1e6 from 1e7 .... HOPE THAT IS ENOUGH  ----- CHECK CHECK CHECK
            options.nb_iter_max = 2.7e3;
            options.constraint_map = ones(size(betaplane_spd));
            options.constraint_map(betaplane_z < 0.5) = -Inf;
            options.dmax = 0.5;
            
           
            [D,S] = perform_fast_marching(double(betaplane_spd), start_point, options);
            
            %  clear options;
            % CHANGED TO 3e3 from 1e6 from 1e7 .... HOPE THAT IS ENOUGH  ----- CHECK CHECK CHECK
            options.nb_iter_max = 2.7e3;
            options.constraint_map = ones(size(Cbetaplane_spd));
            options.constraint_map(Cbetaplane_z < 0.5) = -Inf;
            options.dmax = 0.5;
            
            
            [DC,SC] = perform_fast_marching(double(Cbetaplane_spd), Cstart_point, options);
            %
            
            % keyboard
            D2 = D.*max_real_dist;
            D2C = DC.*max_real_dist_coarse;
            
            D2(D2>max_dist_out) = NaN;
            D2C(D2C>max_dist_out_C) = NaN;
            D2F = D2;
            D2CF = D2C;
            %     end
            % end
            % clf
            % patch(ncst(:,1)+360,ncst(:,2),[0.4 0.4 0.4])
            % hold on
            % patch(ncst(:,1),ncst(:,2),[0.4 0.4 0.4])
            % lm_ang_h = patch(betaplane_lon(:),betaplane_lat(:),D2(:),'Marker','sq','MarkerFaceColor','flat','FaceColor','none','LineStyle','none','markeredgecolor','flat');
            % colorbar
            % set(gca,'xlim',[min(betaplane_lon(:)) max(betaplane_lon(:))],'ylim',[min(betaplane_lat(:)) max(betaplane_lat(:))])
            %
            Dnisnan = ~isnan(D2(:));
            SDnisnin = sum(Dnisnan);
            % Dist_CELL{ystep,xstep} = int16(round([betaplane_lon(~isnan(D2(:))), betaplane_lat(~isnan(D2(:))), (D(~isnan(D2(:))))]));
            
            % now the hard part: get direction.
            Dsiz = size(D);
            D_angle = NaN(Dsiz);
            ct = 0;
            D_angle(m2length,m2length)=0;
            %keyboard
            if sum(D(:)<1e3) > 9
               tmpsus = isinf(D);
               A1 = double(D); A1(tmpsus) = double(max(D(~tmpsus)));
               grad = compute_grad(A1);
               grad = -perform_vf_normalization(grad);
                    
                while true
                    ct = ct + 1;
                    [tmp In]=max(D2(:));
                    if isnan(tmp)
                        break
                    end
                    [I,J]=ind2sub(Dsiz,In);
                    
                   % gpath = compute_geodesic(single(D),[I; J]);
              
                    path = stream2(grad(:,:,2),grad(:,:,1),J,I, [0.1 1500]);
                    for i=1:length(path)
                        path{i} = path{i}(:,[2 1]);
                    end
                    
                    gpath = path{1}';
                    plot(gpath(1,:)+0.5,gpath(2,:)+0.5,'-k','linewidth',2)
                    if max(size(gpath))>2
                        g2m2 = (gpath(2,:)-m2length);
                        g1m2 = (gpath(1,:) - m2length);
                        gpath_imag = sqrt(g2m2.*g2m2 + g1m2.*g1m2);
                        gef_last = find((gpath_imag) > withinone ,1,'last');
                        if isempty(gef_last)
                            [tmp gef_last] = max((gpath_imag));
                        end
                        %pr_ang=angle((gpath(2,gef_last) - m2length)+ 1i * (gpath(1,gef_last)  - m2length));
                        %  pr_ang = atan2((gpath(1,gef_last)  - m2length), (gpath(2,gef_last) -m2length));
                        pr_ang = atan2(g1m2(gef_last), g2m2(gef_last));
                        
                        %address that angle to all boxes line passes through
                      %	 D2(sub2ind_fast_no_error_checks(Dsiz,round(gpath(1,:)),round(gpath(2,:)))) = NaN;		%SUNKE VERSION
                      %  D_angle(sub2ind_fast_no_error_checks(Dsiz,round(gpath(1,:)),round(gpath(2,:)))) = pr_ang;	%SUNKE VERSION

 			D2(sub2ind(Dsiz,round(gpath(1,:)),round(gpath(2,:)))) = NaN;
                        D_angle(sub2ind(Dsiz,round(gpath(1,:)),round(gpath(2,:)))) = pr_ang;


                    else
                        if max(size(gpath))==2
                            gpath = gpath';
                        end
                        pr_ang = atan2((gpath(1)  - m2length), (gpath(2) -m2length));
                       % D_angle(sub2ind_fast_no_error_checks(Dsiz,round(gpath(1)),round(gpath(2)))) = pr_ang;	%SUNKE VERSION
                        
                       % D2(sub2ind_fast_no_error_checks(Dsiz,round(gpath(1)),round(gpath(2)))) = NaN;		%SUNKE VERSION
                        

 			D_angle(sub2ind(Dsiz,round(gpath(1)),round(gpath(2)))) = pr_ang;
                        
                        D2(sub2ind(Dsiz,round(gpath(1)),round(gpath(2)))) = NaN;
                        
                    end
                end
            end
            
            D_angle = round(D_angle/pi*180)-90;
            D_angle(m2length,m2length)=0;
            D_angle(m2length-1,m2length)=180;
            D_angle(m2length+1,m2length)=0;
            D_angle(m2length,m2length-1)=+90;
            D_angle(m2length,m2length+1)=-90;
            D_angle(m2length-1,m2length-1)=+135;
            D_angle(m2length+1,m2length+1)=-45;
            D_angle(m2length-1,m2length+1)=-135;
            D_angle(m2length+1,m2length-1)=+45;
            D_angle(D_angle<-180) = D_angle(D_angle<-180)+360;
            
            %
            %
            DCnisnan = ~isnan(D2C(:));
            SDCnisnin = sum(DCnisnan);
            % Dist_CELL{ystep,xstep} = int16(round([betaplane_lon(~isnan(D2(:))), betaplane_lat(~isnan(D2(:))), (D(~isnan(D2(:))))]));
            
            % now the hard part: get direction.
            CDsiz = size(DC);
            CD_angle = NaN(CDsiz);
            ct = 0;
            CD_angle(m2length,m2length)=0;
            if sum(DC(:)<1e3) > 90
               tmpsus = isinf(DC);
               A1 = double(DC); A1(tmpsus) = double(max(DC(~tmpsus)));
               grad = compute_grad(A1);
               grad = -perform_vf_normalization(grad);
                    
                while true
                    ct = ct + 1;
                    [tmp In]=max(D2C(:));
                    if isnan(tmp)
                        break
                    end
                    [I,J]=ind2sub(CDsiz,In);
                    %  keyboard
                   % gpath = compute_geodesic(double(DC),[I; J]);
              
                    path = stream2(grad(:,:,2),grad(:,:,1),J,I, [0.1 1500]);
                    for i=1:length(path)
                        path{i} = path{i}(:,[2 1]);
                    end
                    
                    gpath = path{1}';

                    if max(size(gpath))>2
                        
                        g2m2 = (gpath(2,:)-m2length_c);
                        g1m2 = (gpath(1,:) - m2length_c);
                        gpath_imag = sqrt(g2m2.*g2m2 + g1m2.*g1m2);
                        gef_last = find((gpath_imag) > withinone ,1,'last');
                        if isempty(gef_last)
                            [tmp gef_last] = max((gpath_imag));
                        end
                        %pr_ang=angle((gpath(2,gef_last) - m2length)+ 1i * (gpath(1,gef_last)  - m2length));
                        %  pr_ang = atan2((gpath(1,gef_last)  - m2length), (gpath(2,gef_last) -m2length));
                        pr_ang = atan2(g1m2(gef_last), g2m2(gef_last));
                        
                        %address that angle to all boxes line passes through	
                     %   D2C(sub2ind_fast_no_error_checks(Dsiz,round(gpath(1,:)),round(gpath(2,:)))) = NaN;		%SUNKE VERSION
                     %   CD_angle(sub2ind_fast_no_error_checks(Dsiz,round(gpath(1,:)),round(gpath(2,:)))) = pr_ang;	%SUNKE VERSION

 			D2C(sub2ind(Dsiz,round(gpath(1,:)),round(gpath(2,:)))) = NaN;
                        CD_angle(sub2ind(Dsiz,round(gpath(1,:)),round(gpath(2,:)))) = pr_ang;


                    else
                        if max(size(gpath))==2
                            gpath = gpath';
                        end
                        pr_ang = atan2((gpath(1)  - m2length_c), (gpath(2) -m2length_c));
                     %   CD_angle(sub2ind_fast_no_error_checks(Dsiz,round(gpath(1)),round(gpath(2)))) = pr_ang;		%SUNKE VERSION
                        
                     %   D2C(sub2ind_fast_no_error_checks(Dsiz,round(gpath(1)),round(gpath(2)))) = NaN; 		%SUNKE VERSION


			CD_angle(sub2ind(Dsiz,round(gpath(1)),round(gpath(2)))) = pr_ang;
                        
                        D2C(sub2ind(Dsiz,round(gpath(1)),round(gpath(2)))) = NaN;
                        
                        
                    end
                end
            end
            
            CD_angle = round(CD_angle/pi*180)-90;
            CD_angle(m2length_c,m2length_c)=0;
            CD_angle(m2length_c-1,m2length_c)=180;
            CD_angle(m2length_c+1,m2length_c)=0;
            CD_angle(m2length_c,m2length_c-1)=+90;
            CD_angle(m2length_c,m2length_c+1)=-90;
            CD_angle(m2length_c-1,m2length_c-1)=+135;
            CD_angle(m2length_c+1,m2length_c+1)=-45;
            CD_angle(m2length_c-1,m2length_c+1)=-135;
            CD_angle(m2length_c+1,m2length_c-1)=+45;
            CD_angle(CD_angle<-180) = CD_angle(CD_angle<-180)+360;
            %
            %keyboard
            % safe data to netcdf file:
            
            Dnisnan = find(Dnisnan);
            DCnisnan = find(DCnisnan);
            if numel(Dnisnan)>2800
                Dnisnan = Dnisnan(1:2800);
                SDnisnin = 2800;
            end
            if numel(DCnisnan)>2500
                DCnisnan = DCnisnan(1:2500);
                SDCnisnin = 2500;
            end
            
            % write MAIN GRID
            netcdf.putVar(ncid,latID,nstep,latgrid(ystep));
            netcdf.putVar(ncid,lonID,nstep,longrid(xstep));
            netcdf.putVar(ncid,dptID,nstep,int16(testcase));
            netcdf.putVar(ncid,SUBdisID,[0 nstep ],[SDnisnin 1],int16(D2F(Dnisnan)));
            netcdf.putVar(ncid,SUBangID,[0 nstep ],[SDnisnin 1],int16(D_angle(Dnisnan)));
            netcdf.putVar(ncid,SUBlatID,[0 nstep ],[SDnisnin 1],betaplane_lat(Dnisnan));
            netcdf.putVar(ncid,SUBlonID,[0 nstep ],[SDnisnin 1],betaplane_lon(Dnisnan));
            netcdf.putVar(ncid,SUB2disID,[0 nstep ],[SDCnisnin 1],int16(D2CF(DCnisnan)));
            netcdf.putVar(ncid,SUB2angID,[0 nstep ],[SDCnisnin 1],int16(CD_angle(DCnisnan)));
            netcdf.putVar(ncid,SUB2latID,[0 nstep ],[SDCnisnin 1],Cbetaplane_lat(DCnisnan));
            netcdf.putVar(ncid,SUB2lonID,[0 nstep ],[SDCnisnin 1],Cbetaplane_lon(DCnisnan));
            nstep = nstep + 1;
        end
    end
end
% instead of close - just sync here

%%
netcdf.close(ncid)

%
disp('##>> it took:')
disp(['##>> ' datestr(now- sttime,'HH:MM:SS')])
disp('##>> to compute global Fast Marching ')
disp('##>> it is now: ')
disp(['##>> ' datestr(now,'HH:MM')])
disp('##>> have a nice day. ')
%save Dist_CELL Dist_CELL LATgrid LONgrid Lat_CELL Lon_CELL
%save Dist_CELL Dist_CELL LATgrid LONgrid
disp('THE END')
