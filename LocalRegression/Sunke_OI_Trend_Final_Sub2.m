






function out=Sunke_OI_Trend_Final_Sub2(INparaXX,varINlist,jstep,istep)
out.istep=NaN;
out.jstep=NaN;
	for ii=1:length(varINlist)
		eval([ varINlist{ii} '=INparaXX.' varINlist{ii}  ';']);
	end
	clear INparaXX

            if x(istep,jstep) >= continue_lon_at && x(istep,jstep) <= continue_lon_upto
                % keeping track of current executions:
                %fprintf(['\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b', num2str(x(istep,jstep),'%05.1f') ' after ' datestr(now-tstart,'HH:MM:SS')]) %sunke version
                
                % Find Fast marching grid point
                FM_pos = find(FM_LAT == y(istep,jstep) & ( FM_LON == x(istep,jstep) | FM_LON == x(istep,jstep)-360));
                %disp(['FM: ' num2str(y(istep,jstep)) ';' num2str(x(istep,jstep))])
                if ~isempty(FM_pos)
                    FM_pos = FM_pos(1);
                    %disp('CA PASSE')
                    if upperocean
                        local_depth = min([upperocean_limit, abs(FM_DPT(FM_pos))]);
                    else
                        local_depth = abs(FM_DPT(FM_pos)); %#ok<UNRCH>
                    end

%%Trying to paralelize%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reopen because I don't know how to pass netcdf-id to a function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncid_FM = netcdf.open(FM_fname,'NOWRITE');
FM_id_LON = netcdf.inqVarID(ncid_FM,'lon');
FM_id_LAT = netcdf.inqVarID(ncid_FM,'lat');
FM_id_DPT = netcdf.inqVarID(ncid_FM,'depth');
%FM_id_NUM = netcdf.inqVarID(ncid_FM,'num_subpoints');
FM_id_SUBDIS = netcdf.inqVarID(ncid_FM,'sub_distance');
FM_id_SUBANG = netcdf.inqVarID(ncid_FM,'sub_angle');
FM_id_SUBLAT = netcdf.inqVarID(ncid_FM,'sub_lat');
FM_id_SUBLON = netcdf.inqVarID(ncid_FM,'sub_lon');

FM_id_SUBDIS_C = netcdf.inqVarID(ncid_FM,'sub_distance_coarse');
FM_id_SUBANG_C = netcdf.inqVarID(ncid_FM,'sub_angle_coarse');
FM_id_SUBLAT_C = netcdf.inqVarID(ncid_FM,'sub_lat_coarse');
FM_id_SUBLON_C = netcdf.inqVarID(ncid_FM,'sub_lon_coarse');
%%Trying to paralelize%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END Reopen because I don't know how to pass netcdf-id to a function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



                    % load fast marching subgrid
                    fm_dis = double(netcdf.getVar(ncid_FM,FM_id_SUBDIS,[0 FM_pos-1],[FM_ALL 1]));
                    
                    % keyboard
                    
                    gefFM = find(fm_dis<2000);
                    
                    if numel(gefFM) > 10
                        
                        % get lon from FM grid
                        fm_lat = double(netcdf.getVar(ncid_FM,FM_id_SUBLAT,[0 FM_pos-1],[FM_ALL 1]));
                        fm_lon = double(netcdf.getVar(ncid_FM,FM_id_SUBLON,[0 FM_pos-1],[FM_ALL 1]));
                        fm_ang = double(netcdf.getVar(ncid_FM,FM_id_SUBANG,[0 FM_pos-1],[FM_ALL 1]));
                        
                        % load coarse grid:
                        fm_lat_c = double(netcdf.getVar(ncid_FM,FM_id_SUBLAT_C,[0 FM_pos-1],[FM_ALL_C 1]));
                        fm_lon_c = double(netcdf.getVar(ncid_FM,FM_id_SUBLON_C,[0 FM_pos-1],[FM_ALL_C 1]));
                        fm_dis_c = double(netcdf.getVar(ncid_FM,FM_id_SUBDIS_C,[0 FM_pos-1],[FM_ALL_C 1]));
                        fm_ang_c = double(netcdf.getVar(ncid_FM,FM_id_SUBANG_C,[0 FM_pos-1],[FM_ALL_C 1]));
netcdf.close(ncid_FM);

                        % use only valid data:
                        fm_dis = fm_dis(gefFM);
                        fm_ang = fm_ang(gefFM);
                        fm_lat = fm_lat(gefFM);
                        fm_lon = fm_lon(gefFM);
                        fm_lon(fm_lon<0) = fm_lon(fm_lon<0) +360;
                        max_small_FM = max(fm_dis);
                        gefFMc = find(fm_dis_c < 2800 & fm_dis_c >= min([max_small_FM 300]));
                        fm_dis_c = fm_dis_c(gefFMc);
                        fm_ang_c = fm_ang_c(gefFMc);
                        fm_lat_c = fm_lat_c(gefFMc);
                        fm_lon_c = fm_lon_c(gefFMc);
                        fm_lon_c(fm_lon_c<0) = fm_lon_c(fm_lon_c<0) +360;
                        %keyboard
                        
                        % take subset of all float/CTD data within range of FM data points
                        maxfmlat = max([max(fm_lat_c) max(fm_lat)]);
                        minfmlat = min([min(fm_lat_c) min(fm_lat)]);
                        maxfmlon = max([max(fm_lon_c) max(fm_lon)]);
                        minfmlon = min([min(fm_lon_c) min(fm_lon)]);
                        
                        fm_fit = find(abs(lat-(maxfmlat+minfmlat)/2) < (maxfmlat-minfmlat)/2 & ...
                            abs(lon-(maxfmlon+minfmlon)/2) < (maxfmlon-minfmlon)/2);
                        
                        % check if enough data is found
                        if ~isempty(fm_fit) && numel(fm_fit) > min_num_profs
                               %keyboard
                            
                            tmp_fm_lon = fm_lon-x(istep,jstep);
                            tmp_fm_lon(tmp_fm_lon>180) = tmp_fm_lon(tmp_fm_lon>180) - 360;
                            tmp_fm_lon(tmp_fm_lon<-180) = tmp_fm_lon(tmp_fm_lon<-180) + 360;
                            tmp_fm_lon = tmp_fm_lon.*cosd(fm_lat);
                            
                            
                            tmp_fm_lon_c = fm_lon_c-x(istep,jstep);
                            tmp_fm_lon_c(tmp_fm_lon_c>180) = tmp_fm_lon_c(tmp_fm_lon_c>180) - 360;
                            tmp_fm_lon_c(tmp_fm_lon_c<-180) = tmp_fm_lon_c(tmp_fm_lon_c<-180) + 360;
                            tmp_fm_lon_c = tmp_fm_lon_c.*cosd(fm_lat_c);
                            
                            % initialize irregular grid interpolant
                            DT = DelaunayTri(double(tmp_fm_lon),double(fm_lat));
                            
                            [PI,D] = nearestNeighbor(DT,tmp_fm_lon_c,fm_lat_c);
                            to_close = find(D<maxdist_finegrid);
                            tmp_fm_lon_c(to_close) = [];
                            fm_lat_c(to_close) = [];
                            fm_dis_c(to_close) = [];
                            fm_ang_c(to_close) = [];
                            
                            if numel(fm_lat_c) > 0
                                % extend Delaunay Triangulation by coarse gridpoints
                                DT.X(end+1 : end+numel(fm_lat_c) ,:) = [tmp_fm_lon_c fm_lat_c];
                            end
                            
                            %
                            if ~isempty(DT.Triangulation)
                                tmp_lon = lon(fm_fit)-x(istep,jstep);
                                tmp_lon(tmp_lon>180) = tmp_lon(tmp_lon>180) - 360;
                                tmp_lon(tmp_lon<-180) = tmp_lon(tmp_lon<-180) + 360;
                                tmp_lon = tmp_lon.*cosd(lat(fm_fit));
                                
                                % find nearest neighbour indices
                                [PI,D] = nearestNeighbor(DT,tmp_lon,lat(fm_fit));
                                
                                fm_dis_vec_combined = [fm_dis; fm_dis_c];
                                
                                % accept divergences up to half FM grid point
                                D_sub = find((D<=maxdist_finegrid+1 & fm_dis_vec_combined(PI) <= max_small_FM+50) | (D<=maxdist_coarsegrid+1+maxdist_finegrid & fm_dis_vec_combined(PI) >= max_small_FM-50));
                                D_subset = fm_fit(D_sub);
                                DD_sub = PI(D_sub);
                                % GRB CLOSE DATA HERE >>>>
                                if numel(D_subset) > min_num_profs &&  numel(D_sub) >= 10
                                    
                                    
                                    % do a more accurate distance computation with subset: (cannot be done earlier due to interpolation errors at bays, inlets, ridges ...
                                    TSI = TriScatteredInterp(DT,fm_dis_vec_combined);
                                    
                                    % in future get rid of one of these 2 dis vectors ... only keep compex one
                                    % put all data into complex dist_vec
                                    dist_vec = dist_dummy;
                                    dist_vec_complex = complex(dist_vec,0);
                                    
                                    % fill distance vector with distance
                                    tmp =  TSI(double(tmp_lon(D_sub)),lat(fm_fit(D_sub)))./110;
                                    gefu = find(isnan(tmp));
                                    if ~isempty(tmp)
                                        dist_vec(D_subset)  = tmp;
                                        gef = isnan(dist_vec(D_subset));
                                        D_subset(gef) = [];
                                        DD_sub(gef) = [];
                                        
                                        if ~isempty(D_subset)
                                            
                                            complex_angle_sin = sind([fm_ang; fm_ang_c]);
                                            complex_angle_cos = cosd([fm_ang; fm_ang_c]);
                                            TSI_sin = TriScatteredInterp(DT,complex_angle_sin);
                                            TSI_cos = TriScatteredInterp(DT,complex_angle_cos);
                                            tmp_sin = TSI_sin(double(tmp_lon(D_sub)),lat(fm_fit(D_sub)));
                                            tmp_cos = TSI_cos(double(tmp_lon(D_sub)),lat(fm_fit(D_sub)));
                                            fm_ang_vec_combined = atan2(tmp_sin,tmp_cos) *180/pi;
                                            fm_ang_vec_combined(gefu) = [];
                                            if numel(fm_ang_vec_combined) ~= numel(D_subset)
                                                disp('this should not happen!!!')
                                                keyboard
                                            end
                                            
                                            % add some noise to direction to prevent colinear data in bays, inlets ...etc
                                            %fm_ang_vec_combined(DD_sub) = fm_ang_vec_combined(DD_sub) + randi(3,numel(DD_sub),1) -2 ;
                                            
                                            %fill complex distance vector with data
                                            dist_vec_complex(D_subset)  = - sind(fm_ang_vec_combined) .* dist_vec(D_subset) + 1i .* cosd(fm_ang_vec_combined) .* dist_vec(D_subset);
                                            
                                            % find profiles within range
                                            ii_dist = D_subset;
                                            %dist_vec(ii_dist) = dist_vec(ii_dist)./maxrange;
                                            %dist_vec_complex(ii_dist) = dist_vec_complex(ii_dist);%./maxrange;
                                            
                                            
                                            % if more than MAX_PROFS_INMEM profiles are in close range, only use closest XXXX don't load more;
                                            if numel(ii_dist) > max_profs_inmem
                                                
                                                % find the closest profs .... BUT maxe sure in open ocean at least 500 CTD full depth profiles are included ...
                                                [~, kk]=sort(dist_vec(ii_dist));
                                                ii_dist_tmp = ii_dist(kk(1:max_profs_inmem));
                                                
                                                if sum(idv(kk(1:max_profs_inmem))<0) < min_CTD_profs
                                                    gef_CTD = find(idv(ii_dist) < 0);
                                                    numCTD = numel(gef_CTD);
                                                    
                                                    [~, kk2]=sort(dist_vec(ii_dist(gef_CTD)));
                                                    ii_dist2 = ii_dist(gef_CTD(kk2(1:min([min_CTD_profs numCTD]))));
                                                    
                                                    [~,IA,~] =setxor(ii_dist_tmp, ii_dist2);
                                                    ii_dist_tmp = ii_dist_tmp(IA);
                                                    
                                                    ii_dist_tmp = [ii_dist2; ii_dist_tmp(1:max_profs_inmem-numel(ii_dist2))];
                                                end
                                                
                                                ii_dist = ii_dist_tmp;
                                            end
                                            
                                            if drawit
                                                figure(1) %#ok<UNRCH>
                                                clf
                                                m_proj('Azimuthal Equal-area','lat',y(istep,jstep),'long',x(istep,jstep),'rad',20,'rect','on');
                                                m_plot(DT.X(1:4:end,1)./cosd(DT.X(1:4:end,2))+x(istep,jstep),DT.X(1:4:end,2),'.k')
                                                hold on
                                                m_plot(lon(ii_dist(1:10:end)),lat(ii_dist(1:10:end)),'+r')
                                                
                                                m_plot(x(istep,jstep),y(istep,jstep),'og')
                                                title('Every 4th FM point and every 10th data point (out of the up to 2000 closest)')
                                                m_coast
                                                m_grid
                                                drawnow
                                            end
                                            % gather all information and to load profiles not yet in memory...
                                            profs_needed = idv(ii_dist);
                                            
                                            % what to read and what not to read
                                            [~, AI, BI] = setxor(profs_needed,profs_inmem);	%
                                            
                                            % these we do need to read, since not yet in memory
                                            profs_toload = profs_needed(AI);
                                            
                                            % these are in memory but no longer needed
                                            profs_redundant = BI;
                                            
                                            % gather all indices that are right next to each other to load at once and not in multiple steps
                                            multi_load_vec = single(diff(profs_toload) == 1) ;
                                            
                                            % load all profiles from file
%%Trying to paralelize%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reopen because I don't know how to pass netcdf-id to a function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncid_CTD = netcdf.open(wod_ctd_data,'NOWRITE');	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CTD_file.nc
CTD_id_LON = netcdf.inqVarID(ncid_CTD,'lon');
CTD_id_LAT = netcdf.inqVarID(ncid_CTD,'lat');
CTD_id_DYR = netcdf.inqVarID(ncid_CTD,'dyr');
CTD_id_SIGMA = netcdf.inqVarID(ncid_CTD,'sigma');
CTD_id_TEMP = netcdf.inqVarID(ncid_CTD,'temp');
CTD_id_SAL = netcdf.inqVarID(ncid_CTD,'sal');
CTD_id_PRES = netcdf.inqVarID(ncid_CTD,'pres');
CTD_id_ML_PRES = netcdf.inqVarID(ncid_CTD,'mld_pres');
CTD_id_ML_PRESERROR = netcdf.inqVarID(ncid_CTD,'mld_presError');
CTD_id_ML_DENS = netcdf.inqVarID(ncid_CTD,'mld_dens');
CTD_id_ML_SALT = netcdf.inqVarID(ncid_CTD,'mld_salt');
CTD_id_ML_TEMP = netcdf.inqVarID(ncid_CTD,'mld_temp');
for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['CTD_id_ML_' field ' = netcdf.inqVarID(ncid_CTD,''' field ''');']);%=============================MODIF
end
    ncid_AFD = netcdf.open(argo_float_data,'NOWRITE');
    disp('loading Argo float profile locations')
    AFD_id_LON = netcdf.inqVarID(ncid_AFD,'lon');
    AFD_id_LAT = netcdf.inqVarID(ncid_AFD,'lat');
    AFD_id_DYR = netcdf.inqVarID(ncid_AFD,'dyr');
    AFD_id_SIGMA = netcdf.inqVarID(ncid_AFD,'sigma');
    AFD_id_TEMP = netcdf.inqVarID(ncid_AFD,'temp');
    AFD_id_SAL = netcdf.inqVarID(ncid_AFD,'sal');
    AFD_id_PRES = netcdf.inqVarID(ncid_AFD,'pres');
    AFD_id_ML_PRES = netcdf.inqVarID(ncid_AFD,'mld_pres');
    AFD_id_ML_PRESERROR = netcdf.inqVarID(ncid_AFD,'mld_presError');
    AFD_id_ML_DENS = netcdf.inqVarID(ncid_AFD,'mld_dens');
    AFD_id_ML_SALT = netcdf.inqVarID(ncid_AFD,'mld_salt');
    AFD_id_ML_TEMP = netcdf.inqVarID(ncid_AFD,'mld_temp');
    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['AFD_id_ML_' field ' = netcdf.inqVarID(ncid_AFD,''' field ''');']);%=============================MODIF
    end
%%Trying to paralelize%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END Reopen because I don't know how to pass netcdf-id to a function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                            while ~isempty(profs_toload)
                                                % do final summarize of data that is adjecent to each other to load
                                                nc_load_width = sum(cumprod(multi_load_vec))+1;
                                                
                                                if numel(profs_redundant)>=nc_load_width
                                                    % assign id in memory block to overwrite
                                                    ID_overwrite = profs_redundant(1:nc_load_width);
                                                    profs_redundant(1:nc_load_width) = [];
                                                else
                                                    % assign  IDs at end of 'known' memory block
                                                    ID_overwrite = (numel(profs_inmem)+1) : (numel(profs_inmem)+nc_load_width);
                                                end
                                                
                                                
                                                switch sign(profs_toload(1)) + sign(profs_toload(nc_load_width)) % determine if CTD or Argo float data is needed:
                                                    case -2 % load data from CTD data base
                                                        ttmp = netcdf.getVar(ncid_CTD,CTD_id_ML_DENS,CTD_offset+profs_toload(1),nc_load_width);
                                                        if ttmp > 100
                                                            MLdens(ID_overwrite) = ttmp -1000; %================== MODIF
                                                        else
                                                            MLdens(ID_overwrite) = ttmp;
                                                        end
                                                        MLpres(ID_overwrite) = netcdf.getVar(ncid_CTD,CTD_id_ML_PRES,CTD_offset+profs_toload(1),nc_load_width);
                                                        MLpresError(ID_overwrite) = netcdf.getVar(ncid_CTD,CTD_id_ML_PRESERROR,CTD_offset+profs_toload(1),nc_load_width);
                                                        MLsalt(ID_overwrite) = netcdf.getVar(ncid_CTD,CTD_id_ML_SALT,CTD_offset+profs_toload(1),nc_load_width);
                                                        MLtemp(ID_overwrite) = netcdf.getVar(ncid_CTD,CTD_id_ML_TEMP,CTD_offset+profs_toload(1),nc_load_width);
    							for ifield=1:length(AdditionalFields);
								field=AdditionalFields{ifield};
								eval(['ML' field '(ID_overwrite)  = netcdf.getVar(ncid_CTD,CTD_id_ML_' field ',CTD_offset+profs_toload(1),nc_load_width);']); %=======================	MODIF
    							end

                                                    case 2 % load data from Argo float data base
                                                        ttmp = netcdf.getVar(ncid_AFD,AFD_id_ML_DENS,profs_toload(1),nc_load_width);
                                                        if ttmp > 100
                                                            MLdens(ID_overwrite) = ttmp-1000;
                                                        else
                                                            MLdens(ID_overwrite) = ttmp;
                                                        end
                                                        MLpres(ID_overwrite) = netcdf.getVar(ncid_AFD,AFD_id_ML_PRES,profs_toload(1),nc_load_width);
                                                        MLpresError(ID_overwrite) = netcdf.getVar(ncid_AFD,AFD_id_ML_PRESERROR,profs_toload(1),nc_load_width);
                                                        MLsalt(ID_overwrite) = netcdf.getVar(ncid_AFD,AFD_id_ML_SALT,profs_toload(1),nc_load_width);
                                                        MLtemp(ID_overwrite) = netcdf.getVar(ncid_AFD,AFD_id_ML_TEMP,profs_toload(1),nc_load_width);
    							for ifield=1:length(AdditionalFields);
								field=AdditionalFields{ifield};
								eval(['ML' field '(ID_overwrite)  = netcdf.getVar(ncid_AFD,AFD_id_ML_' field ',profs_toload(1),nc_load_width);']); %=======================	MODIF
    							end
                                                    otherwise % load data from both data sets
                                                        for loopstep = 1:nc_load_width % read each profile individually - ignore multi-load
                                                            if profs_toload(loopstep)<0
                                                                if compute_interior
                                                                    sg_th(ID_overwrite(loopstep),:) = netcdf.getVar(ncid_CTD,CTD_id_TEMP,[sg_lev(1)-1 CTD_offset+profs_toload(loopstep)],[sg_l 1])';
                                                                    tmp = ones_sig ;
                                                                    tmp(sg_th(ID_overwrite(loopstep),:)>100) = NaN;
                                                                    sg_th(ID_overwrite(loopstep),:) = sg_th(ID_overwrite(loopstep),:) .* tmp;
                                                                    sg_pr(ID_overwrite(loopstep),:) = netcdf.getVar(ncid_CTD,CTD_id_PRES,[sg_lev(1)-1 CTD_offset+profs_toload(loopstep)],[sg_l 1])'.*tmp;
                                                                    sg_sa(ID_overwrite(loopstep),:) = netcdf.getVar(ncid_CTD,CTD_id_SAL,[sg_lev(1)-1 CTD_offset+profs_toload(loopstep)],[sg_l 1])'.*tmp;
                                                                end
                                                                ttmp =  netcdf.getVar(ncid_CTD,CTD_id_ML_DENS,CTD_offset+profs_toload(loopstep),1);
                                                                if ttmp > 100
                                                                    MLdens(ID_overwrite(loopstep)) = ttmp-1000;
                                                                else
                                                                    MLdens(ID_overwrite(loopstep)) = ttmp;
                                                                end
                                                                MLpres(ID_overwrite(loopstep)) = netcdf.getVar(ncid_CTD,CTD_id_ML_PRES,CTD_offset+profs_toload(loopstep),1);
                                                                MLpresError(ID_overwrite(loopstep)) = netcdf.getVar(ncid_CTD,CTD_id_ML_PRESERROR,CTD_offset+profs_toload(loopstep),1);
                                                                MLsalt(ID_overwrite(loopstep)) = netcdf.getVar(ncid_CTD,CTD_id_ML_SALT,CTD_offset+profs_toload(loopstep),1);
                                                                MLtemp(ID_overwrite(loopstep)) = netcdf.getVar(ncid_CTD,CTD_id_ML_TEMP,CTD_offset+profs_toload(loopstep),1);
    								for ifield=1:length(AdditionalFields);
									field=AdditionalFields{ifield};
									eval(['ML' field '(ID_overwrite(loopstep))  = netcdf.getVar(ncid_CTD,CTD_id_ML_' field ',CTD_offset+profs_toload(loopstep),1);']); %=======================	MODIF
    								end
                                                             else
                                                                ttmp = netcdf.getVar(ncid_AFD,AFD_id_ML_DENS,profs_toload(loopstep),1);
                                                                if ttmp > 100
                                                                    MLdens(ID_overwrite(loopstep)) = ttmp-1000; %==================MODIF
                                                                else
                                                                    MLdens(ID_overwrite(loopstep)) = ttmp; %==================MODIF
                                                                end
                                                                MLpres(ID_overwrite(loopstep)) = netcdf.getVar(ncid_AFD,AFD_id_ML_PRES,profs_toload(loopstep),1);
                                                                MLpresError(ID_overwrite(loopstep)) = netcdf.getVar(ncid_AFD,AFD_id_ML_PRESERROR,profs_toload(loopstep),1);
                                                                MLsalt(ID_overwrite(loopstep)) = netcdf.getVar(ncid_AFD,AFD_id_ML_SALT,profs_toload(loopstep),1);
                                                                MLtemp(ID_overwrite(loopstep)) = netcdf.getVar(ncid_AFD,AFD_id_ML_TEMP,profs_toload(loopstep),1);
    								for ifield=1:length(AdditionalFields);
									field=AdditionalFields{ifield};
									eval(['ML' field '(ID_overwrite(loopstep))  = netcdf.getVar(ncid_AFD,AFD_id_ML_' field ',profs_toload(loopstep),1);']); %=======================	MODIF
    								end
                                                            end
                                                        end
                                                end
                                                
                                                gef = find(MLdens(ID_overwrite) < 0 | MLdens(ID_overwrite) > 10000);
                                                MLdens(ID_overwrite(gef)) = NaN;
                                                MLpres(ID_overwrite(gef)) = NaN;
                                                MLpresError(ID_overwrite(gef)) = NaN;
                                                MLsalt(ID_overwrite(gef)) = NaN;
                                                MLtemp(ID_overwrite(gef)) = NaN;
   						for ifield=1:length(AdditionalFields);
							field=AdditionalFields{ifield};
							eval(['ML' field '(ID_overwrite(gef)) = NaN;']); %=======================	MODIF
    						end
                                                
                                                % require ML for data to be used
                                                sg_pr(ID_overwrite(gef),:) = NaN;
                                                
                                                % which are now newly in mem?
                                                profs_inmem(ID_overwrite) = profs_toload(1:nc_load_width);
                                                
                                                % they no longer need to be loaded
                                                profs_toload(1:nc_load_width) = [];
                                                
                                                % remove multiload by same ammount and take care of overflow at end of vector
                                                if numel(multi_load_vec) > nc_load_width
                                                    multi_load_vec(1:nc_load_width) = [];
                                                else
                                                    multi_load_vec = 0;
                                                end
                                                
                                            end
 
netcdf.close(ncid_CTD)
netcdf.close(ncid_AFD)

if doKuusela==0
    Kuusela=Read_Netcdf_JB(Kuusela_fname,{'LONGITUDE','LATITUDE','thetaLongOpt','thetaLatOpt','thetatOpt','thetasOpt','thetasOpt','sigmaOpt'});
    [tt1 ilat]=min(abs(y(istep,1)-Kuusela.LATITUDE));
    [tt2 ilon]=min(abs(x(1,jstep)-Kuusela.LONGITUDE));
    Deccorelation.Time=time_scale;
    Deccorelation.Space=horiz_scale;
    if tt1<1 & tt2<1
            rEarth = 6371;
            thetaLatOptKm = squeeze(Kuusela.thetaLatOpt(:,ilat,ilon))/360*2*pi*rEarth;
            thetaLonOptKm = (squeeze(Kuusela.thetaLongOpt(:,ilat,ilon))/360).*cos(Kuusela.LATITUDE(ilat)/360*2*pi)*2*pi*rEarth;
            Deccorelation.Space=(sqrt(thetaLonOptKm.^2+thetaLatOptKm.^2)/110).^2; % in kms^2
            Deccorelation.Time=(squeeze(Kuusela.thetatOpt(:,ilat,ilon))/365).^2; %in years^2 
            Deccorelation.Phi=squeeze(Kuusela.thetasOpt(:,ilat,ilon)); 
            Deccorelation.nugget=squeeze(Kuusela.sigmaOpt(:,ilat,ilon)).^2; 
    end
else
    Deccorelation.Space=NaN*ones(1,12); 
    Deccorelation.Time=NaN*ones(1,12); 
    Deccorelation.Phi=NaN*ones(1,12);  
    Deccorelation.nugget=NaN*ones(1,12); 
end

                                            %  disp('done loading') --->>>  now sort! to have same indices for all
                                            i2_dist = find(ismember(profs_inmem,profs_needed));
                                            [~, kk] = sort(profs_needed);
                                            ii_dist = ii_dist(kk);
                                            [~, kk] = sort(profs_inmem(i2_dist));
                                            i2_dist = i2_dist(kk);
                                            
                                            sg_pr = double(sg_pr);
                                            sg_sa = double(sg_sa);
                                            sg_th = double(sg_th);
                                            MLdens = double(MLdens);
                                            MLpres = double(MLpres);
                                            MLpresError = double(MLpresError);
                                            MLsalt = double(MLsalt);
                                            MLtemp = double(MLtemp);
                                            for ifield=1:length(AdditionalFields);
                                                field=AdditionalFields{ifield};
                                                eval(['ML' field '= double(ML' field ');']); %=======================	MODIF
                                            end
                                            
                                            dyrTMP = dyr(ii_dist);
                                            dyrMAP = dyr_map(ii_dist);
                                            lonTMP=lon(ii_dist);
                                            latTMP=lat(ii_dist);
                                            dist_vec_complexTMP = dist_vec_complex(ii_dist);
                                            dec_dyrTMP = dec_dyr(ii_dist);
                                            %compute the huge covariancematrix
                                            %   E_all2 = ((double(((repmat(real(dist_vec_complex(ii_dist)),1,numel(ii_dist)) - repmat(real(dist_vec_complex(ii_dist))',numel(ii_dist),1)).^2) ...
                                            %                   + ((repmat(imag(dist_vec_complex(ii_dist)),1,numel(ii_dist)) - repmat(imag(dist_vec_complex(ii_dist))',numel(ii_dist),1)).^2))));
                                            %               E_all2 =( sqrt(E_all2).^(expweight))./horiz_scale.^expweight;
                                            %     E_all2 = E_all2+ double(abs(repmat(dec_dyr(ii_dist),1,numel(ii_dist)) - repmat(dec_dyr(ii_dist)',numel(ii_dist),1)).^expweight ) ./ dec_scale.^expweight ;
                                            
                                            
                                            % C_all =   (sqrt(abs(real(dist_vec_complexTMP).*real(dist_vec_complexTMP)) ...
                                            %     + abs(imag(dist_vec_complexTMP).*imag(dist_vec_complexTMP))));
                                            % identical to:
                                            C_all =abs(dist_vec_complexTMP);
                                            
                                            % EXP 1.5
                                            %C_all = C_all .* sqrt(C_all) ./ (horiz_scale.*sqrt(horiz_scale));
                                            % EXP 2
                                            C_all = C_all .* C_all ./ horiz_scale;
                                            
                                            %EXP 1.5
                                            %E_all_T = make_dist_matrix_mex_exp_scal(...
                                            %    real(dist_vec_complexTMP), ...
                                            %    imag(dist_vec_complexTMP), ...
                                            %    horiz_scale.*sqrt(horiz_scale)) ...
                                            %    + make_time_matrix_mex_exp_scal(dyrTMP,time_scale.*sqrt(time_scale));
                                            %
                                            % EXP 2
%                                             E_all_T = make_dist_matrix_mex_exp_scal_gauss(...
%                                                 real(dist_vec_complexTMP), ...
%                                                 imag(dist_vec_complexTMP), ...
%                                                 horiz_scale) ...
%                                                 + make_time_matrix_mex_exp_scal_gauss(dyrTMP,time_scale);
                                            
                                            
                                            % -------------------------------------
                                            % loop over months/seasons
                                            % -------------------------------------
                                            
                                            i2_dist_tmp2 = i2_dist;
                                            ii_dist_tmp2 = ii_dist;

varlist={'ML_dens_gr';	'ML_pres_gr';'ML_salt_gr';'ML_temp_gr';'ML_year_gr';'ML_salt_gr_wmean';'ML_temp_gr_wmean';'ML_dens_gr_wmean';'ML_pres_gr_wmean';'ML_salt_gr_wstd';'ML_temp_gr_wstd';'ML_dens_gr_wstd';'ML_pres_gr_wstd';'ML_salt_gr_wSE';'ML_temp_gr_wSE';'ML_dens_gr_wSE';
'ML_pres_gr_wSE';'ML_salt_gr_wTrend';'ML_temp_gr_wTrend';'ML_dens_gr_wTrend';'ML_pres_gr_wTrend';'ML_salt_gr_wX_corr';'ML_temp_gr_wX_corr';'ML_dens_gr_wX_corr';'ML_pres_gr_wX_corr';'ML_gr_YearPer5';'ML_gr_YearPer50';'ML_gr_YearPer95';'thetasOpt_wmean';
'thetadOpt_wmean';'thetatOpt_wmean';'sigmaOpt_wmean';'ML_num_used'};

for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['varlist= [varlist ; {''ML_' field '_gr_wTrend''}];']); %=======================MODIF
	eval(['varlist= [varlist ; {''ML_' field '_gr_wX_corr''}];']); %=======================MODIF
	eval(['varlist= [varlist ; {''ML_' field '_gr_wSE''}];']); %=======================MODIF
	eval(['varlist= [varlist ; {''ML_' field '_gr_wstd''}];']); %=======================MODIF
	eval(['varlist= [varlist ; {''ML_' field '_gr_wmean''}];']); %=======================MODIF
end


					    varINlist2=[varINlist {'optionFit','Deccorelation','dist_vec_complexTMP','dec_dyrTMP','lonTMP','latTMP','dyrMAP','C_all','dyrTMP','i2_dist_tmp2','ii_dist_tmp2','istep','jstep'}];
					    for ii=1:length(varINlist2)
						eval([ 'INpara2.' varINlist2{ii}  '=' varINlist2{ii}  ';']);
					    end

					    INpara2=INpara2;
                        tic
					    %parfor tstep =  1:tim % about 12 times ... or 4 times ...
                        for tstep =  1:tim % about 12 times ... or 4 times ...
                            out2(tstep)=Sunke_OI_Trend_Final_Sub3_new(INpara2,varINlist2,tstep);
                        end
                        toc

					    for tstep =  1:tim
					    	for ivar=1:length(varlist)
                                eval([varlist{ivar} '(tstep)= out2(tstep).'  varlist{ivar} '(tstep);']);
                            end
					    end
					    clear out2
                
					    out.istep=istep; out.jstep=jstep;
					    out.idvar{1}=ML_dens_ID;  		out.var{1}=single(ML_dens_gr);
					    out.idvar{2}=ML_pres_ID;  		out.var{2}=single(ML_pres_gr);

					    out.idvar{3}=thetasOpt_wm_ID;  	out.var{3}=single(thetasOpt_wmean);
					    out.idvar{4}=thetadOpt_wm_ID;  	out.var{4}=single(thetadOpt_wmean);
					    out.idvar{5}=thetadOpt_wm_ID; 	out.var{5}=single(thetadOpt_wmean);
					    out.idvar{6}=sigmaOpt_wm_ID;  	out.var{6}=single(sigmaOpt_wmean);

					    out.idvar{7}=ML_dens_wm_ID;  	out.var{7}=single(ML_dens_gr_wmean);
					    out.idvar{8}=ML_pres_wm_ID;  	out.var{8}=single(ML_pres_gr_wmean);
					    out.idvar{9}=ML_salt_wm_ID;  	out.var{9}=single(ML_salt_gr_wmean);
					    out.idvar{10}=ML_temp_wm_ID;	out.var{10}=single(ML_temp_gr_wmean);

					    out.idvar{11}=ML_dens_wstd_ID;  	out.var{11}=single(ML_dens_gr_wstd);
					    out.idvar{12}=ML_pres_wstd_ID;  	out.var{12}=single(ML_pres_gr_wstd);
					    out.idvar{13}=ML_salt_wstd_ID;  	out.var{13}=single(ML_salt_gr_wstd);
					    out.idvar{14}=ML_temp_wstd_ID;	out.var{14}=single(ML_temp_gr_wstd);

					    out.idvar{15}=ML_dens_wSE_ID;  	out.var{15}=single(ML_dens_gr_wSE);
					    out.idvar{16}=ML_pres_wSE_ID;  	out.var{16}=single(ML_pres_gr_wSE);
					    out.idvar{17}=ML_salt_wSE_ID;  	out.var{17}=single(ML_salt_gr_wSE);
					    out.idvar{18}=ML_temp_wSE_ID;	out.var{18}=single(ML_temp_gr_wSE);

					    out.idvar{19}=ML_dens_wTrend_ID;  	out.var{19}=single(ML_dens_gr_wTrend);
					    out.idvar{20}=ML_pres_wTrend_ID;  	out.var{20}=single(ML_pres_gr_wTrend);
					    out.idvar{21}=ML_salt_wTrend_ID;  	out.var{21}=single(ML_salt_gr_wTrend);
					    out.idvar{22}=ML_temp_wTrend_ID;	out.var{22}=single(ML_temp_gr_wTrend);

					    out.idvar{23}=ML_dens_wX_corr_ID;  	out.var{23}=single(ML_dens_gr_wX_corr);
					    out.idvar{24}=ML_pres_wX_corr_ID;  	out.var{24}=single(ML_pres_gr_wX_corr);
					    out.idvar{25}=ML_salt_wX_corr_ID;  	out.var{25}=single(ML_salt_gr_wX_corr);
					    out.idvar{26}=ML_temp_wX_corr_ID;	out.var{26}=single(ML_temp_gr_wX_corr);

					    out.idvar{27}=ML_yprc5_ID;  	out.var{27}=single(ML_gr_YearPer5);
					    out.idvar{28}=ML_yprc50_ID;  	out.var{28}=single(ML_gr_YearPer50);
					    out.idvar{29}=ML_yprc95_ID;  	out.var{29}=single(ML_gr_YearPer95);

					    out.idvar{30}=ML_salt_ID;  		out.var{30}=single(ML_salt_gr);
					    out.idvar{31}=ML_temp_ID;  		out.var{31}=single(ML_temp_gr);
					    out.idvar{32}=ML_nwgt_ID;  		out.var{32}=single(ML_weight);

					    out.idvar{33}=ML_year_ID;  		out.var{33}=int16((ML_year_gr-2000).*1000);
					    out.idvar{34}=ML_dist_M_ID;  	out.var{34}=int16(imag(ML_r_dist));
					    out.idvar{35}=ML_dist_Z_ID;  	out.var{35}=int16(real(ML_r_dist));
					    out.idvar{36}=ML_error_ID;  	out.var{36}=single(ML_error);
					    out.idvar{37}=ML_area_ID;  		out.var{37}=single(ML_area);
 					    out.idvar{38}=ML_num_data_ID;  	out.var{38}=int16(ML_num_used);
					    out.idvar{39}=ML_drho_ID;  		out.var{39}=single(ML_Drho);

					    out.idvar{40}=thetatOpt_wm_ID;  	out.var{40}=single(thetatOpt_wmean);

					    cmptvar=40;
 					    for ifield=1:length(AdditionalFields);
					    	field=AdditionalFields{ifield};
					        eval(['out.idvar{cmptvar+1}=ML_' field '_ID;		out.var{cmptvar+1}=single(ML_' field '_gr);']);
					        eval(['out.idvar{cmptvar+2}=ML_' field '_wm_ID;		out.var{cmptvar+2}=single(ML_' field '_gr_wmean);']);
					        eval(['out.idvar{cmptvar+3}=ML_' field '_wstd_ID;	out.var{cmptvar+3}=single(ML_' field '_gr_wstd);']);
					        eval(['out.idvar{cmptvar+4}=ML_' field '_wSE_ID;	out.var{cmptvar+4}=single(ML_' field '_gr_wSE);']);
					        eval(['out.idvar{cmptvar+5}=ML_' field '_wTrend_ID;	out.var{cmptvar+5}=single(ML_' field '_gr_wTrend);']);
					        eval(['out.idvar{cmptvar+6}=ML_' field '_wX_corr_ID;	out.var{cmptvar+6}=single(ML_' field '_gr_wX_corr);']);
                            cmptvar=cmptvar+6;
					    end

	                                            
                        ML_dens_gr = dummy_T;
                        ML_pres_gr = dummy_T;
                        ML_year_gr = dummy_T;
                        ML_temp_gr = dummy_T;
                        ML_salt_gr = dummy_T;
                        ML_dens_gr_wmean = dummy_T;
                        ML_salt_gr_wmean = dummy_T;
                        ML_temp_gr_wmean = dummy_T;
                        ML_pres_gr_wmean = dummy_T;
                        ML_weight = dummy_T;
                        ML_r_dist = complex(dummy_T);
                        ML_num_used = dummy_T;
                        ML_area = dummy_T;
                        ML_error= dummy_T;
                        ML_Drho = dummy_T;
 					    for ifield=1:length(AdditionalFields);
					    	field=AdditionalFields{ifield};
					    	eval(['ML_' field '_gr_wmean = dummy_T;']);%=======================	MODIF
    					end                                                        
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end

