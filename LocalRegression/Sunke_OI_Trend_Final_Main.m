
function Sunke_OI_Trend_Final_Main(argo_float_data, wod_ctd_data, FM_fname, Kuusela_fname, AdditionalFields, output_fname,  lon_grid, lat_grid, option, logfile);

eval(['!echo Start @'  datestr(now,'HH:MM:SS') ' HH:MM:SS. > ' logfile])           


% climatology from float data - maybe do this basin wide at a time ... ...
% especially if this data set get's combined with CTD data.

% mex make_dist_matrix_mex_exp_scal.c CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"
% gcc 4-2 or newer
% check mex opts file!!!


% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% VARIABLES YOU MAY TUNE
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

[x,y]=meshgrid(lon_grid, lat_grid);	 % LON, LAT

if length(option)~=36
	disp('Option input has wrong dimensions')
	return
end
eps100			=option(1); 
eps999 			=option(2); 
create_new 		=option(3); 
compute_interior  	=option(4); 
drawit 			=option(5); 
useargo 		=option(6); 
upperocean 		=option(7); 
upperocean_limit 	=option(8); 
wfkt_multiplier  	=option(9); 
max_profs_inmem 	=option(10); 
min_CTD_profs		=option(11); 
max_num_profs 		=option(12); 
min_num_profs 		=option(13); 
num_of_rand_profiles 	=option(14); 
min_allowed_weight 	=option(15); 
rnd_within 		=option(16); 
iqr_mult 		=option(17); 
continue_at 		=option(18); 
continue_upto		=option(19); 
continue_lon_at		=option(20); 
continue_lon_upto	=option(21);  
OMP_NUM_THREADS_DAY 	=option(22); 
OMP_NUM_THREADS_NIGHT 	=option(23); 
signal_to_noise 	=option(24); 
horiz_scale_set		=option(25); 
time_scale 		=option(26); 
tim 			=option(27); 
centeryear		=option(28); 
dec_scale		=option(29); 
dec_noise_offset 	=option(30); 
yeartrendmin 		=option(31); 
doIQR			=option(32);
optionFit       =option(33);
reverselat      =option(34);
Parallel        =option(35);
doKuusela       =option(36);
%%
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% VARIABLE YOU MAY   = NOT =   TUNE
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

%% preliminary constants and initialization
%expweight = 1.50; % FIXED, DO NOT CHANGE % use something that can be expressed by x * sqrt(x)
% or similar, since x.^y is VERY time consuming for rem(y,0.5) ~= 0
% (DOUBLE the time for other exponents!!! (25% more time for the whole
% project!!! expweight = 0.5,1,1.5,2 ... the ones to choose from but to implement manually.)

expweight = 2; %!!!!

numsortedprofs = max_num_profs-num_of_rand_profiles;

% Matrix inverse options:
opts.SYM = true;  % COVARIANCE MATRIX ALWAYS symetrical
opts.POSDEF = true; % covariance matrix always positive defined

[m,n]=size(x)  % grid to compute

horiz_scale = (horiz_scale_set / 110)^2;


% find max dist for nearest neighbor
maxdist_finegrid = 0.15 * sqrt(2);
maxdist_coarsegrid = 0.51 * sqrt(2);


%% INITIALIZE ACCESS TO ISOPYCNAL CTD AND FLOAT DATA:    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Argo_file.nc
if useargo
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
    AFD_id_ML_DENS = netcdf.inqVarID(ncid_AFD,'mld_dens');
    AFD_id_ML_SALT = netcdf.inqVarID(ncid_AFD,'mld_salt');
    AFD_id_ML_TEMP = netcdf.inqVarID(ncid_AFD,'mld_temp');
    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['AFD_id_ML_' field ' = netcdf.inqVarID(ncid_AFD,''' field ''');']);%=============================MODIF
    end
    AFD_LAT = netcdf.getVar(ncid_AFD,AFD_id_LAT);
    AFD_LON = netcdf.getVar(ncid_AFD,AFD_id_LON);
    AFD_DYR = netcdf.getVar(ncid_AFD,AFD_id_DYR);
    sg0_pleth = double(netcdf.getVar(ncid_AFD,AFD_id_SIGMA));
    
    % create index vector to find the right profiles in netcdf file.
    AFD_IDVEC = (0:numel(AFD_LAT)-1)';% negative values for CTD file, positive for Argo
    disp(' finished reading Argo float profile locations')
    
else
    AFD_LAT = 190;  %#ok<UNRCH> % arbitary 'off globe' location for a single non-existant profile
    AFD_LON = 0;
    AFD_DYR = -100;
    sg0_pleth = [];
    % create index vector to find the right profiles in netcdf file.
    AFD_IDVEC = 100000;
end



ncid_CTD = netcdf.open(wod_ctd_data,'NOWRITE');	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CTD_file.nc
CTD_id_LON = netcdf.inqVarID(ncid_CTD,'lon');
CTD_id_LAT = netcdf.inqVarID(ncid_CTD,'lat');
CTD_id_DYR = netcdf.inqVarID(ncid_CTD,'dyr');
CTD_id_SIGMA = netcdf.inqVarID(ncid_CTD,'sigma');
CTD_id_TEMP = netcdf.inqVarID(ncid_CTD,'temp');
CTD_id_SAL = netcdf.inqVarID(ncid_CTD,'sal');
CTD_id_PRES = netcdf.inqVarID(ncid_CTD,'pres');
CTD_id_ML_PRES = netcdf.inqVarID(ncid_CTD,'mld_pres');
CTD_id_ML_DENS = netcdf.inqVarID(ncid_CTD,'mld_dens');
CTD_id_ML_SALT = netcdf.inqVarID(ncid_CTD,'mld_salt');
CTD_id_ML_TEMP = netcdf.inqVarID(ncid_CTD,'mld_temp');
for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['CTD_id_ML_' field ' = netcdf.inqVarID(ncid_CTD,''' field ''');']);%=============================MODIF
end

CTD_LAT = netcdf.getVar(ncid_CTD,CTD_id_LAT);
CTD_LON = netcdf.getVar(ncid_CTD,CTD_id_LON);
CTD_DYR = netcdf.getVar(ncid_CTD,CTD_id_DYR);

if isempty(sg0_pleth)
    sg0_pleth = netcdf.getVar(ncid_CTD,CTD_id_SIGMA);
end


% create index vector to find the right profiles in netcdf file.
CTD_IDVEC = (-numel(CTD_LAT):-1)'; % negative values for CTD file, positive for Argo
CTD_offset =numel(CTD_LAT);
disp(' finished reading WOC CTD profile locations')

%%

disp(' pre allocating matrices')


% pick out levels on which to map things
sg_lev=1:1:length(sg0_pleth); % << full length as set in interpolated profiles.
% or subset: ONLY CONTINIOUS SUBSETS POSSIBLE SO FAR!!!
%sg_lev = find(sg0_pleth>24.0,1,'first'):find(sg0_pleth<28.9,1,'last');
%sg_lev = find(sg0_pleth>24.3,1,'first'):find(sg0_pleth<26.65,1,'last');

%sg_l_orig=numel(sg_lev);

% just use those levels - remove others
sg0_pleth=sg0_pleth(sg_lev);


sg_l=numel(sg0_pleth); % number of sigma levels


% initialize time - sigma matrices
E_Z = NaN(max_num_profs);

dummy_T_sig = NaN([tim sg_l],'double');
dummy_T = NaN([tim 1],'double');

ones_sig = ones(1,sg_l);

for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['ML_' field '_gr = dummy_T;']);%=============================MODIF
	eval(['ML_' field '_gr_wmean = dummy_T;']);
	eval(['ML_' field '_gr_wstd = dummy_T;']);
	eval(['ML_' field '_gr_wSE = dummy_T;']);
	eval(['ML_' field '_gr_wTrend = dummy_T;']);
	eval(['ML_' field '_gr_wX_corr = dummy_T;']);
end

ML_area= dummy_T;
ML_error= dummy_T;
ML_weight=dummy_T;
ML_num_used=dummy_T;
ML_r_dist = dummy_T;
ML_Drho = dummy_T;

sg_pr = NaN([max_profs_inmem,sg_l],'single');
sg_sa = sg_pr;
sg_th = sg_pr;

MLpres = NaN([max_profs_inmem,1],'single');
MLpresError = MLpres;
MLdens = MLpres;
MLsalt = MLpres;
MLtemp = MLpres;
for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['ML' field ' = MLpres;']);%=============================MODIF
	eval(['ML' field ' = MLpres;']);
end

profs_inmem = (-100000:-100000+max_profs_inmem-1)-9e12; % set to 10,000 unique values not within IDs from floats or CTD ... ...

%% INITIALIZE ACCESS TO FAST MARCHING DATA:
disp(' loading Fast Marching gridpoints')

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

dimid = netcdf.inqDimID(ncid_FM,'subgrid');
[~,FM_ALL] = netcdf.inqDim(ncid_FM,dimid);
dimid = netcdf.inqDimID(ncid_FM,'subgrid_coarse');
[~,FM_ALL_C] = netcdf.inqDim(ncid_FM,dimid);

FM_LAT = double(netcdf.getVar(ncid_FM,FM_id_LAT));
FM_LON = double(netcdf.getVar(ncid_FM,FM_id_LON));
FM_DPT = double(netcdf.getVar(ncid_FM,FM_id_DPT));
%FM_NUM = double(netcdf.getVar(ncid_FM,FM_id_NUM));
% cut down data set in memory to only data that might get needed:

disp(' wrapping data around globe')

% remove data that won't be used from LAT and LON vectors
% doing this PRIOR merging CTD and Argo data since both vectors already can be VERY large.
found_it = find(FM_LON < 0);

FM_LON(found_it) = FM_LON(found_it)+360;
%%
found_it = find(AFD_LON < 0);

AFD_LON(found_it) = AFD_LON(found_it)+360;

found_it = find(CTD_LON < 0);

CTD_LON(found_it) = CTD_LON(found_it)+360;


% Grid to compute adapted to FM grid: 
FM_LAT=round(FM_LAT.*100)/100;
FM_LON=round(FM_LON.*100)/100;

yy=unique(FM_LAT);
xx=unique(FM_LON);
xx_ok=find(xx>=min(lon_grid) & xx<=max(lon_grid));
yy_ok=find(yy>=min(lat_grid) & yy<=max(lat_grid));
lon_grid=xx(xx_ok); 
lat_grid=yy(yy_ok);
[x,y]=meshgrid(lon_grid, lat_grid);
[m,n]=size(x)  % grid to compute

%% join all data sets

lon = double([AFD_LON; CTD_LON]);
lat = double([AFD_LAT; CTD_LAT]);
dyr = double([AFD_DYR; CTD_DYR]);

dyr_map = dyr;

dec_dyr = double(dyr - centeryear);
%dec_dyr(dec_dyr > 0) = 0;
dec_dyr = (dec_noise_offset + signal_to_noise) - dec_noise_offset*exp(-(abs(dec_dyr./dec_scale)).^expweight) ;

dyr = double(rem(dyr,1));
idv = [AFD_IDVEC; CTD_IDVEC];

dist_dummy = NaN(size(lat));

%free some memory
clear AFD_LAT AFD_LON AFD_DYR CTD_LAT CTD_LON CTD_DYR CTD_IDVEC AFD_IDVEC found_it

%% INITIALIZE FINAL OUTPUT NETCDF FILE

% write to one LARGE netcdf file:::::
if ~exist(output_fname,'file') || create_new
    disp('si il nexiste pas: create new file Clim.nc')
    ncid_CLIM = netcdf.create(output_fname,'NC_64BIT_OFFSET');
    
    % Define the dimensions of the variables.
    dimid_gridLON    = netcdf.defDim(ncid_CLIM,'LON_VEC',n);
    dimid_gridLAT    = netcdf.defDim(ncid_CLIM,'LAT_VEC',m);
    dimid_months     = netcdf.defDim(ncid_CLIM,'MONTH_VEC',tim);
    dimid_isopycnals = netcdf.defDim(ncid_CLIM,'PDENS_VEC',sg_l);
    dimid_scalar     = netcdf.defDim(ncid_CLIM,'SCALAR',1);
    
    
    % Define new variables in file.
    latID  = netcdf.defVar(ncid_CLIM,'LATITUDE','float',dimid_gridLAT);
    lonID  = netcdf.defVar(ncid_CLIM,'LONGITUDE','float',dimid_gridLON);
    sigID  = netcdf.defVar(ncid_CLIM,'NEUTRAL_DENSITY','float',dimid_isopycnals);
    timID  = netcdf.defVar(ncid_CLIM,'MONTH','float',dimid_months);
    
    iqrID  = netcdf.defVar(ncid_CLIM,'IQR_RANGE','float',dimid_scalar);
    max_numID    = netcdf.defVar(ncid_CLIM,'MAX_NUM_PROFS','short',dimid_scalar);
    max_num2ID    = netcdf.defVar(ncid_CLIM,'MAX_NUM_PROFS_IN_MEMORY_GRIDPOINT','short',dimid_scalar);
    min_numID    = netcdf.defVar(ncid_CLIM,'MIN_NUM_PROFS','short',dimid_scalar);
    rnd_numID    = netcdf.defVar(ncid_CLIM,'RND_NUM_PROFS','short',dimid_scalar);
    ctd_numID    = netcdf.defVar(ncid_CLIM,'REQ_CTD_PROFS','short',dimid_scalar);
    horiz_scaleID    = netcdf.defVar(ncid_CLIM,'HORIZONTAL_SCALE','float',dimid_scalar);
    temporal_scaleID    = netcdf.defVar(ncid_CLIM,'TEMPORAL_SCALE','float',dimid_scalar);
    noise_signalID    = netcdf.defVar(ncid_CLIM,'NOISE_SIGNAL_RATIO','float',dimid_scalar);
    decadal_scaleID    = netcdf.defVar(ncid_CLIM,'DECADAL_SCALE','float',dimid_scalar);
    decadal_noiseID    = netcdf.defVar(ncid_CLIM,'DECADAL_NOISE','float',dimid_scalar);
    disp('end: create new file Clim.nc')
    
    if compute_interior
        
        data_pres_ID = netcdf.defVar(ncid_CLIM,'PRESSURE','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_pres_ID,'_FillValue',NaN('single'))
        data_salt_ID = netcdf.defVar(ncid_CLIM,'SALINITY','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_salt_ID,'_FillValue',NaN('single'))
        data_temp_ID = netcdf.defVar(ncid_CLIM,'POTENTIAL_TEMPERATURE','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_temp_ID,'_FillValue',NaN('single'))
        data_pres_wm_ID = netcdf.defVar(ncid_CLIM,'PRESSURE_WEIGHTED_MEAN','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_pres_wm_ID,'_FillValue',NaN('single'))
        data_salt_wm_ID = netcdf.defVar(ncid_CLIM,'SALINITY_WEIGHTED_MEAN','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_salt_wm_ID,'_FillValue',NaN('single'))
        data_temp_wm_ID = netcdf.defVar(ncid_CLIM,'POTENTIAL_TEMPERATURE_WEIGHTED_MEAN','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_temp_wm_ID,'_FillValue',NaN('single'))
        data_error_ID= netcdf.defVar(ncid_CLIM,'ERROR_OBJECTIVE_ANALYSIS','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_error_ID,'_FillValue',NaN('single'))
        data_area_ID = netcdf.defVar(ncid_CLIM,'AREA_OBJECTIVE_ANALYSIS','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_area_ID,'_FillValue',NaN('single'))
        data_nwgt_ID = netcdf.defVar(ncid_CLIM,'WEIGHT_SUMMED','float',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_nwgt_ID,'_FillValue',NaN('single'))
        data_dist_Z_ID = netcdf.defVar(ncid_CLIM,'DISTANCE_WEIGHTED_ZONAL','short',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_dist_Z_ID,'_FillValue',int16(-32000))
        data_dist_M_ID = netcdf.defVar(ncid_CLIM,'DISTANCE_WEIGHTED_MERIDIONAL','short',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_dist_M_ID,'_FillValue',int16(-32000))
        data_num_data_ID = netcdf.defVar(ncid_CLIM,'NUMBER_DATA_POINTS','short',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_num_data_ID,'_FillValue',int16(-32000))
        data_z_weight_ID = netcdf.defVar(ncid_CLIM,'DELTA_Z_SCALE','short',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_z_weight_ID,'_FillValue',int16(-32000))
        
        data_year_ID= netcdf.defVar(ncid_CLIM,'YEAR_OF_DATA','short',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        netcdf.putAtt(ncid_CLIM,data_year_ID,'_FillValue',int16(-32000))
        %  data_year_wm_ID = netcdf.defVar(ncid_CLIM,'YEAR_OF_DATA_WEIGHTED_MEAN','short',[dimid_isopycnals dimid_months dimid_gridLAT dimid_gridLON]);
        %  netcdf.putAtt(ncid_CLIM,data_year_wm_ID,'_FillValue',int16(-32000))
        
    end
    
    disp('Define new variables for ML')
    ML_pres_ID  = netcdf.defVar(ncid_CLIM,'ML_MAX_PRESSURE','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_pres_ID,'_FillValue',NaN('single'))
    ML_dens_ID  = netcdf.defVar(ncid_CLIM,'ML_DENSITY','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_dens_ID,'_FillValue',NaN('single'))
    ML_salt_ID  = netcdf.defVar(ncid_CLIM,'ML_SALINITY','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_salt_ID,'_FillValue',NaN('single'))
    ML_temp_ID  = netcdf.defVar(ncid_CLIM,'ML_TEMPERATURE','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_temp_ID,'_FillValue',NaN('single'))

    thetasOpt_wm_ID  = netcdf.defVar(ncid_CLIM,'thetasOpt','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,thetasOpt_wm_ID,'_FillValue',NaN('single'))
    thetatOpt_wm_ID  = netcdf.defVar(ncid_CLIM,'thetatOpt','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,thetatOpt_wm_ID,'_FillValue',NaN('single'))
    sigmaOpt_wm_ID  = netcdf.defVar(ncid_CLIM,'sigmaOpt','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,sigmaOpt_wm_ID,'_FillValue',NaN('single'))
    thetadOpt_wm_ID  = netcdf.defVar(ncid_CLIM,'thetadOpt','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,thetadOpt_wm_ID,'_FillValue',NaN('single'))


    ML_pres_wm_ID  = netcdf.defVar(ncid_CLIM,'ML_PRESSURE_WEIGHTED_MEAN','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_pres_wm_ID,'_FillValue',NaN('single'))
    ML_dens_wm_ID  = netcdf.defVar(ncid_CLIM,'ML_DENSITY_WEIGHTED_MEAN','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_dens_wm_ID,'_FillValue',NaN('single'))
    ML_salt_wm_ID  = netcdf.defVar(ncid_CLIM,'ML_SALINITY_WEIGHTED_MEAN','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_salt_wm_ID,'_FillValue',NaN('single'))
    ML_temp_wm_ID  = netcdf.defVar(ncid_CLIM,'ML_TEMPERATURE_WEIGHTED_MEAN','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_temp_wm_ID,'_FillValue',NaN('single'))

    ML_pres_wstd_ID  = netcdf.defVar(ncid_CLIM,'ML_PRESSURE_WEIGHTED_STD','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_pres_wstd_ID,'_FillValue',NaN('single'))
    ML_dens_wstd_ID  = netcdf.defVar(ncid_CLIM,'ML_DENSITY_WEIGHTED_STD','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_dens_wstd_ID,'_FillValue',NaN('single'))
    ML_salt_wstd_ID  = netcdf.defVar(ncid_CLIM,'ML_SALINITY_WEIGHTED_STD','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_salt_wstd_ID,'_FillValue',NaN('single'))
    ML_temp_wstd_ID  = netcdf.defVar(ncid_CLIM,'ML_TEMPERATURE_WEIGHTED_STD','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_temp_wstd_ID,'_FillValue',NaN('single'))

    ML_pres_wSE_ID  = netcdf.defVar(ncid_CLIM,'ML_PRESSURE_WEIGHTED_SE','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_pres_wSE_ID,'_FillValue',NaN('single'))
    ML_dens_wSE_ID  = netcdf.defVar(ncid_CLIM,'ML_DENSITY_WEIGHTED_SE','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_dens_wSE_ID,'_FillValue',NaN('single'))
    ML_salt_wSE_ID  = netcdf.defVar(ncid_CLIM,'ML_SALINITY_WEIGHTED_SE','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_salt_wSE_ID,'_FillValue',NaN('single'))
    ML_temp_wSE_ID  = netcdf.defVar(ncid_CLIM,'ML_TEMPERATURE_WEIGHTED_SE','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_temp_wSE_ID,'_FillValue',NaN('single'))

    ML_pres_wTrend_ID  = netcdf.defVar(ncid_CLIM,'ML_PRESSURE_WEIGHTED_Trend','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_pres_wTrend_ID,'_FillValue',NaN('single'))
    ML_dens_wTrend_ID  = netcdf.defVar(ncid_CLIM,'ML_DENSITY_WEIGHTED_Trend','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_dens_wTrend_ID,'_FillValue',NaN('single'))
    ML_salt_wTrend_ID  = netcdf.defVar(ncid_CLIM,'ML_SALINITY_WEIGHTED_Trend','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_salt_wTrend_ID,'_FillValue',NaN('single'))
    ML_temp_wTrend_ID  = netcdf.defVar(ncid_CLIM,'ML_TEMPERATURE_WEIGHTED_Trend','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_temp_wTrend_ID,'_FillValue',NaN('single'))

    ML_pres_wX_corr_ID  = netcdf.defVar(ncid_CLIM,'ML_PRESSURE_WEIGHTED_X_corr','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_pres_wX_corr_ID,'_FillValue',NaN('single'))
    ML_dens_wX_corr_ID  = netcdf.defVar(ncid_CLIM,'ML_DENSITY_WEIGHTED_X_corr','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_dens_wX_corr_ID,'_FillValue',NaN('single'))
    ML_salt_wX_corr_ID  = netcdf.defVar(ncid_CLIM,'ML_SALINITY_WEIGHTED_X_corr','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_salt_wX_corr_ID,'_FillValue',NaN('single'))
    ML_temp_wX_corr_ID  = netcdf.defVar(ncid_CLIM,'ML_TEMPERATURE_WEIGHTED_X_corr','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_temp_wX_corr_ID,'_FillValue',NaN('single'))

    ML_yprc5_ID  = netcdf.defVar(ncid_CLIM,'ML_YearPrc5','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_yprc5_ID,'_FillValue',NaN('single'))
    ML_yprc95_ID  = netcdf.defVar(ncid_CLIM,'ML_YearPrc95','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_yprc95_ID,'_FillValue',NaN('single'))
    ML_yprc50_ID  = netcdf.defVar(ncid_CLIM,'ML_YearPrc50','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_yprc50_ID,'_FillValue',NaN('single'))

    ML_error_ID= netcdf.defVar(ncid_CLIM,'ML_ERROR','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_error_ID,'_FillValue',NaN('single'))
    ML_area_ID = netcdf.defVar(ncid_CLIM,'ML_AREA','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_area_ID,'_FillValue',NaN('single'))
    ML_nwgt_ID = netcdf.defVar(ncid_CLIM,'ML_WEIGHT_SUMMED','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_nwgt_ID,'_FillValue',NaN('single'))
    ML_dist_Z_ID = netcdf.defVar(ncid_CLIM,'ML_DISTANCE_AVG_ZONAL','short',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_dist_Z_ID,'_FillValue',int16(-32000))
    ML_dist_M_ID = netcdf.defVar(ncid_CLIM,'ML_DISTANCE_AVG_MERIDIONAL','short',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_dist_M_ID,'_FillValue',int16(-32000))
    ML_num_data_ID = netcdf.defVar(ncid_CLIM,'ML_NUM_PROFS','short',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_num_data_ID,'_FillValue',int16(0))
    ML_drho_ID = netcdf.defVar(ncid_CLIM,'ML_DELTA_RHO_SCALE','float',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_drho_ID,'_FillValue',NaN('single'))
    ML_year_ID= netcdf.defVar(ncid_CLIM,'ML_YEAR_OF_DATA','short',[dimid_months dimid_gridLAT dimid_gridLON]);
    netcdf.putAtt(ncid_CLIM,ML_year_ID,'_FillValue',int16(-32000))
    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['ML_' field '_ID  = netcdf.defVar(ncid_CLIM,''ML_' field ''',''float'',[dimid_months dimid_gridLAT dimid_gridLON]);']); %=======================	MODIF
	eval(['netcdf.putAtt(ncid_CLIM,ML_' field '_ID,''_FillValue'',NaN(''single''))']);
	eval(['ML_' field '_wm_ID  = netcdf.defVar(ncid_CLIM,''ML_' field '_WEIGHTED_MEAN'',''float'',[dimid_months dimid_gridLAT dimid_gridLON]);']); %=======================	MODIF
	eval(['netcdf.putAtt(ncid_CLIM,ML_' field '_wm_ID,''_FillValue'',NaN(''single''))']);

	eval(['ML_' field '_wstd_ID  = netcdf.defVar(ncid_CLIM,''ML_' field '_WEIGHTED_STD'',''float'',[dimid_months dimid_gridLAT dimid_gridLON]);']); %=======================	MODIF
	eval(['netcdf.putAtt(ncid_CLIM,ML_' field '_wstd_ID,''_FillValue'',NaN(''single''))']);

	eval(['ML_' field '_wSE_ID  = netcdf.defVar(ncid_CLIM,''ML_' field '_WEIGHTED_SE'',''float'',[dimid_months dimid_gridLAT dimid_gridLON]);']); %=======================	MODIF
	eval(['netcdf.putAtt(ncid_CLIM,ML_' field '_wSE_ID,''_FillValue'',NaN(''single''))']);

	eval(['ML_' field '_wTrend_ID  = netcdf.defVar(ncid_CLIM,''ML_' field '_WEIGHTED_TREND'',''float'',[dimid_months dimid_gridLAT dimid_gridLON]);']); %=======================	MODIF
	eval(['netcdf.putAtt(ncid_CLIM,ML_' field '_wTrend_ID,''_FillValue'',NaN(''single''))']);

	eval(['ML_' field '_wX_corr_ID  = netcdf.defVar(ncid_CLIM,''ML_' field '_WEIGHTED_X_corr'',''float'',[dimid_months dimid_gridLAT dimid_gridLON]);']); %=======================	MODIF
	eval(['netcdf.putAtt(ncid_CLIM,ML_' field '_wX_corr_ID,''_FillValue'',NaN(''single''))']);
    end
    % ML_year_wm_ID = netcdf.defVar(ncid_CLIM,'ML_YEAR_OF_DATA_WEIGHTED_MEAN','short',[dimid_months dimid_gridLAT dimid_gridLON]);
    % netcdf.putAtt(ncid_CLIM,ML_year_wm_ID,'_FillValue',int16(-32000))
    disp('end: Define new variables for ML')
    
    % Leave define mode and enter data mode to write data.
    netcdf.endDef(ncid_CLIM)
    
    % Write constants and grid:
    disp('Write constants and grid')
    netcdf.putVar(ncid_CLIM,latID,y(:,1));
    netcdf.putVar(ncid_CLIM,lonID,x(1,:));
    netcdf.putVar(ncid_CLIM,sigID,sg0_pleth);
    netcdf.putVar(ncid_CLIM,timID,(1:tim)/tim*12);
    netcdf.putVar(ncid_CLIM,iqrID,iqr_mult);
    netcdf.putVar(ncid_CLIM,max_numID,max_num_profs);
    
    netcdf.putVar(ncid_CLIM,max_num2ID,max_profs_inmem);
    netcdf.putVar(ncid_CLIM,min_numID,min_num_profs);
    netcdf.putVar(ncid_CLIM,rnd_numID,num_of_rand_profiles);
    netcdf.putVar(ncid_CLIM,ctd_numID,min_CTD_profs);
    netcdf.putVar(ncid_CLIM,horiz_scaleID,horiz_scale_set);
    netcdf.putVar(ncid_CLIM,temporal_scaleID,time_scale);
    
    
    netcdf.putVar(ncid_CLIM,noise_signalID,signal_to_noise-1);
    netcdf.putVar(ncid_CLIM,decadal_scaleID,dec_scale);
    netcdf.putVar(ncid_CLIM,decadal_noiseID,dec_noise_offset);
    disp('end: Write constants and grid')
    
else
    disp('si Clim.nc existe: ouverture')
    ncid_CLIM = netcdf.open(output_fname,'WRITE');
    % get variable IDs in file.
    latID  = netcdf.inqVarID(ncid_CLIM,'LATITUDE');
    lonID  = netcdf.inqVarID(ncid_CLIM,'LONGITUDE');
    sigID  = netcdf.inqVarID(ncid_CLIM,'NEUTRAL_DENSITY');
    timID  = netcdf.inqVarID(ncid_CLIM,'MONTH');
    iqrID  = netcdf.inqVarID(ncid_CLIM,'IQR_RANGE');
    max_numID    = netcdf.inqVarID(ncid_CLIM,'MAX_NUM_PROFS');
    if compute_interior
        data_pres_ID = netcdf.inqVarID(ncid_CLIM,'PRESSURE');
        data_salt_ID = netcdf.inqVarID(ncid_CLIM,'SALINITY');
        data_temp_ID = netcdf.inqVarID(ncid_CLIM,'TEMPERATURE');
        data_pres_wm_ID = netcdf.inqVarID(ncid_CLIM,'PRESSURE_WEIGTHED_MEAN');
        data_salt_wm_ID = netcdf.inqVarID(ncid_CLIM,'SALINITY_WEIGTHED_MEAN');
        data_temp_wm_ID = netcdf.inqVarID(ncid_CLIM,'TEMPERATURE_WEIGTHED_MEAN');
        data_error_ID= netcdf.inqVarID(ncid_CLIM,'ERROR');
        data_area_ID = netcdf.inqVarID(ncid_CLIM,'AREA');
        data_nwgt_ID = netcdf.inqVarID(ncid_CLIM,'WEIGHT_SUMMED');
        data_dist_Z_ID = netcdf.inqVarID(ncid_CLIM,'DISTANCE_AVG_ZONAL');
        data_dist_M_ID = netcdf.inqVarID(ncid_CLIM,'DISTANCE_AVG_MERIDIONAL');
        data_num_data_ID = netcdf.inqVarID(ncid_CLIM,'NUM_PROFS');
    end
    ML_pres_ID  = netcdf.inqVarID(ncid_CLIM,'ML_MAX_PRESSURE');
    ML_dens_ID  = netcdf.inqVarID(ncid_CLIM,'ML_DENSITY');
    ML_salt_ID  = netcdf.inqVarID(ncid_CLIM,'ML_SALINITY');
    ML_temp_ID  = netcdf.inqVarID(ncid_CLIM,'ML_TEMPERATURE');
    ML_pres_wm_ID  = netcdf.inqVarID(ncid_CLIM,'ML_PRESSURE_WEIGHTED_MEAN');
    ML_dens_wm_ID  = netcdf.inqVarID(ncid_CLIM,'ML_DENSITY_WEIGHTED_MEAN');
    ML_error_ID= netcdf.inqVarID(ncid_CLIM,'ML_ERROR');
    ML_area_ID = netcdf.inqVarID(ncid_CLIM,'ML_AREA');
    ML_nwgt_ID = netcdf.inqVarID(ncid_CLIM,'ML_WEIGHT_SUMMED');
    ML_dist_Z_ID = netcdf.inqVarID(ncid_CLIM,'ML_DISTANCE_AVG_ZONAL');
    ML_dist_M_ID = netcdf.inqVarID(ncid_CLIM,'ML_DISTANCE_AVG_MERIDIONAL');
    ML_num_data_ID = netcdf.inqVarID(ncid_CLIM,'ML_NUM_PROFS');
    ML_drho_ID  = netcdf.inqVarID(ncid_CLIM,'ML_DELTA_RHO_SCALE');
    ML_salt_wm_ID  = netcdf.inqVarID(ncid_CLIM,'ML_SALINITY_WEIGHTED_MEAN');
    ML_temp_wm_ID  = netcdf.inqVarID(ncid_CLIM,'ML_TEMPERATURE_WEIGHTED_MEAN');
    ML_year_ID  = netcdf.inqVarID(ncid_CLIM,'ML_YEAR_OF_DATA');

    for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['ML_' field '_ID     = netcdf.inqVarID(ncid_CLIM, ''ML_' field ''');']); 
	eval(['ML_' field '_wm_ID  = netcdf.inqVarID(ncid_CLIM, ''ML_' field '_WEIGHTED_MEAN'');']); 

    end
    disp('end de l ouverture de Clim.nc')
end


[~, seed_rndprofs] = sort(rand(1,max_profs_inmem - (numsortedprofs)));
seed_rndprofs = seed_rndprofs + numsortedprofs;
time_scale = time_scale^2;

%% ALL INITIALIZED FOR WARP ... hmmm COMPUTATION  -   engage

% take time of start

tstart = now;
%disp([' Initializing complete at ' datestr(tstart,'HH:MM:SS') ', starting computation.'])%sunke version
warning('off','MATLAB:lscov:RankDefDesignMat')

if Parallel==1
 nbphyscore=feature('numcores');
 parpool(nbphyscore);
end

profs_inmem = (-100000:-100000+max_profs_inmem-1)-9e12; % set to 10,000 unique values not within IDs from floats or CTD ... ...

varINlist={'doKuusela','signal_to_noise','optionFit','Kuusela_fname', 'logfile', 'tstart','n','continue_upto','continue_at','dummy_T','ML_error','ML_weight','ML_num_used','ML_r_dist','ML_Drho','ML_area','ML_dens_ID' , 'ML_pres_ID', 'thetasOpt_wm_ID', 'thetatOpt_wm_ID','thetadOpt_wm_ID','sigmaOpt_wm_ID',  'ML_dens_wm_ID', 'ML_pres_wm_ID', 'ML_salt_wm_ID',  'ML_temp_wm_ID','ML_dens_wstd_ID',  'ML_pres_wstd_ID',  'ML_salt_wstd_ID','ML_temp_wstd_ID','ML_dens_wSE_ID',  'ML_pres_wSE_ID',  'ML_salt_wSE_ID', 'ML_temp_wSE_ID','ML_dens_wTrend_ID', 'ML_pres_wTrend_ID','ML_salt_wTrend_ID', 'ML_temp_wTrend_ID','ML_dens_wX_corr_ID',  'ML_pres_wX_corr_ID','ML_salt_wX_corr_ID', 'ML_temp_wX_corr_ID','ML_yprc5_ID','ML_yprc50_ID',  'ML_yprc95_ID',  'ML_salt_ID',  'ML_temp_ID',  	'ML_nwgt_ID',  	'ML_year_ID','ML_dist_M_ID',  'ML_dist_Z_ID', 'ML_error_ID',  'ML_area_ID',  'ML_num_data_ID',  'ML_drho_ID','sg0_pleth','opts','wfkt_multiplier','MLpres','MLpresError','MLdens','MLsalt','MLtemp','doIQR','num_of_rand_profiles','numsortedprofs','seed_rndprofs','rnd_within','max_num_profs','yeartrendmin','eps999','eps100','tim','dec_dyr','dyr','dyr_map','sg_th','sg_sa','CTD_offset','profs_inmem','drawit','min_CTD_profs','dist_dummy','maxdist_coarsegrid','maxdist_finegrid','min_num_profs','FM_ALL','FM_ALL_C','FM_LAT','FM_LON','x','y','continue_lon_at','continue_lon_upto','upperocean_limit','FM_DPT','FM_DPT','lat','lon','max_profs_inmem','idv','wod_ctd_data','FM_fname','argo_float_data','AdditionalFields','horiz_scale','time_scale','upperocean'};
for ifield=1:length(AdditionalFields);
		field=AdditionalFields{ifield};
		eval(['varINlist=[varINlist {''ML' field '''}];'])
		eval(['varINlist=[varINlist {''ML_' field '_ID''}];'])
		eval(['varINlist=[varINlist {''ML_' field '_wm_ID''}];'])
		eval(['varINlist=[varINlist {''ML_' field '_wstd_ID''}];'])
		eval(['varINlist=[varINlist {''ML_' field '_wSE_ID''}];'])
		eval(['varINlist=[varINlist {''ML_' field '_wX_corr_ID''}];'])
		eval(['varINlist=[varINlist {''ML_' field '_wTrend_ID''}];'])
		eval(['varINlist=[varINlist {''ML_' field '_gr''}];'])
		eval(['varINlist=[varINlist {''ML_' field '_gr_wmean''}];'])
		eval(['varINlist=[varINlist {''ML_' field '_gr_wstd''}];'])
		eval(['varINlist=[varINlist {''ML_' field '_gr_wSE''}];'])
		eval(['varINlist=[varINlist {''ML_' field '_gr_wTrend''}];'])
		eval(['varINlist=[varINlist {''ML_' field '_gr_wX_corr''}];'])
end

for ii=1:length(varINlist)
		eval([ 'INpara.' varINlist{ii}  '=' varINlist{ii}  ';']);
end


disp('LATITUDE range')
y(:,1)

if Parallel==1
    if reverselat==1
        %%%%%% 
        %reverse latitude 
        i_vec=m:-1:1;
        di=12;
        for iii=1:di:length(i_vec);
            if iii==length(i_vec); i_vec_current=i_vec(iii);
            else 
                i_vec_current=max(i_vec(iii)-(di-1),1):i_vec(iii);
            end
        %%%%%% 


            INpara=INpara;
            parfor istep = i_vec_current 	%Latitude loop 140 times / 280 times if 0.5? ...
            %for istep = 1:m 	%Latitude loop 140 times / 280 times if 0.5? ...
                out{istep}=Sunke_OI_Trend_Final_Sub1(INpara,varINlist,istep)
            end

            disp('Write all now...')
            tic
            for istep = i_vec_current
                for jstep=1:n %longitude loop 360 times / 720 times if 0.5? ...	
                    if ~isnan(out{istep}{jstep}.istep)
                        for ivar=1:length(out{istep}{jstep}.idvar)
                            netcdf.putVar(ncid_CLIM,out{istep}{jstep}.idvar{ivar},[0 out{istep}{jstep}.istep-1 out{istep}{jstep}.jstep-1],[12 1 1],out{istep}{jstep}.var{ivar});
                        end
                    end
                end
            end
            toc
        end
    else
        %%%%%% 
        %keep latitude:
        i_vec=1:m;
        di=12;
        for iii=1:di:length(i_vec);
            if iii==length(i_vec); i_vec_current=i_vec(iii);
            else 
                i_vec_current=i_vec(iii):min(i_vec(iii)+(di-1),length(i_vec));
            end
        %%%%%%% 

            INpara=INpara;
            parfor istep = i_vec_current 	%Latitude loop 140 times / 280 times if 0.5? ...
            %for istep = 1:m 	%Latitude loop 140 times / 280 times if 0.5? ...
                out{istep}=Sunke_OI_Trend_Final_Sub1(INpara,varINlist,istep)
            end

            disp('Write all now...')
            tic
            for istep = i_vec_current
                for jstep=1:n %longitude loop 360 times / 720 times if 0.5? ...	
                    if ~isnan(out{istep}{jstep}.istep)
                        for ivar=1:length(out{istep}{jstep}.idvar)
                            netcdf.putVar(ncid_CLIM,out{istep}{jstep}.idvar{ivar},[0 out{istep}{jstep}.istep-1 out{istep}{jstep}.jstep-1],[12 1 1],out{istep}{jstep}.var{ivar});
                        end
                    end
                end
            end
            toc
        end
    end    
else
    if reverselat==1
        %%%%%% 
        %reverse latitude 
        i_vec=m:-1:1;
        di=12;
        for iii=1:di:length(i_vec);
            if iii==length(i_vec); i_vec_current=i_vec(iii);
            else 
                i_vec_current=max(i_vec(iii)-(di-1),1):i_vec(iii);
            end
        %%%%%% 


            INpara=INpara;
            %parfor istep = i_vec_current 	%Latitude loop 140 times / 280 times if 0.5? ...
            for istep = 1:m 	%Latitude loop 140 times / 280 times if 0.5? ...
                out{istep}=Sunke_OI_Trend_Final_Sub1(INpara,varINlist,istep)
            end

            disp('Write all now...')
            tic
            for istep = i_vec_current
                for jstep=1:n %longitude loop 360 times / 720 times if 0.5? ...	
                    if ~isnan(out{istep}{jstep}.istep)
                        for ivar=1:length(out{istep}{jstep}.idvar)
                            netcdf.putVar(ncid_CLIM,out{istep}{jstep}.idvar{ivar},[0 out{istep}{jstep}.istep-1 out{istep}{jstep}.jstep-1],[12 1 1],out{istep}{jstep}.var{ivar});
                        end
                    end
                end
            end
            toc
        end
    else
        %%%%%% 
        %keep latitude:
        i_vec=1:m;
        di=12;
        for iii=1:di:length(i_vec);
            if iii==length(i_vec); i_vec_current=i_vec(iii);
            else 
                i_vec_current=i_vec(iii):min(i_vec(iii)+(di-1),length(i_vec));
            end
        %%%%%%% 

            INpara=INpara;
            %parfor istep = i_vec_current 	%Latitude loop 140 times / 280 times if 0.5? ...
            for istep = 1:m 	%Latitude loop 140 times / 280 times if 0.5? ...
                out{istep}=Sunke_OI_Trend_Final_Sub1(INpara,varINlist,istep)
            end

            disp('Write all now...')
            tic
            for istep = i_vec_current
                for jstep=1:n %longitude loop 360 times / 720 times if 0.5? ...	
                    if ~isnan(out{istep}{jstep}.istep)
                        for ivar=1:length(out{istep}{jstep}.idvar)
                            netcdf.putVar(ncid_CLIM,out{istep}{jstep}.idvar{ivar},[0 out{istep}{jstep}.istep-1 out{istep}{jstep}.jstep-1],[12 1 1],out{istep}{jstep}.var{ivar});
                        end
                    end
                end
            end
            toc
        end
    end    
end


%%
delete(gcp('nocreate'))
disp(' Closing all open file handles ... ')
% close all handles
netcdf.close(ncid_FM);
netcdf.close(ncid_CLIM);
if useargo
    netcdf.close(ncid_AFD);
end
netcdf.close(ncid_CTD);

disp(' All done! Have a nice day ... ')




% =-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-
% =-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-
%
%           SUB FUNCTIONS
%
% =-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-
% =-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-==-=-



