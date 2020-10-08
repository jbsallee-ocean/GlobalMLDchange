
optionFit=4
reverselat=1
parallel=1
doKuusela=0
yeartrendmin=1970;

% define files to use and output:
argo_float_data='/net/ether/data/proteo1/jbslod/Taf/LOCEAN/MLD/Database/Update2018/MLD003_Gobal_Argo_Seals_withError.nc';
wod_ctd_data = '/net/ether/data/proteo1/jbslod/Taf/LOCEAN/MLD/Database/Update2018/MLD003_Gobal_Ship_withError.nc';	% file arfo_seals avec calcul de Bt & Bs


FM_fname ='/net/ether/data/proteo1/jbslod/Data/Routine-Matlab/Sunke_OI/SO_Gebco_FM_05.nc';
if doKuusela==0
	if reverselat==1
		output_fname = ['/net/ether/data/proteo1/jbslod/Taf/LOCEAN/MLD/Climatology/Update2018/MLD003trend_from' int2str(yeartrendmin) '_Finalnewstd2_rev_Fit' int2str(optionFit) '.nc'];
	else
		output_fname = ['/net/ether/data/proteo1/jbslod/Taf/LOCEAN/MLD/Climatology/Update2018/MLD003trend_from' int2str(yeartrendmin) '_Finalnewstd2_Fit' int2str(optionFit) '.nc'];
	end
else
        if reverselat==1
                output_fname = ['/net/ether/data/proteo1/jbslod/Taf/LOCEAN/MLD/Climatology/Update2018/MLD003trend_from' int2str(yeartrendmin) '_doKuusela_Finalnewstd2_rev_Fit' int2str(optionFit) '.nc'];
        else 
                output_fname = ['/net/ether/data/proteo1/jbslod/Taf/LOCEAN/MLD/Climatology/Update2018/MLD003trend_from' int2str(yeartrendmin) '_doKuusela_Finalnewstd2_Fit' int2str(optionFit) '.nc'];
        end 
end
Kuusela_fname = ['/net/ether/data/proteo1/jbslod/Taf/LOCEAN/MLD/Climatology/Update2018/MLD003trend_from1970NoIQR_climOnly_KuuselaOA_withError.nc'];



%Plot usual fields and additional fields?
AdditionalFields={'NT15','NS15','NS200','NT200','N15','NMLD','N215'};
AdditionalFields={'NT15','NS15','NS200','NT200','N2200','N215'};

%Define grid on which to grid:
% This step is more about defining min/max lon/lat the actual grid will be the FM grid	 
lon_grid=0:0.5:360; % needs to 0-360 
lat_grid=-90:0.5:90; 

if doKuusela==0
	if reverselat==1
		logfile= ['MLDlog_Finalnewstd2_rev_optionFit' int2str(optionFit) '.tmp'];
	else
		logfile= ['MLDlog_Finalnewstd2_optionFit' int2str(optionFit) '.tmp'];
	end
else
        if reverselat==1
                logfile= ['MLDlog_Finalnewstd2_rev_doKuusela_optionFit' int2str(optionFit) '.tmp'];
        else
                logfile= ['MLDlog_Finalnewstd2_doKuusela_optionFit' int2str(optionFit) '.tmp'];
        end
end
% Misc options: 
eps100 = 0.0001;
eps999 = 1 - eps100;
create_new = true; % set to (true) if you start a new project, use (false) to continue writing data into existing nc-file.
compute_interior = false; % set to (false) in case you only want to compute ML, otherwise (true)
drawit = false; % while computing ... show current FM map of gridpoint in figure(1) (true)
useargo = true; % use argo floats + CTD data (true) or ONLY CTD data (false - not recommended for seasonal analysis)
upperocean = true; % if true, only compute down to following limit
upperocean_limit = 1850; %if depth of computed isopycnal exeeds this, move on to next month and/or grid point
wfkt_multiplier = 1.20;
max_profs_inmem = 3000; % ~maxnumprofs * 5
min_CTD_profs = 0; %At least use this many CTD profiles no matter how far away
max_num_profs = 500;% how many data points should be used per gridded point
min_num_profs =  10; % only compute isopycnal if 10 data points available - this is far to low, but can be stripped later
num_of_rand_profiles = 150; % number of random profiles within MAX_NUM_PROFS one fifth of max profs? what ever!
min_allowed_weight = 10^(-20); % has to be >0 otherwise we produce NaNs ...
rnd_within = 900; % take random profiles within closest 900 profiles ...
iqr_mult = 2; % inter-quartile range (IQR) filter width to be applied to last 350 points to be mapped

% only compute subsection of whole data file with following limits:
continue_at = -90;
continue_upto =90;
continue_lon_at = 0;
continue_lon_upto = 360;

% how many cores might the C-programs use?  ... (matlab itself always will use ALL)
OMP_NUM_THREADS_DAY = 4;  % enter number of cores your computer has
OMP_NUM_THREADS_NIGHT = 4;  % enter number of cores your computer has

% added to diagonal of covariance matrix:
signal_to_noise = 1.5; % 1.5 % added to diagonal ONE .... decadal weighting is added on top of this
%signal_to_noise = signal_to_noise + 1;

% covariance weighting
horiz_scale_set = 330; %330; % 550 % in km
time_scale  = 1.5 /12; % % 2  in years

% for mapping each month  set to 12:
tim = 12; % 12 number of steps per year

% for historical data add this noise to diagonal value:
centeryear = 2010;
dec_scale = 12; % yrs
dec_noise_offset = 8.5; 

doIQR=0;

%--------
options= [eps100 eps999 create_new compute_interior  drawit useargo upperocean upperocean_limit wfkt_multiplier  max_profs_inmem min_CTD_profs max_num_profs  min_num_profs num_of_rand_profiles  min_allowed_weight  rnd_within  iqr_mult  continue_at  continue_upto continue_lon_at continue_lon_upto  OMP_NUM_THREADS_DAY  OMP_NUM_THREADS_NIGHT  signal_to_noise horiz_scale_set time_scale  tim  centeryear dec_scale dec_noise_offset yeartrendmin doIQR optionFit reverselat parallel doKuusela];
%--------

Sunke_OI_Trend_Final_Main(argo_float_data, wod_ctd_data, FM_fname, Kuusela_fname, AdditionalFields, output_fname,  lon_grid, lat_grid, options,logfile);


