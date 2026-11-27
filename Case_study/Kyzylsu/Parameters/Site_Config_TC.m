%% SITE CONFIGURATION - Kyzylsu
% Read by Launcher_Dev_Distributed.m with run(); defines the structure CONFIG_vals.
% Empty values ([]) mean "not set".
% Found by the launcher at <Root_Dir>/<site_name>/Parameters/Site_Config_TC.m.
% On the HPC, another config of this site can be used without editing the launcher:
%   export TC_CONFIG=<other_config>.m   (in the job script, before matlab; file in this Parameters/ folder)

CONFIG_vals = struct();

%% T&C code
% Folder of the T&C repository (Functions/, T&C_Code/, Forcing/). Optional: '' = the folder
% of the launcher. Ignored if TC_REPO_DIR is set in the job script.
% (The working directory and the site are set in the launcher or the job script, not here:
% they are needed to find this file.)
CONFIG_vals.REPO_Dir = 'C:\Users\ajoubert\Documents\GitHub\TeC_Source_Code';

%% Site
CONFIG_vals.Site_num  = 4;          % Site number
CONFIG_vals.Site_name = 'Kyzylsu';  % Site name
CONFIG_vals.DeltaGMT  = 5;          % Change in hours from GMT
CONFIG_vals.UTM_zone      = 42;     % UTM zone for catchment
CONFIG_vals.UTM_zone_text = '42 N'; % UTM zone for catchment (text)

%% Run period 
CONFIG_vals.Run_StartDate = datetime(2021,10,01, 1,0,0); % Run start date and time (year,month,day, hour,min,sec)
CONFIG_vals.Run_EndDate   = datetime(2021,10,30, 23,0,0); % Run end date and time (last time step)

% Run name and output directory
CONFIG_vals.Run_Folder    = 'Kyzylsu_test';           % Name of output run folder
CONFIG_vals.Output_Dir    = 'C:\Users\ajoubert\Documents\GitHub\TeC_Source_Code\Case_study\Kyzylsu\Outputs';                          % Folder where the run folder is created; '' = <site folder>/Outputs/Distributed. Ignored if TC_OUTPUT_DIR is set in the job script

%% Restart
CONFIG_vals.Restart_ID   = 0;                % 0 = new run, 1 = continue an unfinished run from its checkpoint
CONFIG_vals.Restart_Date = [];               % Date of the checkpoint to restart from, e.g. datetime(2021,10,25); [] = latest checkpoint
CONFIG_vals.Restart_Run  = 'Kyzylsu_Q_fixed'; % Name of run to restart from

%% Model options
CONFIG_vals.forc_function = 'forcing_gridded'; % Forcing function (maps built once per month):
                               % 'forcing_station_lapse' = one station + lapse rates (settings below),
                               % 'forcing_gridded' = monthly gridded files, 'forcing_downscaling' = on-the-fly
                               % downscaling of ERA5-Land, or a user function, e.g. 'Kyzylsu_forcing' (below)
CONFIG_vals.OPT_Veg_Param = 1; % Vegetation parameters (Ccrown, Cwat, Crock, Curb, Cbare):
                               % 1 = one value per landcover type (change the values in VPAR_T directly),
                               % 2 = gridded and variable per landcover (grids in the dtm_file; VPAR_T still
                               %     needed for the Veg_type names, but its values are not used)
CONFIG_vals.OPT_Decouple  = 1; % Glacier cooling effect on air temperature, over clean ice only: 0 = off,
                               % 1 = Shaw et al. (2025) (needs ModelCoefficients_Decoupling.mat and the FPL grid),
                               % 2 = bias TmodB where Ta > 0,
                               % 3 = bias TmodB, then multiplier TmodM where Ta > 0
                               %     (TmodB, TmodM in Options_and_NoVegParams.xlsx)
							   
CONFIG_vals.OPT_MPType = "POI"; %If using "multipoint' type of Multipoint run (not "point') can choose type of points to run, POI = points of interest, 
                               %"CLUSTER" = cluster points, "ALL" = all points. These should be give in the
                               %Multipoints table under TYPE.
							   
%CONFIG_vals.OPT_Forcing = 1; % Multipoint forcing option (temporarly)							   

%% Meteorological forcing (used by forcing_station_lapse)
CONFIG_vals.Meteo_el        = 3371;                      % Elevation of meteorological station [m a.s.l.]
CONFIG_vals.forc_meteo_file = 'Kyzylsu_Pluvio_meteo.mat'; % File of the meteorological station data (.mat)
CONFIG_vals.forc_meteo_name = 'Meteo_data';               % Variable name of meteorological data in that file
CONFIG_vals.ta_lapse_file   = 'Kyzylsu_hm_Ta_lapse.mat';  % Air temperature lapse rates per hour and month: table with
                                                          % columns hour, month and lapse rate (positive values) (.mat)
CONFIG_vals.ta_lapse_name   = 'hm_Ta_lapse';              % Variable name within the air temperature lapse rate file
CONFIG_vals.pr_lapse_file   = 'Kyzylsu_Pr_ratio_v_el_const.mat'; % Precipitation lapse rates: structure Pr_lapse(m).month,
                                                                 % one linear model per month (.mat)
CONFIG_vals.pr_lapse_name   = 'Pr_lapse';                 % Variable name within the precipitation lapse rate file

%% Preprocessed grids
CONFIG_vals.PreProc_folder = 'Kyzylsu_100m_STEP1'; % Folder containing preprocessing of grids (dtm + POI_dtm file from STEP1_Kyzylsu_100m.m)
CONFIG_vals.dtm_file       = 'dtm_Kyzylsu_100m.mat';                  % File name of gridded data from preprocessing

%% User forcing function Kyzylsu_forcing (forc_function = 'Kyzylsu_forcing')
% Kyzylsu/Forcing/Kyzylsu_forcing.m builds the forcing as in the former launcher
% (Case_study/Kyzylsu_distributed): Ta lapsed between Pluvio and AWS_3900m, other
% variables from Pluvio. Station data: 01-Jun-2021 to 30-Sep-2023. No settings needed.

%% Solar integration interval (used by all forcing functions except forcing_station_lapse)
CONFIG_vals.t_bef = 0.5; % Integration interval for solar variables: hours or fractions before (former launcher: 0.5)
CONFIG_vals.t_aft = 0.5; % Integration interval for solar variables: hours or fractions after (former launcher: 0.5)
% Also needed for forcing_gridded or forcing_downscaling (not set for Kyzylsu):
% CONFIG_vals.forc_meteo_name = '';    % forcing_gridded: 14 variable names; forcing_downscaling: <meteo_name> of the files
% CONFIG_vals.forc_meteo_folder = '';  % forcing_gridded: sub-folder of Forcing/Bias_corrected and Forcing/Output_Downscaling
% CONFIG_vals.OPT_RunBC = 0;           % forcing_downscaling: 1 = bias correction
