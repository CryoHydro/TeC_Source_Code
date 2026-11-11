%% SITE CONFIGURATION - Kyzylsu
% Read by Launcher_Dev_Distributed.m with run(); defines the structure CONFIG_vals.
% Empty values ([]) mean "not set".
% On the HPC, another config file can be used without editing the launcher:
%   export TC_CONFIG=Kyzylsu/Parameters/<other_config>.m   (in the job script, before matlab)

CONFIG_vals = struct();

%% Site
CONFIG_vals.Site_num  = 4;          % Site number
CONFIG_vals.Site_name = 'Kyzylsu';  % Site name
CONFIG_vals.DeltaGMT  = 5;          % Change in hours from GMT
CONFIG_vals.UTM_zone      = 42;     % UTM zone for catchment
CONFIG_vals.UTM_zone_text = '42 N'; % UTM zone for catchment (text)

%% Run period and output folder
CONFIG_vals.Run_StartDate = datetime(2021,10,01, 0,0,0); % Run start date and time (year,month,day, hour,min,sec)
CONFIG_vals.Run_EndDate   = datetime(2022,09,30, 23,0,0); % Run end date and time (last time step)
CONFIG_vals.Run_Folder    = 'Kyzylsu_Run_02';           % Name of output run folder

%% Restart
CONFIG_vals.Restart_ID   = 0;                % 0 = new run, 1 = continue an unfinished run from its checkpoint
CONFIG_vals.Restart_Date = [];               % Date of the checkpoint to restart from, e.g. datetime(2021,10,25); [] = latest checkpoint
CONFIG_vals.Restart_Run  = 'Kyzylsu_Run_02'; % Name of run to restart from

%% Model options
CONFIG_vals.OPT_Forcing   = 4; % Type of forcing: 1 = from AWS which is distributed, 2 = fully gridded, 3 = on-the-fly downscaling,
                               % 4 = user forcing function (maps built at each time step)
CONFIG_vals.OPT_Veg_Param = 1; % Vegetation parameters (Ccrown, Cwat, Crock, Curb, Cbare):
                               % 1 = one value per landcover type (change the values in VPAR_T directly),
                               % 2 = gridded and variable per landcover (grids in the dtm_file; VPAR_T still
                               %     needed for the Veg_type names, but its values are not used)

%% Meteorological forcing (used if OPT_Forcing = 1)
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

%% User forcing function (used if OPT_Forcing = 4)
% Kyzylsu/Forcing/Kyzylsu_forcing.m builds the forcing as in the former launcher
% (Case_study/Kyzylsu_distributed): Ta lapsed between Pluvio and AWS_3900m, other
% variables from Pluvio. Station data: 01-Jun-2021 to 30-Sep-2023.
CONFIG_vals.forc_function = 'Kyzylsu_forcing'; % Name of the forcing function (on the path, in Kyzylsu/Forcing)

%% Solar integration interval (used if OPT_Forcing = 2 or 4)
CONFIG_vals.t_bef = 0.5; % Integration interval for solar variables: hours or fractions before (former launcher: 0.5)
CONFIG_vals.t_aft = 0.5; % Integration interval for solar variables: hours or fractions after (former launcher: 0.5)
% Also needed for OPT_Forcing = 2 or 3 (not set for Kyzylsu):
% CONFIG_vals.forc_meteo_folder = '';  % Sub-folder of Forcing/Bias_corrected and Forcing/Output_Downscaling
% CONFIG_vals.OPT_RunBC = 0;           % OPT_Forcing = 3 only
