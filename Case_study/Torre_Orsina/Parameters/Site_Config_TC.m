%% SITE CONFIGURATION 
% Read by Launcher_Dev_Distributed.m with run(); defines the structure CONFIG_vals.
% Empty values ([]) mean "not set".
% On the HPC, another config file can be used without editing the launcher:
%   export TC_CONFIG=Kyzylsu/Parameters/<other_config>.m   (in the job script, before matlab)

CONFIG_vals = struct();

%% Site
CONFIG_vals.Site_num  = 2;          % Site number
CONFIG_vals.Site_name = 'Torre_Orsina';  % Site name
CONFIG_vals.DeltaGMT  = 1;          % Change in hours from GMT
CONFIG_vals.UTM_zone      = 33;     % UTM zone for catchment
CONFIG_vals.UTM_zone_text = '33 T'; % UTM zone for catchment (text)

%% Run period and output folder
CONFIG_vals.Run_StartDate = datetime(1999,2,1, 0,0,0); % Run start date and time (year,month,day, hour,min,sec)
CONFIG_vals.Run_EndDate   = datetime(1999,3,31, 23,0,0); % Run end date and time (last time step)
CONFIG_vals.Run_Folder    = 'Run_33';           % Name of output run folder

%% Restart
CONFIG_vals.Restart_ID   = 0;                % 0 = new run, 1 = continue an unfinished run from its checkpoint
CONFIG_vals.Restart_Date = [];               % Date of the checkpoint to restart from, e.g. datetime(2021,10,25); [] = latest checkpoint
CONFIG_vals.Restart_Run  = 'Run_33'; % Name of run to restart from

%% Model options
CONFIG_vals.forc_function   = 'forcing_gridded'; % Type of forcing: 'forcing_station_lapse' = from AWS which is distributed, 'forcing_gridded' = fully gridded, 'forcing_downscaling' = on-the-fly downscaling,
                               % 4 = user forcing function (maps built at each time step)
CONFIG_vals.OPT_Veg_Param = 1; % Vegetation parameters (Ccrown, Cwat, Crock, Curb, Cbare):
                               % 1 = one value per landcover type (change the values in VPAR_T directly),
                               % 2 = gridded and variable per landcover (grids in the dtm_file; VPAR_T still
                               %     needed for the Veg_type names, but its values are not used)
CONFIG_vals.OPT_Decouple  = 0; % Glacier cooling effect on air temperature (Shaw et al., 2025):
                               % 1 = on (needs
                               % ModelCoefficients_Decoupling.mat and the
                               % FPL grid), 0 = off, 2= Simple Tmod (removes bias), 3=
                               % Complex Tmod (removes bias and multiplier)
CONFIG_vals.OPT_MPType = "ALL"; %If using "multipoint' type of Multipoint run (not "point') can choose type of points to run, POI = points of interest, 
                                %"CLUSTER" = cluster points, "ALL" = all points. These should be give in the
                                %Multipoints table under TYPE.

%% Meteorological forcing (used if OPT_Forcing = 1)
CONFIG_vals.Meteo_el        = 0;                      % Elevation of meteorological station [m a.s.l.]
CONFIG_vals.forc_meteo_file = ''; % File of the meteorological station data (.mat)
CONFIG_vals.forc_meteo_name = 't2m_3D_bc_c d2m_3D_bc_c tp_3D_bc_c RH_3D_c sp_3D_c ws10_3D_c strd_3D_c N_3D_c PARB_3D_c PARD_3D_c SAB1_3D_c SAB2_3D_c SAD1_3D_c SAD2_3D_c';               % Variable name of meteorological data in that file
CONFIG_vals.forc_meteo_folder = '2026_02_18';   %Folder containing gridded meteo data (only for OPT_Forcing = 2)
CONFIG_vals.ta_lapse_file   = '';  % Air temperature lapse rates per hour and month: table with
                                                          % columns hour, month and lapse rate (positive values) (.mat)
CONFIG_vals.ta_lapse_name   = '';              % Variable name within the air temperature lapse rate file
CONFIG_vals.pr_lapse_file   = ''; % Precipitation lapse rates: structure Pr_lapse(m).month,
                                                                 % one linear model per month (.mat)
CONFIG_vals.pr_lapse_name   = '';                 % Variable name within the precipitation lapse rate file

%% Preprocessed grids
CONFIG_vals.PreProc_folder = '2026_09_12'; % Folder containing preprocessing of grids (dtm + POI_dtm file from STEP1_Kyzylsu_100m.m)
CONFIG_vals.dtm_file       = 'dtm_Torre_Orsina_250m.mat';                  % File name of gridded data from preprocessing

%% User forcing function (used if OPT_Forcing = 4)
% Kyzylsu/Forcing/Kyzylsu_forcing.m builds the forcing as in the former launcher
% (Case_study/Kyzylsu_distributed): Ta lapsed between Pluvio and AWS_3900m, other
% variables from Pluvio. Station data: 01-Jun-2021 to 30-Sep-2023.
% CONFIG_vals.forc_function = 'Kyzylsu_forcing'; % Name of the forcing function (on the path, in Kyzylsu/Forcing)

%% Solar integration interval (used if OPT_Forcing = 2 or 4)
CONFIG_vals.t_bef = 1; % Integration interval for solar variables: hours or fractions before (former launcher: 0.5)
CONFIG_vals.t_aft = 0; % Integration interval for solar variables: hours or fractions after (former launcher: 0.5)
% Also needed for OPT_Forcing = 2 or 3 (not set for Shallap):
% CONFIG_vals.forc_meteo_folder = '';  % Sub-folder of Forcing/Bias_corrected and Forcing/Output_Downscaling
% CONFIG_vals.OPT_RunBC = 0;           % OPT_Forcing = 3 only
