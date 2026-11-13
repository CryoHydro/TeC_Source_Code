%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C) - ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AUTHOR INFO AND STUDY SITE
%==========================================================================
% Updated 15/09/2026 to make into Development version for Mountain Water
% Author: Cat Fyffe
% Code originally from: Max Rodriguez and Achille Jouberton
% Code explanation: This code launches TC model.
%==========================================================================

%% CLEAR ALL
clc; clear;
% clear all   % Don't clear all such that it can be run on the HPC cluster without any issues 
% delete(gcp('nocreate'))

%=========================================================================
%% SITE SETUP
%==========================================================================

site_num_list = [1 2 3 4]; %Please number sequentially
site_name_list = {'Shallap' 'Torre_Orsina' 'Hybas_Lvl7_Catchment_2070497830','Kyzylsu'};
num_sites = size(site_num_list,2);
site_num = 4; %Choose site to run, this then selects the correct site directories
site_name = site_name_list{site_num};

%===================================================
%% Load config info
%==================================================

% Sub-path for Config file
% Config file: plain MATLAB script that defines the structure CONFIG_vals.
% Another file can be chosen without editing this launcher, e.g. in an HPC job script:
%   export TC_CONFIG=Kyzylsu/Parameters/Site_Config_Run03.m
Directories.config = getenv('TC_CONFIG');
if isempty(Directories.config)
    Directories.config = [site_name '/Parameters/Site_Config_TC.m'];
end
if ~isfile(Directories.config); error('Config file not found: %s', Directories.config); end
fprintf('Config file: %s\n', Directories.config);
run(Directories.config); % defines CONFIG_vals

IniCond.DeltaGMT = CONFIG_vals.DeltaGMT; %Define change from GMT
IniCond.run_folder = CONFIG_vals.Run_Folder; %Define folder to save outputs

%%======================================================================
% Model structure choices
%========================================================================

OPT_Forcing = CONFIG_vals.OPT_Forcing; %Choose type of forcing, 1= from AWS which is distributed, 2= fully gridded, 3= on-the-fly downscaling, 4= user forcing function (maps built at each time step)
OPT_Veg_Param = CONFIG_vals.OPT_Veg_Param; %Choose if vegetation parameters (for Ccrown, Cwat, Crock, Curb, Cbare) are one value per landcover type (1), or gridded and variable per landcover (2). 
%For 1 change the values in VPAR_T directly, for 2 provide the grids in the
%dtm_file. Note for 2 you still need the table to read the Veg_type names
%but the values are not used.
OPT_Decouple = CONFIG_vals.OPT_Decouple; %Glacier cooling effect on air temperature (Shaw et al., 2025): 1 = on, 0 = off

%% DIRECTORIES - ALL RELATIVE - NO NEED TO UPDATE 
%==========================================================================
% All paths are here
%==========================================================================

% Sub-path for outputs
Directories.save = [site_name '/Outputs/Distributed/' IniCond.run_folder '/']; %Do we make this site specific?
Directories.restart = Directories.save;

% Sub-path for forcings
Directories.forc = [site_name '/Forcing']; 
if OPT_Forcing == 1
Directories.forc_meteo = [Directories.forc '/' CONFIG_vals.forc_meteo_file]; meteo_name = CONFIG_vals.forc_meteo_name;
ta_lapse_file = CONFIG_vals.ta_lapse_file;  ta_lapse_name = CONFIG_vals.ta_lapse_name; %These are lapse rates per hour and month, in a table with a column with hour, a column with month, and a column with lapse rates (positive values)
pr_lapse_file = CONFIG_vals.pr_lapse_file;  pr_lapse_name = CONFIG_vals.pr_lapse_name; %This is a structure Pr_lapse(m).month where each month has a linear model
Meteo_el = CONFIG_vals.Meteo_el; %Elevation of station, used if OPT_forcing = 1
elseif OPT_Forcing == 2
Directories.forc_meteo = [Directories.forc '/Bias_corrected/' CONFIG_vals.forc_meteo_folder '/'];
Directories.forc_meteoOD = [Directories.forc '/Output_Downscaling/' CONFIG_vals.forc_meteo_folder '/'];
elseif OPT_Forcing == 3
meteo_name = CONFIG_vals.forc_meteo_name; %e.g. 'Hybas_Lvl7_Catchment_2070497830_250m'
elseif OPT_Forcing == 4
forc_function = CONFIG_vals.forc_function; %Name of the user forcing function, e.g. 'Kyzylsu_forcing' (file <site>/Forcing/Kyzylsu_forcing.m)
end

% Parameters
Directories.Vegpar = [site_name '/Parameters/Parameters_TC.xlsx']; 
Directories.Optpar = [site_name '/Parameters/Options_and_NoVegParams.xlsx'];
Directories.VPAR = [site_name,'/Parameters/VPAR_T.xlsx'];

% Preprocessing information
Directories.PreProc = [site_name '/Grids/' CONFIG_vals.PreProc_folder '/']; % Directories.PreProc = [site_name '/Grids/'];
dtm_file = CONFIG_vals.dtm_file;

% Sub-path points of interest for discharge and multipoints
Directories.POIs = [Directories.PreProc 'POI_' erase(dtm_file,'.mat') '.txt']; % POI_dtm_<Site>_<res>m.txt from the preprocessing (next to the dtm)
Directories.Multi = [site_name '/Multipoint/' site_name '_MultiPoints.txt'];   %

% Dependencies
addpath(genpath('Functions')); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([site_name '/Grids'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([site_name '/Forcing'])); % Add path to Ca_Data
addpath(genpath('T&C_Code')); %T&C source code 

% Modelling period and time
%--------------------------------------------------------------------------
dateRun.start = CONFIG_vals.Run_StartDate;  % Starting point of the simulation - start on the first of the month
dateRun.end = CONFIG_vals.Run_EndDate; % Last timestep of the simulation 

x1=datetime(dateRun.start); x2=datetime(dateRun.end);
Date = x1:hours(1):x2;
N_time_step=length(Date);
Nd_time_step = ceil(N_time_step/24)+1;

% Steps
dt=3600; %% [s]
dth=1; %% [h]
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE

% Restart option
%--------------------------------------------------------------------------
restart.id = CONFIG_vals.Restart_ID; % Set to 1 to continue an un-completed T&C run
restart.date = CONFIG_vals.Restart_Date; %Define iter to restart
restart.run = CONFIG_vals.Restart_Run;

%=========================================================================
%% LOCATION OF OUTPUTS - CREATION OF FOLDERS
%==========================================================================
if restart.id ~=1 %If not restarting
    if ~exist(Directories.save, 'dir') 
        disp('Folders do not exist for outputs. Creating new folders')
        mkdir(Directories.save);
        mkdir([Directories.save 'Initial']);
        mkdir([Directories.save 'Store']);
        mkdir([Directories.save 'Spatial_data']); %Create main subfolder where to store spatial results
        mkdir([Directories.save 'Parameters']); %Save outputs of parameters tables run
    addpath(genpath(Directories.save)); 
    end
end

outlocation = [Directories.save];

%=========================================================
%% LOAD SPATIAL DATA
%========================================================
load([Directories.PreProc dtm_file]) %This file is in 2_nputs

[m_cell,n_cell]=size(DTM);
num_cell=numel(DTM);
x_cell=xllcorner:cellsize:(xllcorner+cellsize*(n_cell-1));
y_cell=yllcorner:cellsize:(yllcorner+cellsize*(m_cell-1));

% Central lat and lon Lat
%--------------------------------------------------------------------------
UTM_Y = y(floor(length(y)/2));
UTM_X = x(floor(length(x)/2));

[Lat, Lon] = utm2deg(UTM_X, UTM_Y, CONFIG_vals.UTM_zone_text);

% Load point data for saving
%--------------------------------------------------------------------
% POI file written by the preprocessing together with the dtm
% (POI_dtm_<Site>_<res>m.txt, with at least the columns Name, col, row, idx)
if ~isfile(Directories.POIs)
    error('POI file not found: %s', Directories.POIs)
end
Points = readtable(Directories.POIs,'Delimiter',',');
% Column names differ between sites (e.g. name/Name): use Name, col, row, idx
vn = Points.Properties.VariableNames;
for c = {'Name','col','row','idx'}
    k = find(strcmpi(vn,c{1}),1);
    if isempty(k); error('Column "%s" missing in %s', c{1}, Directories.POIs); end
    Points.Properties.VariableNames{k} = c{1};
end
Points.Name = strtrim(string(Points.Name));
clear vn c k
disp(['Points: ' num2str(height(Points)) ' POIs read from ' Directories.POIs])
% Consistency with the loaded dtm
if any(Points.idx ~= sub2ind([m_cell n_cell], Points.row, Points.col))
    error('POI file: idx does not match row/col for the grid of %s', dtm_file)
end

% Remove points outside the catchment mask (no simulated cell there)
out_mask = MASK(Points.idx) ~= 1;
if any(out_mask)
    warning('Points outside MASK removed: %s', strjoin(Points.Name(out_mask),', '))
    Points(out_mask,:) = [];
end
clear out_mask
UTM_zone = CONFIG_vals.UTM_zone;

% names of points
Points_names = string(Points.Name);

%=========================================================================
%% MODEL PARAMETERS
%==========================================================================

%--------------------------------------
% 1. General options
%---------------------------------------

%Set numerical tolerances
Opt_CR=                     optimset('TolFun',1);%,'UseParallel','always'); %Numerical tolerance for internal CO2 computation
OPT_PH=                     odeset('AbsTol',0.01); %Numerical tolerance for internal Plant Hydraulic Volumes
OPT_SM=                     odeset('AbsTol',0.05,'MaxStep',dth); %Numerical tolerance for soil moisture differential equations
Opt_ST=                     optimset('TolFun',0.1);%,'UseParallel','always'); %Numerical tolerance for surface temperature computation
Opt_ST2 =                   optimset('TolFun',0.1,'Display','off');
OPT_STh=                    odeset('AbsTol',5e+3); %Numerical tolerance for heat transfer differential equations
OPT_VD=                     odeset('AbsTol',0.05); %Numerical tolerance for carbon budget differential equations

%Load general options and parameters
OPT_PARAM = readtable(Directories.Optpar);
OP_vals = OPT_PARAM.Properties.VariableNames == "Value";
OP_vals_id = find(OP_vals);
%Turn values into a structure
OPT_PARAM_vals = OPT_PARAM(:,OP_vals_id);
OPT_PARAM_vals = rows2vars(OPT_PARAM_vals); %Swap around
OPT_PARAM_vals = OPT_PARAM_vals(:,2:end); %Just removed 'values' column
OPT_PARAM_vals.Properties.VariableNames = OPT_PARAM{:,1}; %Give variable names the parameter names
OPT_PARAM_vals = table2struct(OPT_PARAM_vals); %Ok can now get to values using a structure notation :)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v2struct(OPT_PARAM_vals); % Only need this. No need to restate any other parameters from the loaded Excel configs/params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------------------------------
% 2. Vegetation parameters
%--------------------------------------------------------------------------
% Here are the parameters of the model for vegetation. opts check the format of the columns. Use opts to force all the columns with values to be numeric. %}

opts = detectImportOptions(Directories.Vegpar);
opts = setvartype(opts, [7:length(opts.VariableTypes)], 'double');
TT_par = readtable(Directories.Vegpar, opts);

% Vegetation parameters look up table
%--------------------------------------------------------------------------
% Codes from VEG_CODE based on predefined classification of vegetation
% Classes are represented by the vector II. 

VEG_CODE(isnan(VEG_CODE))=0; %If any NaNs replace with 0
ksv=reshape(VEG_CODE,num_cell,1);

%% Veg/land parameters (was POI but renamed to reduce confusion with Points
%of interest)
%Load table (needed for both options)
VPAR_T = readtable(Directories.VPAR);
vpar_list = unique(VPAR_T.Class); %Vegetation classes
No_vpar = size(vpar_list,1);  %Number of vegetation classes

%Check for C values = 1
Headers = VPAR_T.Properties.VariableNames;
Cvals = find(contains(Headers,{'Ccrowns','Cwat','Curb','Crock','Cbare'}));
Ctot = sum(VPAR_T{:,Cvals},2,'omitnan');
if any(Ctot~=1) 
    disp("Warning C values do not equal to 1");
end

% For both types
VPAR.Class = (vpar_list); %Just a list of the numbers
VPAR.Kbot = VPAR_T.Kbot;
VPAR.Krock = VPAR_T.Krock;

%Sort Crowns
Cro_n = contains(Headers,'Ccrowns');
Cro_n_num = sum(Cro_n);
IntT_cc = [VPAR_T{:,Cro_n}]; %array of values
Cro_cell = cell(No_vpar,1);
for v=1:No_vpar
    Cro_cell(v) = {IntT_cc(v,:)}; %This works
end
VPAR_T.Ccrowns = Cro_cell;

%Sort Veg_type
Veg_n = contains(Headers,'Veg_type');
Veg_n_num = sum(Veg_n);

VPAR_T = convertvars(VPAR_T,Veg_n,"cell");
inT_Vo = [VPAR_T{:,Veg_n}]; %array of values

Veg_st = strings(No_vpar,Veg_n_num);
for v=1:No_vpar
    Veg_st(v,:) = string(inT_Vo(v,:)); %This works
    VPAR_T.Veg_type{v} = Veg_st(v,:); %NOTE use of curly braces to take cell array into table
end

%Clean to remove empty data for Veg_type and Ccrowns
for v=1:No_vpar
    vt_temp = VPAR_T.Veg_type{v};
    cc_temp = VPAR_T.Ccrowns{v};
    em = strcmp(vt_temp,"");
    vt_temp(em)=[];
    cc_temp(em)=[];
    VPAR_T.Veg_type{v}=vt_temp;
    VPAR_T.Ccrowns{v}=cc_temp;
end
VPAR.Veg_type = VPAR_T.Veg_type;
VPAR.Ccrowns = VPAR_T.Ccrowns;

% ---- Vegetation parameters of each class, built once (used by PARAMETERS_SOIL_DEV) ----
% TT_vals{iv} holds the TT_par values of the vegetation types of class VPAR.Class(iv),
% one struct field per parameter (same result as the old per-cell code in PARAMETERS_SOIL_DEV)
TT_vals = cell(numel(VPAR.Class),1);
for iv = 1:numel(VPAR.Class)
    try
        [tv_tf,tv_loc] = ismember(TT_par.Properties.VariableNames,VPAR.Veg_type{iv});
        [~,tv_p] = sort(tv_loc(tv_tf)); % keep the order of the veg names
        tv_idx = find(tv_tf);
        tv_T = rows2vars(TT_par(:,tv_idx(tv_p)));
        tv_T = tv_T(:,2:end); % remove the Veg name column
        tv_T.Properties.VariableNames = TT_par{:,TT_par.Properties.VariableNames == "Parameters"};
        TT_vals{iv} = table2struct(tv_T,ToScalar=true);
    catch tv_ME
        warning('No vegetation parameters for class %d (%s): %s', VPAR.Class(iv), strjoin(VPAR.Veg_type{iv},', '), tv_ME.message);
    end
end
clear iv tv_tf tv_loc tv_p tv_idx tv_T tv_ME

if OPT_Veg_Param==1
    %In this case there are single values per land cover
    %Turn into a grid
    VPAR.Ccrown = cell(m_cell,n_cell);
    VPAR.Cwat = zeros(m_cell,n_cell);
    VPAR.Curb = zeros(m_cell,n_cell);
    VPAR.Crock = zeros(m_cell,n_cell);
    VPAR.Cbare = zeros(m_cell,n_cell);
    for v=1:No_vpar
        id_land = VEG_CODE==VPAR_T.Class(v);
        VPAR.Ccrown(id_land)=VPAR_T.Ccrowns(v);
        VPAR.Cwat(id_land)=VPAR_T.Cwat(v);
        VPAR.Curb(id_land)=VPAR_T.Curb(v);
        VPAR.Crock(id_land)=VPAR_T.Crock(v);
        VPAR.Cbare(id_land)=VPAR_T.Cbare(v);
    end
elseif OPT_Veg_Param==2
    %Using grids now (which are created in pre-processing and loaded already), VPAR will be a structure
    VPAR.Ccrown = CCROWN;
    VPAR.Cwat = CWATER;
    if exist('CURB','var') == 1
        VPAR.Curb = CURB;
    else
    VPAR.Curb = CWATER.*0; %As urban 0 everywhere
    end
    VPAR.Crock = CROCK;
    VPAR.Cbare = CBARE;
end  

%Reshape
VPAR.Ccrownr = reshape(VPAR.Ccrown,num_cell,1);
VPAR.Cwatr = reshape(VPAR.Cwat,num_cell,1);
VPAR.Curbr = reshape(VPAR.Curb,num_cell,1);
VPAR.Crockr = reshape(VPAR.Crock,num_cell,1);
VPAR.Cbarer = reshape(VPAR.Cbare,num_cell,1);

%Set cc_max
cc_max = 1;
for z = 1:No_vpar
    if cc_max < length(VPAR_T.Ccrowns{z})
    cc_max = length(VPAR_T.Ccrowns{z});
    end
end

% Crown fractions as a numeric matrix [num_cell x cc_max] for the outputs
% (with OPT_Veg_Param = 1, VPAR.Ccrownr is a cell array of crown vectors)
if iscell(VPAR.Ccrownr)
    Ccrown_OUTv = zeros(num_cell,cc_max);
    for i = 1:num_cell
        Ccrown_OUTv(i,1:numel(VPAR.Ccrownr{i})) = VPAR.Ccrownr{i};
    end
else
    Ccrown_OUTv = VPAR.Ccrownr;
end

% SPATIAL INDICES PER LAND COVER CLASS
for v=1:No_vpar
    fieldName = ['Veg', num2str(v)];
    idxCode.(fieldName) = find(VEG_CODE == v & MASK == 1);
end

%Save the parameter files with the model outputs

writetable(TT_par,[Directories.save 'Parameters/' 'TT_par.xlsx']);
writetable(OPT_PARAM,[Directories.save 'Parameters/' 'OPT_PARAM.xlsx']);
writetable(VPAR_T,[Directories.save 'Parameters/' 'VPAR_T.xlsx']);
copyfile(Directories.config,[Directories.save 'Parameters/' 'Site_Config_TC.m']); % exact config used for this run

%-----------------------------------------------------------------------
% 3. Glacier and snow parameters
% --------------------------------------------------------------------------

% %Tmod parameters
% OPT_Tmod = OPT_PARAM_vals.OPT_Tmod;
% if OPT_Tmod == 0 %No Tmod applied
%     TmodB = 0;
% elseif OPT_Tmod == 1 %Simple bias correction
%     TmodB = OPT_PARAM_vals.TmodB; %Tmod bias (positive reduces temperature over ice)
% elseif OPT_Tmod == 2 %Uses a relationship with optional bias correction 
%     TmodB = OPT_PARAM_vals.TmodB; %Tmod bias (positive reduces temperature over ice)
%     TmodM = OPT_PARAM_vals.TmodM; %Tmod multiplier
%     TmodC = OPT_PARAM_vals.TmodC; %Tmod coefficent
% end

% Glacier dynamics and avalanching
OPT_Idyn = OPT_PARAM_vals.OPT_Idyn; % Switch for glacier dynamics
OPT_Aval = OPT_PARAM_vals.OPT_Aval; % 1 to turn on avalanching, 0 to turn off avalanching
a_aval = OPT_PARAM_vals.a_aval; %avalanche parameters a (Bernhart & Schulz 2010) - from TOPKAPI 0.17245 - 0.12/145 from Jouberton et al.(2025) - 99.05/0.1012 from Buri et al. (2023)
C_aval = OPT_PARAM_vals.C_aval; % avalanche parameters C 

% Precipitation phase partitioning
%1 = 2-threshold, 2 = Ding 2017, 3 = single-threshold, 4 = Pomeroy 2013, 5
%= Wang 2019, 6 = Jennings 2018
parameterize_phase.OPT_Pr_Part = OPT_PARAM_vals.OPT_Pr_Part; % Choice of the precipitation phase scheme
parameterize_phase.Tmax = OPT_PARAM_vals.Tmax; % Maximum air temperature for precipitation phase scheme (2-dual thresholds)
parameterize_phase.Tmin = OPT_PARAM_vals.Tmin; % Minimum air temperature for precipitation phase scheme (2-dual thresholds)
parameterize_phase.Tconst = OPT_PARAM_vals.Tconst; % Air temperature for constant thresholds
parameterize_phase_labels = {'2-Ta','Ding','1-Ta','Pomeroy','Wang','Jennings'};
parameterize_phase_label = parameterize_phase_labels(parameterize_phase.OPT_Pr_Part);

%Multilayer snowpack parameters
% hSTL = OPT_PARAM_vals.hSTL; % Skin layer thickness for the 2-layer snowpack module
% min_SPD = OPT_PARAM_vals.min_SPD;  %% [m] minimum snow pack depth to have a multilayer snow 

% Choice of the snow albedo scheme Evars_load
% OPT_Albsno = OPT_PARAM_vals.OPT_Albsno; % 3 doesn't work, 4 is Brock 2000, 5 is Ding 2017

% Glacier set up
if OPT_Idyn == 0
    GLH(GLH>0) = GLH(GLH>0) +400;  % Only use when ice dynamics are off, to avoid glacier disappearance
end
GLH=GLH.*(GLA_MAP2>0); %Prevent ice thickness outside glaciers
GLHn = reshape(GLH,num_cell,1);
DEB_MAPn = reshape(DEB_MAP,num_cell,1);

% PATCH, check pre-processing!
GLA_ID(MASK==0) = NaN;
GLA_ID(GLH==0) = NaN;

% Compute bedrock DEM for ice flow
DTM_Bedrock = DTM_orig-GLH; 

% Firn albedo - Albedo vs Elevation
% OPT_Alb_v_el = OPT_PARAM_vals.OPT_Alb_v_el; %Switch for using firn albedo versus elevation
% Aice = OPT_PARAM_vals.Aice; %Aice, Artesonraju % 0.2314 Shallap point value
% Amax = OPT_PARAM_vals.Amax; %Max 'firn' albedo value

if OPT_Alb_v_el==1

    % Alb_v_el_MM = OPT_PARAM_vals.Alb_v_el_MM; %Yota remote sensing slope for 25% lowest scenes
    % Alb_v_el_CC = OPT_PARAM_vals.Alb_v_el_CC; %Yota remote sensing intercept for 25% lowest scenes

    Afirn = (DTM.*Alb_v_el_MM) + Alb_v_el_CC; %Construct initial relationship
    Afirn(Afirn<Aice)=Aice; %Prevent values lower than Aice
    Afirn(Afirn>Amax)=Amax;
      
   disp('Albedo varies with elevation')
else
    Afirn = DTM.*0 + Aice; %Artesonraju % 0.2314 Shallap point value
    disp(strcat('Constant bare-ice albedo of',{' '},num2str(Aice),{' '},'used'))
end

% Check if initial snow depth is given, if not assign it to 0
if ~exist('SNOWD','var')
    SNOWD = DTM.*0;
end 

% Check if initial snow albedo is given - this goes into initial conditions as SNOWALB
%SNOWD is taken directly and used in intial conditions as initial snow depth
if ~exist('SNOWALB','var')
    % Ini_Asnow = OPT_PARAM_vals.Ini_Asnow;
    SNOWALB = SNOWD;
    SNOWALB(SNOWD>0) = Ini_Asnow;
end 

%-----------------------------------------------------------------
% 4. Carbon
%-------------------------------------------------------------------------

load('Ca_Data.mat');
% load([meteo_name '_Ca_Data.mat']); % TOM's FILE NAMING HAS SITE IN IT - PROVIDED BY THE PRE_PROCESSING SCRIPTS AT STEP 3

d1 = find(abs(Date_CO2-datenum(Date(1)))<1/36);d2 = find(abs(Date_CO2-datenum(Date(end)))<1/36);
Ca=Ca(d1:d2);
clear d1 d2 Date_CO2 id_df
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -

% Variables unknown for now
a_dis=NaN; pow_dis=NaN;
clear a0 gam1 pow0 k2 DTii

%------------------------------------------------------------------------
% 5. Topography and channel params
%-----------------------------------------------------------------------

MASKn=reshape(MASK,num_cell,1);
DTMn = reshape(DTM,num_cell,1);
Kinde = find(MASK==1);

%%% Slo_top [Fraction] %%% Aspect [rad] from N
[Slo_top,Aspect]=Slope_Aspect_indexes(DTM_orig,cellsize,'mste');
Aspect(isnan(Aspect))=0;
Slo_top(Slo_top<0.001)=0.001;  %to avoid flow routing issues
%%%
Asur=(1./cos(atan(Slo_top))); %% Effective Area / Projected Area
Asur=reshape(Asur,num_cell,1);
aTop= 1000*ones(m_cell,n_cell)*(cellsize^2)/cellsize; %% [mm] Area/Contour length ratio
Ared=ones(num_cell,1);

% Flow Boundary Condition
%--------------------------------------------------------------------------
%"Xoutlet" & "Youtlet": outlet point, predefined in dtm_XXX.mat-file
% Cells where the routing computes the discharge (QpointC, QpointH) and the
% flow velocity (UpointC, UpointH): the outlet, then every point of the POI
% file (e.g. the stream gauges), except the one(s) on the outlet cell.
% With a Point period in Output_Config.m they are written in the point files
% (Points/), with a Map period as time series (Series/, one column per cell).
% If no point is on the outlet cell, the outlet is added as a point
% (<site>_Outlet), so that its outputs are written too.
onOutlet = Points.row == Youtlet & Points.col == Xoutlet;
if ~any(onOutlet)
    newPt = Points(1,:);
    for c = 1:width(newPt)
        if isnumeric(newPt{1,c}); newPt{1,c} = NaN;
        elseif isstring(newPt{1,c}); newPt{1,c} = missing; end
    end
    newPt.Name = string([site_name '_Outlet']);
    newPt.col = Xoutlet; newPt.row = Youtlet; newPt.idx = sub2ind([m_cell n_cell], Youtlet, Xoutlet);
    Points = [Points; newPt];
    Points_names = string(Points.Name);
    onOutlet = [onOutlet; true];
    clear newPt c
end
Xout = [Xoutlet; Points.col(~onOutlet)]; % columns
Yout = [Youtlet; Points.row(~onOutlet)]; % rows
NAMEout = [string([site_name '_Outlet']); Points_names(~onOutlet)];
clear onOutlet

Slo_top(Youtlet,Xoutlet)=0.05; %NOTE on maps and in scatter use (X,Y) but to index use Y (row) X (col)
npoint = length(Xout);
Area= (cellsize^2)*sum(sum(MASK)); %% Projected area [m^2]

% Width channel
%--------------------------------------------------------------------------
WC = cellsize*ones(m_cell,n_cell); %% [m]  Width channel
WC(SN==1)=0.0018*sqrt((cellsize^2)*Aacc(SN==1)); %% [m]
WC=WC.*MASK;
SN(isnan(SN))=0; %% [Stream Identifier]
SNn=reshape(SN,num_cell,1);

NMAN_C=SN*0.040; NMAN_H=0.1; %%[s/(m^1/3)] manning coefficient
MRough = 0.01*(1-SN); %%[m] Microroughness
NMAN_C=NMAN_C.*MASK;
NMAN_H=NMAN_H.*MASK;
MRough=MRough.*MASK;
%%%
Kres_Rock =8760; %%[h] Bedrock aquifer constant
SPRINGn =SNn; %% Spring Location

%----------------------------------------------------------------------
% 6. SOIL PARAMETERS 
%--------------------------------------------------------------------------

%Soil layers and albedo
% rho_g = OPT_PARAM_vals.rho_g; %%% Spatial Albedo - this is quite high??
% md_max = OPT_PARAM_vals.md_max; % Number of debris layers

%Soil parameters spatial
Pss = [800]; Pwp = [3500]; %% [kPa]
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
[Osat,L,Pe,Ks,O33]=Soil_parameters_spatial(PSAN/100,PCLA/100,PORG/100);%NOTE was PORG/1000, but I already have inputs as %,presume inputs as fractions
[Ofc,Oss,Owp,Ohy]=Soil_parametersII_spatial(Osat,L,Pe,Ks,O33,Kfc,Pss,Pwp,Phy);

%Soil outputs
clear Pss Pwp Kfc Phy L Pe Ks O33 Ofc Oss Owp
Osat_OUT = Osat.*MASK; clear Osat
Osat_OUT = reshape(Osat_OUT,num_cell,1);
Ohy_OUT  = Ohy.*MASK; clear Ohy
Ohy_OUT  = reshape(Ohy_OUT,num_cell,1);

%Reshape soil variables
PSANr=reshape(PSAN/100,num_cell,1); %PSAN comes in here as a %,converted to a fraction (0 to 1)
PCLAr=reshape(PCLA/100,num_cell,1); %PCLA comes in here as a %,converted to a fraction (0 to 1)
if max(max(PORG))>40 %Presume PORG in g/kg not %
    PORGr=reshape(PORG/1000,num_cell,1);
else %PORG brought in as % 
    PORGr=reshape(PORG/100,num_cell,1);
end

% WHAT IS THIS?
%Zs_OUT=800*ones(num_cell,1);

%  ---  Set Soil thickness and related params ---
% OPT_soil_th = OPT_PARAM_vals.OPT_soil_th;
%First Zs and the number of soil layers is based on the maximum soil depth
Zs_all = [0 10 20 50 100 150 200 250 300 350 500 1000 1500 2000 2500 3000 3500 4000 4500]; %All potential soil depths (ms_max +1)
vi_all = [0 10 20 50 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100]; %All potential initial soil moistures (ms_max +1)
if OPT_soil_th == 1
    % soil_th_max = OPT_PARAM_vals.soil_th_max; %Max soil thickness, mm
    Zs = Zs_all(Zs_all<=soil_th_max);
    vi = vi_all(Zs_all<=soil_th_max);
    ms_max = size(Zs,2) - 1;
elseif OPT_soil_th == 2
    SOIL_TH = SOIL_TH.*10; %Convert cm to mm
    soil_th_max =  max(max(SOIL_TH)); %Max soil thickness, mm
    Zs = Zs_all(Zs_all<=soil_th_max); %Is it ok that the actual max thickness might not be included? Probably as these are at the top of the layer
    vi = vi_all(Zs_all<=soil_th_max);
    ms_max = size(Zs,2) - 1;
end

%Check
if max(max(SOIL_TH)) < 50 || max(max(SOIL_TH)) > 5000
    disp('Check soil units')
end

% Flow potential
%--------------------------------------------------------------------------
T_pot=cell(1,ms_max);
for jk=1:ms_max
    T_pot{jk}= T_flow;
end

%--------------------------------------------------------------------
% 7.Solar parameters
%------------------------------------------------------------------

% Computation Horizon Angle
[HZ,Zasp] = Horizon_Angle(DTM_orig,cellsize); %%% HZ Horizon angle array [angular degree], %%% Z Azimuth directions  [angular degree] from N
%%% Sky View Factor and Terrain Configuration Factor
[SvF,Ct] = Sky_View_Factor(DTM_orig,atan(Slo_top)*180/pi,Aspect,HZ,Zasp);

%==================================================================
%% INITIAL CONDITIONS
%================================================

if restart.id ~=1
out = [Directories.save 'Initial/INITIAL_CONDITIONS_' site_name '.mat'];
INIT_COND_v6(num_cell,m_cell,n_cell,...
   cc_max,ms_max,md_max,...
   MASKn,GLH,Ca,SNOWD,SNOWALB,out, ...
   VPAR, TT_par, idxCode, Slo_top,vi);
load(out);
end

tic ;
%profile on
%bau = waitbar(0,'Waiting...');

%==============================================================
%% Workers - Only for personal computer (Windows)
%==========================================================================
% 
% %if I am working on my personal computer define this - Think has 10 cores
% if ~contains(Directories.root,"nfs") 
%     numWorkers = 8;
%     %Create a parallel pool with the specified number of workers
%     poolobj = gcp('nocreate'); % Check if a pool already exists
%     if isempty(poolobj)
%         parpool(numWorkers); % Create a new pool if one doesn't exist
%     elseif poolobj.NumWorkers ~= numWorkers
%         %If a pool exists with a different size, delete it and create a new one
%         delete(poolobj);
%         parpool(numWorkers);
%     end
% end

% Restart: handled by Output_Manager ('init' phase, just before the time loop)

%% Display setting of the incoming T&C model runs:
disp(['Site selected: ' site_name])
disp(['Simulation period: ' datestr(x1,'dd-mmm-yyyy HH:MM') ' to ' datestr(x2,'dd-mmm-yyyy HH:MM')])
disp(['Precipitation phase scheme: ' parameterize_phase_label{:}])

if OPT_Idyn == 0; disp('Ice dynamics: off'); else; disp('Ice dynamics: on'); end
if OPT_Aval == 0; disp('Avalanching: off'); else; disp('Avalanching: on'); end


%% ========================================================
% Load meteo data
% ============================================================

%Consider sensor heights - you should add these onto the plant height in
%the zatm parameters

if OPT_Forcing == 1
Meteo_data = load(Directories.forc_meteo,meteo_name);
Meteo_data = Meteo_data.(meteo_name);
Forcing_Date = Meteo_data.DateTime;
Ta = Meteo_data.Ta; %Air temperature, degrees C
U = Meteo_data.RH; %Relative humidity, %
Ws = Meteo_data.WS; %Wind speed, ms-1
Pr = Meteo_data.Pr; %Precipitation, mm 
SWin = Meteo_data.SWin; %Incoming shortwave, Wm-2  
SWout = Meteo_data.SWout; %Outgoing shortwave, Wm-2
Nin = Meteo_data.LWin; %Incoming longwave, Wm-2
LWout = Meteo_data.LWout; %Outgoing longwave, Wm-2

%Lapse rates
hm_Ta_lapse = load(ta_lapse_file); %These are the Ta lapse rates from TOPKAPI, per hour and month for main station 144, catchment 406
hm_Ta_lapse = hm_Ta_lapse.(ta_lapse_name); %Pull out of structure. The lapse rate column should be named 'Ta_lapse'
Ldown_lapse = -0.031; %from Marty et al. (2002) p145
%Pr lapse - based on ratio per month
Pr_lapse = load(pr_lapse_file);
Pr_lapse = Pr_lapse.(pr_lapse_name); %Take out of structure, this is the fit of ratio = (p1 * elevation) + p2. Just multiply ratio by Pr. 
%Radiation lapse rates
SAB1_lapse = 0.0049; %Wm-2 per m (from Simone, for the Alps, may be a little high for Andes), radiation increases with elevation
SAB2_lapse = 0.0080; %Wm-2 per m (from Simone, for the Alps)
PARB_lapse = 0.0047; %Wm-2 per m (from Simone, for the Alps)

%Calculate radiation partition once to get t_bef and t_aft 
a=17.27; b=237.3;
Uf = U/100; %Change to fraction; (this is done again later for _P var)
esat=611.*exp(a.*Ta./(b+Ta)); %Vapour pressure at saturation (Pa)
ea=Uf.*esat;                 %Vapour pressure (Pa)
Ds= esat - ea;              %Vapor Pressure Deficit (Pa)
Ds(Ds<0)=0; 
xr=a.*Ta./(b+Ta)+log10(Uf); 
Tdew=b.*xr./(a-xr);               %Presumed dewpoint temperature (�C)
clear a b xr;

%Run radiation partition once to get t_bef and t_aft, and the radiation for
%distribution
GRAPH_VAR=0;
Forcing_Date_num =datenum(Forcing_Date);
[SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,Rsws,t_bef,t_aft]=Automatic_Radiation_Partition_I(Forcing_Date_num,Lat,Lon,Meteo_el,IniCond.DeltaGMT,Pr,Tdew,SWin,GRAPH_VAR);
clear N Rsws

% =======  Set for fully distributed
elseif OPT_Forcing == 2
t_bef = CONFIG_vals.t_bef;
t_aft = CONFIG_vals.t_aft;

%========== For On The Fly distributed forcing ======
elseif OPT_Forcing == 3
    % -------------------------------------------------------------------------
    % LOAD AND ADJUST ERA5-LAND INPUTS
    % -------------------------------------------------------------------------
    OPT_Run_BC = CONFIG_vals.OPT_RunBC;
    dsERA = load([meteo_name '_ERAdata.mat']);
    % Ensure alignment of ERA data to model run period
    [~,ia,ib] = intersect(Date,dsERA.dateTimes);
    FN = fieldnames(dsERA);
    for vv = 1:length(FN)
        if length(size(dsERA.(FN{vv})))==3
            dsERA.(FN{vv}) = dsERA.(FN{vv})(:,:,ib);
        end
    end
    dsERA.dateTimes = dsERA.dateTimes(ib);
    dsERA.TLR = dsERA.TLR(ib);

    % -------------------------------------------------------------------------
    % LOAD AND ADJUST DEM GRID DATA FOR DOWNSCALING USE
    % -------------------------------------------------------------------------    
    dsGRID = load([meteo_name '_GRIDdata.mat']);
    dsGRID = rmfield(dsGRID,{'dsLats','dsLons','dsZs','dsNames'});  % will remove later

    nYr = years(Date(end)-Date(1)); 
    [~,ia,ib] = intersect(dsERA.dateTimes,Date);
    dsERA.tpMean = nansum(dsERA.tp(:,:,ia),3)/nYr;
    
    dsERA.tpLr = 0.0003; % NEWER VERSIONS HAVE THIS ALREADY INSIDE dsERA.   SO COULD REMOVE IN FUTURE
    
    n = 1500;
    dsGRID.demCurv = getcurvature(dsGRID.demZs,n);
    dsGRID.pixelSize = abs(dsGRID.demLats(1,1)-dsGRID.demLats(2,1))*111.1e3; % (m)
    [dsGRID.demSlope,dsGRID.demAspect] = slopeaspect(dsGRID.demZs,dsGRID.pixelSize);

    load([meteo_name '_BC_Coefficients.mat']);

    % -------------------------------------------------------------------------
    % CONFIGURATION FIELDS FOR ON-THE-FLY GRID DOWNSCALING
    % -------------------------------------------------------------------------
    cfg.runMode = 'OTF';   %  'OTF' or 'save'
    cfg.dim = 'grid';   %  'grid', 'cluster', or 'points'
    cfg.DoubleOutput = 1;  % output both downscaled AND bias-corrected data (if relevant)
    cfg.RunRadiationPartition = 1; % If desire to run radiation partition. 1 = yes
    cfg.Date = Date;  % datetime vector at hourly resolution
    cfg.UTC = IniCond.DeltaGMT;
    cfg.VarNames = {'TA','RH','PP','FF','SWIN','LWIN','PRESS'}; % cell array, e.g. {'TA','RH','PP','SWIN','FF'}
    cfg.VarDS = {'dsT2m','dsRH','dsTp','dsWs10','dsSsrd','dsStrd','dsSp'};  % output field names corresponding to VarNames
    cfg.VarOut =  {'Ta','RH','Pr','Ws','Swin','N','Pre'};  % output field names corresponding to VarNames 
    cfg.FullNames = {'Air Temperature','Relative Humidity','Precipitation','Wind Speed','Shortwave Radiation','Longwave Radiation','Air Pressure'}; % display names corresponding to VarNames
    
    cfg.dsLats = [];         %: grid latitudes
    cfg.dsLons = [];       %: grid longitudes
    cfg.dsZs =  [];      %: grid elevations
    
    cfg.RunBC = OPT_Run_BC;              % Run Bias Correction? 1 = yes
    cfg.BC = COEFFS;             % structure containing monthly EQM information.
    
    cfg.bcLats = COEFFS.bcLats;         %: POI latitudes
    cfg.bcLons = COEFFS.bcLons;       %: POI longitudes
    cfg.bcZs = COEFFS.bcZs;          %: POI elevations

    cfg.foOut = [];   % output directory (required if cfg.runMode = 'save')
    
    cfg.dsNames = [];        % point names, required for dim = 'points' and save
    cfg.clusterNumber = [];  % cluster IDs, required for dim = 'cluster' and save
    
    cfg.precipLimit = 0.1; % precipLimit;    % precipitation threshold (mm hr); default 0.1   - Provide as parameter OPT?  
    cfg.windMinimum = 0.1; % windMinimum;   % wind-speed minimum (m s-1); default 0.3;                - Provide as parameter OPT?
    cfg.writeBiasCorr = true;  % optional logical; default true

% 
% 
%     % Create radiation partition on ERA grid for initial t_bef/t_aft
%     MidCell = [floor(size(dsERA.zs,1)/2),floor(size(dsERA.zs,2)/2)];
%     GRAPH_VAR=0;
%     [SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,Rsws,t_bef,t_aft]=Automatic_Radiation_Partition_I(datenum(Date),Lat,Lon,...
%         dsERA.zs(MidCell(1),MidCell(2)),IniCond.DeltaGMT,squeeze(dsERA.tp(MidCell(1),MidCell(2))),squeeze(dsERA.d2m(MidCell(1),MidCell(2),:)),...
%         squeeze(dsERA.ssrd(MidCell(1),MidCell(2),:)),GRAPH_VAR);
% clear N Rsws


%========== User forcing function (maps built at each time step) ======
% The function (on the path, e.g. in <site>/Forcing) is called with:
%   Forc_data = fun('load', Directories.forc)          once, loads the user data
%   F = fun('maps', Forc_data, Date_t, Forc_grid)       at each time step
% F must contain Ta [C], Pr [mm per time step], U [RH, fraction or %], N [W m-2],
% Ws [m s-1], Pre [hPa or Pa], SAB1, SAB2, SAD1, SAD2, PARB, PARD [W m-2, flat
% surface: the terrain correction is applied here], and optionally Tdew [C].
% Each is a scalar (same value in all cells) or a vector of num_cell values.
elseif OPT_Forcing == 4
    t_bef = CONFIG_vals.t_bef;
    t_aft = CONFIG_vals.t_aft;
    if exist(forc_function,'file') ~= 2
        error('Forcing function not found on the path: %s', forc_function)
    end
    forc_fun = str2func(forc_function);
    Forc_data = forc_fun('load', Directories.forc);
    % Grids given to the forcing function (cell vectors, same order as the model)
    Forc_grid = struct('DTM',DTMn, 'MASK',MASKn, 'GLH',GLHn, 'DEB_MAP',DEB_MAPn, ...
        'num_cell',num_cell, 'cellsize',cellsize, 'x',x_cell, 'y',y_cell);
    forc_vars_4 = {'Ta','Pr','U','N','Ws','Pre','SAB1','SAB2','SAD1','SAD2','PARB','PARD'};
    disp(['Forcing: user function ' forc_function ' (maps built at each time step)'])
end

% Solar variables (just for Lmax_day)
%--------------------------------------------------------------------------
L_day=zeros(length(Datam),1);
for j=2:24:length(Datam)
    [~,~,~,~,~,L_day(j)]= SetSunVariables(Datam(j,:),IniCond.DeltaGMT,Lon,Lat,t_bef,t_aft);
end

Lmax_day = max(L_day);
clear('L_day')

%=========================================================================
%% Iterating on time
%==========================================================================

fts = 2; %First time step
OM_phase = 'init'; Output_Manager; % Reads Output_Config.m; if Restart_ID = 1, loads the last checkpoint. Sets t_start

for t=t_start:N_time_step   
    
    %waitbar(t/N_time_step,bau)
    disp(['Iter: ' char(num2str(t))]);
    
    % CHECK WELL WHY THE MODEL USES t-1 instead of t
    Datam_S=Datam(t-1,:);
    
    % Year and month to load the forcing 
    %----------------------------------------------------------------------
    yy = char(num2str(Datam_S(1,1)));
    mth = char(num2str(Datam_S(1,2)));
    day = char(num2str(Datam_S(1,3)));
    hhlast = char(num2str(Datam_S(1,4)));
    Date_run = datetime(Datam_S(1,1), Datam_S(1,2), Datam_S(1,3));
    
    %% DISTRIBUTED FORCING
    %======================================================================
    % Forcing from stations and distributed across the catchment.
    % Create spatial input data on first day and hour of month 
    %======================================================================      

    %So create new data every month
    %For OPT_Forcing = 4 the forcing is built at every time step (one hour only)
    if OPT_Forcing == 4 || ~exist('Ta_Sm', 'var') || str2num(yy) ~= year_loaded || str2num(mth) ~= month_loaded
    if OPT_Forcing == 4
    date_forMonth = Date(t-1); % Only the hour of this time step
    else
    date_forMonth = Date_generator(str2double(yy),str2double(mth)); % Hours of this month, used to select the forcing

    disp(['New forcing loaded for period: ' char(num2str(yy)) '-' char(num2str(mth)) ])
    end

    %  ===== For weather station point =========
    if OPT_Forcing == 1
    % Pull out meteo data for the month for each variable  (before spatial distribution)
    id_date_for_month = ismember(Forcing_Date,date_forMonth); %So this is the index of the forcing data for this month
    
    %Temperature
    Ta_P = Ta(id_date_for_month); %For the hour
    %Sort later.....Ta_P_day = Ta(pdind); %For the 24 hours before
    %Longwave incoming
    N_P = Nin(id_date_for_month);
    %Relative humidity
    U_P = U(id_date_for_month)/100; %Change to fraction
    %Wind speed
    Ws_P = Ws(id_date_for_month);
    Ws_P(Ws_P < 0.01) = 0.01; %Should not be zero, but double check
    %Precipitation
    Pr_P = Pr(id_date_for_month);
    %Incoming shortwave radiation
    SWin_P = SWin(id_date_for_month);
    %Dew point temperature
    Tdew_P = Tdew(id_date_for_month);
    %Radiation
    SAD1_P = SAD1(id_date_for_month);
    SAD2_P = SAD2(id_date_for_month);
    SAB1_P = SAB1(id_date_for_month);
    SAB2_P = SAB2(id_date_for_month);
    PARB_P = PARB(id_date_for_month);
    PARD_P = PARD(id_date_for_month);
    %Air pressure will be calculated based on elevation
    end

    %=============================================
    %% Distributed meteorological variables
    % ===========================================
    %Then create grids for every cell so time in rows, space in columns
    %Preallocate where necessary

    %==== Load point data =====
    if OPT_Forcing == 1
        size_time_in_month = sum(id_date_for_month);
        Ta_Sm = NaN(size_time_in_month,num_cell);
        Pr_Sm = NaN(size_time_in_month,num_cell);
        N_Sm = NaN(size_time_in_month,num_cell);
        Tdew_Sm = NaN(size_time_in_month,num_cell);
        SAB1_Sm = NaN(size_time_in_month,num_cell);
        SAB2_Sm = NaN(size_time_in_month,num_cell);
        PARB_Sm = NaN(size_time_in_month,num_cell);
        Pr_ratio = NaN(size_time_in_month,num_cell);
    
        for tinm = 1:size_time_in_month % run over each hour
            m_temp = month(date_forMonth(tinm));
            h_temp = hour(date_forMonth(tinm));
            mh_id = hm_Ta_lapse.hour==h_temp & hm_Ta_lapse.month==m_temp;
            Ta_lapse = -hm_Ta_lapse.Ta_lapse(mh_id); % the lapse rate per hour and month - note converting to negative rate
            Tdew_lapse = Ta_lapse; %Its ok to use the same lapse rate
            %Air temperature, degrees C
            Ta_Sm(tinm,:) = Ta_P(tinm) + (Ta_lapse.*(DTMn-Meteo_el)); %Lapse rates should be negative
            %Precipitation, mm - This uses the WRF ratio v el relationship
            Pr_ratio(tinm,:) = (Pr_lapse(m_temp).month.p1.*DTMn) + Pr_lapse(m_temp).month.p2;  %val(x) = p1*x + p2
            Pr_Sm(tinm,:) = Pr_ratio(tinm,:).*Pr_P(tinm); %Apply ratio
            %Longwave radiation, W m-2
            N_Sm(tinm,:) = N_P(tinm) + (Ldown_lapse.*(DTMn-Meteo_el));
            %Dew point temperature, degrees C
            Tdew_Sm(tinm,:) = Tdew_P(tinm) + (Tdew_lapse.*(DTMn-Meteo_el)); %Lapse rates should be negative
            %Direct radiation with elevation
            if SAB1_P(tinm)>0 %Only if radiation values
                SAB1_Sm(tinm,:) = SAB1_P(tinm) + (SAB1_lapse.*(DTMn-Meteo_el)); %Radiation should increase with elevation
            else
                SAB1_Sm(tinm,:) = 0;
            end
            if SAB2_P(tinm)>0
                SAB2_Sm(tinm,:) = SAB2_P(tinm) + (SAB2_lapse.*(DTMn-Meteo_el)); %Radiation should increase with elevation
            else
                SAB2_Sm(tinm,:) = 0;
            end
            if PARB_P(tinm)>0
            PARB_Sm(tinm,:) = PARB_P(tinm) + (PARB_lapse.*(DTMn-Meteo_el)); %Radiation should increase with elevation
            else 
                PARB_Sm(tinm,:) =0;
            end
            %Prevent negative values
            SAB1_Sm(SAB1_Sm<0)=0; SAB2_Sm(SAB2_Sm<0)=0; PARB_Sm(PARB_Sm<0)=0; 
        end %Going over time steps within month
    
        %Air pressure (calculated directly from elevation)
        Pre_temp = (101325*((1-(0.0065.*DTMn./288.15)).^(9.18*0.0289644./(8.31447*0.0065))))./100; %/100 converts to mbar
        Pre_Sm = repmat(Pre_temp',size_time_in_month,1); %Copy to all time steps
    
        %Use the lapsed dew point temprature to derive relative humidity
        c=237.3; b=17.27;
        U_Sm = 100*exp((c*b.*(Tdew_Sm - Ta_Sm))./((c+Ta_Sm).*(c+Tdew_Sm)));
        clear c b 
        if max(max(U_Sm))>=1
        U_Sm = U_Sm./100; %Turn into fraction
        U_Sm(U_Sm>1) = 1; %Cannot be more than 100%
        end
    
        %Here we use the same wind speed everywhere
        %At the moment turn into S vector grids
        Ws_Sm = repmat(Ws_P,1,num_cell); %So the same in all grid cells
        %Diffuse radiation can be the same everywhere
        SAD1_Sm = repmat(SAD1_P',1,num_cell); %So the same in all grid cells
        SAD2_Sm = repmat(SAD2_P',1,num_cell); %So the same in all grid cells
        PARD_Sm = repmat(PARD_P',1,num_cell); %So the same in all grid cell
    end

    %========= Load gridded data ========
    if OPT_Forcing == 2
        %Sort time
        StartMonth = datetime(str2double(yy),str2double(mth),1,0,0,0,0);
        EndMonth = dateshift(StartMonth,'end','month'); EndMonth = EndMonth + hours(23);
        Forcing_Date = StartMonth:hours(1):EndMonth; Forcing_Date = Forcing_Date';
        id_date_for_month = ismember(Forcing_Date,date_forMonth); %So this is the index of the forcing data for this month. date_forMonth created in Initialising outputs
        size_time_in_month = sum(id_date_for_month);
        %Load variables
        forc_vars = ["t2m" "d2m" "tp" "RH" "sp" "ws10" "strd" "N" "PARB" "PARD" "SAB1" "SAB2" "SAD1" "SAD2"]; 
        forc_names = split(CONFIG_vals.forc_meteo_name); 
        OD_vars = contains(forc_vars,["RH","strd","sp","ws10"]);
        num_forc_vars = size(forc_vars,2);
        for v = 1:num_forc_vars
            if OD_vars(v)==0
            fn = char(strcat(Directories.forc_meteo,yy,'/',forc_vars(v), '_',yy,'_',mth,'.mat'));
            Forc_S.forc_vars{v} = load(fn);   
            else
            fn = char(strcat(Directories.forc_meteoOD,yy,'/',forc_vars(v), '_',yy,'_',mth,'.mat'));
            Forc_S.forc_vars{v} = load(fn);
            end
        end

        %Then decompress, flip and reshape
        for v=1:num_forc_vars
            %** Warning the decompress code changes the magnitude of the
            %variables so they should have been compressed before. **
            Forc_S_dc.(forc_vars{v}) = Data_compressor_grid(Forc_S.forc_vars{v}.(forc_names{v}),forc_vars(v),"decompress");
            Forc_S_dc.(forc_vars{v}) =  flipud(Forc_S_dc.(forc_vars{v}));
            for ht = 1:size(Forcing_Date,1) %Probably this could be faster 
            Forc_S_dc_r.(forc_vars{v})(ht,:) = reshape(Forc_S_dc.(forc_vars{v})(:,:,ht),[1,num_cell]); %Reshape to be time (rows) by grid cells (columns)
            end
        end
        %Bring into correct names and apply small corrections
        %Also crops in case the run time is shorter
        Ta_Sm = double(Forc_S_dc_r.t2m(id_date_for_month,:)); %Air temperature, degrees C
        Pr_Sm = double(Forc_S_dc_r.tp(id_date_for_month,:)); %Precipitation, mm
        U_Sm = double(Forc_S_dc_r.RH(id_date_for_month,:)); %Relative humidity, %
        N_Sm = double(Forc_S_dc_r.strd(id_date_for_month,:)); %Incoming Longwave (note N in the ERA data is cloudiness), Wm^-2
        Tdew_Sm = double(Forc_S_dc_r.d2m(id_date_for_month,:)); %Dew point temperature, degrees C
        Pre_Sm = double(Forc_S_dc_r.sp(id_date_for_month,:)); %Surface pressure
        Ws_Sm = double(Forc_S_dc_r.ws10(id_date_for_month,:)); %Wind speed, ms-1
        SAB1_Sm = double(Forc_S_dc_r.SAB1(id_date_for_month,:)); %First band direct radiation, Wm^-2
        SAB2_Sm = double(Forc_S_dc_r.SAB2(id_date_for_month,:)); %Second band direct radiation, Wm^-2
        SAD1_Sm = double(Forc_S_dc_r.SAD1(id_date_for_month,:)); %First band diffuse radiation, Wm^-2
        SAD2_Sm = double(Forc_S_dc_r.SAD2(id_date_for_month,:)); %Second band diffuse radiation, Wm^-2
        PARB_Sm = double(Forc_S_dc_r.PARB(id_date_for_month,:)); %PAR radiation direct, Wm^-2
        PARD_Sm = double(Forc_S_dc_r.PARD(id_date_for_month,:)); %PAR radiation diffuse Wm^-2

        %Corrections
        if max(max(U_Sm))>=1
            U_Sm = U_Sm./100; %Turn into fraction
            U_Sm(U_Sm>1) = 1; %Cannot be more than 100%
        end
        Ws_Sm(Ws_Sm < 0.01) = 0.01; %To prevent 0 turbulent fluxes
        if max(max(Pre_Sm))>10000
        Pre_Sm = Pre_Sm./100; %Convert from Pascals to mbar/hectopascal
        end
    end

     %========= On-The_Fly Downscaling ========
    if OPT_Forcing == 3

    StartMonth = datetime(str2double(yy),str2double(mth),1,0,0,0,0); % Derive the monthly timesteps
        EndMonth = dateshift(StartMonth,'end','month'); EndMonth = EndMonth + hours(23);
        Forcing_Date = StartMonth:hours(1):EndMonth; Forcing_Date = Forcing_Date';

    dsERAsub = dsERA; % Make a copy of the ERA5 data
    [~,ia,ib] = intersect(Forcing_Date,dsERA.dateTimes); % Find the overlapping timesteps
     FN = fieldnames(dsERAsub); % get field names from the ERA5 data structure

     % %%%%% A SUB-TIMING TEST DUE TO A SLOW RADIATION PARTITION
     % ib = 1:48; 
     % cfg.Date = Date(1:48);
     % cfg.RunRadiationPartition = 0;
     % %%%%

    for vv = 1:length(FN) % Adjust the timesteps on the 3D matrices
        if length(size(dsERAsub.(FN{vv})))==3
            dsERAsub.(FN{vv}) = dsERAsub.(FN{vv})(:,:,ib);
        end
    end
    % adjust the timing of the date itself and the lapse rates
    dsERAsub.dateTimes = dsERAsub.dateTimes(ib);
    dsERAsub.TLR = dsERAsub.TLR(ib);
    
    % MAIN DOWNSCALING CODE
    %--------------------
    OTF = downscale_4_T_C(cfg,dsERAsub,dsGRID);
    %--------------------

    try % Open structure containing bias-corrected fields. If no bias correction is found /or applied, the non-bias-corrected files are opened from the structure
        FN = fieldnames(OTF.biasCorrected);
        for ff = 1:length(FN)
            [m,n,q] = size(OTF.biasCorrected.(FN{ff}));
            OTF.biasCorrected.(FN{ff}) = transpose(reshape(flipud(OTF.biasCorrected.(FN{ff})),[m*n,q]));
        end
        v2struct(OTF.biasCorrected) % this v2struct should stay ideally
    catch
        FN = fieldnames(OTF.downscaled);
        for ff = 1:length(FN)
            [m,n,q] = size(OTF.downscaled.(FN{ff}));
            OTF.downscaled.(FN{ff}) = transpose(reshape(flipud(OTF.downscaled.(FN{ff})),[m*n,q]));
        end
        v2struct(OTF.downscaled) % this v2struct should stay ideally
    end

    U_Sm = RH / 100;
    % ADJUST OTHER NAMES HERE (FOR THE CURRENT WORKING CONVENTION) - (30/10/26)
    Ta_Sm = Ta;
    Pr_Sm = Pr;
    Ws_Sm = Ws;
    N_Sm = N;
    Pre_Sm = Pre;
    SAB1_Sm = SAB1;
    SAB2_Sm = SAB2;
    SAD1_Sm = SAD1;
    SAD2_Sm = SAD2;
    PARB_Sm = PARB;
    PARD_Sm = PARD;
   
    Tdew_Sm = (243.12 .* (log(U_Sm) + (17.62 .* Ta_Sm)./(243.12 + Ta_Sm))) ./ ...
     (17.62 - (log(U_Sm) + (17.62 .* Ta_Sm)./(243.12 + Ta_Sm)));

    end

    %========= User forcing function (one time step) ========
    if OPT_Forcing == 4
        size_time_in_month = 1;
        F = forc_fun('maps', Forc_data, date_forMonth, Forc_grid);
        if ~isfield(F,'Tdew'); forc_list = forc_vars_4; else; forc_list = [forc_vars_4 {'Tdew'}]; end
        for v = 1:numel(forc_list)
            if ~isfield(F,forc_list{v})
                error('Forcing function %s: variable %s missing', forc_function, forc_list{v})
            end
            val = double(F.(forc_list{v})(:))'; % row vector (time x cells)
            if isscalar(val)
                val = repmat(val,1,num_cell); % same value in all cells
            elseif numel(val) ~= num_cell
                error('Forcing function %s: %s has %d values, expected 1 or %d', forc_function, forc_list{v}, numel(val), num_cell)
            end
            F.(forc_list{v}) = val;
        end
        Ta_Sm = F.Ta; Pr_Sm = F.Pr; U_Sm = F.U; N_Sm = F.N; Ws_Sm = F.Ws; Pre_Sm = F.Pre;
        SAB1_Sm = F.SAB1; SAB2_Sm = F.SAB2; SAD1_Sm = F.SAD1; SAD2_Sm = F.SAD2; PARB_Sm = F.PARB; PARD_Sm = F.PARD;

        %Corrections (as for OPT_Forcing = 2)
        if max(U_Sm(MASKn==1))>1
            U_Sm = U_Sm./100; %Turn into fraction
        end
        U_Sm(U_Sm>1) = 1; U_Sm(U_Sm<0) = 0; %Between 0 and 100%
        Ws_Sm(Ws_Sm < 0.01) = 0.01; %To prevent 0 turbulent fluxes
        Pr_Sm(Pr_Sm < 0) = 0;
        if max(Pre_Sm(MASKn==1))>10000
            Pre_Sm = Pre_Sm./100; %Convert from Pascals to mbar/hectopascal
        end

        %Dew point temperature: from the function, or from Ta and U (same formula as OPT_Forcing = 1)
        if isfield(F,'Tdew')
            Tdew_Sm = F.Tdew;
        else
            a=17.27; b=237.3;
            xr = a.*Ta_Sm./(b+Ta_Sm) + log10(U_Sm);
            Tdew_Sm = b.*xr./(a-xr);
            clear a b xr
        end
        clear F val v forc_list
    end

    %========= Further forcing steps for all types ==========
    if OPT_Decouple == 1

    % Glacier Cooling Effect - Shaw et al. (2025)
    load('ModelCoefficients_Decoupling.mat'); % Should be in 'Functions' path that is already generated
    
    SpecHum = convert_humidity(Pre_Sm*100,Ta_Sm+273.15,U_Sm*100,'relative humidity','specific humidity');
    
    D = DTM(:)'; F = FPL(:)'; % NEED TO CHANGE ORIENTATION OF GRIDS?

    for tt = 1:size(Ta_Sm,1)
        
        Swin = SAD1_Sm + SAB1_Sm + SAD2_Sm + SAB2_Sm;
        if min(Swin(tt,:),[],'all') > 30 % day time only (without very overcast conditions)
        
        T = Ta_Sm(tt,:); S = SpecHum(tt,:); FF = Ws_Sm(tt,:); Td = Tdew_Sm(tt,:);
        K = predict(DECOUPLE_COEFF,[T',S',FF',F',D'])'; % Add into prediction as vector (nx1)
        K(K<0.2) = 0.2; K(K>1) = 1; %limits of the observations

        ValidCooling = ~isnan(F) & ~isnan(D) & T > 0;

        Td(ValidCooling) = Td(ValidCooling) .* K(ValidCooling);      
        T(ValidCooling) = T(ValidCooling) .* K(ValidCooling);
                  
        Ta_Sm(tt,:) = T;
        Tdew_Sm(tt,:) = Td;

        end
    end

    % Reconvert to relative humidity (and fraction)
    U_Sm = exp((237.3*17.227.*(Tdew_Sm - Ta_Sm))./((237.3+Ta_Sm).*(237.3+Tdew_Sm))); 
    U_Sm(U_Sm>1) = 1; %Cannot be more than 100%

    end


    % Vapor pressure - calculate based on distributed Ta and U
    % esat/ea/Ds/Tdew grids of month of time x cell
    a=17.27; b=237.3;
    esat_Sm=611.*exp(a.*Ta_Sm./(b+Ta_Sm)); %Vapour pressure at saturation (Pa)
    ea_Sm=U_Sm.*esat_Sm;                 %Vapour pressure (Pa)
    Ds_Sm= esat_Sm - ea_Sm;              %Vapor Pressure Deficit (Pa)
    Ds_Sm(Ds_Sm<0)=0; 
    %Don't recalculate Tdew - its above
    %NOTE SHOULD WE ADJUST DEW POINT WITH LOWER AIR TEMP FROM TMOD?
    clear a b xr;

    %---- Check variables -----
    %We only check inside the mask in case there are NaNs outside.
    check_var_3D(Ta_Sm(:,MASKn==1),Tdew_Sm(:,MASKn==1),Pr_Sm(:,MASKn==1),N_Sm(:,MASKn==1),Ws_Sm(:,MASKn==1),Pre_Sm(:,MASKn==1),esat_Sm(:,MASKn==1),ea_Sm(:,MASKn==1),U_Sm(:,MASKn==1),SAD1_Sm(:,MASKn==1),SAD2_Sm(:,MASKn==1),SAB1_Sm(:,MASKn==1),SAB2_Sm(:,MASKn==1),PARB_Sm(:,MASKn==1),PARD_Sm(:,MASKn==1));
    
%   % Store year and month loaded
%   %----------------------------------------------------------------------
    year_loaded = str2num(yy);
    month_loaded = str2num(mth);
% 
    end %Loading month of data

    %% All now on single timestep
    % Finding the row in the forcing of the modeling date      
    t_forc = find(Date(t-1) == date_forMonth); %Pulls out the index for the hour - note this means that it actually starts on t=1

    %Put here downscaling
    %And the partition of SW
    
    %Export timestep (a row for all the cells)
    Ta_S = Ta_Sm(t_forc,:);
    Pr_S = Pr_Sm(t_forc,:);
    N_S = N_Sm(t_forc,:);
    Pre_S = Pre_Sm(t_forc,:);
    Ws_S = Ws_Sm(t_forc,:);
    U_S = U_Sm(t_forc,:);
    ea_S = ea_Sm(t_forc,:);
    Ds_S = Ds_Sm(t_forc,:);
    Tdew_S = Tdew_Sm(t_forc,:);
    SAD1_S = SAD1_Sm(t_forc,:);
    SAD2_S = SAD2_Sm(t_forc,:);
    SAB1_S = SAB1_Sm(t_forc,:);
    SAB2_S = SAB2_Sm(t_forc,:);
    PARB_S = PARB_Sm(t_forc,:);
    PARD_S = PARD_Sm(t_forc,:);

    %% Sort Ta for day before
    if t == fts
        Ta_Spdind = Ta_S;
    elseif t > fts && t < fts+24
        Ta_Spdind = [Ta_Spdind;Ta_S];
    elseif t >= fts+24
        Ta_Spdind = [Ta_Spdind(2:24,:);Ta_S];
    end

    %% 

    % Reshape - This because of the problem in the Hydrological module
    %----------------------------------------------------------------------
    Slo_top2 = reshape(Slo_top,num_cell,1); % Creating aux variable for the hydrological module
    aTop = reshape(aTop,num_cell,1);

    %% Radiation Part B 
    %Run per timestep now
    [jDay]= JULIAN_DAY(Datam_S);
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day]= SetSunVariables(Datam_S,IniCond.DeltaGMT,Lon,Lat,t_bef,t_aft); %These are one value per time step
    [ShF] = Shadow_Effect(DTM,h_S,zeta_S,HZ,Zasp); %ShF is x,y grid

    %needed, if terrain effects have not been considered during pre-processing
    cos_fst = cos(atan(Slo_top))*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);
    cos_fst(cos_fst<0)=0; %cos_fst is x,y grid   
    
    %Reshape radiation grids for calculations
    %Note shifting to rows to match radiation inputs
    SvFn = reshape(SvF,num_cell,1)';
    Ctn = reshape(Ct,num_cell,1)';
    ShFn = reshape(ShF,num_cell,1)';
    cos_fstn = reshape(cos_fst,num_cell,1)';

    %rho_g and h_S are scalars
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SAD1_S = SAD1_S.*SvFn + Ctn.*rho_g.*(SAD1_S./sin(h_S).*cos_fstn + (1-SvFn).*SAD1_S);
    SAD2_S = SAD2_S.*SvFn + Ctn.*rho_g.*(SAD2_S./sin(h_S).*cos_fstn + (1-SvFn).*SAD2_S);
    PARD_S = PARD_S.*SvFn + Ctn.*rho_g.*(PARD_S./sin(h_S).*cos_fstn + (1-SvFn).*PARD_S);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SAB1_S =(SAB1_S./sin(h_S)).*cos_fstn.*ShFn;
    SAB2_S =(SAB2_S./sin(h_S)).*cos_fstn.*ShFn;
    PARB_S = (PARB_S./sin(h_S)).*cos_fstn.*ShFn;        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SAB1_S(SAB1_S<0)=0;
    SAB2_S(SAB2_S<0)=0;
    PARB_S(PARB_S<0)=0;
    PARD_S(PARD_S<0)=0;
    SAD1_S(SAD1_S<0)=0;
    SAD2_S(SAD2_S<0)=0;
    SAB1_S(isnan(SAB1_S)) = 0;
    SAB2_S(isnan(SAB2_S)) = 0;
    SAD1_S(isnan(SAD1_S)) = 0;
    SAD2_S(isnan(SAD2_S)) = 0;
    PARB_S(isnan(PARB_S)) = 0;
    PARD_S(isnan(PARD_S)) = 0;

    % Sun below ~5.7 deg: the slope correction (division by sin(h_S)) is unstable;
    % set all shortwave to 0, as in the former launcher
    if sin(h_S) <= 0.10
        SAB1_S(:) = 0; SAB2_S(:) = 0; SAD1_S(:) = 0; SAD2_S(:) = 0; PARB_S(:) = 0; PARD_S(:) = 0;
    end
       
    %% Carbon
    Ca_S = Ca(t)*MASKn; %Carbon is already set on the correct time series
    IrD_S =  MASKn*0; %Setting as 0
    Salt_S =  MASKn*0; %Setting as 0

    %% Swapping other grids
    Afirnn = reshape(Afirn,num_cell,1);
    SOIL_THn = reshape(SOIL_TH,num_cell,1);
    
    %% SPATIAL INITIALIZATION VECTOR PREDEFINING
    %======================================================================
   
    if t == 2
    Qi_in_Ro = Qi_out;  Ccrown_t = Ccrown_t_tm1;
    v = [split("alp_soil b_soil Bam Bem BLit Cice Cicew CK1 Csno Csnow dQ_S DQ_S dQVEG DT_S dw_SNO e_sno " + ...
        "EG EICE EIn_rock EIn_urb EK ELitter er ESN_In SSN_In ESN SSN EWAT FROCK f Gfin G H HV ICE_D ICE " + ...
        "Imelt In_Litter In_rock In_SWE In_urb IP_wc Lk_rock Lk_wat Lk Lpho NavlI NDVI NIce NIn_SWE OF " + ...
        "Oice OS O POT Pr_liq Pr_sno QE QEV Qfm Qi_out Qv r_litter r_soil ra Rd Rh Rn ros SE_rock SE_urb " + ...
        "Smelt SND snow_albedo soil_albedo SP_wc surface_albedo SWE SWE_avalanched t_sls tau_sno Tdamp Tdeb " + ...
        "Tdp Tdp_snow Tice Ts Ts_under TsVEG U_SWE Vice V WAT WIS WR_IP WR_SP Ws_under ZWT"); ...
         split("AgeDL_H AgeDL_L AgeL_H AgeL_L AgePl_H AgePl_L An_H An_L ANPP_H ANPP_L B_H B_L BA_H BA_L " + ...
         "Bfac_dayH Bfac_dayL Bfac_weekH Bfac_weekL dflo_H dflo_L Dr_H Dr_L e_rel_H e_rel_L e_relN_H e_relN_L " + ...
         "EIn_H EIn_L fapar_H fapar_L FNC_H FNC_L gsr_H gsr_L hc_H hc_L In_H In_L ISOIL_H ISOIL_L Jsx_H Jsx_L " + ...
         "Jxl_H Jxl_L Kleaf_H Kleaf_L Kreserve_H Kreserve_L Kuptake_H Kuptake_L Kx_H Kx_L LAI_H LAI_L " + ...
         "LAIdead_H LAIdead_L ManIH ManIL NBLeaf_H NBLeaf_L NBLI_H NBLI_L NPP_H NPP_L NPPI_H NPPI_L " + ...
         "Nreserve_H Nreserve_L NuLit_H NuLit_L NupI_H NupI_L Nuptake_H Nuptake_L OH OL PARI_H PARI_L " + ...
         "PHE_S_H PHE_S_L Preserve_H Preserve_L Psi_l_H Psi_l_L Psi_s_H Psi_s_L Psi_x_H Psi_x_L Puptake_H " + ...
         "Puptake_L RA_H RA_L rap_H rap_L RB_H rb_H RB_L rb_L Rdark_H Rdark_L Rexmy_H Rexmy_L Rg_H Rg_L " + ...
         "rKc_H rKc_L Rmc_H Rmc_L Rmr_H Rmr_L Rms_H Rms_L rNc_H rNc_L rPc_H rPc_L Rrootl_H Rrootl_L " + ...
         "rs_shdH rs_shdL rs_sunH rs_sunL SAI_H SAI_L Sfr_H Sfr_L SIF_H SIF_L Slf_H Slf_L Sll_H Sll_" + ...
         "L Sr_H Sr_L SupK_H SupK_L SupN_H SupN_L SupP_H SupP_L Swm_H Swm_L T_H T_L TBio_H TBio_L " + ...
         "Tden_H Tden_L Tdp_H Tdp_L TdpI_H TdpI_L TexC_H TexC_L TexK_H TexK_L TexN_H TexN_L TexP_H " + ...
         "TexP_L TNIT_H TNIT_L TPHO_H TPHO_L TPOT_H TPOT_L Vl_H Vl_L Vx_H Vx_L")];
        for x = v', eval(sprintf('%s = %stm1;', x, x)); end
        for x = ["shdH","shdL","sunH","sunL"], eval(sprintf('Ci_%s = Citm1_%s;', x, x)); end
    end
  

    %% LOOP
    %======================================================================
    % LOOP OVER CELLS
    %======================================================================
    %P = (Xout - 1) * 133 + Yout; AJ: What's this? To remove?

    %Follow=MASK; %Debugging
   parfor ij= 1:num_cell % this is a parfor
        % Good practice to use a simple for loop for debugging/testing
        % ij is the index to go pixel by pixel through the mask
        % ij=1:num_cell
        %disp(strcat('in the loop', ij))
   
        %% Debugging for a for loop
%         if ij == 3211 
%             disp('wait') %261 first 6 cell
%         end
        % =================================================================
        
        if MASKn(ij)== 1
	    try % report which cell fails (see catch at the end of the cell block)

            %disp(['Cell: ' char(num2str(ij)) ', Veg Type: ' char(num2str(ksv(ij)))])
            Elev=DTMn(ij);
            %%
            % tic % timing of PARAMETERS_SOIL_DEV (prints one line per cell and step)
            % BOUNDARY CONDITION
            % INTRODUCED SOIL AND VEG. for ij
            %--------------------------------------------------------------
            [aR,             EvL_Zs,       Inf_Zs,     Bio_Zs,      Zinf, ...
             RfH_Zs,         RfL_Zs,         dz,           Ks_Zs,      Dz,          ms, ...
             Kbot,           Krock,          zatm,         Color_Class,  OM_H,       OM_L,        PFT_opt_H, ...
             PFT_opt_L,      d_leaf_H,       d_leaf_L,     SPAR,       Phy,         Soil_Param, ...
             Interc_Param,   SnowIce_Param,  VegH_Param,   VegL_Param, fpr,         VegH_Param_Dyn, ...
             VegL_Param_Dyn, Stoich_H,       aSE_H,        Stoich_L,   aSE_L,       fab_H,...
             fbe_H,          fab_L,          fbe_L,        ZR95_H,     ZR95_L,      ZRmax_H,      ZRmax_L,...
             In_max_urb,     In_max_rock,    K_usle,         Urb_Par,      Deb_Par,    Zs_deb,      Sllit,...
             Kct,            ExEM,           ParEx_H,      Mpar_H,     ParEx_L,     Mpar_L] = ....
                                    PARAMETERS_SOIL_DEV( ...
             ksv(ij),        PSANr(ij),       PCLAr(ij),     PORGr(ij),   DEB_MAPn(ij),  md_max,...
             Afirnn(ij),      SOIL_THn(ij),    VPAR,        TT_vals,      OPT_PARAM_vals, Zs, ms_max);
            % toc
%% 
            %Extract C values directly from grids
            Ccrown = VPAR.Ccrownr(ij);
            if iscell(Ccrown)==1
                Ccrown = Ccrown{:};
            end
            Cbare = VPAR.Cbarer(ij); 
            Crock = VPAR.Crockr(ij); 
            Curb = VPAR.Curbr(ij); 
            Cwat = VPAR.Cwatr(ij); 
            
            % If hour (Datam_S(4)) is equal to 1
            if (Datam_S(4)==1)                                
                %% SOIL BIOGEOCHEMISTRY MODULE
                [Se_bio,Se_fc,Psi_bio,Tdp_bio,VSUM,VTSUM]=Biogeo_environment([squeeze(Tdp_t(ij,:,:))]',[squeeze(O_t(ij,:,:))]',[squeeze(V_t(ij,:,:))]',...
                    Soil_Param,Phy,SPAR,Bio_Zs);%
                
                % Biogeochemistry Unit
                Nuptake_H(ij,:)= 0.0;
                Puptake_H(ij,:)= 0.0;
                Kuptake_H(ij,:)= 0.0; %% [gK/m^2 day]
                %%%
                Nuptake_L(ij,:)= 0.0; %% [gN/m^2 day]
                Puptake_L(ij,:)= 0.0;
                Kuptake_L(ij,:)= 0.0;
                %%%
                NavlI(ij,:)=[1 1 1];
                Bam(ij)=0; Bem(ij)=0;
                %%%              
              
               
        %% VEGETATION MODULE
        %==================================================================
        % FUNCTION: VEGETATION_MODULE_PAR
        %==================================================================

       [LAI_H(ij,:),          B_H(ij,:,:),        NPP_H(ij,:),          ANPP_H(ij,:),        Rg_H(ij,:), ...
        RA_H(ij,:),           Rms_H(ij,:),        Rmr_H(ij,:),          Rmc_H(ij,:),         PHE_S_H(ij,:),...
        dflo_H(ij,:),         AgeL_H(ij,:),       e_rel_H(ij,:),        e_relN_H(ij,:),      LAI_L(ij,:),...
        B_L(ij,:,:),          NPP_L(ij,:),        ANPP_L(ij,:),         Rg_L(ij,:),          RA_L(ij,:),...
        Rms_L(ij,:),          Rmr_L(ij,:),        Rmc_L(ij,:),          PHE_S_L(ij,:),       dflo_L(ij,:), ...
        AgeL_L(ij,:),         e_rel_L(ij,:),      e_relN_L(ij,:),       SAI_H(ij,:),         hc_H(ij,:), ...
        SAI_L(ij,:),          hc_L(ij,:),         LAIdead_H(ij,:),      NBLeaf_H(ij,:),      Sr_H(ij,:), ...
        Slf_H(ij,:),          Sfr_H(ij,:),        Sll_H(ij,:),          Swm_H(ij,:),         Rexmy_H(ij,:,:), ...
        NupI_H(ij,:,:),       NuLit_H(ij,:,:),    LAIdead_L(ij,:),      NBLeaf_L(ij,:),      Sr_L(ij,:), ...
        Slf_L(ij,:),          Sfr_L(ij,:),        Sll_L(ij,:),          Swm_L(ij,:),         Rexmy_L(ij,:,:),... 
        NupI_L(ij,:,:),       NuLit_L(ij,:,:),    Rrootl_H(ij,:),       AgeDL_H(ij,:),       Bfac_dayH(ij,:), ...
        Bfac_weekH(ij,:),     NPPI_H(ij,:),       TdpI_H(ij,:),         PARI_H(ij,:,:),      NBLI_H(ij,:), ...
        RB_H(ij,:,:),         FNC_H(ij,:),        Nreserve_H(ij,:),     Preserve_H(ij,:),    Kreserve_H(ij,:), ...
        rNc_H(ij,:),          rPc_H(ij,:),        rKc_H(ij,:),          ManIH(ij,:),         Rrootl_L(ij,:), ...
        AgeDL_L(ij,:),        Bfac_dayL(ij,:),    Bfac_weekL(ij,:),     NPPI_L(ij,:),        TdpI_L(ij,:), ...
        PARI_L(ij,:,:),       NBLI_L(ij,:),       RB_L(ij,:,:),         FNC_L(ij,:),         Nreserve_L(ij,:), ...
        Preserve_L(ij,:),     Kreserve_L(ij,:),   rNc_L(ij,:),          rPc_L(ij,:),         rKc_L(ij,:), ...
        ManIL(ij,:),          TexC_H(ij,:),       TexN_H(ij,:),         TexP_H(ij,:),        TexK_H(ij,:), ...
        TNIT_H(ij,:),         TPHO_H(ij,:),       TPOT_H(ij,:),         SupN_H(ij,:),        SupP_H(ij,:), ...
        SupK_H(ij,:),         ISOIL_H(ij,:,:),    TexC_L(ij,:),         TexN_L(ij,:),        TexP_L(ij,:), ...
        TexK_L(ij,:),         TNIT_L(ij,:),       TPHO_L(ij,:),         TPOT_L(ij,:),        SupN_L(ij,:), ...
        SupP_L(ij,:),         SupK_L(ij,:),       ISOIL_L(ij,:,:),      BA_H(ij,:),          Tden_H(ij,:), ...
        AgePl_H(ij,:),        BA_L(ij,:),         Tden_L(ij,:),         AgePl_L(ij,:),       Ccrown_t(ij,:)]= ...
                          VEGETATION_MODULE_PAR( ...
        cc_max,               Ccrown,              ZR95_H,              ZR95_L,              B_Htm1(ij,:,:),...
        PHE_S_Htm1(ij,:),     dflo_Htm1(ij,:),     AgeL_Htm1(ij,:),     AgeDL_Htm1(ij,:),    Ta_t(ij,:), ...
        PAR_t(ij,:),          Tdp_H_t(ij,:,:),     Psi_x_H_t(ij,:,:),   Psi_l_H_t(ij,:,:),   An_H_t(ij,:,:), ...
        Rdark_H_t(ij,:,:),    NPP_Htm1(ij,:),      jDay,                Datam_S,             NPPI_Htm1(ij,:), ...
        TdpI_Htm1(ij,:),      Bfac_weekHtm1(ij,:), Stoich_H,            aSE_H,               VegH_Param_Dyn,...
        Nreserve_Htm1(ij,:),  Preserve_Htm1(ij,:), Kreserve_Htm1(ij,:), Nuptake_H(ij,:),     Puptake_H(ij,:), ...
        Kuptake_H(ij,:),      FNC_Htm1(ij,:),      Tden_Htm1(ij,:),     AgePl_Htm1(ij,:),    fab_H, ...
        fbe_H,                ParEx_H,             Mpar_H,              TBio_H(ij,:),        SAI_Htm1(ij,:), ...
        hc_Htm1(ij,:),        B_Ltm1(ij,:,:),      PHE_S_Ltm1(ij,:),    dflo_Ltm1(ij,:),     AgeL_Ltm1(ij,:), ...
        AgeDL_Ltm1(ij,:),     Tdp_L_t(ij,:,:),     Psi_x_L_t(ij,:,:),   Psi_l_L_t(ij,:,:),   An_L_t(ij,:,:), ...
        Rdark_L_t(ij,:,:),    NPP_Ltm1(ij,:),      NPPI_Ltm1(ij,:),     TdpI_Ltm1(ij,:),     Bfac_weekLtm1(ij,:),...
        NupI_Htm1(ij,:,:),    NupI_Ltm1(ij,:,:),   NuLit_Htm1(ij,:,:),  NuLit_Ltm1(ij,:,:),  NBLeaf_Htm1(ij,:), ...
        NBLeaf_Ltm1(ij,:),    PARI_Htm1(ij,:,:),   NBLI_Htm1(ij,:),     PARI_Ltm1(ij,:,:),   NBLI_Ltm1(ij,:),...
        Stoich_L,             aSE_L,               VegL_Param_Dyn,      NavlI(ij,:),         Bam(ij), ...
        Bem(ij),              Ccrown_t_tm1(ij,:),  Nreserve_Ltm1(ij,:), Preserve_Ltm1(ij,:), Kreserve_Ltm1(ij,:), ...
        Nuptake_L(ij,:),      Puptake_L(ij,:),     Kuptake_L(ij,:),     FNC_Ltm1(ij,:),      Tden_Ltm1(ij,:), ...
        AgePl_Ltm1(ij,:),     fab_L,               fbe_L,               ParEx_L,             Mpar_L, ...
        TBio_L(ij,:),         SAI_Ltm1(ij,:),      hc_Ltm1(ij,:),       ExEM,                Lmax_day, ...
        L_day,                Se_bio,              Tdp_bio,             OPT_EnvLimitGrowth,  OPT_VD, ...
        OPT_VCA,              OPT_ALLOME,          OPT_SoilBiogeochemistry);
    
            BLit(ij,:)= 0.0 ; % %% %%[kg DM / m2] - unsure what this is?
            end        

            %% HYDROLOGY MODULE
            %==============================================================
            % FUNTION: HYDROLOGY_MODULE_PAR
            %==============================================================                                                 

    [V(ij,:),           O(ij,:),          Vice(ij,:),       Oice(ij,:),          ZWT(ij), ...         % 1
     OF(ij),            OS(ij),           OH(ij,:),         OL(ij,:),            Psi_s_H(ij,:),...    % 2
     Psi_s_L(ij,:),     Rd(ij),           Qi_out(ij,:),     Rh(ij),              Lk(ij), ...          % 3
     f(ij),             WIS(ij),          Ts(ij),           Csno(ij),            Cice(ij), ...        % 4 
     NDVI(ij),          Pr_sno(ij),       Pr_liq(ij),       rb_H(ij,:),          rb_L(ij,:), ...      % 5 
     rs_sunH(ij,:),     rs_sunL(ij,:),    rs_shdH(ij,:),    rs_shdL(ij,:),       rap_H(ij,:), ...     % 6
     rap_L(ij,:),       r_soil(ij),       b_soil(ij),       alp_soil(ij),        ra(ij), ...          % 7
     r_litter(ij,:),    WR_SP(ij),        U_SWE(ij),        NIn_SWE(ij),         dQ_S(ij), ...        % 8 
     DQ_S(ij),          DT_S(ij),         WAT(ij),          ICE(ij),             ICE_D(ij), ...       % 9
     IP_wc(ij),         WR_IP(ij),        NIce(ij),         Cicew(ij),           Csnow(ij), ...       % 10
     FROCK(ij),         Dr_H(ij,:),       Dr_L(ij,:),       SE_rock(ij),         SE_urb(ij), ...      % 11
     Lk_wat(ij),        Lk_rock(ij),      An_L(ij,:),       An_H(ij,:),          Rdark_L(ij,:), ...   % 12
     Rdark_H(ij,:),     Ci_sunH(ij,:),    Ci_sunL(ij,:),    Ci_shdH(ij,:),       Ci_shdL(ij,:), ...   % 13
     Rn(ij),            H(ij),            QE(ij),           Qv(ij),              Lpho(ij), ...        % 14
     T_H(ij,:),         T_L(ij,:),        EIn_H(ij,:),      EIn_L(ij,:),         EG(ij), ...          % 15
     ELitter(ij),       ESN(ij),          ESN_In(ij),       EWAT(ij),            EICE(ij), ...        % 16
     EIn_urb(ij),       EIn_rock(ij),     dw_SNO(ij),       Imelt(ij),           Smelt(ij), ...       % 17
     G(ij),             Gfin(ij),         Tdp(ij,:),        Tdp_snow(ij,:),      Tdeb(ij,:), ...      % 18
     Tice(ij),          Tdamp(ij),        Tdp_H(ij,:),      Tdp_L(ij,:),         SWE(ij), ...         % 19
     SND(ij),           ros(ij),          In_SWE(ij),       SP_wc(ij),           Qfm(ij), ...         % 20 
     t_sls(ij),         In_H(ij,:),       In_L(ij,:),       In_Litter(ij),       In_urb(ij), ...      % 21
     In_rock(ij),       gsr_H(ij,:),      Psi_x_H(ij,:),    Psi_l_H(ij,:),       Jsx_H(ij,:), ...     % 22
     Jxl_H(ij,:),       Kleaf_H(ij,:),    Kx_H(ij,:),       Vx_H(ij,:),          Vl_H(ij,:),...       % 23
     gsr_L(ij,:),       Psi_x_L(ij,:),    Psi_l_L(ij,:),    Jsx_L(ij,:),         Jxl_L(ij,:), ...     % 24
     Kleaf_L(ij,:),     Kx_L(ij,:),       Vx_L(ij,:),       Vl_L(ij,:),          fapar_H(ij,:), ...   % 25
     fapar_L(ij,:),     SIF_H(ij,:),      SIF_L(ij,:),      Ws_under(ij),        er(ij), ...          % 26
     snow_albedo(ij,:), tau_sno(ij),      e_sno(ij),        HV(ij),              QEV(ij), ...         % 27
     dQVEG(ij),         TsVEG(ij),        Ts_under(ij),     EK(ij),              POT(ij,:), ...       % 28
     CK1(ij)] = ...
                HYDROLOGY_MODULE_PAR(...
     Vtm1(ij,:),            Oicetm1(ij,:),       aR,               Zs,                  EvL_Zs, ...
     Inf_Zs,                Zinf,                RfH_Zs,           RfL_Zs,              dz,  ...
     Dz,                    ms,                  Kbot,             Pr_S(ij),           Ta_S(ij), ...
     Ds_S(ij),             Ws_S(ij),           zatm,             Tstm1(ij),           dt, ...
     dth,                   ea_S(ij),           N_S(ij),         Pre_S(ij),          Tstm0(ij), ...
     LAI_H(ij,:),           SAI_H(ij,:),         LAI_L(ij,:),      SAI_L(ij,:),         LAIdead_H(ij,:), ...
     LAIdead_L(ij,:),       Rrootl_H(ij,:),      Rrootl_L(ij,:),   BLit(ij,:),          Sllit,   ...
     Kct,                   Datam_S,             IniCond.DeltaGMT, Lon,                 Lat,     ...
     t_bef,                 t_aft,               Ccrown,        Cbare,               Crock,  ...
     Curb,                  Cwat,                SAB1_S(ij),      SAB2_S(ij),         SAD1_S(ij), ...
     SAD2_S(ij),           PARB_S(ij),         PARD_S(ij),      SvFn(ij),            SNDtm1(ij),  ...
     snow_albedotm1(ij,:),  Color_Class,         OM_H,             OM_L,                PFT_opt_H,     ...
     PFT_opt_L,             hc_H(ij,:),          hc_L(ij,:),       d_leaf_H,            d_leaf_L,     ...
     Soil_Param,            Interc_Param,        SnowIce_Param,    VegH_Param,          VegL_Param,  ...
     Ca_S(ij),              Oa,                  Citm1_sunH(ij,:), Citm1_shdH(ij,:),    Citm1_sunL(ij,:), ...
     Citm1_shdL(ij,:),      e_rel_H(ij,:),       e_relN_H(ij,:),   e_rel_L(ij,:),       e_relN_L(ij,:), ...
     e_snotm1(ij),          In_Htm1(ij,:),       In_Ltm1(ij,:),    In_Littertm1(ij),    In_urbtm1(ij), ...
     In_rocktm1(ij),        SWEtm1(ij),          In_SWEtm1(ij),    Tdebtm1(ij,:),       Ticetm1(ij), ...
     Tdptm1(ij,:),          Tdp_snowtm1(ij,:),   Tdamptm1(ij),     Ts_undertm1(ij),     WATtm1(ij), ...
     ICEtm1(ij),            IP_wctm1(ij),        ICE_Dtm1(ij),     Cicewtm1(ij),        Vx_Htm1(ij,:), ...
     Vl_Htm1(ij,:),         Vx_Ltm1(ij,:),       Vl_Ltm1(ij,:),    Psi_x_Htm1(ij,:),    Psi_l_Htm1(ij,:), ...
     Psi_x_Ltm1(ij,:),      Psi_l_Ltm1(ij,:),    ZR95_H,           ZR95_L,              FROCKtm1(ij),  ...
     Krock,                 Urb_Par,             Deb_Par,          Zs_deb,              Tdew_S(ij), ...
     t_slstm1(ij),          rostm1(ij),          SP_wctm1(ij),     fpr,                 IrD_S(ij), ...
     In_max_urb,            In_max_rock,         K_usle,           tau_snotm1(ij),      Ta_Spdind(:,ij), ...
     Slo_top2(ij),          Slo_head(ij,:),      Asur(ij),         Ared(ij),            aTop(ij), ...
     EKtm1(ij),             q_runon(ij),         Qi_in(ij,:),      Ws_undertm1(ij),     Pr_sno_t(ij,:), ...
     pow_dis,               a_dis,               Salt_S(ij),       SPAR,                SNn(ij), ...
     min_SPD,               OPT_VegSnow,         OPT_SoilTemp,     OPT_PlantHydr,       Opt_CR, ...
     Opt_ST,                Opt_ST2,             OPT_SM,           OPT_STh,             OPT_FR_SOIL, ...
     OPT_PH,                parameterize_phase,  hSTL,             OPT_Albsno);
    %catch ME
    %disp(['Error in HYDROLOGY_MODULE_PAR in ij:' char(num2str(ij))])
    
    % Get the full error identifier (type/ID)
    %errorID = ME.identifier;
    %disp(['Error Type (ID): ', errorID]);

    %end
    
          catch ME
            [r_err, c_err] = ind2sub([m_cell n_cell], ij);
            if isempty(ME.stack)
                where_err = 'unknown location';
            else
                where_err = sprintf('%s, line %d', ME.stack(1).name, ME.stack(1).line);
            end
            msg_err = sprintf(['Cell ij = %d (row %d, col %d), elevation %.0f m, class %d, ' ...
                'model date %04d-%02d-%02d %02d:00\n%s\n(in %s)'], ...
                ij, r_err, c_err, DTMn(ij), ksv(ij), Datam_S(1:4), ME.message, where_err);
            error('TC:cellFailed', '%s', msg_err);
          end
    end
end
    
%% END OF PARFOR
%======================================================================
    %%post-compute sublimation from ESN, do this inside hydrology module
    %%once it turns out to be useful
    SSN = ESN.*(Ts<0);
    SSN_In = ESN_In.*(Ts<0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii=1:ms_max
        Qi_out_Rout(:,:,ii)=   reshape(Qi_out(:,ii),m_cell,n_cell); %[mm/h]
    end
    Rd= reshape(Rd,m_cell,n_cell); %%[mm]
    Rh= reshape(Rh,m_cell,n_cell); %%[mm]
    
    %% ROUTING MODULE
    %======================================================================
    %Should values be full grid on vector grid? Here especially DTM?
    
    [q_runon,    Q_channel,    Qi_in_Rout,   Slo_pot,     Q_exit, ...
     Qsub_exit,  T_pot,        QpointH,      QpointC,     UpointH, ...
     UpointC]= ...
                     ROUTING_MODULE( ...
     dt,         dth,          Rd,           Rh,          Qi_out_Rout, ...
     Q_channel,  cellsize,     Area,         DTM,         NMAN_H, ...
     NMAN_C,     MRough,       WC,           SN,          T_flow, ...
     T_pot,      Slo_top,      ms_max,       POT,         ZWT, ...
     OPT_HEAD,   Xout,         Yout); 
    
    for ii=1:ms_max
        Qi_in(:,ii)=	reshape(Qi_in_Rout(:,:,ii),num_cell,1); %[mm]
        Slo_head(:,ii)=	reshape(Slo_pot(:,:,ii),num_cell,1);
    end

    Rd = reshape(Rd,num_cell,1); %%[mm]
    Rh = reshape(Rh,num_cell,1); %%[mm]
    q_runon = reshape(q_runon,num_cell,1); %%[mm]
    Qi_in=Qi_in/dth;%% [mm/h]
    q_runon = q_runon/dth; %%% [mm/h]
    %%% Q_exit Qsub_exit [mm] over the entire domain
    if not(isreal(sum(sum(q_runon))))
        disp('The Program fails because of runon numerical instability')
        break
    end

    
    %% Glacier volume redistribution 

if (Datam_S(2)==1) && (Datam_S(3)==1) && (Datam_S(4)==1) && (OPT_Idyn > 0)
        
        GLA = (GLH > 0) & (MASK == 1); % Glacier mask

        % sum up total water equivalent in [m w.e.], store ratios
        WEbef = ICE + SWE; %total w.e. before redistribution
        WEym1 = ICEym1 + SWEym1; %total w.e. from the year before

        raICE = ICE./WEbef; raICE(isnan(raICE))=0; %ice ratio of total
        raSWE = SWE./WEbef; raSWE(isnan(raSWE))=0; %snow ratio of total mm w.e.

        THbef = reshape(WEbef,m_cell,n_cell)./916; %conversion to thickness in m, rhi parameter needed, and reshape
        THym1 = reshape(WEym1,m_cell,n_cell)./916; %conversion to thickness in m, rhi parameter needed, and reshape

        THbef(GLA ~= 1) = 0;
        THym1(GLA ~= 1) = 0;
        
        SMB = (WEbef - WEym1)./916; % Distributed surface mass balance of the previous year in m/y ice-eq (IGM)
        SMB(GLH == 0 | MASK~=1 | GLA ~= 1) = -10;
        SMB(isnan(SMB)) = -10;
        SMB = reshape(SMB,m_cell,n_cell);
 
        %%%% Create or update the .nc file %%%%%%%%%%%%

        x(2) = x(1) + (y(2)-y(1)); % Temporary fix
        IGM_netcdf([outlocation 'igm_inputs3.nc'],x,y,rot90(flipud(double(GLA)),3),rot90(flipud(THym1),3),rot90(flipud(DTM_Bedrock),3)...
            , rot90(flipud(DTM_orig),3),rot90(flipud(SMB),3)) % Create .nc file

        %%%% Run IGM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pyenv('Version',path_igmEnv, 'ExecutionMode','OutOfProcess');   % setup python virtual env
        igm_run_path =  [path_igm '/igm/igm_run.py']; % path of IGM .py script
        addpath(genpath([path_igm '/igm/']))
 
        igm_arg = [' --lncd_input_file ' [outlocation 'igm_inputs3.nc']...
            ' --working_dir ' outlocation...
            ' --time_start ' num2str(1)...
            ' --time_end ' num2str(2) ... 
            ' --time_save ' num2str(1)];

        pyrunfile([igm_run_path igm_arg]) % Run IGM

        thk = ncread([outlocation 'output.nc'],'thk'); % Load updated glacier thickness from IGM

        THnew = double(rot90(thk(:,:,2))); % New thickness after IGM

        %Check IGM mass conservation
        dMB_igm = nansum(THnew,'all') - nansum(THbef,'all');

        % Add the icethickness which leaks outside of glacier mask back on existing glacier
        IGM_smb = nansum(THnew,'all') - nansum(THym1,'all');
        TC_smb = nansum(SMB(GLA>0),'all');

        THnew(flipud(GLA == 1) & THnew > 0) = THnew(flipud(GLA == 1) & THnew > 0) + (TC_smb-IGM_smb)./nansum(GLA>0,'all');
        THnew(flipud(GLA ~= 1)) = 0; %No glacier growing in area

        % reshape total water equivalent back to array map
        WEnew = reshape(flipud(THnew),num_cell,1).*916; % conversion back to mm w.e., rhi parameter needed, and reshape
        
        % recalculate snow/ice w.e. maps based on old ratios
        ICE = raICE.*WEnew; % re-calculate ice w.e.
        SWE(WEnew>0) = raSWE(WEnew>0).*WEnew(WEnew>0); % re-calculate SWE, only where glacier

        %Mass balance check
        dMass = nansum(ICE(WEnew >0) + SWE(WEnew >0)) - nansum(WEbef(WEnew >0));

        disp(['Total SMB: ' num2str(TC_smb)])
        disp(['SMB after IGM redistribution: ' num2str(dMass/916)])
        
        % store for next year redistribution
        ICEym1 = ICE;     
        SWEym1 = SWE;

        % calculate new snow/ice depths
        ICE_D = ICE/916; ICE_D(isnan(ICE_D))=0;
        SND = SWE./ros; SND(isnan(SND))=0;
end


    %% AVALANCHES COMPONENT

    SND=    reshape(SND,m_cell,n_cell); SND(isnan(SND)) = 0;
    SWE=    reshape(SWE,m_cell,n_cell); SWE(isnan(SWE)) = 0;
    ros=    reshape(ros,m_cell,n_cell); ros(isnan(ros)) = 0;

    if OPT_Aval == 1
        SWEpreava = SWE; 

        [SND,SWE,ros,Swe_exit]= AVALANCHES(DTM,cellsize,Area,...
        reshape(Asur,m_cell,n_cell),Slo_top,SND,SWE,ros,a_aval,C_aval);
        SWE_avalanched = SWE-SWEpreava;   
        clear SWEpreava;
    else 
        SWE_avalanched = SWE.*0;
        Swe_exit = SWE.*0;
    end 

    %% Catchment average snowline    
    SND=	reshape(SND,num_cell,1); SND(isnan(SND)) = 0;
    SWE=    reshape(SWE,num_cell,1); SWE(isnan(SWE)) = 0;
    ros=    reshape(ros,num_cell,1); ros(isnan(ros)) = 0;
    SWE_avalanched = reshape(SWE_avalanched,num_cell,1);

    %% FRACTURED ROCK COMPONENT   
    [Q_channel,FROCK,Qflow_rock]= FRACTURED_ROCK(Q_channel,FROCK,...
        SPRINGn,dth,m_cell,n_cell,num_cell,Kres_Rock);

    if t==2
        V_tgtm1=    sum(Ared.*Asur.*sum(Vtm1,2))*(cellsize^2)/Area;
        Vice_tgtm1= sum(Ared.*Asur.*sum(Vicetm1,2))*(cellsize^2)/Area;
        SWE_tgtm1=  sum(SWEtm1)*(cellsize^2)/Area;
        In_tgtm1=   (sum(sum(In_Htm1)) + sum(sum(In_Ltm1)) + ...
            sum(sum(In_Littertm1)) + sum(SP_wctm1) + ...
            sum(In_SWEtm1) + sum(In_urbtm1) + ...
            sum(In_rocktm1) + sum(IP_wctm1) ) * (cellsize^2)/Area;
        ICE_tgtm1=  sum(ICEtm1)*(cellsize^2)/Area; %%
        WAT_tgtm1=  sum(WATtm1)*(cellsize^2)/Area; %%
        FROCK_tgtm1=sum(FROCKtm1)*(cellsize^2)/Area; %%
    end
       
    %% ALBEDO MAP

    % if t==2 
    %     snow_albedo_out = snow_albedotm1;
    %     surface_albedo_out = surface_albedotm1;
    % elseif Datam_S(4)==12
    % snow_albedo_out = snow_albedo;
    % surface_albedo_out = surface_albedotm1;
    % surface_albedo_out(ICE_D>0)=0.28; %THESE ARE NOT CORRECT?
    % surface_albedo_out(DEB_MAP>0)=0.13; %THESE ARE NOT CORRECT?
    % surface_albedo_out(SND>0)=snow_albedo(SND>0);
    % end

%% INITIAL CONDITION FOR THE NEXT STEP
%==========================================================================
    Cicewtm1=          Cicew;
    Citm1_shdH=        Ci_shdH;
    Citm1_shdL=        Ci_shdL;
    Citm1_sunH=        Ci_sunH;
    Citm1_sunL=        Ci_sunL;
    e_snotm1=          e_sno;
    EKtm1=             EK;
    FROCKtm1=          FROCK;
    ICE_Dtm1=          ICE_D;
    ICEtm1=            ICE;
    In_Htm1=           In_H ;
    In_Littertm1=      In_Litter;
    In_Ltm1=           In_L;
    In_rocktm1=        In_rock;
    In_SWEtm1=         In_SWE;
    In_urbtm1=         In_urb;
    IP_wctm1=          IP_wc;
    Oicetm1=           Oice;
    Psi_l_Htm1=	       Psi_l_H;
    Psi_l_Ltm1=	       Psi_l_L;
    Psi_x_Htm1=	       Psi_x_H;
    Psi_x_Ltm1=	       Psi_x_L;
    rostm1=            ros;
    SNDtm1=            SND ;
    snow_albedotm1=    snow_albedo ;
    soil_albedotm1=    soil_albedo ;
    SP_wctm1=          SP_wc;
    surface_albedotm1= surface_albedo ;
    SWEtm1=            SWE;
    SWE_avalanchedtm1= SWE_avalanched ;
    t_slstm1=          t_sls;
    tau_snotm1=        tau_sno;
    Tdamptm1=          Tdamp;
    Tdebtm1=           Tdeb;
    Tdp_snowtm1 =      Tdp_snow;
    Tdptm1=            Tdp;
    Ticetm1=           Tice;
    %%% order!
    Tstm1=             Ts;
    Tstm0=             2*Ts-Tstm1;
    Ts_undertm1=       Ts_under;
    %%%
    Vl_Htm1=           Vl_H;
    Vl_Ltm1=           Vl_L;
    Vtm1=              V;
    Vx_Htm1=           Vx_H;
    Vx_Ltm1=           Vx_L;
    WATtm1=            WAT;
    Ws_undertm1=Ws_under;
    
    if (Datam_S(4)==1)
    AgeDL_Htm1=     AgeDL_H;
    AgeDL_Ltm1=     AgeDL_L;
    AgeL_Htm1=      AgeL_H;
    AgeL_Ltm1=      AgeL_L ;
    AgePl_Htm1=     AgePl_H;
    AgePl_Ltm1=     AgePl_L;
    B_Htm1=         B_H;
    B_Ltm1=         B_L;
    Bfac_weekHtm1=  Bfac_weekH;
    Bfac_weekLtm1=  Bfac_weekL;
    Ccrown_t_tm1=   Ccrown_t;
    dflo_Htm1=      dflo_H;
    dflo_Ltm1=      dflo_L ;
    FNC_Htm1=       FNC_H;
    FNC_Ltm1=       FNC_L;
    hc_Htm1=        hc_H;
    hc_Ltm1=        hc_L;
    Kreserve_Htm1=  Kreserve_H;
    Kreserve_Ltm1=  Kreserve_L;
    NBLeaf_Htm1=    NBLeaf_H;
    NBLeaf_Ltm1=    NBLeaf_L;
    NBLI_Htm1=      NBLI_H;
    NBLI_Ltm1=      NBLI_L;
    NPP_Htm1=       NPP_H;
    NPP_Ltm1=       NPP_L;
    NPPI_Htm1=      NPPI_H;
    NPPI_Ltm1=      NPPI_L;
    Nreserve_Htm1=  Nreserve_H;
    Nreserve_Ltm1=  Nreserve_L;
    NuLit_Htm1=     NuLit_H;
    NuLit_Ltm1=     NuLit_L;
    NupI_Htm1=      NupI_H;
    NupI_Ltm1=      NupI_L;
    PARI_Htm1=      PARI_H;
    PARI_Ltm1=      PARI_L;
    PHE_S_Htm1=     PHE_S_H;
    PHE_S_Ltm1=     PHE_S_L;
    Preserve_Htm1=  Preserve_H;
    Preserve_Ltm1=  Preserve_L;
    SAI_Htm1=       SAI_H;
    SAI_Ltm1=       SAI_L;
    Tden_Htm1=      Tden_H;
    Tden_Ltm1=      Tden_L;
    TdpI_Htm1=      TdpI_H;
    TdpI_Ltm1=      TdpI_L;
    end
    
    
    % MEMORY CONDITION FOR VEGETATION  MODEL
    %----------------------------------------------------------------------
    % 1 day
    if t > 24
        An_H_t(:,:,1:23)=       An_H_t(:,:,2:24);
        An_H_t(:,:,24)=         An_H;
        An_L_t(:,:,1:23)=       An_L_t(:,:,2:24);
        An_L_t(:,:,24)=         An_L;
        O_t(:,:,1:23)=          O_t(:,:,2:24);
        O_t(:,:,24)=            O;
        PAR_t(:,1:23)=          PAR_t(:,2:24);
        PAR_t(:,24)=            PARB_S + PARD_S;
        Pr_sno_t(:,1:23)=       Pr_sno_t(:,2:24);
        Pr_sno_t(:,24)=         Pr_sno;
        Psi_l_H_t(:,:,1:23)=    Psi_l_H_t(:,:,2:24);
        Psi_l_H_t(:,:,24)=      Psi_l_H;
        Psi_l_L_t(:,:,1:23)=    Psi_l_L_t(:,:,2:24);
        Psi_l_L_t(:,:,24)=      Psi_l_L;
        Psi_x_H_t(:,:,1:23)=    Psi_x_H_t(:,:,2:24);
        Psi_x_H_t(:,:,24)=      Psi_x_H;
        Psi_x_L_t(:,:,1:23)=    Psi_x_L_t(:,:,2:24);
        Psi_x_L_t(:,:,24)=      Psi_x_L;
        Rdark_H_t(:,:,1:23)=    Rdark_H_t(:,:,2:24);
        Rdark_H_t(:,:,24)=      Rdark_H;
        Rdark_L_t(:,:,1:23)=    Rdark_L_t(:,:,2:24);
        Rdark_L_t(:,:,24)=      Rdark_L;
        Ta_t(:,1:23)=           Ta_t(:,2:24);
        Ta_t(:,24)=             Ta_S;
        Tdp_H_t(:,:,1:23)=      Tdp_H_t(:,:,2:24);
        Tdp_H_t(:,:,24)=        Tdp_H;
        Tdp_L_t(:,:,1:23)=      Tdp_L_t(:,:,2:24);
        Tdp_L_t(:,:,24)=        Tdp_L;
        Tdp_t(:,:,1:23)=        Tdp_t(:,:,2:24);
        Tdp_t(:,:,24)=          Tdp;
        V_t(:,:,1:23)=          V_t(:,:,2:24);
        V_t(:,:,24)=            V;
    else
        An_H_t(:,:,t)=          An_H;
        An_L_t(:,:,t)=          An_L;
        O_t(:,:,t)=             O;
        PAR_t(:,24)=            PARB_S + PARD_S;
        Pr_sno_t(:,t)=          Pr_sno;
        Psi_l_H_t(:,:,t)=       Psi_l_H;
        Psi_l_L_t(:,:,t)=       Psi_l_L;
        Psi_x_H_t(:,:,t)=       Psi_x_H;
        Psi_x_L_t(:,:,t)=       Psi_x_L;
        Rdark_H_t(:,:,t)=       Rdark_H;
        Rdark_L_t(:,:,t)=       Rdark_L;
        Ta_t(:,t)=              Ta_S;
        Tdp_H_t(:,:,t)=         Tdp_H;
        Tdp_L_t(:,:,t)=         Tdp_L;
        Tdp_t(:,:,t)=           Tdp;
        V_t(:,:,t)=             V;
    end



    %% OUTPUTS AND CHECKPOINTS
    %======================================================================
    % Aggregation and writing of the outputs chosen in
    % <site>/Parameters/Output_Config.m, restart checkpoints (T&C_Code/Output_Manager.m)
    %======================================================================
    OM_phase = 'step'; Output_Manager;

end
%close(bau)
OM_phase = 'close'; Output_Manager;

%% Display
Computational_Time =toc;
disp('COMPUTATIONAL TIME [h] ')
disp(Computational_Time/3600)

%Q_channel