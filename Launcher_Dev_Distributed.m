%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C) - ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AUTHOR INFO AND STUDY SITE
%==========================================================================
% Updated 15/09/2026 to makeinto Development version for Mountain Water
% Author: Cat Fyffe
% Code originally from: Max Rodriguez and ACHILLE JOUBERTON
% Area of Study: Rio Santa
% Region: Shallap
% Code explanation: This code launches TC model.
%==========================================================================

%% THINGS TO CHECK/WATCH LATER

%UPDATE PARAMS WITH SIMONES
%Initial albedo map
%rho_g
%Soil paramters issues Max mentioned
%CHECK times as starts Date-1 but saves against Date?

%% CLEAR ALL
clc; clear;
% clear all   % Don't clear all such that it can be run on the HPC cluster without any issues 
% delete(gcp('nocreate'))

%=========================================================================
%% SITE SETUP
%==========================================================================

site_num_list = [1 2]; %Please number sequentially
site_name_list = {'Shallap' 'Artesonraju'};
num_sites = size(site_num_list,2);
site_num = 1; %Choose site to run, this then selects the correct site directories
site_name = site_name_list{site_num};

%===================================================
%% Load config info
%==================================================

% Sub-path for Config file
Directories.config = [site_name '/Parameters/Site_config_TC.xlsx'];

CONFIG = readcell(Directories.config);
CON_vals_id = find(CONFIG(1,:)=="Value");
CON_label_id = find(CONFIG(1,:)=="Variable");
CONFIG_vals = cell2struct(CONFIG(2:end,CON_vals_id),CONFIG(2:end,CON_label_id)'); %Now a structure of the config info

IniCond.DeltaGMT = CONFIG_vals.DeltaGMT; %Define change from GMT
IniCond.run_folder = CONFIG_vals.Run_Folder; %Define folder to save outputs

%%======================================================================
% Model structure choices
%========================================================================

OPT_Forcing = CONFIG_vals.OPT_Forcing; %Choose type of forcing, 1= from AWS which is distributed
OPT_Veg_Param = CONFIG_vals.OPT_Veg_Param; %Choose if vegetation parameters (for Ccrown, Cwat, Crock, Curb, Cbare) are one value per landcover type (1), or gridded and variable per landcover (2). 
%For 1 change the values in VPAR_T directly, for 2 provide the grids in the
%dtm_file. Note for 2 you still need the table to read the Veg_type names
%but the values are not used. 

%% DIRECTORIES - ALL RELATIVE - NO NEED TO UPDATE 
%==========================================================================
% All paths are here
%==========================================================================

% Sub-path for outputs
Directories.save = ['Outputs/Distributed/' IniCond.run_folder '/']; %Do we make this site specific?
Directories.restart = Directories.save;

% Sub-path for forcings
Directories.forc = [site_name '/Forcing']; 
if OPT_Forcing == 1
Directories.forc_meteo = [Directories.forc '/' CONFIG_vals.forc_meteo_file]; meteo_name = CONFIG_vals.forc_meteo_name;
ta_lapse_file = CONFIG_vals.ta_lapse_file;  ta_lapse_name = CONFIG_vals.ta_lapse_name; %These are lapse rates per hour and month, in a table with a column with hour, a column with month, and a column with lapse rates (positive values)
pr_lapse_file = CONFIG_vals.pr_lapse_file;  pr_lapse_name = CONFIG_vals.pr_lapse_name; %This is a structure Pr_lapse(m).month where each month has a linear model
Meteo_el = CONFIG_vals.Meteo_el; %Elevation of station, used if OPT_forcing = 1
end

% Parameters
Directories.Vegpar = [site_name '/Parameters/Parameters_TC.xlsx']; 
Directories.Optpar = [site_name '/Parameters/Options_and_NoVegParams.xlsx'];
Directories.VPAR = [site_name,'/Parameters/VPAR_T.xlsx'];

% Preprocessing information
Directories.PreProc = [site_name '/Preprocessing/OUTPUTS/' CONFIG_vals.PreProc_folder '/']; 
dtm_file = CONFIG_vals.dtm_file;

% Sub-path points of interest for discharge and multipoints
Directories.POIs = [site_name '/Preprocessing/OUTPUTS/']; %POI_dtm_Shallap_50m.txt
Directories.Multi = ['Multipoint/OUTPUTS/' site_name '_MultiPoints.txt'];

% Dependencies
addpath(genpath('Functions')); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([site_name '/Preprocessing/OUTPUTS'])); % Where are distributed model set-up files (needed ? yes to load dtm)
%addpath(genpath([Directories.model,'/T_and_C/TC_setups/' site_name '/RUNS/INPUTS'])); % Where is located the meteorological forcing and Shading matrix 
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
Points = readtable(Directories.Multi); %import table with points info
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

OPT_ALLOME = OPT_PARAM_vals.OPT_ALLOME; %% Option for vegetation - structural attributes?
OPT_EnvLimitGrowth = OPT_PARAM_vals.OPT_EnvLimitGrowth; %Option for introducing Environmental Limitation of Growth
OPT_FR_SOIL = OPT_PARAM_vals.OPT_FR_SOIL; % Option for freezing soil
OPT_HEAD = OPT_PARAM_vals.OPT_HEAD; %% Option for determining routing method
OPT_PlantHydr = OPT_PARAM_vals.OPT_PlantHydr; %Option for including Plant Hydrualic
OPT_SoilBiogeochemistry = OPT_PARAM_vals.OPT_SoilBiogeochemistry; % Option for including Soil Biochemistry
OPT_SoilTemp = OPT_PARAM_vals.OPT_SoilTemp; % Option for computing soil temperature or not
OPT_VCA = OPT_PARAM_vals.OPT_VCA; %% Option for vegetation 
OPT_VegSnow = OPT_PARAM_vals.OPT_VegSnow; % Option for computing energy budget of vegetation when there is snow at the ground

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

ksv=reshape(VEG_CODE,num_cell,1);
vpar_list = unique(ksv(ksv~=0)); %Vegetation classes
No_vpar = size(vpar_list,1);  %Number of vegetation classes

%Veg/land parameters (was POI but renamed to reduce confusion with Points
%of interest)
%Load table (needed for both options)
VPAR_T = readtable(Directories.VPAR);

if OPT_Veg_Param==1
    %In this case there are single values per land cover
    %Turn into a grid
    VPAR.Ccrown = zeros(m_cell,n_cell);
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
    VPAR.Curb = CWATER.*0; %As urban 0 everywhere
    VPAR.Crock = CROCK;
    VPAR.Cbare = CBARE;
end  

%For both types
VPAR.Class = (vpar_list); %Just a list of the numbers
VPAR.Veg_type = string(VPAR_T.Veg_type);

%Reshape
VPAR.Ccrownr = reshape(VPAR.Ccrown,num_cell,1);
VPAR.Cwatr = reshape(VPAR.Cwat,num_cell,1);
VPAR.Curbr = reshape(VPAR.Curb,num_cell,1);
VPAR.Crockr = reshape(VPAR.Crock,num_cell,1);
VPAR.Cbarer = reshape(VPAR.Cbare,num_cell,1);
cc_max = 1;

% SPATIAL INDICES PER LAND COVER CLASS
for v=1:No_vpar
    fieldName = ['Veg', num2str(v)];
    idxCode.(fieldName) = find(VEG_CODE == v & MASK == 1);
end

%Save the parameter files with the model outputs
writetable(TT_par,[Directories.save 'Parameters/' 'TT_par.xlsx']);
writetable(OPT_PARAM,[Directories.save 'Parameters/' 'OPT_PARAM.xlsx']);
writetable(VPAR_T,[Directories.save 'Parameters/' 'VPAR_T.xlsx']);
writecell(CONFIG,[Directories.save 'Parameters/' 'CONFIG.xlsx']);

%-----------------------------------------------------------------------
% 3. Glacier and snow parameters
% --------------------------------------------------------------------------

%Tmod parameters
OPT_Tmod = OPT_PARAM_vals.OPT_Tmod;
if OPT_Tmod == 0 %No Tmod applied
    TmodB = 0;
elseif OPT_Tmod == 1 %Simple bias correction
    TmodB = OPT_PARAM_vals.TmodB; %Tmod bias (positive reduces temperature over ice)
elseif OPT_Tmod == 2 %Uses a relationship with optional bias correction 
    TmodB = OPT_PARAM_vals.TmodB; %Tmod bias (positive reduces temperature over ice)
    TmodM = OPT_PARAM_vals.TmodM; %Tmod multiplier
    TmodC = OPT_PARAM_vals.TmodC; %Tmod coefficent
end

% Glacier dynamics and avalanching
OPT_Idyn = OPT_PARAM_vals.OPT_Idyn; % Switch for glacier dynamics
OPT_Aval = OPT_PARAM_vals.OPT_Aval; % 1 to turn on avalanching, 0 to turn off avalanching
a_aval = OPT_PARAM_vals.a_aval; %avalanche parameters a (Bernhart & Schulz 2010) - from TOPKAPI 0.17245 - 0.12/145 from Jouberton et al.(2025) - 99.05/0.1012 from Buri et al. (2023)
C_aval = OPT_PARAM_vals.C_aval; % avalanche parameters C 

% Initial snow depth and albedo
%diff_IniSND = 1; %NOT USED
%fn_IniSnowDepth = SNOWD; %NOT USED!! - SNOWD used directly in initial conditions
%fn_IniSnowAlbedo = 'Cinca_Init_Snow_Albedo_virtual.mat'; %TO CREATE

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
hSTL = OPT_PARAM_vals.hSTL; % Skin layer thickness for the 2-layer snowpack module
min_SPD = OPT_PARAM_vals.min_SPD;  %% [m] minimum snow pack depth to have a multilayer snow 

% Choice of the snow albedo scheme Evars_load
OPT_Albsno = OPT_PARAM_vals.OPT_Albsno; % 3 doesn't work, 4 is Brock 2000, 5 is Ding 2017

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
OPT_Alb_v_el = OPT_PARAM_vals.OPT_Alb_v_el; %Switch for using firn albedo versus elevation
Aice = OPT_PARAM_vals.Aice; %Aice, Artesonraju % 0.2314 Shallap point value
Amax = OPT_PARAM_vals.Amax; %Max 'firn' albedo value

if OPT_Alb_v_el==1

    Alb_v_el_MM = OPT_PARAM_vals.Alb_v_el_MM; %Yota remote sensing slope for 25% lowest scenes
    Alb_v_el_CC = OPT_PARAM_vals.Alb_v_el_CC; %Yota remote sensing intercept for 25% lowest scenes

    Afirn = (DTM.*Alb_v_el_MM) + Alb_v_el_CC; %Construct initial relationship
    Afirn(Afirn<Aice)=Aice; %Prevent values lower than Aice
    Afirn(Afirn>Amax)=Amax;
      
   disp('Albedo varies with elevation')
else
    Afirn = DTM.*0 + Aice; %Artesonraju % 0.2314 Shallap point value
    disp(strcat('Constant bare-ice albedo of',{' '},num2str(Aice),{' '},'used'))
end

% Check if initial snow albedo is given - 
% this goes into initial conditions as SNOWALB
%SNOWD is taken directly and used in intial conditions as initial snow depth
if ~exist('SNOWALB','var')
    Ini_Asnow = OPT_PARAM_vals.Ini_Asnow;
    SNOWALB = SNOWD;
    SNOWALB(SNOWD>0) = Ini_Asnow;
end 

%-----------------------------------------------------------------
% 4. Carbon
%-------------------------------------------------------------------------

load(['Ca_Data.mat']);
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

% MASK=ones(m_cell,n_cell); MASK(isnan(DTM))=0;
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
Xout_long = Xout; %This is the list of POIs (loaded with spatial data)
Yout_long = Yout; %This is the list of POIs (loaded with spatial data)
Xout = Xoutlet; % Location of outlet discharge (column)
Yout = Youtlet; % Location of outlet discharge (this is the row)
NAMEout = POI_names(1);

Slo_top(Youtlet,Xoutlet)=0.05; %NOTE on maps and in scatter use (X,Y) but to index use Y (row) X (col)
npoint = length(Xout);
Area= (cellsize^2)*sum(sum(MASK)); %% Projected area [m^2]


% Flow potential
%--------------------------------------------------------------------------
ms_max = OPT_PARAM_vals.ms_max; %% Number of soil layers
T_pot=cell(1,ms_max);
for jk=1:ms_max
    T_pot{jk}= T_flow;
end

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
rho_g = OPT_PARAM_vals.rho_g; %%% Spatial Albedo - this is quite high??
md_max = OPT_PARAM_vals.md_max; % Number of debris layers

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
PSANr=reshape(PSAN/100,num_cell,1);
PCLAr=reshape(PCLA/100,num_cell,1);
PORGr=reshape(PORG/100,num_cell,1);

% WHAT IS THIS?
Zs_OUT=800*ones(num_cell,1);

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

%Note slight change to VPAR/POI in INI_COND, using standard not curly
%brackets
if restart.id ~=1
out = [Directories.save 'Initial/INITIAL_CONDITIONS_' site_name '.mat'];
INIT_COND_v6(num_cell,m_cell,n_cell,...
   cc_max,ms_max,md_max,...
   MASKn,GLH,Ca,SNOWD,SNOWALB,out, ...
   VPAR, TT_par, idxCode, Slo_top);
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

%--------------------------------------------------------------------------
%% Restart condition here restart the simulation from a specific month
%--------------------------------------------------------------------------
if restart.id == 1  
    load([Directories.save 'Store/Final_reached_step_' site_name '.mat'])
end

% Label for the creation of outputs
%--------------------------------------------------------------------------
output_creation = 0; 

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

Meteo_data = load(Directories.forc_meteo,meteo_name);
Meteo_data = Meteo_data.(meteo_name);
%** site specific
idForc = isbetween(Meteo_data.DateTime,Meteo_data.DateTime(1),"2024-09-30 23:00"); %This is to prevent Pr issues where NaNs
%**
Forcing_Date = Meteo_data.DateTime(idForc);
Ta = Meteo_data.Ta(idForc); %Air temperature, degrees C
U = Meteo_data.RH(idForc); %Relative humidity, %
Ws = Meteo_data.u(idForc); %Wind speed, ms-1
Pr = Meteo_data.Pr_CorRH(idForc); %Precipitation, mm -  with undercatch and RH correction
SWin = Meteo_data.SWinCor(idForc); %Incoming shortwave, Wm-2  SWinCor with correction for maxSWin using 20260202
SWout = Meteo_data.SWout(idForc); %Outgoing shortwave, Wm-2
Nin = Meteo_data.LWin(idForc); %Incoming longwave, Wm-2
LWout = Meteo_data.LWout(idForc); %Outgoing longwave, Wm-2

%Lapse rates
hm_Ta_lapse = load(ta_lapse_file); %These are the Ta lapse rates from TOPKAPI, per hour and month for main station 144, catchment 406
hm_Ta_lapse = hm_Ta_lapse.(ta_lapse_name); %Pull out of structure
%** site specific
hm_Ta_lapse.Properties.VariableNames{5}='Ta_lapse'; 
%**
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

% Solar variables (just for Lmax_day)
%--------------------------------------------------------------------------
L_day=zeros(length(Datam),1);
for j=2:24:length(Datam)
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),IniCond.DeltaGMT,Lon,Lat,t_bef,t_aft);
end

Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')

%=========================================================================
%% Iterating on time
%==========================================================================
fts = 2; %First time step
for t=fts:N_time_step   
    
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
    
    %% Matrix for storing - Initializing outputs
    %==========================================================================
    % Made by month
    % Only for the first day of the month
    % Because of this, the modelling must start on day 1 for now
    %==========================================================================
    
    if output_creation == 0
    disp("Creating matrices for storing")
    
    %Choose in here the outputs to save
    Initialising_Outputs; %Run initilising outputs script
    
    % Series
    QpointC_series         = single(zeros(zdim, length(Xout)));
    
    % Changing the label to not create it again at the next hour
    %--------------------------------------------------------------------------
    output_creation = 1;
    end


    %% DISTRIBUTED FORCING
    %======================================================================
    % Forcing from stations and distributed across the catchment.
    % Create spatial input data on first day and hour of month 
    %======================================================================      

    %So create new data every month
    if ~exist('Ta_P', 'var') || str2num(yy) ~= year_loaded || str2num(mth) ~= month_loaded

    disp(['New forcing loaded for period: ' char(num2str(yy)) '-' char(num2str(mth)) ])    

    %% Pull out meteo data for the month for each variable  (before spatial distribution)
    id_date_for_month = ismember(Forcing_Date,date_forMonth); %So this is the index of the forcing data for this month
    
    % So for weather station point
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

    %=============================================
    %% Disribute meteorological variables
    % ===========================================
    %Then create grids for every cell so time in rows, space in columns
    %Preallocate where necessary
    size_time_in_month = sum(id_date_for_month);
    Ta_S = NaN(size_time_in_month,num_cell);
    Pr_S = NaN(size_time_in_month,num_cell);
    N_S = NaN(size_time_in_month,num_cell);
    Tdew_S = NaN(size_time_in_month,num_cell);
    SAB1_S = NaN(size_time_in_month,num_cell);
    SAB2_S = NaN(size_time_in_month,num_cell);
    PARB_S = NaN(size_time_in_month,num_cell);
    Pr_ratio = NaN(size_time_in_month,num_cell);

    for tinm = 1:size_time_in_month %So run over each hour
        m_temp = month(date_forMonth(tinm));
        h_temp = hour(date_forMonth(tinm));
        mh_id = hm_Ta_lapse.hour==h_temp & hm_Ta_lapse.month==m_temp;
        Ta_lapse = -hm_Ta_lapse.Ta_lapse(mh_id); %So this is the lapse rate per hour and month - note converting to negative rate
        Tdew_lapse = Ta_lapse; %Its ok to use the same lapse rate
        %Air temperature, degrees C
        Ta_S(tinm,:) = Ta_P(tinm) + (Ta_lapse.*(DTMn-Meteo_el)); %Lapse rates should be negative
        %Precipitation, mm - This uses the WRF ratio v el relationship
        Pr_ratio(tinm,:) = (Pr_lapse(m_temp).month.p1.*DTMn) + Pr_lapse(m_temp).month.p2;  %val(x) = p1*x + p2
        Pr_S(tinm,:) = Pr_ratio(tinm,:).*Pr_P(tinm); %Apply ratio
        %Longwave radiation, W m-2
        N_S(tinm,:) = N_P(tinm) + (Ldown_lapse.*(DTMn-Meteo_el));
        %Dew point temperature, degrees C
        Tdew_S(tinm,:) = Tdew_P(tinm) + (Tdew_lapse.*(DTMn-Meteo_el)); %Lapse rates should be negative
        %Direct radiation with elevation
        if SAB1_P(tinm)>0 %Only if radiation values
            SAB1_S(tinm,:) = SAB1_P(tinm) + (SAB1_lapse.*(DTMn-Meteo_el)); %Radiation should increase with elevation
        else
            SAB1_S(tinm,:) = 0;
        end
        if SAB2_P(tinm)>0
            SAB2_S(tinm,:) = SAB2_P(tinm) + (SAB2_lapse.*(DTMn-Meteo_el)); %Radiation should increase with elevation
        else
            SAB2_S(tinm,:) = 0;
        end
        if PARB_P(tinm)>0
        PARB_S(tinm,:) = PARB_P(tinm) + (PARB_lapse.*(DTMn-Meteo_el)); %Radiation should increase with elevation
        else 
            PARB_S(tinm,:) =0;
        end
    end %Going over time steps within month

    %Apply Tmod to whole month - note this is based on initial ice thickness
    idcli = GLHn>0 & DEB_MAPn == 0; %Clean ice only
    idcli_r = repmat(idcli',size_time_in_month,1);
    Ta_S(idcli_r) = Ta_S(idcli_r) - TmodB; %So remove bias first (over clean ice)
    idTa = Ta_S>0;  %Apply multiplier when >0
    idapply = idTa & idcli_r; %So where Temp threshold exceeded and over clean ice
    Ta_S(idapply) = Ta_S(idapply).*TmodM; %Apply Tmod multiplier

    %Air pressure (calculated directly from elevation)
    Pre_temp = (101325*((1-(0.0065.*DTMn./288.15)).^(9.18*0.0289644./(8.31447*0.0065))))./100; %/100 converts to mbar
    Pre_S = repmat(Pre_temp',size_time_in_month,1); %Copy to all time steps

    %Use the lapsed dew point temprature to derive relative humidity
    c=237.3; b=17.27;
    U_S = 100*exp((c*b.*(Tdew_S - Ta_S))./((c+Ta_S).*(c+Tdew_S)));
    clear c b 
    U_S = U_S./100; %Turn into fraction
    U_S(U_S>1) = 1; %Cannot be more than 100%

    %Do we apply a wind speed lapse rate?
    %At the moment turn into S vector grids
    Ws_S = repmat(Ws_P,1,num_cell); %So the same in all grid cells
    %Diffuse radiation can be the same everywhere
    SAD1_S = repmat(SAD1_P',1,num_cell); %So the same in all grid cells
    SAD2_S = repmat(SAD2_P',1,num_cell); %So the same in all grid cells
    PARD_S = repmat(PARD_P',1,num_cell); %So the same in all grid cells

    % Vapor pressure - calculate based on distributed Ta and U
    % esat/ea/Ds/Tdew grids of month of time x cell
    a=17.27; b=237.3;
    esat_S=611.*exp(a.*Ta_S./(b+Ta_S)); %Vapour pressure at saturation (Pa)
    ea_S=U_S.*esat_S;                 %Vapour pressure (Pa)
    Ds_S= esat_S - ea_S;              %Vapor Pressure Deficit (Pa)
    Ds_S(Ds_S<0)=0; 
    %Don't recalculate Tdew - its above
    clear a b xr;
    
%     % Store year and month loaded
%     %----------------------------------------------------------------------
    year_loaded = str2num(yy);
    month_loaded = str2num(mth);
% 
     end

    %% All now on single timestep
    % Finding the row in the forcing of the modeling date      
    t_forc = find(Date(t-1) == date_forMonth); %Pulls out the index for the hour - note this means that it actually starts on t=1
    
    %Export timestep (a row for all the cells)
    Ta_St = Ta_S(t_forc,:);
    Pr_St = Pr_S(t_forc,:);
    N_St = N_S(t_forc,:);
    Pre_St = Pre_S(t_forc,:);
    Ws_St = Ws_S(t_forc,:);
    U_St = U_S(t_forc,:);
    ea_St = ea_S(t_forc,:);
    Ds_St = Ds_S(t_forc,:);
    Tdew_St = Tdew_S(t_forc,:);
    SAD1_St = SAD1_S(t_forc,:);
    SAD2_St = SAD2_S(t_forc,:);
    SAB1_St = SAB1_S(t_forc,:);
    SAB2_St = SAB2_S(t_forc,:);
    PARB_St = PARB_S(t_forc,:);
    PARD_St = PARD_S(t_forc,:);

    %% Sort Ta for day before
    if t == fts
        Ta_Spdind = Ta_St;
    elseif t > fts && t < fts+24
        Ta_Spdind = [Ta_Spdind;Ta_St];
    elseif t >= fts+24
        Ta_Spdind = [Ta_Spdind(2:24,:);Ta_St];
    end

    %% 

    % Other parameters
    %----------------------------------------------------------------------
    %t_bef=1; t_aft=0; % otherwise problems when loading SWPART

    % Reshape - This because of the problem in the Hydrological module
    %----------------------------------------------------------------------
    Slo_top2 = reshape(Slo_top,num_cell,1); % Creating aux variable for the hydrological module
    aTop = reshape(aTop,num_cell,1);
    %Slo_top = reshape( Slo_top,num_cell,1);

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
    SAD1_St = SAD1_St.*SvFn + Ctn.*rho_g.*(SAD1_St./sin(h_S).*cos_fstn + (1-SvFn).*SAD1_St);
    SAD2_St = SAD2_St.*SvFn + Ctn.*rho_g.*(SAD2_St./sin(h_S).*cos_fstn + (1-SvFn).*SAD2_St);
    PARD_St = PARD_St.*SvFn + Ctn.*rho_g.*(PARD_St./sin(h_S).*cos_fstn + (1-SvFn).*PARD_St);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SAB1_St =(SAB1_St./sin(h_S)).*cos_fstn.*ShFn;
    SAB2_St =(SAB2_St./sin(h_S)).*cos_fstn.*ShFn;
    PARB_St = (PARB_St./sin(h_S)).*cos_fstn.*ShFn;        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    SAB1_St(SAB1_St<0)=0;
    SAB2_St(SAB2_St<0)=0;
    PARB_St(PARB_St<0)=0;
    PARD_St(PARD_St<0)=0;
    SAD1_St(SAD1_St<0)=0;
    SAD2_St(SAD2_St<0)=0;
    SAB1_St(isnan(SAB1_St)) = 0;
    SAB2_St(isnan(SAB2_St)) = 0;
    SAD1_St(isnan(SAD1_St)) = 0;
    SAD2_St(isnan(SAD2_St)) = 0;
    PARB_St(isnan(PARB_St)) = 0;
    PARD_St(isnan(PARD_St)) = 0;
       
    %% Carbon
    Ca_S = Ca(t)*MASKn; %Carbon is already set on the correct time series
    IrD_S =  MASKn*0; %Setting as 0
    Salt_S =  MASKn*0; %Setting as 0
    %N_S = N(t)*MASKn;

    %% Swapping other grids
    Afirnn = reshape(Afirn,num_cell,1);
    SOIL_THn = reshape(SOIL_TH,num_cell,1);

    %{
    forcing.t2m = t2m;
    forcing.d2m = d2m;
    forcing.tp = tp;
    forcing.ws10 = ws10;
    forcing.sp = sp; 

    check_var(forcing, ...
    ["t2m" ... % Temperature
    "d2m" ...  % Dew Point temperature
    "tp" ...   % Precipitation
    "ssrd" ... % Downward short wave radiation
    "strd" ... % Downdward Long wave radiation
    "ws10" ... % Wind speed
    "sp" ...   % Air pressure
    "es" ...   % saturation vapor pressure
    "ea" ...   % actual vapor pressure
    "RH" ...   % Relative humidity
    "SAD1" ... % SAD1
    "SAD2" ... % SAD2
    "SAB1" ... % SAB1
    "SAB2" ... % SAB2
    "PARB" ... % PARB
    "PARD" ... % PARD
    "N" ...    % Cloudiness
    ],Point)
    %}
    
    %% SPATIAL INITIALIZATION VECTOR PREDEFINING
    %======================================================================
   
    if t == 2
        % General vegetation/hydrology
        %------------------------------------------------------------------
        alp_soil=	     alp_soiltm1;
        b_soil=          b_soiltm1;
        Bam=             Bamtm1;
        Bem=             Bemtm1;
        BLit=            BLittm1;
        Ccrown_t=	     Ccrown_t_tm1;
        Cice=            Cicetm1;
        Cicew=           Cicewtm1;
        CK1=             CK1tm1;
        Csno=            Csnotm1;
        Csnow=           Csnowtm1;
        dQ_S=            dQ_Stm1;
        DQ_S=            DQ_Stm1;
        dQVEG=           dQVEGtm1;
        DT_S=            DT_Stm1;
        dw_SNO=          dw_SNOtm1;
        e_sno=           e_snotm1;
        EG=              EGtm1;
        EICE=            EICEtm1;
        EIn_rock=	     EIn_rocktm1;
        EIn_urb=	     EIn_urbtm1;
        EK=              EKtm1;
        ELitter=	     ELittertm1;
        er=              ertm1;
        ESN_In=          ESN_Intm1;
        SSN_In =         SSN_Intm1;
        ESN=             ESNtm1; 
        SSN=             SSNtm1;  
        EWAT=            EWATtm1;
        FROCK=           FROCKtm1;
        f=               ftm1;
        Gfin=            Gfintm1;
        G=               Gtm1;
        H=               Htm1;
        HV=              HVtm1;
        ICE_D=           ICE_Dtm1;
        ICE=             ICEtm1;
        Imelt=           Imelttm1;
        In_H=            In_Htm1;
        In_Litter=	     In_Littertm1;
        In_L=            In_Ltm1;
        In_rock=	     In_rocktm1;
        In_SWE=          In_SWEtm1;
        In_urb=          In_urbtm1;
        IP_wc=           IP_wctm1;
        Lk_rock=	     Lk_rocktm1;
        Lk_wat=          Lk_wattm1;
        Lk=              Lktm1;
        Lpho=            Lphotm1;
        NavlI=           NavlItm1;
        NDVI=            NDVItm1;
        NIce=            NIcetm1;
        NIn_SWE=	     NIn_SWEtm1;
        OF=              OFtm1;
        Oice=            Oicetm1;
        OS=              OStm1;
        O=               Otm1;
        POT=             POTtm1;
        Pr_liq=          Pr_liqtm1;
        Pr_sno=          Pr_snotm1;
        Q_channel=	     Q_channel;
        Q_exit=          Q_exit;
        q_runon=         q_runon;
        QE=              QEtm1;
        QEV=             QEVtm1;
        Qfm=             Qfmtm1;
        Qi_in=           Qi_in;
        Qi_in_Ro=        Qi_out;
        Qi_out_Rout=     Qi_out_Rout;
        Qi_out=          Qi_outtm1;
        Qsub_exit=	     Qsub_exit;
        Qv=              Qvtm1;
        r_litter=	     r_littertm1;
        r_soil=          r_soiltm1;
        ra=              ratm1;
        Rd=              Rdtm1;
        Rh=              Rhtm1;
        Rn=              Rntm1;
        ros=             rostm1;
        SE_rock=	     SE_rocktm1;
        SE_urb=          SE_urbtm1;
        %Slo_head=	     Slo_head;
        Smelt=           Smelttm1;
        SND=             SNDtm1;
        snow_albedo=     snow_albedotm1;
        soil_albedo=     soil_albedotm1;
        SP_wc=           SP_wctm1;
        surface_albedo=  surface_albedotm1;
        SWE=             SWEtm1;
        SWE_avalanched=  SWE_avalanchedtm1;
        t_sls=           t_slstm1;
        tau_sno=	     tau_snotm1;
        Tdamp=           Tdamptm1;
        Tdeb=            Tdebtm1;
        Tdp=             Tdptm1;
        Tdp_snow =       Tdp_snowtm1;
        Tice=            Ticetm1;
        Tstm0=           Tstm0;
        Ts=              Tstm1;
        Ts_under =       Ts_undertm1;
        TsVEG=           TsVEGtm1;
        U_SWE=           U_SWEtm1;
        Vice=            Vicetm1;
        V=               Vtm1;
        WAT=             WATtm1;
        WIS=             WIStm1;
        WR_IP=           WR_IPtm1;
        WR_SP=           WR_SPtm1;
        Ws_under=	     Ws_undertm1;
        ZWT=             ZWTtm1;
        
        % Specifications for high & low vegetation
        %------------------------------------------------------------------
        AgeDL_H=	AgeDL_Htm1;
        AgeDL_L=	AgeDL_Ltm1;
        AgeL_H=     AgeL_Htm1;
        AgeL_L=     AgeL_Ltm1;
        AgePl_H=	AgePl_Htm1;
        AgePl_L=	AgePl_Ltm1;
        An_H=       An_Htm1;
        An_L=       An_Ltm1;
        ANPP_H=     ANPP_Htm1;
        ANPP_L=     ANPP_Ltm1;
        B_H=        B_Htm1;
        B_L=        B_Ltm1;
        BA_H=       BA_Htm1;
        BA_L=       BA_Ltm1;
        Bfac_dayH=	Bfac_dayHtm1;
        Bfac_dayL=	Bfac_dayLtm1;
        Bfac_weekH=	Bfac_weekHtm1;
        Bfac_weekL=	Bfac_weekLtm1;
        Ci_shdH=	Citm1_shdH;
        Ci_shdL=	Citm1_shdL;
        Ci_sunH=	Citm1_sunH;
        Ci_sunL=	Citm1_sunL;
        dflo_H=     dflo_Htm1;
        dflo_L=     dflo_Ltm1;
        Dr_H=       Dr_Htm1;
        Dr_L=       Dr_Ltm1;
        e_rel_H=	e_rel_Htm1;
        e_rel_L=	e_rel_Ltm1;
        e_relN_H=	e_relN_Htm1;
        e_relN_L=	e_relN_Ltm1;
        EIn_H=      EIn_Htm1;
        EIn_L=      EIn_Ltm1;
        fapar_H=	fapar_Htm1;
        fapar_L=	fapar_Ltm1;
        FNC_H=      FNC_Htm1;
        FNC_L=      FNC_Ltm1;
        gsr_H=      gsr_Htm1;
        gsr_L=      gsr_Ltm1;
        hc_H=       hc_Htm1;
        hc_L=       hc_Ltm1;
        In_H=       In_Htm1;
        In_L=       In_Ltm1;
        ISOIL_H=	ISOIL_Htm1;
        ISOIL_L=	ISOIL_Ltm1;
        Jsx_H=      Jsx_Htm1;
        Jsx_L=      Jsx_Ltm1;
        Jxl_H=      Jxl_Htm1;
        Jxl_L=      Jxl_Ltm1;
        Kleaf_H=	Kleaf_Htm1;
        Kleaf_L=	Kleaf_Ltm1;
        Kreserve_H=	Kreserve_Htm1;
        Kreserve_L=	Kreserve_Ltm1;
        Kuptake_H=	Kuptake_Htm1;
        Kuptake_L=	Kuptake_Ltm1;
        Kx_H=       Kx_Htm1;
        Kx_L=       Kx_Ltm1;
        LAI_H=      LAI_Htm1;
        LAI_L=      LAI_Ltm1;
        LAIdead_H=	LAIdead_Htm1;
        LAIdead_L=	LAIdead_Ltm1;
        ManIH=      ManIHtm1;
        ManIL=      ManILtm1;
        NBLeaf_H=	NBLeaf_Htm1;
        NBLeaf_L=	NBLeaf_Ltm1;
        NBLI_H=     NBLI_Htm1;
        NBLI_L=     NBLI_Ltm1;
        NPP_H=      NPP_Htm1;
        NPP_L=      NPP_Ltm1;
        NPPI_H=     NPPI_Htm1;
        NPPI_L=     NPPI_Ltm1;
        Nreserve_H=	Nreserve_Htm1;
        Nreserve_L=	Nreserve_Ltm1;
        NuLit_H=	NuLit_Htm1;
        NuLit_L=	NuLit_Ltm1;
        NupI_H=     NupI_Htm1;
        NupI_L=     NupI_Ltm1;
        Nuptake_H=	Nuptake_Htm1;
        Nuptake_L=	Nuptake_Ltm1;
        OH=         OHtm1;
        OL=         OLtm1;
        PARI_H=     PARI_Htm1;
        PARI_L=     PARI_Ltm1;
        PHE_S_H=	PHE_S_Htm1;
        PHE_S_L=	PHE_S_Ltm1;
        Preserve_H=	Preserve_Htm1;
        Preserve_L=	Preserve_Ltm1;
        Psi_l_H=	Psi_l_Htm1;
        Psi_l_L=	Psi_l_Ltm1;
        Psi_s_H=	Psi_s_Htm1;
        Psi_s_L=	Psi_s_Ltm1;
        Psi_x_H=	Psi_x_Htm1;
        Psi_x_L=	Psi_x_Ltm1;
        Puptake_H=	Puptake_Htm1;
        Puptake_L=	Puptake_Ltm1;
        RA_H=       RA_Htm1;
        RA_L=       RA_Ltm1;
        rap_H=      rap_Htm1;
        rap_L=      rap_Ltm1;
        RB_H=       RB_Htm1;
        rb_H=       rb_Htm1;
        RB_L=       RB_Ltm1;
        rb_L=       rb_Ltm1;
        Rdark_H=	Rdark_Htm1;
        Rdark_L=	Rdark_Ltm1;
        Rexmy_H=	Rexmy_Htm1;
        Rexmy_L=	Rexmy_Ltm1;
        Rg_H=       Rg_Htm1;
        Rg_L=       Rg_Ltm1;
        rKc_H=      rKc_Htm1;
        rKc_L=      rKc_Ltm1;
        Rmc_H=      Rmc_Htm1;
        Rmc_L=      Rmc_Ltm1;
        Rmr_H=      Rmr_Htm1;
        Rmr_L=      Rmr_Ltm1;
        Rms_H=      Rms_Htm1;
        Rms_L=      Rms_Ltm1;
        rNc_H=      rNc_Htm1;
        rNc_L=      rNc_Ltm1;
        rPc_H=      rPc_Htm1;
        rPc_L=      rPc_Ltm1;
        Rrootl_H=	Rrootl_Htm1;
        Rrootl_L=	Rrootl_Ltm1;
        rs_shdH=	rs_shdHtm1;
        rs_shdL=	rs_shdLtm1;
        rs_sunH=	rs_sunHtm1;
        rs_sunL=	rs_sunLtm1;
        SAI_H=      SAI_Htm1;
        SAI_L=      SAI_Ltm1;
        Sfr_H=      Sfr_Htm1;
        Sfr_L=      Sfr_Ltm1;
        SIF_H=      SIF_Htm1;
        SIF_L=      SIF_Ltm1;
        Slf_H=      Slf_Htm1;
        Slf_L=      Slf_Ltm1;
        Sll_H=      Sll_Htm1;
        Sll_L=      Sll_Ltm1;
        Sr_H=       Sr_Htm1;
        Sr_L=       Sr_Ltm1;
        SupK_H=     SupK_Htm1;
        SupK_L=     SupK_Ltm1;
        SupN_H=     SupN_Htm1;
        SupN_L=     SupN_Ltm1;
        SupP_H=     SupP_Htm1;
        SupP_L=     SupP_Ltm1;
        Swm_H=      Swm_Htm1;
        Swm_L=      Swm_Ltm1;
        T_H=        T_Htm1;
        T_L=        T_Ltm1;
        TBio_H=     TBio_Htm1;
        TBio_L=     TBio_Ltm1;
        Tden_H=     Tden_Htm1;
        Tden_L=     Tden_Ltm1;
        Tdp_H=      Tdp_Htm1;
        Tdp_L=      Tdp_Ltm1;
        TdpI_H=     TdpI_Htm1;
        TdpI_L=     TdpI_Ltm1;
        TexC_H=     TexC_Htm1;
        TexC_L=     TexC_Ltm1;
        TexK_H=     TexK_Htm1;
        TexK_L=     TexK_Ltm1;
        TexN_H=     TexN_Htm1;
        TexN_L=     TexN_Ltm1;
        TexP_H=     TexP_Htm1;
        TexP_L=     TexP_Ltm1;
        TNIT_H=     TNIT_Htm1;
        TNIT_L=     TNIT_Ltm1;
        TPHO_H=     TPHO_Htm1;
        TPHO_L=     TPHO_Ltm1;
        TPOT_H=     TPOT_Htm1;
        TPOT_L=     TPOT_Ltm1;
        Vl_H=       Vl_Htm1;
        Vl_L=       Vl_Ltm1;
        Vx_H=       Vx_Htm1;
        Vx_L=       Vx_Ltm1;
    end
  

    %% LOOP
    %======================================================================
    % LOOP OVER CELLS
    %======================================================================
    P = (Xout - 1) * 133 + Yout;

    %Follow=MASK; %Debugging
   for ij= 1:num_cell % this is a parfor
        % Good practice to use a simple for loop for debugging/testing
        % ij is the index to go pixel by pixel through the mask
        % ij=1:num_cell
        %disp(strcat('in the loop', ij))
        
        %% Debugging for a for loop
        %Follow(ij) = 2222; %22307 - 22308 (row 107 - col 151)
        %disp('bye')
        %{
        if ismember(ij, [num_cell/8, num_cell/4, num_cell/2, 3*num_cell/4])   
        disp(ij)
        end        
        %Ta,Ts,Pre,zatm,disp_h,zom,zoh,Ws,ea
        %} 
        % =================================================================
        
        if MASKn(ij)== 1
            %disp(['Cell: ' char(num2str(ij)) ', Veg Type: ' char(num2str(ksv(ij)))])
            Elev=DTMn(ij);
            %[i,j] = ind2sub([m_cell,n_cell],ij);

            % BOUNDARY CONDITION  
            % INTRODUCED SOIL AND VEG. for ij
            %--------------------------------------------------------------
            [aR,             Zs,             EvL_Zs,       Inf_Zs,     Bio_Zs,      Zinf, ...
             RfH_Zs,         RfL_Zs,         dz,           Ks_Zs,      Dz,          ms, ...
             Kbot,           Krock,          zatm,         Color_Class,  OM_H,       OM_L,        PFT_opt_H, ...
             PFT_opt_L,      d_leaf_H,       d_leaf_L,     SPAR,       Phy,         Soil_Param, ...
             Interc_Param,   SnowIce_Param,  VegH_Param,   VegL_Param, fpr,         VegH_Param_Dyn, ...
             VegL_Param_Dyn, Stoich_H,       aSE_H,        Stoich_L,   aSE_L,       fab_H,...
             fbe_H,          fab_L,          fbe_L,        ZR95_H,     ZR95_L,      In_max_urb,...
             In_max_rock,    K_usle,         Urb_Par,      Deb_Par,    Zs_deb,      Sllit,...
             Kct,            ExEM,           ParEx_H,      Mpar_H,     ParEx_L,     Mpar_L,] = ....
                                    PARAMETERS_SOIL_DEV( ...
             ksv(ij),        PSANr(ij),       PCLAr(ij),     PORGr(ij),   DEB_MAPn(ij),  md_max,...
             Afirnn(ij),      SOIL_THn(ij),    VPAR,        TT_par,        OPT_PARAM_vals);

            %Extract C values directly from grids
            Ccrown = VPAR.Ccrownr(ij);
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
    
                BLit(ij,:)= 0.0 ; % %% %%[kg DM / m2]
            end        

            %% HYDROLOGY MODULE
            %==============================================================
            % FUNTION: HYDROLOGY_MODULE_PAR
            %==============================================================                                                 
    %try
 %{   
mm = {Vtm1,        Oicetm1,       aR,           Zs,                  EvL_Zs,        ...
     Inf_Zs, ...
     Zinf,         RfH_Zs,        RfL_Zs,       dz,                  Dz, ...
     ms,           Kbot,          Pr_S,         Ta_S,                Ds_S, ...
     Ws_S,         zatm,          Tstm1,        dt,                  dth, ...
     ea_S,         N_S,           Pre_S,        Tstm0,               LAI_H, ...
     SAI_H,        LAI_L,         SAI_L,        LAIdead_H,           LAIdead_L,...
     Rrootl_H,     Rrootl_L,      BLit,         Sllit,               Kct,...
     Datam_S,      IniCond.DeltaGMT,      Lon,          Lat,                 t_bef, ...
     t_aft,        Ccrown,        Cbare,        Crock,               Curb, ...
     Cwat,         SAB1_S,        SAB2_S,       SAD1_S,              SAD2_S, ...
     PARB_S,       PARD_S,        SvF,          SNDtm1,              snow_albedotm1, ...
     Color_Class,  OM_H,          OM_L,         PFT_opt_H,           PFT_opt_L, ...
     hc_H,         hc_L,          d_leaf_H,     d_leaf_L,            Soil_Param, ...
     Interc_Param, SnowIce_Param, VegH_Param,   VegL_Param,          Ca_S, ...
     Oa,           Citm1_sunH,    Citm1_shdH,   Citm1_sunL,          Citm1_shdL,...
     e_rel_H,      e_relN_H,      e_rel_L,      e_relN_L,            e_snotm1, ...
     In_Htm1,      In_Ltm1,       In_Littertm1, In_urbtm1,           In_rocktm1(ij), ...
     SWEtm1,       In_SWEtm1,     Tdebtm1,      Ticetm1,             Tdptm1(ij,:), ...
     Tdp_snowtm1,  Tdamptm1,      Ts_undertm1,  WATtm1,              ICEtm1(ij), ...
     IP_wctm1,     ICE_Dtm1,      Cicewtm1,     Vx_Htm1,             Vl_Htm1(ij,:), ...
     Vx_Ltm1,      Vl_Ltm1,       Psi_x_Htm1,   Psi_l_Htm1,          Psi_x_Ltm1(ij,:), ...
     Psi_l_Ltm1,   ZR95_H,        ZR95_L,       FROCKtm1,            Krock,...
     Urb_Par, ...
     Deb_Par,      Zs_deb,        Tdew_S,       t_slstm1,            rostm1, ...
     SP_wctm1,     fpr,           IrD_S,        In_max_urb,          In_max_rock, ...
     K_usle,       tau_snotm1,    Ta_day,       Slo_top2,            Slo_head, ...
     Asur,         Ared,          aTop,         EKtm1,               q_runon, ...
     Qi_in,        Ws_undertm1,   Pr_sno_t,     pow_dis,             a_dis, ...
     Salt_S,       SPAR,SNn,      min_SPD,  OPT_VegSnow,         OPT_SoilTemp, ...
     OPT_PlantHydr,Opt_CR,        Opt_ST,       Opt_ST2,             OPT_SM, ...
     OPT_STh,      OPT_FR_SOIL,   OPT_PH,       parameterize_phase,  hSTL, ...
     OPT_Albsno};

sizes = cellfun(@size, mm, 'UniformOutput', false);
stringCell = cellfun(@(x) sprintf('[%.0f,%.0f]', x(1), x(2)), sizes, 'UniformOutput', false);

T = table('Size', [31,5], ...
    'VariableTypes', {'string', 'string', 'string', 'string', 'string'});

%i = 5;
z = 1;
y = 1;
        for i = 1:length(sizes)
        T(z,y) = stringCell(i);
        y = y+1;
            if mod(i, 5) == 0ij
            z = z+1;
            y=1;
            end
        
        end
 
    if t == 2    
    writetable(T, [Directories.save ' output_data.csv']);
    end
    
 %}

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
     Dz,                    ms,                  Kbot,             Pr_St(ij),           Ta_St(ij), ...
     Ds_St(ij),             Ws_St(ij),           zatm,             Tstm1(ij),           dt, ...
     dth,                   ea_St(ij),           N_St(ij),         Pre_St(ij),          Tstm0(ij), ...
     LAI_H(ij,:),           SAI_H(ij,:),         LAI_L(ij,:),      SAI_L(ij,:),         LAIdead_H(ij,:), ...
     LAIdead_L(ij,:),       Rrootl_H(ij,:),      Rrootl_L(ij,:),   BLit(ij,:),          Sllit,   ...
     Kct,                   Datam_S,             IniCond.DeltaGMT, Lon,                 Lat,     ...
     t_bef,                 t_aft,               Ccrown,           Cbare,               Crock,  ...
     Curb,                  Cwat,                SAB1_St(ij),      SAB2_St(ij),         SAD1_St(ij), ...
     SAD2_St(ij),           PARB_St(ij),         PARD_St(ij),      SvFn(ij),            SNDtm1(ij),  ...
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
     Krock,                 Urb_Par,             Deb_Par,          Zs_deb,              Tdew_St(ij), ...
     t_slstm1(ij),          rostm1(ij),          SP_wctm1(ij),     fpr,                 IrD_S(ij), ...
     In_max_urb,            In_max_rock,         K_usle,           tau_snotm1(ij),      Ta_Spdind(:,ij), ...
     Slo_top2(ij),          Slo_head(ij,:),      Asur(ij),         Ared(ij),            aTop(ij), ...
     EKtm1(ij),             q_runon(ij),         Qi_in(ij,:),      Ws_undertm1(ij),     Pr_sno_t(ij,:), ...
     pow_dis,               a_dis,               Salt_S(ij),       SPAR,                SNn(ij), ...
     min_SPD,           OPT_VegSnow,         OPT_SoilTemp,     OPT_PlantHydr,       Opt_CR, ...
     Opt_ST,                Opt_ST2,             OPT_SM,           OPT_STh,             OPT_FR_SOIL, ...
     OPT_PH,                parameterize_phase,  hSTL,             OPT_Albsno);
    %catch ME
    %disp(['Error in HYDROLOGY_MODULE_PAR in ij:' char(num2str(ij))])
    
    % Get the full error identifier (type/ID)
    %errorID = ME.identifier;
    %disp(['Error Type (ID): ', errorID]);

    %end
    
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

    if t==2 
        snow_albedo_out = snow_albedotm1;
        surface_albedo_out = surface_albedotm1;
    elseif Datam_S(4)==12
    snow_albedo_out = snow_albedo;
    surface_albedo_out = surface_albedotm1;
    surface_albedo_out(ICE_D>0)=0.28; %THESE ARE NOT CORRECT?
    surface_albedo_out(DEB_MAP>0)=0.13; %THESE ARE NOT CORRECT?
    surface_albedo_out(SND>0)=snow_albedo(SND>0);
    end

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
        PAR_t(:,24)=            PARB_St + PARD_St;
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
        Ta_t(:,24)=             Ta_St;
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
        PAR_t(:,24)=            PARB_St + PARD_St;
        Pr_sno_t(:,t)=          Pr_sno;
        Psi_l_H_t(:,:,t)=       Psi_l_H;
        Psi_l_L_t(:,:,t)=       Psi_l_L;
        Psi_x_H_t(:,:,t)=       Psi_x_H;
        Psi_x_L_t(:,:,t)=       Psi_x_L;
        Rdark_H_t(:,:,t)=       Rdark_H;
        Rdark_L_t(:,:,t)=       Rdark_L;
        Ta_t(:,t)=              Ta_St;
        Tdp_H_t(:,:,t)=         Tdp_H;
        Tdp_L_t(:,:,t)=         Tdp_L;
        Tdp_t(:,:,t)=           Tdp;
        V_t(:,:,t)=             V;
    end



    %% OUTPUT WRITING
    %======================================================================
    %
    %======================================================================
       
    % Ccrown_OUT and EVcode are used in the OUTPUT_MANAGER_DIST_LABEL.m
    %--------------------------------------------------------------------------
    % This only for the OUTPUT_MANAGER_DIST_LABEL
    %Ccrown_OUT =[ 1 ; 1; 1 ; 0.9 ; 1; 1; 0];  %% Ccrown fraction for Plant Functional Types (PFT)
    %EVcode = [ 1 2 3 4 5 6 7 8 9 10];             %% code of each PFTs
    %Ccrown_OUT2 = VPAR.Ccrowns;
    %{
    %Mike
    Ccrown_OUT = [ 1   0 ; 0.5  0.5 ; 1   0 ; 1   0];  %% Ccrown fraction for PFT
    EVcode     = [ 1      2           3       4    ];  %% Code of each PFT
    cc_max = size(Ccrown_OUT,2);                       %% Number of vegetation types per cell (max.)
    %}

    % Generating the vector for the calculation of the outputs by crowns
    %----------------------------------------------------------------------
    
%     VPARc.Ccrowns = num2cell(VPAR.Ccrowns); %Turn into cell, as was originally, but mine is a simple table
%     max_len = max(cellfun(@numel, VPARc.Ccrowns)); %So this is the number of crown values per land cover
%     padded_cells = cellfun(@(v) [v, zeros(1, max_len - numel(v))], VPARc.Ccrowns, 'UniformOutput', false);
%     Ccrown_OUT = cell2mat(padded_cells);
%     
%     %I don't know if this is right!! PLEASE CHECK
%     Ccrown_OUTv=zeros(num_cell,1);
%     for c=1:size(Ccrown_OUT,1)
%         v_id = ksv==c; %So this is the land class, we start from 1
%         Ccrown_OUTv(v_id) = Ccrown_OUT(c);
%     end
    %EVcode     = VPAR.Class;

    %Can directly use Ccrown vector of the grid
    Ccrown_OUTv = VPAR.Ccrownr; 
  
    % Assignation to temporal matrix
    %----------------------------------------------------------------------
    
    t_date = datetime(Datam_S(1), Datam_S(2), Datam_S(3), Datam_S(4), 0, 0); %Note Datam_S is set to t-1
    t_assign =  find(t_date  ==  date_forMonth);

    %Transpose the meteo variables
    Pr_St = Pr_St'; Ta_St = Ta_St'; Tdew_St = Tdew_St'; Ws_St = Ws_St';
    U_St = U_St'; N_St = N_St'; ea_St = ea_St'; Pre_St = Pre_St';
    SAD1_St = SAD1_St'; SAD2_St = SAD2_St'; SAB1_St = SAB1_St'; SAB2_St = SAB2_St';

    if t==fts %Sort variable dimensions of chosen variables on first timestep
    %NOTE if you can new variables in Initialing variable check if they are
    %multiplied by Ccrown
    VAR_CROWNS = ["ANPP_H"; "ANPP_L"; "LAI_H"; "LAI_L"; "NPP_H"; "NPP_L"; 
                  "hc_H";   "hc_L";   "EIn_H"; "EIn_L"; "T_H";   "T_L" ;
                  "An_H"; "AN_L"; "B_H"; "B_L"; "Bfac_dayH"; "Bfac_dayL";
                  "RA_H"; "RA_L"; "Rg_H"; "Rg_L"; "SAI_H";"SAI_L"; 
                  "PHE_S_H"; "PHE_S_L"; "In_L"; "In_H"; "Ci_shdH"; "Ci_shdL";
                  "Ci_sunH";'Ci_sunL';"LAIdead_H"; "LAIdead_L"; "Rdark_H";"Rdark_L";
                  "Rmc_H";"Rmc_L";"Rmr_H";"Rmr_L";"Rms_H";"Rms_L"];

    %Check dimensions
    WSinfo = whos; %Extracts all names and sizes
    Allnames = {WSinfo.name}; Allnames=Allnames';
    inVAR = ismember(Allnames,VAR_OUT);
    num_vars = size(VAR_OUT,1);
    VARSinfo = WSinfo(inVAR); %Great this is a structure with just the dimensions of the output variables
    Varsinfonames = {VARSinfo.name}; Varsinfonames=Varsinfonames';
    all_sizes1 = zeros(num_vars,1); %Needs to be three in case of array %First dimension
    all_sizes2 = zeros(num_vars,1); %Needs to be three in case of array %Second dimension
    all_sizes3 = zeros(num_vars,1); %Needs to be three in case of array %Third dimension
    for v=1:num_vars
            all_sizes1(v,1) = VARSinfo(v).size(1); %1st dimension
            all_sizes2(v,1) = VARSinfo(v).size(2); %2nd dimension
            if size(VARSinfo(v).size,3)>1
            all_sizes3(v,1) = VARSinfo(v).size(3); %3rd dimension only if exists
            end
    end
    %Note all sizes at the moment refer to the variable order in VARSinfo,
    %not VAR_OUT so make a table
    VAR_OUT_T = array2table(VAR_OUT,'VariableNames',{'Name'});
    for v=1:num_vars
        idvar = matches(Varsinfonames,VAR_OUT_T.Name(v));
        VAR_OUT_T.Dim1(v) = all_sizes1(idvar);
        VAR_OUT_T.Dim2(v) = all_sizes2(idvar);
        VAR_OUT_T.Dim3(v) = all_sizes3(idvar);
    end
    end
    
    %----------------------------------------------------------------
    %Assign the gridded values to an hourly time in the initial grids
    %-----------------------------------------------------------------

    for k=1:length(VAR_OUT)
        if VAR_OUT(k) == "QpointC" % case of discharge
        QpointC_hourly(t_assign,:) = QpointC';                                  % runoff/discharge    
    
        elseif VAR_OUT(k) == "O"
        O1_hourly(:,:,t_assign)       = reshape(O(:,1), xdim,ydim);            % Soil mositure layer 1
        O2_hourly(:,:,t_assign)       = reshape(O(:,2), xdim,ydim);            % Soil mositure layer 2
        O3_hourly(:,:,t_assign)       = reshape(O(:,3), xdim,ydim);            % Soil mositure layer 3
        O4_hourly(:,:,t_assign)       = reshape(O(:,4), xdim,ydim);            % Soil mositure layer 4
        O5_hourly(:,:,t_assign)       = reshape(O(:,5), xdim,ydim);            % Soil mositure layer 5
        O6_hourly(:,:,t_assign)       = reshape(O(:,6), xdim,ydim);            % Soil mositure layer 6
        O7_hourly(:,:,t_assign)       = reshape(O(:,7), xdim,ydim);            % Soil mositure layer 7
        O8_hourly(:,:,t_assign)       = reshape(O(:,8), xdim,ydim);            % Soil mositure layer 8
        O9_hourly(:,:,t_assign)       = reshape(O(:,9), xdim,ydim);            % Soil mositure layer 9
        O10_hourly(:,:,t_assign)      = reshape(O(:,10), xdim,ydim);           % Soil mositure layer 10
    
        elseif VAR_OUT(k) == "snow_albedo"
        snow_albedo1_hourly(:,:,t_assign)       = reshape(snow_albedo(:,1), xdim,ydim);            % albedo layer 1
        snow_albedo2_hourly(:,:,t_assign)       = reshape(snow_albedo(:,2), xdim,ydim);            % albedo layer 2
        snow_albedo3_hourly(:,:,t_assign)       = reshape(snow_albedo(:,3), xdim,ydim);            % albedo layer 3
        snow_albedo4_hourly(:,:,t_assign)       = reshape(snow_albedo(:,4), xdim,ydim);            % albedo layer 4

        elseif VAR_OUT(k) == "Tdp"
        Tdp1_hourly(:,:,t_assign)       = reshape(Tdp(:,1), xdim,ydim);            % Tdp layer 1
        Tdp2_hourly(:,:,t_assign)       = reshape(Tdp(:,2), xdim,ydim);            % Tdp layer 2
        Tdp3_hourly(:,:,t_assign)       = reshape(Tdp(:,3), xdim,ydim);            % Tdp layer 3
        Tdp4_hourly(:,:,t_assign)       = reshape(Tdp(:,4), xdim,ydim);            % Tdp layer 4
        Tdp5_hourly(:,:,t_assign)       = reshape(Tdp(:,5), xdim,ydim);            % Tdp layer 5
        Tdp6_hourly(:,:,t_assign)       = reshape(Tdp(:,6), xdim,ydim);            % Tdp layer 6
        Tdp7_hourly(:,:,t_assign)       = reshape(Tdp(:,7), xdim,ydim);            % Tdp layer 7
        Tdp8_hourly(:,:,t_assign)       = reshape(Tdp(:,8), xdim,ydim);            % Tdp layer 8
        Tdp9_hourly(:,:,t_assign)       = reshape(Tdp(:,9), xdim,ydim);            % Tdp layer 9
        Tdp10_hourly(:,:,t_assign)      = reshape(Tdp(:,10), xdim,ydim);           % Tdp layer 10

        elseif ismember(VAR_OUT(k), VAR_CROWNS) % Case of variables with crowns
            commandStr = [char(VAR_OUT_HOUR(k)) '(:,:,t_assign) = reshape(sum(' char(VAR_OUT(k)) '.*Ccrown_OUTv,2),xdim,ydim)' ';']; %Other matrices
            eval(commandStr);
    
        elseif VAR_OUT_T.Dim1(v)==(xdim*ydim) && VAR_OUT_T.Dim2(v)==1 %So correct dimensions
                commandStr = [char(VAR_OUT_HOUR(k)) '(:,:,t_assign) = reshape(' char(VAR_OUT(k)) ',xdim,ydim)' ';']; %Other matrices
                eval(commandStr);
        
        else % Other variables with different sizes
                disp(['Variable ' char(VAR_OUT(k)) ' has a different dimension to Ccrown and other dimensions considered'])
        end   
    end

    %------------------------------------
    % Extracting data for points
    %--------------------------------------------------

    %NOTE VEG DATA NOT CORRECTED FOR CCROWN

    % Here pull out from the grids the hourly variables from points, add
    % Date info into table
    %This needs to match the list in Timestep_T!
     if t==fts %Only run on first time step or it will over-write
        vars_to_keep = {'EICE','ESN','SND','SWE','Ta_St','Ws_St','U_St','N_St','ea_St','Rn',...
        'SAD1_St','SAD2_St','SAB1_St','SAB2_St','Pre_St','Pr_St','Pr_sno','Pr_liq', 'snow_albedo','Smelt','Imelt', ...
        'ICE','ICE_D','QE', 'ros', 'NDVI', 'T_H', 'T_L', 'O', 'FROCK', ...
        'LAI_H', 'LAI_L', 'NPP_H', 'NPP_L','Csno','An_L','An_H',...
        'EG','EIn_L','EIn_H','EIn_rock','H','HV','G','In_L','In_H','In_SWE','In_rock',...
        'OF','OH','OL','OS','PHE_S_L','PHE_S_H','Qfm','Qv','RA_L','RA_H','Rg_L','Rg_H','Rh','SAI_L','SAI_H',...
        'Tdp','Tdp_L','Tdp_H','Ts','ZWT','dQ_S','dQVEG','IP_wc'};
        num_vk = size(vars_to_keep,2);
        run_len = size(Date,2);
        num_pts = size(Points,1);
    
        %Due to snow_albedo and O need to work harder to make the table blank
        BC = zeros(run_len,1);
        BC_sa = [BC BC BC BC];
        BC_O = [BC BC BC BC BC BC BC BC BC BC]; %Also for Tdp
        
        for p=1:num_pts
            Points_Data(p).Table = table(BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC_sa,...
                BC,BC,BC,BC,BC,BC,BC,BC,BC,BC_O,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,...
                BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC,BC_O,BC,BC,BC,BC,BC,BC,BC,'VariableNames',vars_to_keep);
            Points_Data(p).Table.Date = Date';
            Points_Data(p).Table = movevars(Points_Data(p).Table,'Date','Before',1);
        end 
        %Note you can use te idx from the Points table to get values in the
        %column variables, or use the ROW and COL value to get the grid values.
     end

    %Create temp table for this time step of all variables
    Timestep_T = table(EICE,ESN,SND,SWE,Ta_St,Ws_St,U_St,N_St,ea_St,Rn,...
    SAD1_St,SAD2_St,SAB1_St,SAB2_St,Pre_St,Pr_St,Pr_sno,Pr_liq, snow_albedo,Smelt,Imelt, ...
    ICE,ICE_D,QE, ros, NDVI, T_H, T_L, O, FROCK, ...
    LAI_H, LAI_L, NPP_H, NPP_L,Csno,An_L,An_H,...
    EG,EIn_L,EIn_H,EIn_rock,H,HV,G,In_L,In_H,In_SWE,In_rock,...
    OF,OH,OL,OS,PHE_S_L,PHE_S_H,Qfm,Qv,RA_L,RA_H,Rg_L,Rg_H,Rh,SAI_L,SAI_H,...
    Tdp,Tdp_L,Tdp_H,Ts,ZWT,dQ_S,dQVEG,IP_wc,'VariableNames',vars_to_keep);

    %Extract only the cells necessary
    Timestep_T_Pts = Timestep_T(Points.idx,:); %Each row here corresponds to a point location
   
    for p=1:num_pts
        for v=1:num_vk
        tt_col = matches(Timestep_T_Pts.Properties.VariableNames,vars_to_keep(v));
        pd_col = matches(Points_Data(p).Table.Properties.VariableNames,vars_to_keep(v));
        Points_Data(p).Table(t-1,pd_col) = Timestep_T_Pts(p,tt_col); %So transfer point data to the correct timestep %Note t-1, so in correct time!
        end
    end

    clear Timestep_T Timestep_T_Pts %Clear these tables as they are remade next timestep

    %Old for output manager
% %     % Outputs for vegetation
% %     ANPP_H_spatial(:,:,t_assign)   = reshape(sum(ANPP_H.*Ccrown_OUTv,2),xdim,ydim); 
% %     ANPP_L_spatial(:,:,t_assign)   = reshape(sum(ANPP_L.*Ccrown_OUTv,2),xdim,ydim);
% %     LAI_H_spatial(:,:,t_assign)    = reshape(sum(LAI_H.*Ccrown_OUTv,2),xdim,ydim);
% %     LAI_L_spatial(:,:,t_assign)    = reshape(sum(LAI_L.*Ccrown_OUTv,2),xdim,ydim);
% %     NDVI_spatial(:,:,t_assign)     = reshape(NDVI,xdim,ydim);
% % 
% %     % Output for ET
% %     EG_spatial(:,:,t_assign)       = reshape(EG,xdim,ydim);
% %     EICE_spatial(:,:,t_assign)     = reshape(EICE,xdim,ydim);    
% %     EIn_H_spatial(:,:,t_assign)    = reshape(sum(EIn_H.*Ccrown_OUTv,2),xdim,ydim);
% %     EIn_L_spatial(:,:,t_assign)    = reshape(sum(EIn_L.*Ccrown_OUTv,2),xdim,ydim);
% %     EIn_rock_spatial(:,:,t_assign) = reshape(EIn_rock,xdim,ydim); 
% %     EIn_urb_spatial(:,:,t_assign)  = reshape(EIn_urb,xdim,ydim); 
% %     EWAT_spatial(:,:,t_assign)     = reshape(EWAT,xdim,ydim); 
% %     T_H_spatial(:,:,t_assign)      = reshape(sum(T_H.*Ccrown_OUTv,2),xdim,ydim); 
% %     T_L_spatial(:,:,t_assign)      = reshape(sum(T_L.*Ccrown_OUTv,2),xdim,ydim); 
% %     SSN_spatial(:,:,t_assign)      = reshape(SSN,xdim,ydim);  
% %     ESN_spatial(:,:,t_assign)      = reshape(ESN,xdim,ydim);
% %     ELitter_spatial(:,:,t_assign)  = reshape(ELitter,xdim,ydim);
% %     ESN_In_spatial(:,:,t_assign)  = reshape(ESN_In,xdim,ydim);
% % 
% %     % Outputs for snow
% %     SND_spatial(:,:,t_assign)      = reshape(SND,xdim,ydim);
% %     SWE_spatial(:,:,t_assign)      = reshape(SWE,xdim,ydim);    
% %         
% %     % Outputs for groundwater
% %     FROCK_spatial(:,:,t_assign)    = reshape(FROCK,xdim,ydim); 
% %     f_spatial(:,:,t_assign)        = reshape(f,xdim,ydim); 
% %     Gfin_spatial(:,:,t_assign)     = reshape(Gfin,xdim,ydim); 
% %     %GPP_H_spatial(:,:,t_assign) 
% %     %GPP_L_spatial(:,:,t_assign) 
% %     G_spatial(:,:,t_assign)        = reshape(G,xdim,ydim);
% %     hc_H_spatial(:,:,t_assign)     = reshape(sum(hc_H.*Ccrown_OUTv,2),xdim,ydim);
% %     hc_L_spatial(:,:,t_assign)     = reshape(sum(hc_H.*Ccrown_OUTv,2),xdim,ydim);
% % 
% %     % Outputs runoff
% %     %Q_channel_spatial(:,:,t_assign)       = reshape(Q_channel,xdim,ydim);
% %     %q_runon_spatial(:,:,t_assign)         = reshape(q_runon,xdim,ydim);
% %     QpointC_series(t_assign,:)     = QpointC;
% % 
% %     Rd_spatial(:,:,t_assign)       = reshape(Rd,xdim,ydim);

    %% SAVING
    %======================================================================
    % Only save when it is the last day of the month or the last iteration
    %======================================================================
    if t_date  == t_store | t == N_time_step

    % Saving
    %----------------------------------------------------------------------
    disp(['Storing results in .mat file for: ' char(string(Date_run, 'MMM')) '_' char(string(Date_run, 'yyyy'))])

    % Daily calculation before saving
    %----------------------------------------------------------------------
    all_hours = length(date_forMonth); % all hours within the month
    div = 1:24:(all_hours+1); %divisions for calculation
    
    % Variables that must be summed - check that there aren't more!
    var_summ = ["Pr_sno_hourly";   "Pr_liq_hourly";    "Pr_St_hourly";  "EG_hourly";       "EICE_hourly";    "EIn_H_hourly"; 
                "EIn_L_hourly";    "EIn_rock_hourly";   "EIn_urb_hourly";  "EWAT_hourly";    "T_H_hourly";
                "T_L_hourly";      "ESN_hourly";        "ELitter_hourly";  "ESN_In_hourly";  "Smelt_hourly";
                "Imelt_hourly";    "Rd_hourly";         "Rh_hourly";       "f_hourly";       "Lk_rock_hourly";
                "Lk_hourly";       "Lk_wat_hourly"]; % variables that must be summed
    
    for k=1:length(VAR_OUT_HOUR)
        for y = 1:(length(div)-1)
            if VAR_OUT_HOUR(k) == "QpointC_hourly" 
                QpointC_daily(y,:) = mean(QpointC_hourly(div(y):(div(y+1)-1),:)); % Discharge
            else
                if ismember(VAR_OUT_HOUR(k), var_summ)
                commandStr = [char(VAR_OUT_DAY(k)) '(:,:,y) = sum(' char(VAR_OUT_HOUR(k)) '(:,:,div(y):(div(y+1)-1)),3)' ';']; %Other matrices
                eval(commandStr);                
                else
                commandStr = [char(VAR_OUT_DAY(k)) '(:,:,y) = mean(' char(VAR_OUT_HOUR(k)) '(:,:,div(y):(div(y+1)-1)),3)' ';']; %Other matrices
                eval(commandStr);
                end
            end
        end
    end

    % Compressing daily results  
    %----------------------------------------------------------------------   
%     %i=1; %Works only on Max' meteo inputs
%     for i=1:length(VAR_OUT_DAY)     
%         if VAR_OUT_DAY(i) ~= "QpointC_daily" 
%             commandStr = [char(VAR_OUT_DAY(i)) ' = Data_compressor_grid(' char(VAR_OUT_DAY(i)) ',"mm","compress")' ';']; %Other matrices
%             eval(commandStr);
%         end
%     end
	
    % Converting discharge to table to preserve names
    QpointC_daily = array2table(QpointC_daily, 'VariableNames', NAMEout);
    
    % Saving daily results
    %----------------------------------------------------------------------    
    VAR_CELL = cellstr(VAR_OUT_DAY);
    save([Directories.save '3D_Outputs_'  char(string(Date_run, 'MMM')) '_' char(string(Date_run, 'yyyy')) '.mat'], ...
            VAR_CELL{:});

    % Saving point data
    %------------------------------------------------------------------------
    %Note this doesn't need to happen every month, but it means data is
    %saved in case of crash

    save(strcat([Directories.save,'Points_Data.mat']),'Points_Data','Points','-v7.3','-nocompression'); %Saving the data and point info for completeness

    % Saving
    %----------------------------------------------------------------------
%     disp('Storing results in .mat file')
% 
%     save([Directories.save '/3D_Outputs_' yy '_' mth '.mat'], ...
%          'ANPP_H_spatial',   'ANPP_L_spatial',    'EG_spatial',      'EIn_H_spatial',  'EIn_L_spatial', ...
%          'EIn_rock_spatial', 'EIn_urb_spatial',   'EWAT_spatial',    'T_H_spatial',    'T_L_spatial', ...
%          'SSN_spatial',      'ESN_spatial',       'ELitter_spatial', 'ESN_In_spatial', 'FROCK_spatial', ...
%          'f_spatial',        'Gfin_spatial',      'G_spatial',       'hc_H_spatial',   'hc_L_spatial', ...
%          'EICE_spatial',     'LAI_H_spatial',     'LAI_L_spatial', ...
%          'NDVI_spatial',     'SND_spatial' ,      'SWE_spatial', ...
%          'Rd_spatial',       'QpointC_series' ...
%           );

    % Changing the label to create again the matrices for storing
    %--------------------------------------------------------------------------
    output_creation = 0;
    
    %% Memory 
    %----------------------------------------------------------------------
    % Information about the memory used by variables when saving occurs
    %----------------------------------------------------------------------
    
    % Get a structure array of all variables in the current workspace
    %--------------------------------------------------------------------------
    vars = whos;
    
    % Iterate through the variables and sum their sizes
    %--------------------------------------------------------------------------
    totalSize = 0;
    for i = 1:length(vars)
        totalSize = totalSize + vars(i).bytes;
    end
    
    % Display variables size
    %--------------------------------------------------------------------------
    disp(['Memory used by variables at the end of MAIN_FRAME: ', num2str(totalSize/1e6), ' MB']);
    
    % Saving actual conditions in the model to restart if needed
    %----------------------------------------------------------------------
    %run(['INIT_COND_MID.m'])

    end

    %run(['OUTPUT_MANAGER_DIST_LABEL.m']);

    % Save the workspace at frequent interval. Very useful in case it crashes 
    %if  mod(t,25)==0    
    %    save([outlocation, Fstep], '-regexp', '^(?!(FF_BC|LWIN_BC|PARB|PARD|PP_BC|PRESS_BC|RH_BC|SAB1|SAB2|SAD1|SAD2|TA_BC|WS)$).');
    %end   

    %if  mod(t,8760)==0  ||  t==N_time_step
    %    Fstep2= strcat(Fstep,'_',num2str(t));
    %    save([outlocation, Fstep2], '-regexp', '^(?!(FF_BC|LWIN_BC|PARB|PARD|PP_BC|PRESS_BC|RH_BC|SAB1|SAB2|SAD1|SAD2|TA_BC|WS)$).');
    %end


end
%close(bau)

%% Display
Computational_Time =toc;
disp('COMPUTATIONAL TIME [h] ')
disp(Computational_Time/3600)

%Q_channel