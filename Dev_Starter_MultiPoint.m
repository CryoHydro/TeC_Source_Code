%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C) - ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%
%%%%%%%%%%%%%%         MULTIPOINT SCALE MODEL                   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AUTHOR INFO AND STUDY SITE
%==========================================================================
%{
Created on 1st December 2025
Author: Catriona Fyffe
Original code: MAXIMILIANO RODRIGUEZ and Achille Jouberton

Code explanation: 
This code launches the Point Scale version of TC model.
%}

%Note 23...... is x/Lon in UTM and 89....... is y/Lat in UTM
% Latitute is -9......... and Longitute is -77.......
%==========================================================================
%% SITE SETUP
%==========================================================================

site_num_list = [1 2 3]; %Please number sequentially
site_name_list = {'Shallap' 'Torre_Orsina' 'Hybas_Lvl7_Catchment_2070497830'};
num_sites = size(site_num_list,2);
site_num = 1; %Choose site to run, this then selects the correct site directories
site_name = site_name_list{site_num};

% Choose point or multipoint
%--------------------------------------------------------------------------
IniCond.mode = "point";

% if point mode, then select the point
selected_ID = 152; %(only if running on laptop, SLURM can be used to state when using HPC) 

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
IniCond.SITE = site_name;

%%======================================================================
% Model structure choices
%========================================================================

OPT_Forcing = CONFIG_vals.OPT_Forcing; %Choose type of forcing, 1= from AWS which is distributed, 2= fully gridded
OPT_Veg_Param = CONFIG_vals.OPT_Veg_Param; %Choose if vegetation parameters (for Ccrown, Cwat, Crock, Curb, Cbare) are one value per landcover type (1), or gridded and variable per landcover (2). 
%For 1 change the values in VPAR_T directly, for 2 provide the grids in the
%dtm_file. Note for 2 you still need the table to read the Veg_type names
%but the values are not used. 

%% DIRECTORIES - ALL RELATIVE - NO NEED TO UPDATE 
%==========================================================================
% All paths are here
%=========================================================================

% Sub-path for outputs
Directories.save = ['Outputs/PointScale/' IniCond.run_folder '/'];  %Do we make this site specific?

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
end 

% Parameters
Directories.Vegpar = [site_name '/Parameters/Parameters_TC.xlsx']; 
Directories.Optpar = [site_name '/Parameters/Options_and_NoVegParams.xlsx'];
Directories.VPAR = [site_name,'/Parameters/VPAR_T.xlsx'];

% Preprocessing information
Directories.PreProc = [site_name '/Preprocessing/OUTPUTS/' CONFIG_vals.PreProc_folder '/']; 
dtm_file = CONFIG_vals.dtm_file;

% Sub-path points of interest for discharge and multipoints
Directories.POIs = [site_name '/Preprocessing/OUTPUTS/']; %POI_dtm_Site_50m.txt
Directories.Multi = [site_name '/Multipoint/OUTPUTS/' site_name '_MultiPoints.txt'];

% Dependencies
addpath(genpath('Functions')); %Location of additional functions
addpath(genpath([site_name '/Forcing'])); % Add path to Ca_Data
addpath(genpath([site_name '/Preprocessing/OUTPUTS'])); % Where are distributed model set-up files 
addpath(genpath('T&C_Code')); %T&C source code 
addpath(genpath('Outputs'));


% Modelling period and time
%--------------------------------------------------------------------------
dateRun.start = CONFIG_vals.Run_StartDate;  % Starting point of the simulation - start on the first of the month
dateRun.end = CONFIG_vals.Run_EndDate; % Last timestep of the simulation 
Date = datetime(dateRun.start):hours(1):datetime(dateRun.end);
Date = Date'; dateRun.list = Date; 

dt=3600; %% [s]
dth=1; %% [h]
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE

%% LOCATION OF OUTPUTS - CREATION OF FOLDERS
%==========================================================================

%% CHECK ALL ARE NEEDED!

if ~exist(Directories.save, 'dir') 
    disp('Folders do not exist for outputs. Creating new folders')
    mkdir(Directories.save); 
    mkdir([Directories.save '1_Hourly/']); 
    mkdir([Directories.save '2_Daily/']); 
    mkdir([Directories.save '3_Parameters/']); %Use this
    mkdir([Directories.save '4_INIT_COND/']); 
    mkdir([Directories.save '5_Env/']); 
    mkdir([Directories.save '6_POI_table/']); %Use this
addpath(genpath(Directories.save)); 
end

disp('Directories and folders done')

%% MODEL INPUTS
%==========================================================================

% Points of interest
%--------------------------------------------------------------------------
POI = readtable([site_name,'/Multipoint/OUTPUTS/' site_name '_MultiPoints.txt']); %import table with points info
UTM_zone = CONFIG_vals.UTM_zone;

% names of points
names = string(POI.Name);
selected_point = string(POI.Name(POI.ID==selected_ID)); %so choose selected point based on selected ID

% Parameters for vegetation
%--------------------------------------------------------------------------
%{
Here are the parameters of the model for vegetation.
opts check the format of the columns. 
Use opts to force all the columns with values to be numeric.
%}
%--------------------------------------------------------------------------
opts = detectImportOptions(Directories.Vegpar);
opts = setvartype(opts, [7:length(opts.VariableTypes)], 'double');
TT_par = readtable(Directories.Vegpar, opts);

% ==================================================================
%% LOAD SPATIAL DATA
%==========================================================================

% Load preprocessed data
%--------------------------------------------------------------------------
mm = load([Directories.PreProc dtm_file]); % Distributed maps pre-processing. Useful here to get the DTM and initial snow depth

% Central lat and lon Lat
%--------------------------------------------------------------------------
UTM_Y = mm.y(floor(length(mm.y)/2));
UTM_X = mm.x(floor(length(mm.x)/2));

[Lat, Lon] = utm2deg(UTM_X, UTM_Y, CONFIG_vals.UTM_zone_text);

% Extraction of DEM and features
%--------------------------------------------------------------------------
DTM = mm.DTM_orig; % Use the full DEM in case running POI outside of mask
DTM(isnan(DTM)) = 0; %%??

[m_cell,n_cell]=size(DTM);
num_cell=numel(DTM);

xllcorner = mm.xllcorner; %bottom corner
yllcorner = mm.yllcorner; %bottom corner
cellsize = mm.cellsize;
VEG_CODE = mm.VEG_CODE; 

%=========================================================================
%% MODEL PARAMETERS
%==========================================================================

%--------------------------------------
% 1. General options
%---------------------------------------

%Set numerical tolerances
OPT_Num.Opt_CR=                     optimset('TolFun',1);%,'UseParallel','always'); %Numerical tolerance for internal CO2 computation
OPT_Num.OPT_PH=                     odeset('AbsTol',0.01); %Numerical tolerance for internal Plant Hydraulic Volumes
OPT_Num.OPT_SM=                     odeset('AbsTol',0.05,'MaxStep',dth); %Numerical tolerance for soil moisture differential equations
OPT_Num.Opt_ST=                     optimset('TolFun',0.1);%,'UseParallel','always'); %Numerical tolerance for surface temperature computation
OPT_Num.Opt_ST2 =                   optimset('TolFun',0.1,'Display','off');
OPT_Num.OPT_STh=                    odeset('AbsTol',5e+3); %Numerical tolerance for heat transfer differential equations
OPT_Num.OPT_VD=                     odeset('AbsTol',0.05); %Numerical tolerance for carbon budget differential equations

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


%% PIXELS
%==========================================================================
% Construction of POI matrix with details of points.
%==========================================================================
% k = 10;
for k = 1:height(POI)
id_location = char(string(POI.Name(k))); %id

y_coord = POI.LAT_UTM(k); %Y
x_coord = POI.LON_UTM(k); %X

pos_col = floor((x_coord - xllcorner) / cellsize) + 1;
pos_row = floor((y_coord - yllcorner) / cellsize) + 1;

%ij = POI.idx(loc);
ij = sub2ind(size(DTM),pos_row, pos_col); % Location
[j, i] = ind2sub(size(DTM), ij); % Location

POI.ij(k) = ij;
POI.i(k) = i; %Same as column
POI.j(k) = j; %Same as row

POI.Zbas(k) = DTM(j,i); % Altitude
end

%---------------------------------------------------------------------
% 2. Vegetation parameters
%--------------------------------------------------------------------------

%% Surface elevation for vegetation
%--------------------------------------------------------------------------
% Here, note that it is needed to define the starting point of the data in
% the excel file. In this case it is from column 5 to the end. 
%--------------------------------------------------------------------------
%NOTE zatm_surface needs adjusted
zatm_surface = TT_par(strcmp(TT_par.Parameters,'zatm_surface'),7:size(TT_par,2));
zatm_hourly_on = 0;

%% Sort VPAR info
%Veg/land parameters (was POI but renamed to reduce confusion with Points
%of interest)
ksv=reshape(VEG_CODE,num_cell,1);

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
    VPAR.Ccrown = mm.CCROWN;
    VPAR.Cwat = mm.CWATER;
    if exist('mm.CURB','var') == 1
        VPAR.Curb = mm.CURB;
    else
    VPAR.Curb = mm.CWATER.*0; %As urban 0 everywhere
    end
    VPAR.Crock = mm.CROCK;
    VPAR.Cbare = mm.CBARE;
end  

%Reshape
VPAR.Ccrownr = reshape(VPAR.Ccrown,num_cell,1);
VPAR.Cwatr = reshape(VPAR.Cwat,num_cell,1);
VPAR.Curbr = reshape(VPAR.Curb,num_cell,1);
VPAR.Crockr = reshape(VPAR.Crock,num_cell,1);
VPAR.Cbarer = reshape(VPAR.Cbare,num_cell,1); 

%% Loop for vegetation classification

%Now take from grids directly
for k = 1:height(POI)
    loc_p = POI.ij(k); %Index of the point
    POI.VEG_CLASS(k) = ksv(loc_p); %Landcover
    POI.Cwat(k) = VPAR.Cwatr(loc_p); % Cwat
    POI.Curb(k) = VPAR.Curbr(loc_p); %Curb
    POI.Crock(k) = VPAR.Crockr(loc_p); %Crock
    POI.Cbare(k) = VPAR.Cbarer(loc_p); %Cbare
    POI.Ccrown(k) = VPAR.Ccrownr(loc_p); %Ccrown
    %% CHECK - what if classes are not sequential?
    POI.II(k) = VPAR.Veg_type(POI.VEG_CLASS(k)); %Veg type name, use the veg number to find
    POI.NCrown(k) = length(POI.Ccrown(k)); %Number of crown values, in my set up always 1
    POI.zatm(k) = zatm_surface{:,POI.VEG_CLASS(k)}; %choose correct atmospheric reference height, note the column of zatm_surface = row of VPAR
end
POI.cc_max = repmat(max(POI.NCrown),height(POI),1);

%% SAVE PARAMS and POI TABLES FOR REFERENCE
%==========================================================================

disp('POI and veg params done')

%Save the parameter files with the model outputs
%Clean CONFIG for saving
CONFIG_s = CONFIG;
Cmask = cellfun(@(CONFIG_s) any(isa(CONFIG_s,'missing')), CONFIG_s);
CONFIG_s(Cmask) = {[]}; % 

writetable(POI, [Directories.save '6_POI_table/POI_' IniCond.run_folder '.xlsx']);
writetable(TT_par,[Directories.save '3_Parameters/' 'TT_par.xlsx']);
writetable(OPT_PARAM,[Directories.save '3_Parameters/' 'OPT_PARAM.xlsx']);
writetable(VPAR_T,[Directories.save '3_Parameters/' 'VPAR_T.xlsx']);
writecell(CONFIG_s,[Directories.save '3_Parameters/' 'CONFIG.xlsx']);

%-----------------------------------------------------------------------
% 3. Glacier and snow parameters - in Dev_run_Point_Pro
% --------------------------------------------------------------------------
%-----------------------------------------------------------------
% 4. Carbon
%-------------------------------------------------------------------------

load('Ca_Data.mat');
d1 = find(abs(Date_CO2-datenum(Date(1)))<1/36);d2 = find(abs(Date_CO2-datenum(Date(end)))<1/36);
Ca=Ca(d1:d2);
clear d1 d2 Date_CO2 id_df
Oa = 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -


%% LOAD AND PART A OF METEOROLOGY
%==============================================================================

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
Meteo_Params.hm_Ta_lapse = load(ta_lapse_file); %These are the Ta lapse rates from TOPKAPI, per hour and month for main station 144, catchment 406
Meteo_Params.hm_Ta_lapse = Meteo_Params.hm_Ta_lapse.(ta_lapse_name); %Pull out of structure. The lapse rate column should be named 'Ta_lapse'
Meteo_Params.Ldown_lapse = -0.031; %from Marty et al. (2002) p145
%Pr lapse - based on ratio per month
Meteo_Params.Pr_lapse = load(pr_lapse_file);
Meteo_Params.Pr_lapse = Meteo_Params.Pr_lapse.(pr_lapse_name); %Take out of structure, this is the fit of ratio = (p1 * elevation) + p2. Just multiply ratio by Pr. 
%Radiation lapse rates
Meteo_Params.SAB1_lapse = 0.0049; %Wm-2 per m (from Simone, for the Alps, may be a little high for Andes), radiation increases with elevation
Meteo_Params.SAB2_lapse = 0.0080; %Wm-2 per m (from Simone, for the Alps)
Meteo_Params.PARB_lapse = 0.0047; %Wm-2 per m (from Simone, for the Alps)
%Station elevations
Meteo_Params.St_el = CONFIG_vals.Meteo_el;

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
    %Load data
    OPT_Run_BC = CONFIG_vals.OPT_RunBC;
    dsERA = load([meteo_name '_ERAdata.mat']);
    dsGRID = load([meteo_name '_GRIDdata.mat']);
    dsCoeff = load([meteo_name '_BC_Coefficients.mat']);
    dim = 'grid';
    n = 1500;
    demCurv = getcurvature(dsGRID.demZs,n);
    dsGRID.pixelSize = abs(dsGRID.demLats(1,1)-dsGRID.demLats(2,1))*111.1e3; % (m)
    [dsGRID.demSlope,dsGRID.demAspect] = slopeaspect(dsGRID.demZs,dsGRID.pixelSize);
    % Precipitation gradient + Mean Annual Precipitation
    dsERA.tpLr = 0.003; % Precipitation gradient (mm yr-1 m-1) - CAREFUL - should this be a param?
    nYrs = years(Date(end)+(hours(1)) - Date(1));
    dsERA.tpMean = nansum(dsERA.tp,3)/nYrs;
    % Set output grids to be the dem grid resolutions (default are provided
    % bias-correction locations if not)
    dsGRID.dsLats = dsGRID.demLats; dsGRID.dsLons = dsGRID.demLons; dsGRID.dsZs = dsGRID.demZs;
end

%Tmod parameters
OPT_Tmod = OPT_PARAM_vals.OPT_Tmod;
if OPT_Tmod == 0 %No Tmod applied
    Meteo_Params.TmodB = 0;
elseif OPT_Tmod == 1 %Simple bias correction
    Meteo_Params.TmodB = OPT_PARAM_vals.TmodB; %Tmod bias (positive reduces temperature over ice)
elseif OPT_Tmod == 2 %Uses a relationship with optional bias correction 
    Meteo_Params.TmodB = OPT_PARAM_vals.TmodB; %Tmod bias (positive reduces temperature over ice)
    Meteo_Params.TmodM = OPT_PARAM_vals.TmodM; %Tmod multiplier
    Meteo_Params.TmodC = OPT_PARAM_vals.TmodC; %Tmod coefficent
end


% Solar variables (just for Lmax_day)- already cropped to Run time - maybe
% not necessary as given in param/init conditions?
%--------------------------------------------------------------------------
L_day=zeros(length(Datam),1); 
for j=2:24:length(Datam)
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),IniCond.DeltaGMT,Lon,Lat,t_bef,t_aft);
end

Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day');

% Radiation variables for  topographical influence
%............................................................
%%% Slo_top [Fraction] %%% Aspect [rad] from N
[Slo_top,Aspect]=Slope_Aspect_indexes(mm.DTM_orig,cellsize,'mste');
Aspect(isnan(Aspect))=0;
Slo_top(Slo_top<0.001)=0.001;  %to avoid flow routing issues
Slo_top(mm.Youtlet,mm.Xoutlet)=0.05; %NOTE on maps and in scatter use (X,Y) but to index use Y (row) X (col)
%Note other variables (Asur,Ared,aTop,SN,Slo_pot,a_dis, pow_dis) given
%directly in Dev_Run_point_Pro as were in MOD_PARAM Multipoint

% Computation Horizon Angle
[HZ,Zasp] = Horizon_Angle(mm.DTM_orig,cellsize);
%%% HZ Horizon angle array [angular degree]
%%% Z Azimuth directions  [angular degree] from N

%%% Sky View Factor and Terrain Configuration Factor
[SvF,Ct] = Sky_View_Factor(mm.DTM_orig,atan(Slo_top)*180/pi,Aspect,HZ,Zasp);
ShF=zeros(m_cell,n_cell,size(Date,1),'single'); %Needs to be single for memory for long runs
for t=1:size(Date,1)
[h_S(t),delta_S(t),zeta_S(t),T_sunrise(t),T_sunset(t),L_day(t)]= SetSunVariables(Datam(t,:),IniCond.DeltaGMT,Lon,Lat,t_bef,t_aft); %These are one value per time step
[ShF(:,:,t)] = Shadow_Effect(DTM,h_S(t),zeta_S(t),HZ,Zasp); %ShF is x,y grid
end

%Collect sun variables
Sun_Vars.ShF = ShF; %grid by timestep
Sun_Vars.zeta_S = zeta_S; %timeseries
Sun_Vars.h_S = h_S; %timeseries
Sun_Vars.SvF = SvF; %single grid
Sun_Vars.Ct = Ct; %single grid
Sun_Vars.Slo_top = Slo_top; %single grid
Sun_Vars.Aspect = Aspect; %Single grid

% Crop to time
%.........................................................................

id_date_for_run = ismember(Forcing_Date,Date); %So this is the index of the forcing data for this month
%  ===== For weather station point =========
if OPT_Forcing == 1
% So for weather station point Meteo_P is the meteorology at the station
Meteo_P=array2table(Date);
Meteo_P.Properties.VariableNames = {'Date'};
%Temperature
Meteo_P.Ta = Ta(id_date_for_run); %For the hour
%Sort later.....Ta_P_day = Ta(pdind); %For the 24 hours before
%Longwave incoming
Meteo_P.N = Nin(id_date_for_run);
%Relative humidity
Meteo_P.U = U(id_date_for_run)./100; %Change to fraction
%Wind speed
Meteo_P.Ws = Ws(id_date_for_run);
Meteo_P.Ws(Meteo_P.Ws < 0.01) = 0.01; %Should not be zero, but double check
%Precipitation
Meteo_P.Pr = Pr(id_date_for_run);
%Incoming shortwave radiation
Meteo_P.SWin = SWin(id_date_for_run);
%Dew point temperature
Meteo_P.Tdew = Tdew(id_date_for_run);
%Radiation
Meteo_P.SAD1 = SAD1(id_date_for_run)'; %Note need to transpose here
Meteo_P.SAD2 = SAD2(id_date_for_run)';
Meteo_P.SAB1 = SAB1(id_date_for_run)';
Meteo_P.SAB2 = SAB2(id_date_for_run)';
Meteo_P.PARB = PARB(id_date_for_run)';
Meteo_P.PARD = PARD(id_date_for_run)';
%Air pressure will be calculated based on elevation
%Vapour pressure will be calculated per point
end

disp('Meteo loaded')

%% SINGLE POINT LAUNCHER FOR PERSONAL COMPUTER
%==========================================================================

if IniCond.mode == "point"

% Time counter
%--------------------------------------------------------------------------
tic;

%Main function
%--------------------------------------------------------------------------
Dev_run_Point_Pro(Directories, IniCond, selected_point, POI, ksv, dateRun, TT_par, zatm_surface, dtm_file,t_bef,t_aft,Meteo_P,Meteo_Params,Datam,Sun_Vars,dt,dth,OPT_PARAM_vals,mm,Ca,Oa,CONFIG_vals,VPAR,OPT_Num);

% Computational time
%--------------------------------------------------------------------------
Computational_Time =toc;
disp('Computation time - Single Point')
disp([num2str(round(Computational_Time/60,1)) ' mins'])

%% Memory use - Windows only
if ~contains(Directories.root,"nfs") 
[user, sys] = memory; % Windows only
disp(['Mem used in the Starter: ', num2str(user.MemUsedMATLAB/1e6), ' MB'])
end

% Memory out
%--------------------------------------------------------------------------
%profile off
%profile viewer


%% MULTIPLE POINT LAUNCHER 
%==========================================================================
% Specify the desired number of workers
%==========================================================================

elseif IniCond.mode == "multipoint"

% Workers - Only for personal computer
%--------------------------------------------------------------------------

%if I am working on my personal computer define this
if ~contains(Directories.root,"nfs") 

    numWorkers = 2; % Example: Set to 4 workers
    
    % Create a parallel pool with the specified number of workers
    poolobj = gcp('nocreate'); % Check if a pool already exists
    if isempty(poolobj)
        parpool(numWorkers); % Create a new pool if one doesn't exist
    elseif poolobj.NumWorkers ~= numWorkers
        % If a pool exists with a different size, delete it and create a new one
        delete(poolobj);
        parpool(numWorkers);
    end

end


% Computational time - Time counter - For personal computer or HPC
%--------------------------------------------------------------------------
tic;

% Specific tasks - Here define the runs - TO CHECK
%--------------------------------------------------------------------------

Torun = names; %{"Cluster_26","Cluster_11"}; %Set names to run

% %Not sure this works
% idxn=zeros(size(names,1),size(Torun,2));
% for r=1:size(Torun,2)
% idxn(:,r) = ismember(names,Torun{r});
% end
% idxn2 = sum(idxn,2); %Index to run

% %specific_indices = 1:length(names); % for all the runs
% specific_indices = 1:20;
% 
% % Create a new cell array of names based on the specific indices
% selected_names = names(specific_indices);

% Parallel computing launch for the Model
%--------------------------------------------------------------------------
parfor k = 1:length(Torun) %length(names)  % 3    
    try
        Dev_run_Point_Pro(Directories, IniCond, Torun{k}, POI, ksv, dateRun, TT_par, zatm_surface, dtm_file,t_bef,t_aft,Meteo_P,Meteo_Params,Datam,Sun_Vars,dt,dth,OPT_PARAM_vals,mm,Ca,Oa,CONFIG_vals,VPAR,OPT_Num);
    catch ME
        warning('Error occurred on worker %d: %s', k, ME.message);       
    end
end

% Computational time - End of processing
%--------------------------------------------------------------------------
Computational_Time =toc;
disp('Computation time - MultiPoint ')
disp([num2str(round(Computational_Time/60,1)) ' mins'])

%% Memory use - Windows only
if ~contains(Directories.root,"nfs") 
[user, sys] = memory; % Windows only
disp(['Mem used by the Starter: ', num2str(user.MemUsedMATLAB/1e6), ' MB'])
end

else %if the label is not point or multipoint, define it well
disp('The label is not right. Choose "point" or "multipoint"')
end