%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C) - ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%
%%%%%%%%%%%%%%              POINT SCALE MODEL                   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AUTHOR INFO AND STUDY SITE
%==========================================================================
%{
Created on 9 Oct 2026
Author: Cat Fyffe
Code originally from: MAXIMILIANO RODRIGUEZ and ACHILLE JOUBERTON

Code explanation: 
     This code launches the Point Scale version of TC model.
     This function runs only one point. There is no for loop.       
%}
%==========================================================================

%% Debugger
%Point = "VelinoCluster95";


%% As a Function 
function result=Dev_run_Point_Pro(Directories, IniCond, Point, POI, ksv, dateRun, TT_par, zatm_surface, dtm_file,t_bef,t_aft,Meteo_P,Meteo_Params,Datam,Sun_Vars,dt,dth,OPT_PARAM_vals,mm,Ca,Oa,CONFIG_vals,VPAR,OPT_Num)

%% Inputs

% DeltaGMT must be defined for MAIN_FRAME_Pro
DeltaGMT = IniCond.DeltaGMT;

%% MODEL PARAMETERS
%==========================================================================

%% Data storing  
%==========================================================================
% How to store outputs
% Recommended for long-term simulations (> 10-20 years, otherwise data volume is too important)
%==========================================================================
output_daily = 1; 

%=========================================================
%% RELEASE SPATIAL DATA
%========================================================
%Already loaded as mm
disp(strcat('Model resolution: ',num2str(mm.cellsize))); % simulation resolution [m]
v2struct(mm);
DTM = DTM_orig; % Use the full DEM in case running POI outside of mask
DTM(isnan(DTM)) = 0; %%??

[m_cell,n_cell]=size(DTM);
num_cell=numel(DTM);

%Swapping somegrids


% ----------- Pull out info for location -------------------

Names = string(POI.Name);
loc = find(Point == Names);  %Finds ID of point

%Definitions of site
% Crowns
cc = POI.NCrown(loc); 
II = POI.II{loc};
Cwat = POI.Cwat(loc); 
Curb = POI.Curb(loc); 
Crock = POI.Crock(loc);
Cbare = POI.Cbare(loc);
Ccrown = POI.Ccrown(loc);
zatm = POI.zatm(loc);
id_location = string(POI.Name(loc));
cc_max = POI.cc_max(loc);

%% Locations
Zbas = POI.Zbas(loc);
Lat = POI.LAT(loc);
Lon = POI.LON(loc);
ij = POI.ij(loc);
ROW = POI.ROW(loc);
COL = POI.COL(loc);

%------------------------------------------------
% 3. Glacier and snow parameters
%---------------------------------------------------

%Precipitation phase parametrization
%1 = 2-threshold, 2 = Ding 2017, 3 = single-threshold, 4 = Pomeroy 2013, 5
%= Wang 2019, 6 = Jennings 2018
parameterize_phase.OPT_Pr_Part = OPT_PARAM_vals.OPT_Pr_Part; % Choice of the precipitation phase scheme
parameterize_phase.Tmax = OPT_PARAM_vals.Tmax; % Maximum air temperature for precipitation phase scheme (2-dual thresholds)
parameterize_phase.Tmin = OPT_PARAM_vals.Tmin; % Minimum air temperature for precipitation phase scheme (2-dual thresholds)
parameterize_phase.Tconst = OPT_PARAM_vals.Tconst; % Air temperature for constant thresholds

%Multilayer snowpack parameters
hSTL = OPT_PARAM_vals.hSTL; % Skin layer thickness for the 2-layer snowpack module
%CHECK - min_SPD not used in MAIN_FRAME
min_SPD = OPT_PARAM_vals.min_SPD;  %% [m] minimum snow pack depth to have a multilayer snow 

% Albedo scheme choice
OPT_Albsno = OPT_PARAM_vals.OPT_Albsno; % 3 doesn't work, 4 is Brock 2000, 5 is Ding 2017
rho_g = OPT_PARAM_vals.rho_g; %%% Spatial Albedo

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

Afirn = reshape(Afirn,num_cell,1);

%SNOWD is taken directly and used in intial conditions as initial snow depth
SNOWD=reshape(SNOWD,num_cell,1);

% Initial snow albedo
if ~exist('SNOWALB','var')
    Ini_Asnow = OPT_PARAM_vals.Ini_Asnow;
    SNOWALB = SNOWD;
    SNOWALB(SNOWD>0) = Ini_Asnow;
else 
    SNOWALB=reshape(SNOWALB,num_cell,1);
end 

%Glacier and Debris info
DEB = DEB_MAP(ROW,COL); 
idg = mm.GLH>0;
GLH = mm.GLH;
GLH(idg) = GLH(idg) + 400; %To avoid glacier disappearance. 
Glh = GLH(ROW,COL); 
md_max = OPT_PARAM_vals.md_max; % Number of debris layers

SNOWD=reshape(SNOWD,num_cell,1);

%----------------------------------------------------------------------
% 6. SOIL PARAMETERS 
%--------------------------------------------------------------------------
%==========================================================================
% Vector is divided by 100 to put the numbers within 0-1
% Original PSAN, PCLA and PORG come with values between 0-100, g/100g
%==========================================================================

%Reshape soil variables
SOIL_TH=reshape(SOIL_TH,num_cell,1);
PSAN=reshape(PSAN/100,num_cell,1); Psan = PSAN(ij); % Soil sand content at pixel ij
PCLA=reshape(PCLA/100,num_cell,1);  Pcla = PCLA(ij); % Soil clay content at pixel ij
if max(max(PORG>40)) %Presume PORG in g/kg not %
    PORGr=reshape(PORG/1000,num_cell,1);
else %PORG brought in as % 
    PORGr=reshape(PORG/100,num_cell,1);
end
Porg = PORG(ij); % Soil organic content at pixel ij

%  ---  Set Soil thickness and related params ---
OPT_soil_th = OPT_PARAM_vals.OPT_soil_th;
%First Zs and the number of soil layers is based on the maximum soil depth
Zs_all = [0 10 20 50 100 150 200 250 300 350 500 1000 1500 2000 2500 3000 3500 4000 4500]; %All potential soil depths (ms_max +1)
vi_all = [0 10 20 50 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100]; %All potential initial soil moistures (ms_max +1)
if OPT_soil_th == 1
    soil_th_max = OPT_PARAM_vals.soil_th_max; %Max soil thickness, mm
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
ms = ms_max; %ms used in MAIN_FRAME

%Check
if max(max(SOIL_TH)) < 50 || max(max(SOIL_TH)) > 5000
    disp('Check soil units')
end

%% Initial condition
% Saving initial conditions of the model
out = strcat(Directories.save, '4_INIT_COND/INIT_COND_', IniCond.SITE ,'_MultiPoint_',Point,'.mat'); % file path initial conditions

%% CHECK THIS
topo = 1;



%% FORCING
%==========================================================================
% Distribute for full period based on lapse rates and cell info
%==========================================================================

% % Years needed - maybe don't need
% years_forc = year(dateRun.start):year(dateRun.end);
% forc_opened = struct(); 

%% Extraction of forcing data
% DEBUGGER
%--------------------------------------------------------------------------
%Already downloaded and cropped
% forc_file= [Directories.forc 'Forcing_ERA5Land_' char(Point) '.mat']; % Put here the path of where you downloaded the repository
% load(forc_file); % Load forcing table for the current POI

%forc(forc.t2m<forc.d2m,:)

% % Decompressing - Only for big data sets
% %--------------------------------------------------------------------------
% forc = Data_compressor(forc,'decompress');

%DEBUGGING
%find(forc.t2m<forc.d2m)


%% Combination of data into a single table
Date = dateRun.list;

% %define period and time zone info
% x1=find(dateRun.start == dateRun.list,1);
% x2=find(dateRun.end == dateRun.list,1);


%% Displaying modelling parameters
disp(strcat("Site selected: ", IniCond.SITE));
disp(['Running T&C for pixel: ' char(id_location)]);
disp(['Simulation period: ' datestr(dateRun.start) ' to ' datestr(dateRun.end)]);

%% Fetch time and do date handling
%Date = dateRun.list;
% [YE,MO,DA,HO,MI,SE] = datevec(Date); %Datam and input
% Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
% clear YE MO DA HO MI SE

%% Meteorological input
%==========================================================================
%{
FORCINGS: 
    Precipitation [mm]
    Air Pressure [Pa]
    Temperature [C]
    Wind Speed [m s-1]    
    Radiation [W m-2]    
    Dew Point temperature [C]

CALCULATIONS:
    Relative Humidity [-]
%}
%==========================================================================

% % Period of forcing data
% forcing = forc(x1:x2,:); %Don't need this

% Check function for variables - check_var()
%--------------------------------------------------------------------------
%{
Define variables names in the string in order based on your input file.
Do no change the order of variables in check_var. Just change the names
based on your file of forcing. So, if your column for temperature is called
Temp instead of t2m, then change t2m to Temp below.
%}

%Define variables names in order based on your input file

% check_var(Meteo_P, ...
%     ["Ta" ... % Temperature
%     "Tdew" ...  % Dew Point temperature
%     "Pr" ...   % Precipitation
%     "SWin" ... % Downward short wave radiation
%     "N" ... % Downdward Long wave radiation
%     "Ws" ... % Wind speed
%     "sp" ...   % Air pressure
%     "es" ...   % saturation vapor pressure
%     "ea" ...   % actual vapor pressure
%     "RH" ...   % Relative humidity
%     "SAD1" ... % SAD1
%     "SAD2" ... % SAD2
%     "SAB1" ... % SAB1
%     "SAB2" ... % SAB2
%     "PARB" ... % PARB
%     "PARD" ... % PARD
%     "N" ...    % Cloudiness
%     ],Point)

NN = size(dateRun.list,1);%%% time Step

% Height of virtual station
zatm_hourly = repmat(2.00,NN,1);
zatm_hourly_on = 0;


%Distribution of Meteo variables

%Sort lapse rates
if CONFIG_vals.OPT_Forcing == 1
    Meteo_P.month = month(Meteo_P.Date);
    Meteo_P.hour = hour(Meteo_P.Date);
    Meteo_P.Ta_Lapse = zeros(NN,1);
    for m=1:12
        for h=0:23
            id_mh = Meteo_P.month==m & Meteo_P.hour==h; %Find index of hour and month
            id_lamh = Meteo_Params.hm_Ta_lapse.month==m & Meteo_Params.hm_Ta_lapse.hour==h; %Find index of lapse rate
            Meteo_P.Ta_Lapse(id_mh)=Meteo_Params.hm_Ta_lapse.Ta_lapse(id_lamh);
        end
    end
    Meteo_P.Ta_Lapse = -Meteo_P.Ta_Lapse; %Ta lapse rates should be negative
    Meteo_P.Tdew_Lapse = Meteo_P.Ta_Lapse; %Its ok to use the same lapse rate
    
    %Air temperature, degrees C
    Ta = Meteo_P.Ta + (Meteo_P.Ta_Lapse.*(Zbas-Meteo_Params.St_el)); %Lapse rates should be negative
    %Apply Tmod
    if DEB==0 && Glh>0 %On clean ice
        %Remove bias for all time first 
        Ta = Ta - Meteo_Params.TmodB;
        idapply = Ta>0; %Apply multiplier when >0
        Ta(idapply) = Ta(idapply).* Meteo_Params.TmodM; %Apply Tmod multiplier
    end
    %Precipitation, mm - This uses the WRF ratio v el relationships
    %First find the ratio for each month related to the elevation of the cell
    Pr_ratio = zeros(NN,1);
    for m=1:12
        mid = Meteo_P.month==m;
        Pr_ratio(mid) = (Meteo_Params.Pr_lapse(m).month.p1.*Zbas) + Meteo_Params.Pr_lapse(m).month.p2;  %val(x) = p1*x + p2
    end
    Pr = Pr_ratio.*Meteo_P.Pr; %Apply ratio
    %Longwave radiation, W m-2
    N = Meteo_P.N + (Meteo_Params.Ldown_lapse.*(Zbas-Meteo_Params.St_el)); %IN LW decreases with elevation
    %Dew point temperature, degrees C
    Tdew = Meteo_P.Tdew + (Meteo_P.Tdew_Lapse.*(Zbas-Meteo_Params.St_el)); %Lapse rates should be negative
    %Direct radiation with elevation
    %SAB1 
    id0B1 = Meteo_P.SAB1==0; 
    SAB1 = zeros(NN,1);
    SAB1(~id0B1)=Meteo_P.SAB1(~id0B1) + (Meteo_Params.SAB1_lapse.*(Zbas-Meteo_Params.St_el)); %Radiation should increase with elevation
    SAB1(SAB1<0)=0;
    %SAB2 
    id0B2 = Meteo_P.SAB2==0; 
    SAB2 = zeros(NN,1);
    SAB2(~id0B2)=Meteo_P.SAB2(~id0B2) + (Meteo_Params.SAB2_lapse.*(Zbas-Meteo_Params.St_el)); %Radiation should increase with elevation
    SAB2(SAB2<0)=0;
    %PARB
    id0RB = Meteo_P.PARB==0; 
    PARB = zeros(NN,1);
    PARB(~id0RB)=Meteo_P.PARB(~id0RB) + (Meteo_Params.PARB_lapse.*(Zbas-Meteo_Params.St_el)); %Radiation should increase with elevation
    PARB(PARB<0)=0;
    %Pre
    Pre = (101325*((1-(0.0065.*Zbas./288.15)).^(9.18*0.0289644./(8.31447*0.0065))))./100; %/100 converts to mbar
    Pre = repmat(Pre,NN,1);

    %Do we apply a wind speed lapse rate?
    Ws = Meteo_P.Ws;
    %Diffuse radiation can be the same everywhere
    SAD1 = Meteo_P.SAD1;
    SAD2 = Meteo_P.SAD2;
    PARD = Meteo_P.PARD;

    %Use the lapsed dew point temprature to derive relative humidity
    c=237.3; b=17.27;
    U = 100*exp((c*b.*(Tdew - Ta))./((c+Ta).*(c+Tdew)));
    clear c b 
    U = U./100; %Turn into fraction
    U(U>1) = 1; %Cannot be more than 100%

elseif OPT_Forcing == 2
    
    %TO CHECK
    %Load table of data for given point
    forc_file= [Directories.forc 'Forcing_ERA5Land_' char(Point) '.mat']; % Put here the path of where you downloaded the repository
    load(forc_file); % Load forcing table for the current POI
    
    % Decompressing - Only for big data sets
    %--------------------------------------------------------------------------
    forc = Data_compressor(forc,'decompress');

    %** Also crop to time
    %Turn table to double in one go

    % Precipitation from ERA5Land
    Pr=double(forcing.tp);
    Pr(isnan(Pr))=0;
    Pr(Pr<0.01)=0;
    
    %Precipitation from stations
    
    % Air pressure
    %--------------------------------------------------------------------------
    % Air Pressure in mbar based on the PDF for variables and parameters of TC
    % Pressure comes in Pa from ERA5
    Pre=double(forcing.sp/100);    
    
    % Temperature
    %--------------------------------------------------------------------------
    % 2m air temperature
    Ta=double(forcing.t2m);
    
    % Wind Speed
    %--------------------------------------------------------------------------
    Ws=double(forcing.ws10);
    Ws(Ws < 0.01) = 0.01;
    
    % Relative humidity
    %--------------------------------------------------------------------------
    % Divided by 100 to set the number in the range 0-1
    U=double(forcing.RH/100);
    
    % Longwave radiation
    %--------------------------------------------------------------------------
    % N can be cloud cover [-] or longwave incoming radiation [W m-2]
    % Latm=forcing.LW_rad_downward_HH; % Latm:Incoming long wave radiation [W m-2]
    % N=ones(NN,1); % cloud cover [-]
    % N = forcing.LW_rad_downward_HH;
    N = double(forcing.N);
    
    % Radiation partition
    %--------------------------------------------------------------------------
    SAD1=double(forcing.SAD1); SAD2=double(forcing.SAD2); 
    SAB1=double(forcing.SAB1); SAB2=double(forcing.SAB2);
    PARB=double(forcing.PARB); PARD=double(forcing.PARD);
    
    % Vapor pressure - Dew Point temperature
    %--------------------------------------------------------------------------
    %esat/ea/Ds/Tdew
    esat = double(forcing.es);   % Vapour pressure at saturation (Pa)
    ea = double(forcing.ea);     % Vapour pressure (Pa)
    Tdew = double(forcing.d2m);

end


% Vapor pressure - calculate based on distributed Ta and U
a=17.27; b=237.3;
esat=611.*exp(a.*Ta./(b+Ta)); %Vapour pressure at saturation (Pa)
ea=U.*esat;                 %Vapour pressure (Pa)
Ds= esat - ea;              %Vapor Pressure Deficit (Pa)
Ds(Ds<0)=0; 
%Don't recalculate Tdew - its above
clear a b xr;

%% Then check_var

% Albedo parameters
%--------------------------------------------------------------------------
%Ameas = ones(NN,1);
%alpha = 0; % switch for albedo
%Ameas_t=0; % albedo
%Aice_meas_on_hourly = ones(height(forcing),1)/2; % albedo
%Asno_meas_on_hourly = ones(height(forcing),1)/2; % albedo

% Clear forcing and forc - Already used - Important for several clusters
%--------------------------------------------------------------------------
clear forc forcing


%% DING PARAMETERIZATION
% Initial daily mean values for Ding parametrization
%Are they used??

% Ta_Ding_d = nanmean(Ta(1:24));
% Pre_Ding_d = nanmean(Pre(1:24));
% ea_Ding_d = nanmean(ea(1:24));

%% RADIATION AND TOPOGRAPHY

%MASK = MASK.*0+1;
MASKn=reshape(MASK,num_cell,1);

if topo == 1
 
    % Topography
    %----------------------------------------------------------------------

    %Pull out sun variables
    ShF = squeeze(Sun_Vars.ShF(ROW,COL,:)); %Now a timeseries per pixel
    ShF = double(ShF);
    zeta_S = Sun_Vars.zeta_S; %timeseries already
    zeta_S = zeta_S';
    h_S = Sun_Vars.h_S; %timeseries already
    h_S = h_S';
    SvF = Sun_Vars.SvF(ROW,COL); % Sky view factor at pixel ij
    Ct = Sun_Vars.Ct(ROW,COL); %For pixel
    Slo_top = Sun_Vars.Slo_top(ROW,COL); % Slope at pixel ij
    Aspect = Sun_Vars.Aspect(ROW,COL); % Aspect at pixel ij
    
    cos_fst = cos(atan(Slo_top)).*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);
    cos_fst(cos_fst<0)=0;
    
    % Radiative parameters
    %----------------------------------------------------------------------
    SAB1(sin(h_S) <= 0.10)  =  0; %The threshold in the distributed version is 0?
    SAB2(sin(h_S) <= 0.10)  =  0;
    SAD1(sin(h_S) <= 0.10)  =  0;
    SAD2(sin(h_S) <= 0.10)  =  0;
    PARB(sin(h_S) <= 0.10)  =  0;
    PARD(sin(h_S) <= 0.10)  =  0;

    SAD1 = SAD1.*SvF + Ct.*rho_g.*((SAB1./sin(h_S)).*cos_fst + (1-SvF).*SAD1);
    SAD2 = SAD2.*SvF + Ct.*rho_g.*((SAB2./sin(h_S)).*cos_fst + (1-SvF).*SAD2);
    PARD = PARD.*SvF + Ct.*rho_g.*((PARB./sin(h_S)).*cos_fst + (1-SvF).*PARD);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SAB1 =(SAB1./sin(h_S)).*cos_fst.*ShF;
    SAB2 =(SAB2./sin(h_S)).*cos_fst.*ShF;
    PARB = (PARB./sin(h_S)).*cos_fst.*ShF;

    %corrections, temporary
    SAB1(SAB1<0)=0; SAB2(SAB2<0)=0;
    PARB(PARB<0)=0; PARD(PARD<0)=0;
    SAD1(SAD1<0)=0; SAD2(SAD2<0)=0;

    SAB1(isnan(SAB1))= 0; SAB2(isnan(SAB2)) = 0;
    SAD1(isnan(SAD1))= 0; SAD2(isnan(SAD2)) = 0;
    PARB(isnan(PARB))= 0; PARD(isnan(PARD)) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   SvF=1; %% Sky View factor = 1 if topography is not considered
   Slo_top_S = MASKn*0;
   Slo_top=Slo_top_S(ij);
end


%% SAVING INITIAL CONDITIONS AND PARAMETERS 
%==========================================================================
% In MultiPoint analysis INIT_COND_Tiber_MultiPoint.mat is created.
% This can cause problems with the initial conditions. Specially with Ca.
%==========================================================================
% (run this only once in MultiPoint model!)
%if exist(out, 'file') == 2
%else

if topo == 1
INIT_COND_v4(num_cell,m_cell,n_cell,...
   cc_max,ms_max,md_max,...
   MASKn,GLH,Ca,SNOWD,SNOWALB,out, II, ij, TT_par,ksv,vi);
load(out);
else 
    %Warning this does not look correct now....
    disp('Check init cond');
INIT_COND_v4(num_cell,m_cell,n_cell,...
   cc_max,ms_max,md_max,...
   MASKn,GLH,Slo_top_S,ksv,Ca,SNOWD,SNOWALB,out, TT_par);
load(out);
end

%% Debugger
%disp(Nreserve)

%% RUN MODEL
%==========================================================================
% PARAM_IC: Define parameter file
% MAIN_FRAME: Contains the model
%==========================================================================
%PARAM_IC = strcat(Directories.model,'/TC_setups/Shallap/Shallap_MOD_PARAM_Multipoint.m');

[aR,             Zs,             EvL_Zs,       Inf_Zs,     Bio_Zs,      Zinf, ...
 RfH_Zs,         RfL_Zs,         dz,           Ks_Zs,      Dz,          ms, ...
 Kbot,           Krock,          zatm,         Color_Class,  OM_H,       OM_L,        PFT_opt_H, ...
 PFT_opt_L,      d_leaf_H,       d_leaf_L,     SPAR,       Phy,         Soil_Param, ...
 Interc_Param,   SnowIce_Param,  VegH_Param,   VegL_Param, fpr,         VegH_Param_Dyn, ...
 VegL_Param_Dyn, Stoich_H,       aSE_H,        Stoich_L,   aSE_L,       fab_H,...
 fbe_H,          fab_L,          fbe_L,        ZR95_H,     ZR95_L,      ZRmax_H,       ZRmax_L,...
 In_max_urb,     In_max_rock,    K_usle,         Urb_Par,      Deb_Par,    Zs_deb,      Sllit,...
 Kct,            ExEM,           ParEx_H,      Mpar_H,     ParEx_L,     Mpar_L] = ....
                        PARAMETERS_SOIL_DEV_SHORT_RP( ...
 ksv(ij),        Psan,       Pcla,     Porg,   DEB,  md_max,...
 Afirn(ij),      SOIL_TH(ij),    VPAR,        TT_par,       OPT_PARAM_vals, Zs, ms_max);

%These variables were in MOD_PARAM_Multipoint, now given here 
Asur=(1./cos(atan(Slo_top))); %% Effective Area / Projected Area
Ared = 1; % Reduction due to soil rock content 
aTop = 1000*cellsize^2./cellsize; %  Ratio betweeen Area/ContourLenght [mm]
SN = 0; % Stream Identifier
Slo_pot=zeros(1,ms); % Slope of hydraulic head [fraction dy/dx]
% Rainfall disaggregation information 
a_dis = NaN ;
pow_dis = NaN;

%Ok so previously, some initial conditions were overrun in
%MOD_PARAM_Multipoint, esp. snow_albedo (=0.6)

%% DEBUGGER
%Deleting initial conditions 
%% DEBUGGER
%------------------------------------------------------------------
%[user, sys] = memory; % Windows only
%disp(['Mem used by the Starter: ', num2str(user.MemUsedMATLAB/1e6), ' MB'])
%clear An_H_t An_L_t O_t Tdp_t V_t Psi_I_H_t Psi_l_L_t Psi_l_H_t ...
%    Psi_x_H_t Psi_x_L_t Rdark_H_t Rdark_L_t Tdp_H_t Tdp_L_t
%------------------------------------------------------------------

MAIN_FRAME_Pro; % Launch the main frame of T&C. Most of the things happen in this line of code %Was MAIN_FRAM_Pro ? What's the difference?

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

%% Post-compute calculations
%==========================================================================
% Negative evapotranspiration (ET) is dew
%==========================================================================
% post-compute sublimation from ESN
SSN = ESN.*(Ts<0);
ET(ET < 0) = 0; 

%% REMOVAL OF PARAMETERS TO FREE SPACE -SIZE OF ELEMENT ARE IMPORTANT
% clear variables not needed any more
%--------------------------------------------------------------------------
all_vars = who;
% 
% vars_to_keep = {'Date', 'EICE','ESN','SND','SWE', 'Ta','Ws','U','N', 'Ts', ...
%     'SAD1','SAD2','SAB1','SAB2','Pre','Pr','Pr_sno', 'ALB','Smelt','Imelt','SSN', ...
%     'ICE','ET', 'ETen' ,'QE', 'ros', 'NDVI', 'T_H', 'T_L', 'O', 'FROCK', ...
%     'Directories', 'IniCond', 'Point', 'POI', 'ksv', 'dateRun', ...
%     'Lat','Lon','Zbas','dbThick', 'SnowIce_Param', 'Deb_Par','id_location','output_daily', ...
%     'LAI_H', 'LAI_L', 'NPP_H', 'NPP_L', 'Ccrown'};

vars_to_keep = {'Date', 'EICE','ESN','SND','SWE', 'Ta','Ws','U','N', 'Ts', 'ea','Rn',...
    'SAD1','SAD2','SAB1','SAB2','Pre','Pr','Pr_sno','Pr_liq', 'ALB','Smelt','Imelt','SSN', ...
    'ICE','ICE_D','ET', 'ETen' ,'QE', 'ros', 'NDVI', 'T_H', 'T_L', 'O', 'FROCK', ...
    'Directories', 'IniCond', 'Point', 'POI', 'ksv', 'dateRun', ...
    'Lat','Lon','Zbas','dbThick', 'SnowIce_Param', 'Deb_Par','id_location','output_daily', ...
    'LAI_H', 'LAI_L', 'NPP_H', 'NPP_L','Ccrown','Csno','B_L','B_H','LAI_L','LAI_H','NPP_L','NPP_H','An_L','An_H',...
    'Bfac_dayL','Bfac_dayH','EG','EIn_L','EIn_H','EIn_rock','H','HV','G','In_L','In_H','In_SWE','In_rock'...
    'OF','OH','OL','OS','PHE_S_L','PHE_S_H','Qfm','Qv','RA_L','RA_H','Rg_L','Rg_H','Rh','SAI_L','SAI_H','T_L','T_H'...
    'Tdp','Tdp_L','Tdp_H','Ts','Tice','TsV','ZWT','dQ','dQVEG','IP_wc'};

%Note don't have Rsw_out, L_out, Rsw, Latm=N, 

%Save as a mat file
save(strcat(Directories.save, '1_Hourly/',id_location,'_all.mat'),vars_to_keep{:});


vars_to_delete = setdiff(all_vars, vars_to_keep);
for i = 1:length(vars_to_delete)
    clear(vars_to_delete{i});
end

%% Output manager
%==========================================================================
%{
As outputs:
    Hourly outputs
    Daily outputs
    Parameters
%}
%==========================================================================

% Make date usable for R
Date_R = char(Date); 

%% Parameter output
Param_t = table(Lat,Lon,Zbas,dbThick,'VariableNames',{'Lat','Lon','Zbas','dbThick'});
Param_t = [Param_t, struct2table(SnowIce_Param), struct2table(Deb_Par)];
Param_t = rows2vars(Param_t);
Param_t = renamevars(Param_t,{'OriginalVariableNames','Var1'},{'Parameter','Value'});

% Exporting as .txt
writetable(Param_t, strcat(Directories.save, '3_Params/',id_location,'_param.txt') )

% Here I manually choose the T&C outputs I want to save at each point.
Outputs_t = table(Date, ...
    EICE,ESN,SND,SWE,...
    Ta,Ws,U,N, ...
    SAD1+SAD2+SAB1+SAB2,Pre,Pr,Pr_sno, ...
    ALB,Smelt,Imelt,SSN, ...
    ICE,ET, ETen ,QE, ...
    ros, NDVI, T_H, T_L, ...
    O, FROCK, ...
'VariableNames',{'Date', ...
    'EICE','ESN','SND','SWE',...
    'Ta','Ws','U','N', ...
    'Rsw', 'Pre','Pr','Pr_sno', ...
    'Albedo','Smelt','Imelt','SSN', ...
    'ICE','ET','ETen','QE', ...
    'ros', 'NDVI', 'T_H', 'T_L', ...
    'O_mois','FROCK'});

%% Hourly output
writetable(Outputs_t, strcat(Directories.save, '1_Hourly/',id_location,'_hourly_results.txt'))

% Save mat file of table
%--------------------------------------------------------------------------
save(strcat(Directories.save, '1_Hourly/',id_location,'_hourly_results.mat'),'Outputs_t');

%% Daily outputs
% If daily outputs are activated


if output_daily == 1

% Timetable creation
%--------------------------------------------------------------------------
Outputs_tt = table2timetable(Outputs_t); 


% Mean of variables
%--------------------------------------------------------------------------
Outputs_dm = retime(Outputs_tt,'daily',@nanmean); % For average variables
% Initial matrix set based on averages
Outputs_d = Outputs_dm; 


% Some variables must be summed
%--------------------------------------------------------------------------
Outputs_ds = retime(Outputs_tt,'daily',@nansum); % For average variables

% Changing variables that must be summed in Outputs_d
Outputs_d.Pr = Outputs_ds.Pr;
Outputs_d.Pr_sno = Outputs_ds.Pr_sno;
Outputs_d.ET = Outputs_ds.ET;
Outputs_d.ETen = Outputs_ds.ETen;
Outputs_d.T_H = Outputs_ds.T_H;
Outputs_d.T_L = Outputs_ds.T_L;



% Average of NDVI must be calculated without zeros during the night
%--------------------------------------------------------------------------
% 1. Create a logical index for rows where Rsw is positive and not NaN
idx_positive_Rsw = Outputs_tt.Rsw > 0 & ~isnan(Outputs_tt.Rsw);

% 2. Filter the timetable based on this index
Filtered_Outputs_tt = Outputs_tt(idx_positive_Rsw, :);

% It filters out zeros and NaNs before calculating the mean.
customMeanNoZeros = @(x) mean(x(~isnan(x)), 'omitnan');

table_ndvi = retime(Filtered_Outputs_tt(:,{'NDVI'}),'daily',customMeanNoZeros); % For sum variables


% Adding daily outputs
%--------------------------------------------------------------------------
Outputs_d.LAI_H = sum(LAI_H.*Ccrown,2); % Leaf area index
Outputs_d.LAI_L = sum(LAI_L.*Ccrown,2);
Outputs_d.LAI = Outputs_d.LAI_H+Outputs_d.LAI_L;

Outputs_d.NPP_H = sum(NPP_H.*Ccrown,2); % Net primary production
Outputs_d.NPP_L = sum(NPP_L.*Ccrown,2);
Outputs_d.NPP = Outputs_d.NPP_H + Outputs_d.NPP_L;

Outputs_d.NDVI = table_ndvi.NDVI;


% Change of format
%--------------------------------------------------------------------------
Outputs_t = timetable2table(Outputs_d);


% Exporting as .txt
%--------------------------------------------------------------------------
writetable(Outputs_t, strcat(Directories.save, '2_Daily/',id_location,'_daily_results.txt'))

% Save mat file of table
%--------------------------------------------------------------------------
save(strcat(Directories.save, '2_Daily/',id_location,'_daily_results.mat'),'Outputs_t');

end 


%% Save of the environment
% Here only for one point - For the snow station. More points can be
% defined
% if strcmp(POI{POI.Name == Point,'Feature'}{1}, 'Snow_station')
% save(strcat(Directories.save, '5_Env/Env_',id_location,'.mat'))
% end

%% Memory
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
disp(['Memory used by variables at the end of run_Point_Pro_BC: ', num2str(totalSize/1e6), ' MB']);

%% OTHER CALCULATIONS OUT OF THE MODEL
%==========================================================================
%
%==========================================================================
%{
% Resample the DEM to a new resolution (e.g., 30 meters)
DEM_resampled_out = resample(DEM, 500); 

% Display the original and resampled DEMs
imagesc(DEM);
figure;
imagesc(DEM_resampled);
%}

end