%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PARAMETERS SOILS AND VEGETATION  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% --> II (define high and low vegetation
%%% Fir (H) / Larch (H) / Grass C3 (L) / Shrub Winter Dec. (L) 
%%% zatm_surface - instrument height above shrub, grass meadow and rock

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[aR,              EvL_Zs,      Inf_Zs,     Bio_Zs, ...
         Zinf,            RfH_Zs,     RfL_Zs,      dz,         Ks_Zs, ...
         Dz,              ms,         Kbot,        Krock,      zatm,...
         Color_Class,     OM_H,       OM_L,        PFT_opt_H,  PFT_opt_L, ...
         d_leaf_H,        d_leaf_L,   SPAR,Phy,    Soil_Param, Interc_Param, ...
         SnowIce_Param,   VegH_Param, VegL_Param,  fpr,        VegH_Param_Dyn, ...
         VegL_Param_Dyn,  Stoich_H,   aSE_H,       Stoich_L,   aSE_L, ...
         fab_H,           fbe_H,      fab_L,       fbe_L,      ZR95_H, ...
         ZR95_L,          ZRmax_H,    ZRmax_L,     In_max_urb, In_max_rock, K_usle,     Urb_Par, ...
         Deb_Par,         Zs_deb,     Sllit,       Kct,        ExEM, ...
         ParEx_H,         Mpar_H,     ParEx_L,     Mpar_L] ...
                     =PARAMETERS_SOIL_DEV( ...
         cell_class,      Psan,       Pcla,        Porg,       dbThick, ...
         md_max,          Afirn,      Soil_th,     VPAR,        TT_vals, ...
         OPT_PARAM_vals,  Zs,         ms)

%% SOIL PARAMETERS
%==========================================================================

% Unpack params from OPT_PARAM_vals (explicit, replaces v2struct which is slow:
% a parameter added to Options_and_NoVegParams.xlsx and needed here must also be added below)
Zdes           = OPT_PARAM_vals.Zdes;
Zinf           = OPT_PARAM_vals.Zinf;
Zbio           = OPT_PARAM_vals.Zbio;
SPAR           = OPT_PARAM_vals.SPAR;
Kfc            = OPT_PARAM_vals.Kfc;
Phy            = OPT_PARAM_vals.Phy;
In_max_urb     = OPT_PARAM_vals.In_max_urb;
In_max_rock    = OPT_PARAM_vals.In_max_rock;
ros_max1       = OPT_PARAM_vals.ros_max1;
ros_max2       = OPT_PARAM_vals.ros_max2;
Th_Pr_sno      = OPT_PARAM_vals.Th_Pr_sno;
Ice_wc_sp      = OPT_PARAM_vals.Ice_wc_sp;
ros_Ice_thr    = OPT_PARAM_vals.ros_Ice_thr;
WatFreez_Th    = OPT_PARAM_vals.WatFreez_Th;
dz_ice         = OPT_PARAM_vals.dz_ice;
CASE_ROOT      = OPT_PARAM_vals.CASE_ROOT;
Kct            = OPT_PARAM_vals.Kct;
gcI            = OPT_PARAM_vals.gcI;
KcI            = OPT_PARAM_vals.KcI;
Sp_SN_In       = OPT_PARAM_vals.Sp_SN_In;
Sllit          = OPT_PARAM_vals.Sllit;
fpr            = OPT_PARAM_vals.fpr;
aR             = OPT_PARAM_vals.aR;
Color_Class    = OPT_PARAM_vals.Color_Class;
Urb_Par_alb    = OPT_PARAM_vals.Urb_Par_alb;
Urb_Par_e_sur  = OPT_PARAM_vals.Urb_Par_e_sur;
Urb_Par_BuildH = OPT_PARAM_vals.Urb_Par_BuildH;

%Where:
%Zdes % Depth of evaporation layer [mm]
%Zinf % Depth of infiltration layer (=first layer) [mm]
%Zbio % Depth of the active Biogeochemistry zone [mm]
%SPAR %%% SOIL PARAMETER TYPE (1-VanGenuchten 2-Saxton-Rawls)
%Kfc %Soil matric potential [mm/h]
%Phy %Soil matric potential [kPa]
%In_max_urb % [mm] maximum interception capacity in urban
%In_max_rock %% [mm] Maximum interception capacity in rocks
%ros_max1  %600; %%% [kg/m^3] - Snow density
%ros_max2  %450; %%% [kg/m^3] - Snow density
%Th_Pr_sno [mm/day] Threshold Intensity of snow to consider a New SnowFall - Shallap point value
%Ice_wc_sp [-] Specific Maximum water content ice
%ros_Ice_thr [kg/m^3] Density Thrshold to transform snow into ice
%WatFreez_Th  [?C] Threshold for freezing lake water
%dz_ice Water Freezing Layer progression without snow-layer
%CASE_ROOT % Type of Root Profile
%Kct  %%% Factor Vegetation Cover --- for throughfall
%gcI  %%% [1/mm] % Interception Parameter
%KcI %%%% [mm] -- Mahfouf and Jacquemin 1989 % Interception Parameter
%Sp_SN_In [mm/LAI] Interception Parameter
%Sllit  % Litter Specific Leaf area [m2 Litter / kg DM]
% fpr Unknown
% aR Anisotrophy ration
% Color_Class Soil Class Parameters from Oleson et al., 2010 (updated color class)
% Urb_Par_alb Urban Albedo
% Urb_Par_e_sur urban surface emissivity
% Urb_Par_BuildH Urban BuildH

% Row number of the class
cell_class = find(cell_class == VPAR.Class,1);

%  One value per class
Kbot  = VPAR.Kbot(cell_class); % Conductivity at the bedrock layer [mm/h] 
Krock =  VPAR.Krock(cell_class); % Conductivity of Fractured Rock [mm/h] 

% Warning for soil
%--------------------------------------------------------------------------
if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    return
end

[EvL_Zs]=Evaporation_layers(Zs,Zdes); % Evaporation Layer fraction
[Inf_Zs]=Evaporation_layers(Zs,Zinf); % Infiltration Depth Layer fraction
[Bio_Zs]=Evaporation_layers(Zs,Zbio);
dz= diff(Zs); % Thickness of the Layers [mm]
Dz=zeros(1,ms);

for ii = 1:ms
    if ii>1
        Dz(ii)= (dz(ii)+ dz(ii-1))/2; %%% Delta Depth Between Middle Layer  [mm]
    else
        Dz(ii)=dz(1)/2; %%% Delta Depth Between First Middle Layer and soil surface [mm]
    end
end

%% VEGETATION

%C values now given as grids after Parameters Soil
cc      = length(VPAR.Ccrowns{cell_class}); %Note this was POI.Crown and is used in case there are multiple Crowns per landcover
II      = VPAR.Veg_type{cell_class}; %Name of cell type

% ===== Vegetation params of this class =========
% TT_vals{k} is built once in the launcher from TT_par, for row k of VPAR.Class
TP = TT_vals{cell_class};
if isempty(TP)
    error('PARAMETERS_SOIL_DEV:noVegParams', 'No vegetation parameters could be built from TT_par for class %d (%s)', VPAR.Class(cell_class), strjoin(II, ', '));
end
% Unpack (explicit, replaces v2struct: a parameter added to Parameters_TC.xlsx
% and needed here must also be added below)
zatm_surface = TP.zatm_surface;
ExEM         = TP.ExEM;
ZR95_H       = TP.ZR95_H;
ZR95_L       = TP.ZR95_L;
ZR50_H       = TP.ZR50_H;
ZR50_L       = TP.ZR50_L;
ZRmax_H      = TP.ZRmax_H;
ZRmax_L      = TP.ZRmax_L;
Sp_LAI_H_In  = TP.Sp_LAI_H_In;
Sp_LAI_L_In  = TP.Sp_LAI_L_In;
d_leaf_H     = TP.d_leaf_H;
d_leaf_L     = TP.d_leaf_L;
KnitH        = TP.KnitH;
KnitL        = TP.KnitL;
mSl_H        = TP.mSl_H;
mSl_L        = TP.mSl_L;
FI_H         = TP.FI_H;
Do_H         = TP.Do_H;
a1_H         = TP.a1_H;
go_H         = TP.go_H;
CT_H         = TP.CT_H;
DSE_H        = TP.DSE_H;
Ha_H         = TP.Ha_H;
gmes_H       = TP.gmes_H;
rjv_H        = TP.rjv_H;
FI_L         = TP.FI_L;
Do_L         = TP.Do_L;
a1_L         = TP.a1_L;
go_L         = TP.go_L;
CT_L         = TP.CT_L;
DSE_L        = TP.DSE_L;
Ha_L         = TP.Ha_L;
gmes_L       = TP.gmes_L;
rjv_L        = TP.rjv_L;
Vmax_H       = TP.Vmax_H;
Vmax_L       = TP.Vmax_L;
Psi_sto_00_H = TP.Psi_sto_00_H;
Psi_sto_50_H = TP.Psi_sto_50_H;
PsiL00_H     = TP.PsiL00_H;
PsiL50_H     = TP.PsiL50_H;
Kleaf_max_H  = TP.Kleaf_max_H;
Cl_H         = TP.Cl_H;
Axyl_H       = TP.Axyl_H;
Kx_max_H     = TP.Kx_max_H;
PsiX50_H     = TP.PsiX50_H;
Cx_H         = TP.Cx_H;
Psi_sto_00_L = TP.Psi_sto_00_L;
Psi_sto_50_L = TP.Psi_sto_50_L;
PsiL00_L     = TP.PsiL00_L;
PsiL50_L     = TP.PsiL50_L;
Kleaf_max_L  = TP.Kleaf_max_L;
Cl_L         = TP.Cl_L;
Axyl_L       = TP.Axyl_L;
Kx_max_L     = TP.Kx_max_L;
PsiX50_L     = TP.PsiX50_L;
Cx_L         = TP.Cx_L;
PsiG50_H     = TP.PsiG50_H;
PsiG99_H     = TP.PsiG99_H;
gcoef_H      = TP.gcoef_H;
PsiG50_L     = TP.PsiG50_L;
PsiG99_L     = TP.PsiG99_L;
gcoef_L      = TP.gcoef_L;
OPT_PROP_H   = TP.OPT_PROP_H;
OPT_PROP_L   = TP.OPT_PROP_L;
OM_H         = TP.OM_H;
OM_L         = TP.OM_L;
aSE_H        = TP.aSE_H;
Sl_H         = TP.Sl_H;
Nl_H         = TP.Nl_H;
r_H          = TP.r_H;
gR_H         = TP.gR_H;
Tcold_H      = TP.Tcold_H;
age_cr_H     = TP.age_cr_H;
Trr_H        = TP.Trr_H;
LtR_H        = TP.LtR_H;
eps_ac_H     = TP.eps_ac_H;
fab_H        = TP.fab_H;
ff_r_H       = TP.ff_r_H;
Wm_H         = TP.Wm_H;
dd_max_H     = TP.dd_max_H;
dc_C_H       = TP.dc_C_H;
drn_H        = TP.drn_H;
dsn_H        = TP.dsn_H;
Mf_H         = TP.Mf_H;
Klf_H        = TP.Klf_H;
Bfac_lo_H    = TP.Bfac_lo_H;
Bfac_ls_H    = TP.Bfac_ls_H;
Tlo_H        = TP.Tlo_H;
Tls_H        = TP.Tls_H;
PAR_th_H     = TP.PAR_th_H;
dmg_H        = TP.dmg_H;
LAI_min_H    = TP.LAI_min_H;
mjDay_H      = TP.mjDay_H;
LDay_min_H   = TP.LDay_min_H;
LDay_cr_H    = TP.LDay_cr_H;
aSE_L        = TP.aSE_L;
Sl_L         = TP.Sl_L;
Nl_L         = TP.Nl_L;
r_L          = TP.r_L;
gR_L         = TP.gR_L;
Tcold_L      = TP.Tcold_L;
dd_max_L     = TP.dd_max_L;
dc_C_L       = TP.dc_C_L;
drn_L        = TP.drn_L;
dsn_L        = TP.dsn_L;
Mf_L         = TP.Mf_L;
Klf_L        = TP.Klf_L;
age_cr_L     = TP.age_cr_L;
Trr_L        = TP.Trr_L;
LtR_L        = TP.LtR_L;
Wm_L         = TP.Wm_L;
eps_ac_L     = TP.eps_ac_L;
fab_L        = TP.fab_L;
ff_r_L       = TP.ff_r_L;
Bfac_lo_L    = TP.Bfac_lo_L;
Bfac_ls_L    = TP.Bfac_ls_L;
Tlo_L        = TP.Tlo_L;
Tls_L        = TP.Tls_L;
PAR_th_L     = TP.PAR_th_L;
dmg_L        = TP.dmg_L;
LAI_min_L    = TP.LAI_min_L;
mjDay_L      = TP.mjDay_L;
LDay_min_L   = TP.LDay_min_L;
LDay_cr_L    = TP.LDay_cr_L;

%Where:
%zatm = reference height
%-- Vegetation parameters
%ExEM, ZR95_H, ZR95_L, ZR50_H, ZR50_L, ZRmax_H, ZRmax_L
% -- Interception Parameter
%Sp_LAI_H_In, Sp_LAI_L_In 
%-- Leaf Dimension
%d_leaf_H, d_leaf_L
% ==== Biochemical parameterS
%-- Veg Biochemical parameter
%KnitH, KnitL, mSl_H ([m2 PFT /gC]  Linear increase in Sla with LAI), ( %
%[m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree)
%-- Photosynthesis Parameter
% High vegetation
%FI_H, Do_H, a1_H, go_H, CT_H, DSE_H, Ha_H, gmes_H, rjv_H 
% Low vegetation
%FI_L, Do_L, a1_L, go_L , CT_L, DSE_L , Ha_L, gmes_L, rjv_L, Vmax_H, Vmax_L 
%-- Hydraulic Parameters
%Psi_sto_00_H, Psi_sto_50_H 
%-- Leaf
%PsiL00_H, PsiL50_H, Kleaf_max_H, Cl_H
%-- Xylem
%Axyl_H, Kx_max_H, PsiX50_H, Cx_H 
%-- Stomata
%Psi_sto_00_L, Psi_sto_50_L
%-- Leaf
%categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A    grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
%PsiL00_L, PsiL50_L, Kleaf_max_L, Cl_L 
%-- Xylem
%Axyl_L, Kx_max_L, PsiX50_L, Cx_L 
%--- Growth Parameters
%-- High vegetation
%categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R    grass_A   grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
%PsiG50_H, PsiG99_H, gcoef_H 
%-- Low vegetation
%PsiG50_L, PsiG99_L, gcoef_L 
% ==== Vegetation Optical Parameter
%OPT_PROP_H, OPT_PROP_L, OM_H, OM_L 
% VEGETATION PART
% High Vegetation
%BLever: broadleaf evergreen vegetation dec;
%BLdec: broadleaf deciduous vegetation dec
%aSE_H, Sl_H, Nl_H, r_H, gR_H, Tcold_H, age_cr_H, Trr_H, LtR_H, eps_ac_H,fab_H,ff_r_H,Wm_H, dd_max_H
%dc_C_H, drn_H, dsn_H, Mf_H, Klf_H 
%-- Phenology 
%Bfac_lo_H, Bfac_ls_H, Tlo_H, Tls_H, PAR_th_H, dmg_H 
%LAI_min_H, mjDay_H, LDay_min_H, LDay_cr_H 
% LOW VEGETATION
%aSE_L, Sl_L, Nl_L, r_L, gR_L, Tcold_L, dd_max_L, dc_C_L, drn_L, dsn_L, Mf_L, Klf_L 
%age_cr_L, Trr_L, LtR_L, Wm_L, eps_ac_L, fab_L, ff_r_L 
%-- Phenology 
%Bfac_lo_L, Bfac_ls_L, Tlo_L, Tls_L, PAR_th_L, dmg_L, LAI_min_L, mjDay_L, LDay_min_L, LDay_cr_L 

%========= Elevation ======
zatm = max(zatm_surface); %Takes the highest if multiple land covers

%=======  Soil parameters ======
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);
%%%%%
rsd=rsd*ones(1,ms);
lan_dry=lan_dry*ones(1,ms);
lan_s =lan_s*ones(1,ms);
cv_s = cv_s*ones(1,ms);
%%%
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%%%%%
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]

%%%%%%%%%%%%%%%% Define soil layer depth by specifying 
%%%%%%%%%%%%%%%% impermeable layer at given layer
%Only use if a %...
%Soil_th=Zs(end)*(Soil_th/100); %% convert from relative to absolute depth
[~, ix]= min(abs(Zs-Soil_th));
if ix==1
    ix=2;
end
Ks_Zs(ix-1:end)=1e-12; 
clear ix

%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
clear Oss_Lp Owp_Lp
% force field capacity (Ofc) to equal Osat for impermeable layers
Ofc(Ks_Zs<0.2)=Osat(Ks_Zs<0.2);

Aice = Afirn; %% [-] Ice albedo Use standard albedo

% Root Parameters (Function)
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);

for i=1:cc
    %%%%%%%% Vegetation Optical Parameter
    [PFT_opt_H(i)]=Veg_Optical_Parameter(OPT_PROP_H(i));
    [PFT_opt_L(i)]=Veg_Optical_Parameter(OPT_PROP_L(i));
end

fbe_H = 1-fab_H; %% fraction below-ground sapwood and reserve

% Phenology 
%--------------------------------------------------------------------------

for i=1:cc
    [Stoich_H(i)]=Veg_Stoichiometric_Parameter(Nl_H(i));
    [ParEx_H(i)]=Exudation_Parameter(0);
    [Mpar_H(i)]=Vegetation_Management_Parameter;
end

fbe_L = 1-fab_L; %% fraction below-ground sapwood and reserve

% Phenology 
%--------------------------------------------------------------------------
for i=1:cc
    [Stoich_L(i)]=Veg_Stoichiometric_Parameter(Nl_L(i));
    [ParEx_L(i)]=Exudation_Parameter(0);
    [Mpar_L(i)]=Vegetation_Management_Parameter;
end

%% ORGANISING PARAMETERS
%Was restating parameters

%%% Soil Structural Parameters
if not(isfield(TP,'Omac'))
    if SPAR == 3
        disp('NOT SPECIFIED SOIL STRUCTURE PARAMETERS')
        return
    end
    nVGM = 3.0*ones(1,ms); %  %% [-]
    alpVGM= 30*alpVG; %% [1/mm]
    Omac = 0.00*Osat; % 0.021*Osat;
    Ks_mac = 100*Ks_Zs; %%[mm/h]
    lVGM =0.5*ones(1,ms); %  %% [-]
else
    nVGM = TP.nVGM; alpVGM = TP.alpVGM; Omac = TP.Omac; Ks_mac = TP.Ks_mac; lVGM = TP.lVGM;
end
%%% Soil VG modifications
if not(isfield(TP,'bVG'))
    if SPAR == 1
        disp('NOT SPECIFIED CORRECTION FOR VG')
        return
    end
    s_SVG =zeros(1,ms);
    bVG = zeros(1,ms);
else
    s_SVG = TP.s_SVG; bVG = TP.bVG;
end
if not(isfield(TP,'lVG'))
    lVG =0.5*ones(1,ms); %  %% [-]
else
    lVG = TP.lVG;
end
%%%%%%% Soil_Parameters
Soil_Param.Ks_Zs=Ks_Zs;
Soil_Param.Osat=Osat;
Soil_Param.Ohy=Ohy;
Soil_Param.L=L;
Soil_Param.Pe=Pe;
Soil_Param.O33=O33;
Soil_Param.Ofc=Ofc;
Soil_Param.alpVG=alpVG;
Soil_Param.nVG=nVG;
Soil_Param.Ks_mac=Ks_mac;
Soil_Param.Omac=Omac;
Soil_Param.alpVGM=alpVGM;
Soil_Param.nVGM=nVGM;
Soil_Param.rsd=rsd;
Soil_Param.lan_dry=lan_dry ;
Soil_Param.lan_s=lan_s;
Soil_Param.cv_s=cv_s;
Soil_Param.s_SVG=s_SVG;
Soil_Param.bVG=bVG;
Soil_Param.lVG=lVG;
Soil_Param.lVGM=lVGM;
Soil_Param.Phy = Phy;

%%%%%%%%%% Interception Parameters    
Interc_Param.Sp_SN_In=Sp_SN_In;
Interc_Param.Sp_LAI_L_In=Sp_LAI_L_In;
Interc_Param.Sp_LAI_H_In=Sp_LAI_H_In;
Interc_Param.gcI=gcI;
Interc_Param.KcI=KcI;

%%%%%% Snow-Ice Parameters
SnowIce_Param.WatFreez_Th=WatFreez_Th;
SnowIce_Param.dz_ice=dz_ice;
SnowIce_Param.Th_Pr_sno=Th_Pr_sno;
SnowIce_Param.ros_max1=ros_max1;
SnowIce_Param.ros_max2=ros_max2;
SnowIce_Param.Ice_wc_sp=Ice_wc_sp;
SnowIce_Param.ros_Ice_thr=ros_Ice_thr;
SnowIce_Param.Aice=Aice;

%%%%%%%%% Debris-Cover properties
Deb_Par.alb = OPT_PARAM_vals.Deb_Par_alb; 
Deb_Par.e_sur =  OPT_PARAM_vals.Deb_Par_e_sur;
Deb_Par.lan = OPT_PARAM_vals.Deb_Par_lan; 
Deb_Par.rho = OPT_PARAM_vals.Deb_Par_rho;  % [kg/m^3]
Deb_Par.cs = OPT_PARAM_vals.Deb_Par_cs;  % [J/kg K]
Deb_Par.zom = OPT_PARAM_vals.Deb_Par_zom; 
Deb_Par.min_th = OPT_PARAM_vals.Deb_Par_min_th; 

Min_layer_thickness = Deb_Par.min_th - 5;
if Min_layer_thickness <= 0
    error('Minimum debris thickness too small!')
end

if dbThick>0
    if dbThick < Deb_Par.min_th
        dbThick = Deb_Par.min_th; % Set debris thickness to be at least equal to the minimum value
    end
    nst=md_max;
    k=(dbThick/(Min_layer_thickness))^(1/(nst-1));
    Zs_deb = [0 (Min_layer_thickness)*k.^(1:nst-1)]; %% [mm]
    %clear nst k
else
    Zs_deb = [];
end

%Urban Parameters
Urb_Par.alb = Urb_Par_alb;
Urb_Par.e_sur= Urb_Par_e_sur;
Urb_Par.BuildH= Urb_Par_BuildH;

%%%%% Vegetation High Parameters
VegH_Param.KnitH=KnitH;
VegH_Param.mSl_H=mSl_H;
VegH_Param.Sl_H=Sl_H;
VegH_Param.CT_H=CT_H ;
VegH_Param.Vmax_H=Vmax_H;
VegH_Param.FI_H=FI_H;
VegH_Param.a1_H=a1_H;
VegH_Param.go_H=go_H;
VegH_Param.DSE_H=DSE_H;
VegH_Param.Ha_H=Ha_H;
VegH_Param.Do_H=Do_H;
VegH_Param.gmes_H=gmes_H;
VegH_Param.rjv_H=rjv_H;
VegH_Param.Psi_sto_50_H=Psi_sto_50_H;
VegH_Param.Psi_sto_00_H=Psi_sto_00_H;
VegH_Param.Axyl_H=Axyl_H;
VegH_Param.PsiL50_H=PsiL50_H;
VegH_Param.PsiL00_H=PsiL00_H;
VegH_Param.Kleaf_max_H=Kleaf_max_H;
VegH_Param.Cl_H=Cl_H;
VegH_Param.Kx_max_H=Kx_max_H;
VegH_Param.PsiX50_H=PsiX50_H;
VegH_Param.Cx_H=Cx_H;

if not(isfield(TP,'Osm_reg_Max_H'))
    VegH_Param.Osm_reg_Max_H = 0*ones(1,cc);
    VegH_Param.eps_root_base_H=  0.9*ones(1,cc);
else
    VegH_Param.Osm_reg_Max_H = TP.Osm_reg_Max_H;   
    VegH_Param.eps_root_base_H=  TP.eps_root_base_H; 
end

%%%%% Vegetation Low Parameters
VegL_Param.KnitL=KnitL;
VegL_Param.mSl_L=mSl_L;
VegL_Param.Sl_L=Sl_L;
VegL_Param.CT_L=CT_L;
VegL_Param.Vmax_L=Vmax_L;
VegL_Param.FI_L=FI_L;
VegL_Param.a1_L=a1_L;
VegL_Param.go_L=go_L;
VegL_Param.DSE_L=DSE_L;
VegL_Param.Ha_L=Ha_L;
VegL_Param.Do_L=Do_L;
VegL_Param.gmes_L=gmes_L;
VegL_Param.rjv_L=rjv_L;
VegL_Param.Psi_sto_50_L = Psi_sto_50_L;
VegL_Param.Psi_sto_00_L=Psi_sto_00_L;
VegL_Param.Axyl_L=Axyl_L;
VegL_Param.PsiL50_L=PsiL50_L;
VegL_Param.PsiL00_L=PsiL00_L;
VegL_Param.Kleaf_max_L= Kleaf_max_L;
VegL_Param.Cl_L=Cl_L;
VegL_Param.Kx_max_L=Kx_max_L;
VegL_Param.PsiX50_L=PsiX50_L;
VegL_Param.Cx_L=Cx_L;

if not(isfield(TP,'Osm_reg_Max_L'))
    VegL_Param.Osm_reg_Max_L = 0*ones(1,cc);
    VegL_Param.eps_root_base_L=  0.9*ones(1,cc);
else
    VegL_Param.Osm_reg_Max_L = TP.Osm_reg_Max_L;   
    VegL_Param.eps_root_base_L=  TP.eps_root_base_L; 
end

%%%%% Vegetation High Parameters for vegetation dynamics
VegH_Param_Dyn.Sl = Sl_H;
VegH_Param_Dyn.mSl = mSl_H;
VegH_Param_Dyn.Stoich = Stoich_H;
VegH_Param_Dyn.r = r_H;
VegH_Param_Dyn.gR = gR_H;
VegH_Param_Dyn.LtR = LtR_H;
VegH_Param_Dyn.eps_ac = eps_ac_H;
VegH_Param_Dyn.aSE = aSE_H;
VegH_Param_Dyn.Trr = Trr_H;
VegH_Param_Dyn.dd_max = dd_max_H;
VegH_Param_Dyn.dc_C = dc_C_H;
VegH_Param_Dyn.Tcold = Tcold_H;
VegH_Param_Dyn.drn = drn_H;
VegH_Param_Dyn.dsn= dsn_H;
VegH_Param_Dyn.age_cr = age_cr_H;
VegH_Param_Dyn.Bfac_lo = Bfac_lo_H;
VegH_Param_Dyn.Bfac_ls = Bfac_ls_H;
VegH_Param_Dyn.Tlo= Tlo_H;
VegH_Param_Dyn.Tls = Tls_H;
VegH_Param_Dyn.mjDay = mjDay_H;
VegH_Param_Dyn.LDay_min =LDay_min_H;
VegH_Param_Dyn.dmg= dmg_H;
VegH_Param_Dyn.Mf = Mf_H;
VegH_Param_Dyn.Wm = Wm_H;
VegH_Param_Dyn.LAI_min =LAI_min_H;
VegH_Param_Dyn.LDay_cr= LDay_cr_H;
VegH_Param_Dyn.PsiG50 = PsiG50_H;
VegH_Param_Dyn.PsiG99 =PsiG99_H;
VegH_Param_Dyn.gcoef =gcoef_H;
VegH_Param_Dyn.fab= fab_H;
VegH_Param_Dyn.fbe = fbe_H;
VegH_Param_Dyn.Klf =Klf_H;
VegH_Param_Dyn.ff_r =ff_r_H;
VegH_Param_Dyn.PAR_th =PAR_th_H;
VegH_Param_Dyn.PsiL50 = PsiL50_H;
VegH_Param_Dyn.PsiL00 = PsiL00_H;

if not(isfield(TP,'soCrop_H'))
    VegH_Param_Dyn.soCrop = NaN*ones(1,cc);
    VegH_Param_Dyn.MHcrop=  NaN*ones(1,cc);
    VegH_Param_Dyn.Sl_emecrop=  NaN*ones(1,cc);
else
    VegH_Param_Dyn.soCrop = TP.soCrop_H;  
    VegH_Param_Dyn.MHcrop=  TP.MHcrop_H; 
    VegH_Param_Dyn.Sl_emecrop= TP.Sl_emecrop_H;
end

%%%%% Vegetation Low Parameters for vegetation dynamics
VegL_Param_Dyn.Sl = Sl_L;
VegL_Param_Dyn.mSl = mSl_L;
VegL_Param_Dyn.Stoich = Stoich_L;
VegL_Param_Dyn.r = r_L;
VegL_Param_Dyn.gR = gR_L;
VegL_Param_Dyn.LtR = LtR_L;
VegL_Param_Dyn.eps_ac = eps_ac_L;
VegL_Param_Dyn.aSE = aSE_L;
VegL_Param_Dyn.Trr = Trr_L;
VegL_Param_Dyn.dd_max = dd_max_L;
VegL_Param_Dyn.dc_C = dc_C_L;
VegL_Param_Dyn.Tcold = Tcold_L;
VegL_Param_Dyn.drn = drn_L;
VegL_Param_Dyn.dsn = dsn_L;
VegL_Param_Dyn.age_cr = age_cr_L;
VegL_Param_Dyn.Bfac_lo = Bfac_lo_L;
VegL_Param_Dyn.Bfac_ls = Bfac_ls_L;
VegL_Param_Dyn.Tlo= Tlo_L;
VegL_Param_Dyn.Tls = Tls_L;
VegL_Param_Dyn.mjDay = mjDay_L;
VegL_Param_Dyn.LDay_min =LDay_min_L;
VegL_Param_Dyn.dmg= dmg_L;
VegL_Param_Dyn.Mf = Mf_L;
VegL_Param_Dyn.Wm = Wm_L;
VegL_Param_Dyn.LAI_min =LAI_min_L;
VegL_Param_Dyn.LDay_cr= LDay_cr_L;
VegL_Param_Dyn.PsiG50 = PsiG50_L;
VegL_Param_Dyn.PsiG99 =PsiG99_L;
VegL_Param_Dyn.gcoef =gcoef_L;
VegL_Param_Dyn.fab= fab_L;
VegL_Param_Dyn.fbe = fbe_L;
VegL_Param_Dyn.Klf =Klf_L;
VegL_Param_Dyn.ff_r =ff_r_L;
VegL_Param_Dyn.PAR_th =PAR_th_L;
VegL_Param_Dyn.PsiL50 = PsiL50_L;
VegL_Param_Dyn.PsiL00 = PsiL00_L;

if not(isfield(TP,'soCrop_L'))
    VegL_Param_Dyn.soCrop = NaN*ones(1,cc);
    VegL_Param_Dyn.MHcrop=  NaN*ones(1,cc);
    VegL_Param_Dyn.Sl_emecrop=  NaN*ones(1,cc);
else
    VegL_Param_Dyn.soCrop = TP.soCrop_L;  
    VegL_Param_Dyn.MHcrop=  TP.MHcrop_L; 
    VegL_Param_Dyn.Sl_emecrop= TP.Sl_emecrop_L;
end

end




 
