clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C)ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parpool('local',10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% METEO-INPUT PARAMETER %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Date,Pr,Ta,Ds,Ws,ea,N,Pre,Tdew,SAB1,SAB2,SAD1,SAD2,PARB,PARD,Ca
%%% 1:N_time_step
TITLE_SAVE = 'SIM_Gletsch';
%%%%%%%%%%%% Data loading
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%
load('INPUTS\Data_Grimsel_RUN.mat')
load('INPUTS\Pr_dist_Gletsch.mat')
load('INPUTS\Virt_Stati_Gletsch.mat')
clear Xvs Yvs
%Output location
outlocation = 'OUTPUTS/';
%%% 26-Apr-1989 12:00:00 --> 14-Apr-2011 09:00:00
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%x1=3782; x2=181308;
x1=3782;
x2=38846;
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%5
Date=Date(x1:x2);
Pr=Pr(x1:x2); Pr_dist=Pr_dist(x1:x2,:);
Ta=Ta(x1:x2); Pre=Pre(x1:x2);
Ws=Ws(x1:x2); ea=ea(x1:x2);  SAD1=SAD1(x1:x2);
SAD2=SAD2(x1:x2); SAB1=SAB1(x1:x2);
SAB2=SAB2(x1:x2); N=N(x1:x2); Tdew=Tdew(x1:x2);esat=esat(x1:x2);
PARB=PARB(x1:x2); PARD = PARD(x1:x2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
qta_Sta=Zbas; %%[m a.s.l.]
clear Zbas;
Ds=esat-ea; %% [Pa] Vapor Pressure Deficit
Ds(Ds<0)=0;
Ws(Ws<=0)=0.01;
%%%%%%%%%%%%%%
t_bef= -0.67; t_aft= 1.67;
%%%%%%%%%%%%%%%%%%%%%%%
load('INPUTS\Ca_Data.mat');
d1 = find(abs(Date_CO2-Date(1))<1/36);d2 = find(abs(Date_CO2-Date(end))<1/36);
Ca=Ca(d1:d2);
clear d1 d2 Date_CO2
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_time_step= 35064;%  8760;
Nd_time_step = ceil(N_time_step/24)+1;
tstore= 35064 ; % 8766:8766:N_time_step;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% GENERAL PARAMETER %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=3600; %%[s] %%%
dth=1; %%[h]
DeltaGMT=1 ;
Lon=8.33 ;
Lat=46.56;
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE
%%%%%%%%
L_day=zeros(length(Datam),1);
for j=2:24:length(Datam)
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_dis=NaN; pow_dis=NaN;
clear a0 gam1 pow0 k2 DTii
%%%%%%%%%%%%%%%%%%%%
rho_g = 0.35; %%% Spatial Albedo
cc_max = 1; %% Number of vegetation
ms_max = 11; %% Number of soil layers
md_max = 10; %% % Number of debris layers

%Set zatm
zatm_hourly_on=0; %Switch for using variable instrument heights, 1=on 
if zatm_hourly_on==1
    zatm_hourly=((rand(N_time_step,1)*10)+250)/100; % Instrument height timeseries (m), this is turned into zatm in Main_Frame [Dummy variable made here for testing]
else
    zatm_surface = 2.5; %% Reference Height single value
end

%Set if using measured albedo
Ameas=rand(N_time_step,1); %Making dummy variable here, but add measured if you have it, it will only be used if one of the switches is on. Any missing data leave as NaN and modelled Albedo will be used.
Aice_meas_on_hourly=NaN(N_time_step,1);
Asno_meas_on_hourly=NaN(N_time_step,1);
idxa=isnan(Ameas)==1;
Aice_meas_on_hourly(idxa)=0;  %Use modelled ice albedo if not measured
Aice_meas_on_hourly(~idxa)=0; %When = 1 use measured ice albedo when available (Change to 0 if you don't want to use measured albedo at all)
Asno_meas_on_hourly(idxa)=0;  %Use modelled snow albedo if not measured
Asno_meas_on_hourly(~idxa)=0; %When = 1 use measured snow albedo when available (Change to 0 if you don't want to use measured albedo at all)


%%% switch Huss parameterization; 0 - no redistribution; 1 - standard Huss parameters; 2 - dh/dt derived parameters
%%% glacier size "small","med" or "large"
% Idyns = 0; 
% resol = 100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% TOPOGRAPHIC PARAMETER %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%   Matrix [m_cell x n_cell]
%DTM ;T_flow ; cellsize; xllcorner ; yllcorner ; SN ; outlet ; Aacc
load('INPUTS\dtm_Gletsch.mat')
GLH=GLH.*(GLA_MAP2>0);
clear m n MASK
[m_cell,n_cell]=size(DTM);
num_cell=numel(DTM);
x_cell=xllcorner:cellsize:(xllcorner+cellsize*(n_cell-1));
y_cell=yllcorner:cellsize:(yllcorner+cellsize*(m_cell-1));
MASK=ones(m_cell,n_cell); MASK(isnan(DTM))=0;
MASKn=reshape(MASK,num_cell,1);
Kinde = find(MASK==1);
Xout=75; Yout = 84; %%%% Outlet Point %%%%%%%%
%%% Slo_top [Fraction] %%% Aspect [rad] from N
[Slo_top,Aspect]=Slope_Aspect_indexes(DTM_orig,cellsize,'mste');
Aspect(isnan(Aspect))=0;
Slo_top(Slo_top<0.0005)=0.0005;
%Slo_top=reshape(Slo_top,num_cell,1);
%%%%%%%%%%%%%%%%%%
Asur=(1./cos(atan(Slo_top))); %% Effective Area / Projected Area
Asur=reshape(Asur,num_cell,1);
aTop= 1000*ones(m_cell,n_cell)*(cellsize^2)/cellsize; %%[mm] Area/Contour length ratio
Ared=ones(num_cell,1);
%%%%%%%%%%%%%%%%%
%%%% Flow Boundary Condition
Slo_top(Yout,Xout)=0.05;
Xout= [114 151 135 165]; Yout = [122 78 230 60]; %%%% Tracked Points %%%%%%%%%% Lysimeter
npoint = length(Xout);
Area= (cellsize^2)*sum(sum(MASK)); %% Projected area [m^2]
%%%%% Flow potential
T_pot=cell(1,ms_max);
for jk=1:ms_max;
    T_pot{jk}= T_flow;
end
%%%%%%%%%%%%%%%%%%%%%
WC = cellsize*ones(m_cell,n_cell); %%% [m]  Width channel
WC(SN==1)=0.0018*sqrt((cellsize^2)*Aacc(SN==1)); %% [m]
WC=WC.*MASK;
SN(isnan(SN))=0; %% [Stream Identifier]
SNn=reshape(SN,num_cell,1);
NMAN_C=SN*0.040; NMAN_H=0.1; %%[s/(m^1/3)] manning coefficient
MRough = 0.01*(1-SN); %%[m] Microroughness
NMAN_C=NMAN_C.*MASK;
NMAN_H=NMAN_H.*MASK;
MRough=MRough.*MASK;
%%%%
Kres_Rock =8760; %%[h] Bedrock aquifer constant
SPRINGn =SNn; %% Spring Location
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% VEGETATION      PARAMETERS LOOK-UP TABLE %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ksv=1*ones(num_cell,1);
ksv=reshape(VEG_CODE,num_cell,1);
%%% 1 Evergreen forest
%%% 2 Dec. Shrub
%%% 3 Grassland
%%% 4-5 Rock
ksv(ksv==5)=4;
Ccrown_OUT =[ 1 ; 0.9 ; 1 ; 0];  %% Ccrown fraction for PFT
EVcode = [ 1 2 3 4 ];  %% code of each PFTs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%% SOIL PARAMETER ---> ???
%%%% 2 3 10
%    90  91
%   113  114 115 116
%   118   119   121
PSOIL=[70,12,3
    85,8,1
    55,8,1];
%%%%
PSAN = 1*MASK;
PCLA = 1*MASK;
PORG = 1*MASK;
%%%%
PSAN(SOIL_ID==2)=PSOIL(1,1); PSAN(SOIL_ID==3)=PSOIL(2,1); PSAN(SOIL_ID==10)=PSOIL(3,1);
PCLA(SOIL_ID==2)=PSOIL(1,2); PCLA(SOIL_ID==3)=PSOIL(2,2); PCLA(SOIL_ID==10)=PSOIL(3,2);
PORG(SOIL_ID==2)=PSOIL(1,3); PORG(SOIL_ID==3)=PSOIL(2,3); PORG(SOIL_ID==10)=PSOIL(3,3);
%%%%%%%%%%%%%%%%%%%%%%
Pss = [800]; Pwp = [3500]; %%% [kPa]
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
[Osat,L,Pe,Ks,O33]=Soil_parameters_spatial(PSAN/100,PCLA/100,PORG/100);
[Ofc,Oss,Owp,Ohy]=Soil_parametersII_spatial(Osat,L,Pe,Ks,O33,Kfc,Pss,Pwp,Phy);
%%%%%%%%%%%%%%%%%%%%
clear Pss Pwp Kfc Phy L Pe Ks O33 Ofc Oss Owp
Osat_OUT =   Osat.*MASK; clear Osat
Osat_OUT=reshape(Osat_OUT,num_cell,1);
Ohy_OUT =   Ohy.*MASK; clear Ohy
Ohy_OUT=reshape(Ohy_OUT,num_cell,1);
%Mdep1 = 50; Mdep2= 100; Mdep3=250;
%%%%%%%%%%%%%%%%
PSAN=reshape(PSAN/100,num_cell,1);
PCLA=reshape(PCLA/100,num_cell,1);
PORG=reshape(PORG/100,num_cell,1);
Zs_OUT = 800*ones(num_cell,1); %%[mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SOLAR PARAMETER     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Computation Horizon Angle
%[HZ,Zasp] = Horizon_angle_polar(DTM,cellsize);
[HZ,Zasp] = Horizon_Angle(DTM_orig,cellsize);
%%% HZ Horizon angle array [angular degree]
%%% Z Azimuth directions  [angualr degree] from N
%%% Sky View Factor and Terrain Configuration Factor
[SvF,Ct] = Sky_View_Factor(DTM_orig,atan(Slo_top)*180/pi,Aspect,HZ,Zasp);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      TEMPORAL INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% TOTAL INITIAL CONDITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PREVIOUS CONDITION %%%%%%%%%%%%%%%%%%%%%
Ta_t=zeros(num_cell,24);
PAR_t=zeros(num_cell,24);
An_H_t=zeros(num_cell,cc_max,24);
An_L_t=zeros(num_cell,cc_max,24);
Rdark_H_t=zeros(num_cell,cc_max,24);
Rdark_L_t=zeros(num_cell,cc_max,24);
Tdp_H_t=zeros(num_cell,cc_max,24);
Tdp_L_t=zeros(num_cell,cc_max,24);
Psi_x_H_t=zeros(num_cell,cc_max,24);
Psi_l_H_t=zeros(num_cell,cc_max,24);
Psi_x_L_t=zeros(num_cell,cc_max,24);
Psi_l_L_t=zeros(num_cell,cc_max,24);
%%%
Tdp_t=zeros(num_cell,ms_max,24);
O_t=zeros(num_cell,ms_max,24);
V_t=zeros(num_cell,ms_max,24);
Pr_sno_t=zeros(num_cell,24);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CONDITION INITIAL 1 STEP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SOIL MOISTURE
vi=[  1.3108    1.3108    3.9323    6.5538    6.5538    6.5538   13.1077   13.1077   13.1077   13.1077   13.1077 ];
Vtm1=zeros(num_cell,ms_max);
for jk=1:ms_max
    Vtm1(:,jk)=vi(jk)*MASKn;
end
%%%%%%%%
Vicetm1=zeros(num_cell,ms_max);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% HYDROLOGY PART %%%%%
Tstm1 =0*MASKn;
Tdptm1 = 0*ones(num_cell,ms_max);
Tdamptm1 = Tstm1;
TdpI_Htm1= Tdamptm1*ones(1,cc_max);
TdpI_Ltm1= Tdamptm1*ones(1,cc_max);
%%%%
SWEtm1 =0*MASKn;
SNDtm1 =0*MASKn;
snow_albedotm1 = 0.2*ones(num_cell,4) ;
e_snotm1 = 0.97*MASKn;
t_slstm1= 0*MASKn;
rostm1=0*MASKn;
SP_wctm1=0*MASKn ;
tau_snotm1 = 0*MASKn;
In_Htm1 = 0*ones(num_cell,cc_max);
In_Ltm1 = 0*ones(num_cell,cc_max);
In_urbtm1 = 0*MASKn;
In_rocktm1 = 0*MASKn;
In_Littertm1 = 0*MASKn;
In_SWEtm1 = 0*MASKn;
EKtm1 = 0*MASKn;
WATtm1= 0*MASKn;
ICEtm1=  1000*0.916*reshape(GLH,num_cell,1);
IP_wctm1= 0*MASKn;
ICE_Dtm1= reshape(GLH,num_cell,1); %
Cicewtm1= 0*MASKn;
FROCKtm1= 0*MASKn;
Ticetm1= 0*MASKn;
Oicetm1 =zeros(num_cell,ms_max);
Tdebtm1 = 0*ones(num_cell,md_max-1);
Ws_undertm1 = 1*MASKn;
%%%
Ccrown_t_tm1 = ones(num_cell,cc_max); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tstm0 =Tstm1;
Slo_head = reshape(Slo_top,num_cell,1)*ones(1,ms_max); %%%
q_runon =  0*MASKn;
Q_channel = zeros(m_cell,n_cell);
Qi_in = zeros(num_cell,ms_max);
Qi_out = zeros(num_cell,ms_max);
Qi_out_Rout= zeros(m_cell,n_cell,ms_max);
Qi_in_Rout= zeros(m_cell,n_cell,ms_max);
Q_exit=0;
Qsub_exit=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('D:\T&C_code\DATA_FOLDER\Ini_condition_Rietho.mat');
%Tstm1=Ts;
%Tstm0=Ts;clear Ts
%V=V(:,1:ms_max); Vtm1=V; clear V
%%%%%%%%%%%%%%%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%
Citm1_sunH = 370*ones(num_cell,cc_max);
Citm1_shdH = 370*ones(num_cell,cc_max);
LAI_Htm1 = zeros(num_cell,cc_max);
B_Htm1 =  zeros(num_cell,cc_max,8);
PHE_S_Htm1=1*ones(num_cell,cc_max);
NPP_Htm1=zeros(num_cell,cc_max);
NPPI_Htm1=zeros(num_cell,cc_max);
Bfac_weekHtm1 = ones(num_cell,cc_max);
dflo_Htm1= zeros(num_cell,cc_max);
AgeL_Htm1= zeros(num_cell,cc_max);
AgeDL_Htm1= zeros(num_cell,cc_max);
e_rel_Htm1=ones(num_cell,cc_max);
e_relN_Htm1=ones(num_cell,cc_max);
SAI_Htm1=zeros(num_cell,cc_max);
hc_Htm1=zeros(num_cell,cc_max);
Nreserve_Htm1= 1000*ones(num_cell,cc_max);
Preserve_Htm1= 1000*ones(num_cell,cc_max);
Kreserve_Htm1= 1000*ones(num_cell,cc_max);
NupI_Htm1 = zeros(num_cell,cc_max,3);
Bfac_weekLtm1 = ones(num_cell,cc_max);
FNC_Htm1=1*ones(num_cell,cc_max);
Vx_Htm1= zeros(num_cell,cc_max);
Vl_Htm1= zeros(num_cell,cc_max);
Psi_x_Htm1= zeros(num_cell,cc_max);
Psi_l_Htm1 = zeros(num_cell,cc_max);
Tden_Htm1=zeros(num_cell,cc_max);
AgePl_Htm1 =zeros(num_cell,cc_max);


%%% Evergreen
LAI_Htm1(ksv==1,1)=2.8;
B_Htm1(ksv==1,1,1)= 280;
B_Htm1(ksv==1,1,2)= 590;
B_Htm1(ksv==1,1,3)= 395;
B_Htm1(ksv==1,1,4)= 395;
AgeL_Htm1(ksv==1,1)=1700;
SAI_Htm1(ksv==1,1)=0.2;
hc_Htm1(ksv==1,1)=13;


%%%%%%%%%%%%%%%%% LOW VEGETATION
Citm1_sunL = 370*ones(num_cell,cc_max);
Citm1_shdL = 370*ones(num_cell,cc_max);
LAI_Ltm1= zeros(num_cell,cc_max);
B_Ltm1= zeros(num_cell,cc_max,8);
PHE_S_Ltm1=1*ones(num_cell,cc_max);
NPP_Ltm1=zeros(num_cell,cc_max);
NPPI_Ltm1=zeros(num_cell,cc_max);
dflo_Ltm1= zeros(num_cell,cc_max);
AgeL_Ltm1= zeros(num_cell,cc_max);
AgeDL_Ltm1= zeros(num_cell,cc_max);
SAI_Ltm1=zeros(num_cell,cc_max);
hc_Ltm1=zeros(num_cell,cc_max);
e_rel_Ltm1=ones(num_cell,cc_max);
e_relN_Ltm1=ones(num_cell,cc_max);
Nreserve_Ltm1= 1000*ones(num_cell,cc_max);
Preserve_Ltm1= 1000*ones(num_cell,cc_max);
Kreserve_Ltm1= 1000*ones(num_cell,cc_max);
NupI_Ltm1 = zeros(num_cell,cc_max,3);
FNC_Ltm1=1*ones(num_cell,cc_max);
Vx_Ltm1= zeros(num_cell,cc_max);
Vl_Ltm1= zeros(num_cell,cc_max);
Psi_x_Ltm1 = zeros(num_cell,cc_max);
Psi_l_Ltm1 = zeros(num_cell,cc_max);
Tden_Ltm1=zeros(num_cell,cc_max);
AgePl_Ltm1=zeros(num_cell,cc_max);
%%%%%%%%%%%%%%%%%%%%%%

%%%% Shrub Dec
LAI_Ltm1(ksv==2,1)=0.0;
B_Ltm1(ksv==2,1,1)= 0;
B_Ltm1(ksv==2,1,2)= 130;
B_Ltm1(ksv==2,1,3)= 105;
B_Ltm1(ksv==2,1,4)= 120;
AgeL_Ltm1(ksv==2,1)=0;
SAI_Ltm1(ksv==2,1)=0.1;
hc_Ltm1(ksv==2,1)=0.5;
%%%% Grass
LAI_Ltm1(ksv==3,1)=0.75;
B_Ltm1(ksv==3,1,1)= 25;
B_Ltm1(ksv==3,1,2)= 0;
B_Ltm1(ksv==3,1,3)= 570;
B_Ltm1(ksv==3,1,4)= 415;
dflo_Ltm1(ksv==3,1)=0;
AgeL_Ltm1(ksv==3,1)=120;
SAI_Ltm1(ksv==3,1)=0.001;
hc_Ltm1(ksv==3,1)=0.15;
%%%%%%
NuLit_Htm1= zeros(num_cell,cc_max,3);
NuLit_Ltm1= zeros(num_cell,cc_max,3);
NBLeaf_Htm1 =zeros(num_cell,cc_max);
NBLeaf_Ltm1 =zeros(num_cell,cc_max);
PARI_Htm1 =zeros(num_cell,cc_max,3);
NBLI_Htm1  =zeros(num_cell,cc_max);
PARI_Ltm1 =zeros(num_cell,cc_max,3);
NBLI_Ltm1  =zeros(num_cell,cc_max);
%%%
TBio_H = zeros(num_cell,cc_max);
TBio_L = zeros(num_cell,cc_max);
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ij  --> Cells [num_cell]
% t, tday --> Time step [hour] [day]
Fstep=strcat('Final_reached_step_',TITLE_SAVE);
ij=0;
t=0; tday=0;
%%%%%%%%%%%%%%%% NUMERICAL METHODS OPTIONS
Opt_CR = optimset('TolFun',1);%,'UseParallel','always');
Opt_ST = optimset('TolFun',0.1);%,'UseParallel','always');
OPT_SM=  odeset('AbsTol',0.05,'MaxStep',dth);
OPT_VD=  odeset('AbsTol',0.05);
OPT_STh = odeset('AbsTol',5e+3);
OPT_PH= odeset('AbsTol',0.01);
OPT_VegSnow = 1;
OPT_SoilTemp = 0;
OPT_FR_SOIL =0;
OPT_PlantHydr = 0;
OPT_EnvLimitGrowth = 0;
OPT_VCA=0;
OPT_ALLOME=0;
OPT_Pr = 2; % Option for precipitation phase calculation, 1: Wigmosta et.al 1994, 2: Ding et.al 2014
%%%
OPT_HEAD = 0;
OPT_SoilBiogeochemistry = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic ;
%profile on
%bau = waitbar(0,'Waiting...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=2:N_time_step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %waitbar(t/N_time_step,bau)
    disp('Iter:'); disp(t);
    
    pdind = [max(1,t-24):t-1]; %% previous day indexes
    %%%%%%%%%%%%%%%% DATE and METEO SPATIAL;
    Datam_S=Datam(t,:);
    if (Datam_S(4)==1)
        tday = tday+1;
    end
    [jDay]= JULIAN_DAY(Datam_S);
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day]= SetSunVariables(Datam_S,DeltaGMT,Lon,Lat,t_bef,t_aft);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [ShF] = Shadow_Effect(DTM,h_S,zeta_S,HZ,Zasp);

    %Set zatm
    if zatm_hourly_on==1
    zatm_surface=zatm_hourly(t); %Set zatm for timestep
    end
    
    %Measured ice/snow albedo
    Aice_meas_on=Aice_meas_on_hourly(t);
    Asno_meas_on=Asno_meas_on_hourly(t);
    Ameas_t=Ameas(t); %Need to create Ameas_t to be brought into Parameters_Soil_Gletsch. Only used if Aice_meas_on=1.
    Asno_meas=Ameas(t); %Need to create Asno_meas. Will only be used if Asno_meas_on=1. 
    
    %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% INPUT  Vector [num_cell]
    %Datam_S Pr_S,Ta_S,Ds_S,Ws_S ,ea_S,N_S,Pre_S, Tdew_S
    % SAB1_S , SAB2_S,SAD1_S,SAD2_S,PARB_S,PARD_S,SvF
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%% SPATIAL DISTRIBUTION OF INPUT
    %%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PARD       Rsw_Rif    SAD1       SD         U_Rif      esat_Rif   qta_T
    % Ds_Rif     Pr         SAB1       SAD2       Ta_Rif     Ws_Rif     qta_P
    % PARB       Pre_Rif    SAB2       SB         Tdew_Rif   ea_Rif
    % qta_Pr
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IrD_S =  MASK*0; IrD_S=reshape(IrD_S,num_cell,1);
    %Pr_S = Pr(t)*MASKn;
    Pr_S = MASK*0;
    for igpr=1:21
        Pr_S = Pr_S + Pr_dist(t,igpr)*MASK.*(PID==igpr);
    end
    Pr_S=reshape(Pr_S,num_cell,1);
    %Ta_S = Ta(t)*MASKn;
    %%%%%
    lapse_rate = -0.0055; %%Thermal Lapse rate [�C/m]
    Ta_S= Ta(t)+ lapse_rate*(DTM-qta_Sta); %%%
    Ta_S=reshape(Ta_S,num_cell,1);
    %%%%%%%
    Ws_S = Ws(t)*MASKn;
    %%%
    %U_S= U(t)*MASKn;
    %Ds_S = Ds(t)*MASKn;
    %ea_S = ea(t)*MASKn;
    Vgrad= 0;% gradient [%/m]
    U_S= U(t) + Vgrad*(DTM-qta_Sta); %%%
    U_S(U_S<0)=0; U_S(U_S>1)=1;
    U_S=reshape(U_S,num_cell,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    esat_S=611*exp(17.27*Ta_S./(237.3+Ta_S)); %% [Pa]
    ea_S=U_S.*esat_S;
    Ds_S= esat_S - ea_S; %%[Pa]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Pre_S = Pre(t)*MASKn;
    Pre_S=Pre(t)*exp((-9.81/287)*(DTM-qta_Sta)/(Ta(t)+273.15));
    Pre_S=reshape(Pre_S,num_cell,1);
    %%%%%%%%%%%%%%%%%
    Tdew_S = Tdew(t)*MASKn;
    %%%%%%%%%%%%%%
    N_S = N(t)*MASKn;
    Ca_S = Ca(t)*MASKn;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cos_fst = cos(atan(Slo_top))*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);
    cos_fst(cos_fst<0)=0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sin(h_S) <= 0.10 %%[5.73�] Numerical problems
        SAB1_S  =  0*MASKn;
        SAB2_S  =  0*MASKn;
        SAD1_S  =  0*MASKn;
        SAD2_S  =  0*MASKn;
        PARB_S  =  0*MASKn;
        PARD_S  =  0*MASKn;
    else
        SAD1_S = SAD1(t)*SvF + Ct.*rho_g.*((SAB1(t)/sin(h_S)).*cos_fst + (1-SvF).*SAD1(t));
        SAD2_S = SAD2(t)*SvF + Ct.*rho_g.*((SAB2(t)/sin(h_S)).*cos_fst + (1-SvF).*SAD2(t));
        PARD_S = PARD(t).*SvF + Ct*rho_g.*((PARB(t)/sin(h_S)).*cos_fst + (1-SvF).*PARD(t));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SAB1_S =(SAB1(t)/sin(h_S)).*cos_fst.*ShF;
        SAB2_S =(SAB2(t)/sin(h_S)).*cos_fst.*ShF;
        PARB_S = (PARB(t)/sin(h_S)).*cos_fst.*ShF;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        grad1= 0.0049; % gradient [[W/m^2]/m]
        SAB1_S= (SAB1_S + grad1*(DTM-qta_Sta)) ; %%%
        SAB1_S(SAB1_S<0)=0;
        grad1= 0.0080; % gradient  [[W/m^2]/m]
        SAB2_S= (SAB2_S + grad1*(DTM-qta_Sta)) ; %%%
        SAB2_S(SAB2_S<0)=0;
        grad1=  0.0047; % gradient [[W/m^2]/m]
        PARB_S= (PARB_S + grad1*(DTM-qta_Sta)) ; %%%
        PARB_S(PARB_S<0)=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SAB1_S  =  reshape(SAB1_S,num_cell,1);
        SAB2_S  =  reshape(SAB2_S,num_cell,1);
        SAD1_S  =  reshape(SAD1_S,num_cell,1);
        SAD2_S  =  reshape(SAD2_S,num_cell,1);
        PARB_S  =  reshape(PARB_S,num_cell,1);
        PARD_S  =  reshape(PARD_S,num_cell,1);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% SPATIAL INITIALIZATION VECTOR PREDEFINING %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if t == 2
        %%%%%%%% VEGETATION
        %%%%%%%% VEGETATION
        NPP_H=zeros(num_cell,cc_max);
        ANPP_H=zeros(num_cell,cc_max);
        Rg_H=zeros(num_cell,cc_max);
        RA_H=zeros(num_cell,cc_max);
        Rms_H=zeros(num_cell,cc_max);
        Rmr_H=zeros(num_cell,cc_max);
        Rmc_H=zeros(num_cell,cc_max);
        LAIdead_H=zeros(num_cell,cc_max);
        Sr_H=zeros(num_cell,cc_max);
        Slf_H=zeros(num_cell,cc_max);
        Sfr_H=zeros(num_cell,cc_max);
        Swm_H=zeros(num_cell,cc_max);
        Sll_H=zeros(num_cell,cc_max);
        Rexmy_H= zeros(num_cell,cc_max,3);
        Rrootl_H=zeros(num_cell,cc_max);
        Bfac_dayH=zeros(num_cell,cc_max);
        Bfac_weekH=Bfac_weekHtm1;
        NPPI_H= NPPI_Htm1;
        RB_H = zeros(num_cell,cc_max,7);
        rNc_H =zeros(num_cell,cc_max);
        rPc_H =zeros(num_cell,cc_max);
        rKc_H = zeros(num_cell,cc_max);
        ManIH = zeros(num_cell,cc_max);
        %%%
        NPP_L=zeros(num_cell,cc_max);
        ANPP_L=zeros(num_cell,cc_max);
        Rg_L=zeros(num_cell,cc_max);
        RA_L=zeros(num_cell,cc_max);
        Rms_L=zeros(num_cell,cc_max);
        Rmr_L=zeros(num_cell,cc_max);
        Rmc_L=zeros(num_cell,cc_max);
        LAIdead_L=zeros(num_cell,cc_max);
        Sr_L=zeros(num_cell,cc_max);
        Slf_L=zeros(num_cell,cc_max);
        Sfr_L=zeros(num_cell,cc_max);
        Swm_L=zeros(num_cell,cc_max);
        Sll_L=zeros(num_cell,cc_max);
        Rexmy_L= zeros(num_cell,cc_max,3);
        Rrootl_L=zeros(num_cell,cc_max);
        Bfac_dayL =zeros(num_cell,cc_max);
        Bfac_weekL=Bfac_weekLtm1;
        NPPI_L=NPPI_Ltm1;
        RB_L  = zeros(num_cell,cc_max,7);
        rNc_L=zeros(num_cell,cc_max);
        rPc_L=zeros(num_cell,cc_max);
        rKc_L=zeros(num_cell,cc_max);
        ManIL=zeros(num_cell,cc_max);
        %%%%
        e_rel_H=e_rel_Htm1;
        e_relN_H = e_relN_Htm1;
        LAI_H=LAI_Htm1;
        B_H= B_Htm1;
        PHE_S_H=PHE_S_Htm1;
        dflo_H= dflo_Htm1;
        AgeL_H= AgeL_Htm1;
        AgeDL_H=AgeDL_Htm1;
        SAI_H= SAI_Htm1;
        hc_H= hc_Htm1;
        NBLeaf_H =NBLeaf_Htm1;
        Nreserve_H = Nreserve_Htm1;
        Preserve_H = Preserve_Htm1;
        Kreserve_H = Kreserve_Htm1;
        FNC_H =  FNC_Htm1;
        NupI_H = NupI_Htm1;
        TdpI_H= TdpI_Htm1;
        PARI_H =PARI_Htm1 ;
        NBLI_H  =NBLI_Htm1;
        %%%
        e_rel_L=e_rel_Ltm1;
        e_relN_L = e_relN_Ltm1;
        LAI_L= LAI_Ltm1;
        B_L= B_Ltm1;
        PHE_S_L= PHE_S_Ltm1;
        dflo_L= dflo_Ltm1 ;
        AgeL_L= AgeL_Ltm1 ;
        AgeDL_L =AgeDL_Ltm1;
        SAI_L= SAI_Ltm1;
        hc_L= hc_Ltm1;
        NBLeaf_L =NBLeaf_Ltm1;
        Nreserve_L = Nreserve_Ltm1;
        Preserve_L = Preserve_Ltm1;
        Kreserve_L = Kreserve_Ltm1;
        FNC_L =  FNC_Ltm1;
        NupI_L = NupI_Ltm1;
        TdpI_L= TdpI_Ltm1;
        PARI_L =PARI_Ltm1 ;
        NBLI_L  =NBLI_Ltm1;
        %%%%%
        Nuptake_H= zeros(num_cell,cc_max);
        Puptake_H= zeros(num_cell,cc_max);
        Kuptake_H= zeros(num_cell,cc_max);
        Nuptake_L= zeros(num_cell,cc_max);
        Puptake_L= zeros(num_cell,cc_max);
        Kuptake_L= zeros(num_cell,cc_max);
        NavlI = zeros(num_cell,3);
        Bam=zeros(num_cell,1);
        Bem=zeros(num_cell,1);
        %%%%
        TexC_H=zeros(num_cell,cc_max);
        TexN_H=zeros(num_cell,cc_max);
        TexP_H=zeros(num_cell,cc_max);
        TexK_H=zeros(num_cell,cc_max);
        TNIT_H=zeros(num_cell,cc_max);
        TPHO_H=zeros(num_cell,cc_max);
        TPOT_H=zeros(num_cell,cc_max);
        SupN_H=zeros(num_cell,cc_max);
        SupP_H=zeros(num_cell,cc_max);
        SupK_H=zeros(num_cell,cc_max);
        ISOIL_H= zeros(num_cell,cc_max,18);
        TexC_L=zeros(num_cell,cc_max);
        TexN_L=zeros(num_cell,cc_max);
        TexP_L=zeros(num_cell,cc_max);
        TexK_L=zeros(num_cell,cc_max);
        TNIT_L=zeros(num_cell,cc_max);
        TPHO_L=zeros(num_cell,cc_max);
        TPOT_L=zeros(num_cell,cc_max);
        SupN_L=zeros(num_cell,cc_max);
        SupP_L=zeros(num_cell,cc_max);
        SupK_L=zeros(num_cell,cc_max);
        ISOIL_L= zeros(num_cell,cc_max,18);
        BA_H =zeros(num_cell,cc_max);
        Tden_H=zeros(num_cell,cc_max);
        AgePl_H=zeros(num_cell,cc_max);
        BA_L =zeros(num_cell,cc_max);
        Tden_L =zeros(num_cell,cc_max);
        AgePl_L =zeros(num_cell,cc_max);
        Ccrown_t =zeros(num_cell,cc_max);
        %%%
        NuLit_H= NuLit_Htm1;
        NuLit_L= NuLit_Ltm1;
        BLit=zeros(num_cell,cc_max);
        %%%%%%% HYDROLOGY
        V=zeros(num_cell,ms_max);
        O=zeros(num_cell,ms_max);
        ZWT=zeros(num_cell,1);
        POT=zeros(num_cell,ms_max);
        OF=zeros(num_cell,1);
        OS=zeros(num_cell,1);
        OH=zeros(num_cell,cc_max);
        OL=zeros(num_cell,cc_max);
        Rd=zeros(num_cell,1);
        Qi_out=zeros(num_cell,ms_max);
        Rh=zeros(num_cell,1);
        Lk=zeros(num_cell,1);
        f=zeros(num_cell,1);
        WIS=zeros(num_cell,1);
        Ts=zeros(num_cell,1);
        Csno=zeros(num_cell,1);
        Cice=zeros(num_cell,1);
        Pr_sno=zeros(num_cell,1);
        Pr_liq=zeros(num_cell,1);
        rb_H=zeros(num_cell,cc_max);
        rb_L=zeros(num_cell,cc_max);
        rs_sunH=zeros(num_cell,cc_max);
        rs_sunL=zeros(num_cell,cc_max);
        rs_shdH=zeros(num_cell,cc_max);
        rs_shdL=zeros(num_cell,cc_max);
        rap_H=zeros(num_cell,cc_max);
        rap_L=zeros(num_cell,cc_max);
        r_soil=zeros(num_cell,1);
        b_soil=zeros(num_cell,1);
        alp_soil=zeros(num_cell,1);
        ra=zeros(num_cell,1);
        r_litter=zeros(num_cell,cc_max);
        WR_SP=zeros(num_cell,1);
        U_SWE=zeros(num_cell,1);
        NIn_SWE=zeros(num_cell,1);
        dQ_S=zeros(num_cell,1);
        DQ_S=zeros(num_cell,1);
        DT_S=zeros(num_cell,1);
        Dr_H=zeros(num_cell,cc_max);
        Dr_L=zeros(num_cell,cc_max);
        SE_rock=zeros(num_cell,1);
        SE_urb=zeros(num_cell,1);
        An_L=zeros(num_cell,cc_max);
        An_H=zeros(num_cell,cc_max);
        Rdark_L=zeros(num_cell,cc_max);
        Rdark_H=zeros(num_cell,cc_max);
        Ci_sunH=zeros(num_cell,cc_max);
        Ci_sunL=zeros(num_cell,cc_max);
        Ci_shdH=zeros(num_cell,cc_max);
        Ci_shdL=zeros(num_cell,cc_max);
        Rn=zeros(num_cell,1);
        H=zeros(num_cell,1);
        QE=zeros(num_cell,1);
        Qv=zeros(num_cell,1);
        Lpho=zeros(num_cell,1);
        T_H=zeros(num_cell,cc_max);
        T_L=zeros(num_cell,cc_max);
        EIn_H=zeros(num_cell,cc_max);
        EIn_L=zeros(num_cell,cc_max);
        EG=zeros(num_cell,1);
        ELitter=zeros(num_cell,1);
        ESN=zeros(num_cell,1);
        ESN_In=zeros(num_cell,1);
        EWAT=zeros(num_cell,1);
        EICE=zeros(num_cell,1);
        EIn_urb=zeros(num_cell,1);
        EIn_rock=zeros(num_cell,1);
        dw_SNO=zeros(num_cell,1);
        G=zeros(num_cell,1);
        Gfin=zeros(num_cell,1);
        Tdp=zeros(num_cell,ms_max);
        Tdeb=zeros(num_cell,md_max-1);
        Vice=zeros(num_cell,ms_max);
        Oice=zeros(num_cell,ms_max);
        Tice=zeros(num_cell,1);
        Imelt=zeros(num_cell,1);
        Smelt=zeros(num_cell,1);
        Tdamp=zeros(num_cell,1);
        Tdp_H=zeros(num_cell,cc_max);
        Tdp_L=zeros(num_cell,cc_max);
        SWE=zeros(num_cell,1);
        SND=zeros(num_cell,1);
        ros=zeros(num_cell,1);
        In_SWE=zeros(num_cell,1);
        SP_wc=zeros(num_cell,1);
        WAT=zeros(num_cell,1);
        ICE=zeros(num_cell,1);
        ICE_D=zeros(num_cell,1);
        IP_wc=zeros(num_cell,1);
        WR_IP=zeros(num_cell,1);
        NIce=zeros(num_cell,1);
        Cicew=zeros(num_cell,1);
        Csnow=zeros(num_cell,1);
        FROCK=zeros(num_cell,1);
        Lk_wat=zeros(num_cell,1);
        Lk_rock=zeros(num_cell,1);
        Qfm=zeros(num_cell,1);
        t_sls=zeros(num_cell,1);
        In_H=In_Htm1;
        In_L=In_Ltm1;
        In_Litter=In_Littertm1;
        In_urb=zeros(num_cell,1);
        In_rock=zeros(num_cell,1);
        er=zeros(num_cell,1);
        snow_albedo=zeros(num_cell,4);
        tau_sno=zeros(num_cell,1);
        e_sno=zeros(num_cell,1);
        EK=zeros(num_cell,1);
        dQVEG=zeros(num_cell,1);
        TsVEG=zeros(num_cell,1);
        Psi_s_H=zeros(num_cell,cc_max);
        Psi_s_L=zeros(num_cell,cc_max);
        gsr_H=zeros(num_cell,cc_max);
        Psi_x_H=Psi_x_Htm1;
        Psi_l_H=Psi_l_Htm1;
        Jsx_H=zeros(num_cell,cc_max);
        Jxl_H=zeros(num_cell,cc_max);
        Kleaf_H=zeros(num_cell,cc_max);
        Kx_H=zeros(num_cell,cc_max);
        Vx_H=Vx_Htm1;
        Vl_H=Vl_Htm1;
        gsr_L=zeros(num_cell,cc_max);
        Psi_x_L=Psi_x_Ltm1;
        Psi_l_L=Psi_l_Ltm1;
        Jsx_L=zeros(num_cell,cc_max);
        Jxl_L=zeros(num_cell,cc_max);
        Kleaf_L=zeros(num_cell,cc_max);
        Kx_L=zeros(num_cell,cc_max);
        Vx_L=Vx_Ltm1;
        Vl_L=Vl_Ltm1;
        fapar_H=zeros(num_cell,cc_max);
        fapar_L=zeros(num_cell,cc_max);
        SIF_H=zeros(num_cell,cc_max);
        SIF_L=zeros(num_cell,cc_max);
        NDVI=zeros(num_cell,1);
        %CK2=zeros(num_cell,1);
        Ws_under=zeros(num_cell,1);
        HV=zeros(num_cell,1);
        QEV=zeros(num_cell,1);
        CK1=zeros(num_cell,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %parfor (ij=1:num_cell)
    for (ij=1:num_cell)
        if MASKn(ij)== 1
            Elev=DTM(ij);
            %[i,j] = ind2sub([m_cell,n_cell],ij);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% BOUNDARY CONDITION  %%% INTRODUCED SOIL AND VEG. for ij
            [aR,Zs,...
                EvL_Zs,Inf_Zs,Bio_Zs,Zinf,RfH_Zs,RfL_Zs,dz,Ks_Zs,Dz,...
                ms,Kbot,Krock,zatm,...
                Ccrown,Cbare,Crock,Curb,Cwat,...
                Color_Class,OM_H,OM_L,PFT_opt_H,PFT_opt_L,d_leaf_H,d_leaf_L,...
                SPAR,Phy,Soil_Param,Interc_Param,SnowIce_Param,VegH_Param,VegL_Param,fpr,...
                VegH_Param_Dyn,VegL_Param_Dyn,...
                Stoich_H,aSE_H,Stoich_L,aSE_L,fab_H,fbe_H,fab_L,fbe_L,...
                ZR95_H,ZR95_L,In_max_urb,In_max_rock,K_usle,Deb_Par,Zs_deb,...
                Sllit,Kct,ExEM,ParEx_H,Mpar_H,ParEx_L,Mpar_L]=PARAMETERS_SOIL_Gletsch(ksv(ij),PSAN(ij),PCLA(ij),PORG(ij), DEB_MAP(ij), md_max,zatm_surface,Aice_meas_on,Ameas_t);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (Datam_S(4)==1)
                %%%%%%%%%%%%%
                %%%%% SOIL BIOGEOCHEMISTRY MODULE
                [Se_bio,Se_fc,Psi_bio,Tdp_bio,VSUM,VTSUM]=Biogeo_environment([squeeze(Tdp_t(ij,:,:))]',[squeeze(O_t(ij,:,:))]',[squeeze(V_t(ij,:,:))]',...
                    Soil_Param,Phy,SPAR,Bio_Zs);%
                %%% Biogeochemistry Unit
                Nuptake_H(ij,:)= 0.0;
                Puptake_H(ij,:)= 0.0;
                Kuptake_H(ij,:)= 0.0; %% [gK/m^2 day]
                %%%%
                Nuptake_L(ij,:)= 0.0;  %% [gN/m^2 day]
                Puptake_L(ij,:)= 0.0;
                Kuptake_L(ij,:)= 0.0;
                %%%% %%%%
                NavlI(ij,:)=[1 1 1];
                Bam(ij)=0; Bem(ij)=0;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %VEGETATION MODULE
                [LAI_H(ij,:),B_H(ij,:,:),NPP_H(ij,:),ANPP_H(ij,:),Rg_H(ij,:),RA_H(ij,:),...
                    Rms_H(ij,:),Rmr_H(ij,:),Rmc_H(ij,:),PHE_S_H(ij,:),...
                    dflo_H(ij,:),AgeL_H(ij,:),e_rel_H(ij,:),e_relN_H(ij,:),...
                    LAI_L(ij,:),B_L(ij,:,:),NPP_L(ij,:),ANPP_L(ij,:),Rg_L(ij,:),RA_L(ij,:),...
                    Rms_L(ij,:),Rmr_L(ij,:),Rmc_L(ij,:),PHE_S_L(ij,:),...
                    dflo_L(ij,:),AgeL_L(ij,:),e_rel_L(ij,:),e_relN_L(ij,:),...
                    SAI_H(ij,:),hc_H(ij,:),SAI_L(ij,:),hc_L(ij,:),...
                    LAIdead_H(ij,:),NBLeaf_H(ij,:),Sr_H(ij,:),Slf_H(ij,:),Sfr_H(ij,:),Sll_H(ij,:),Swm_H(ij,:),Rexmy_H(ij,:,:),NupI_H(ij,:,:),NuLit_H(ij,:,:),...
                    LAIdead_L(ij,:),NBLeaf_L(ij,:),Sr_L(ij,:),Slf_L(ij,:),Sfr_L(ij,:),Sll_L(ij,:),Swm_L(ij,:),Rexmy_L(ij,:,:),NupI_L(ij,:,:),NuLit_L(ij,:,:),...
                    Rrootl_H(ij,:),AgeDL_H(ij,:),Bfac_dayH(ij,:),Bfac_weekH(ij,:),NPPI_H(ij,:),TdpI_H(ij,:),PARI_H(ij,:,:),NBLI_H(ij,:),RB_H(ij,:,:),FNC_H(ij,:),...
                    Nreserve_H(ij,:),Preserve_H(ij,:),Kreserve_H(ij,:),rNc_H(ij,:),rPc_H(ij,:),rKc_H(ij,:),ManIH(ij,:),...
                    Rrootl_L(ij,:),AgeDL_L(ij,:),Bfac_dayL(ij,:),Bfac_weekL(ij,:),NPPI_L(ij,:),TdpI_L(ij,:),PARI_L(ij,:,:),NBLI_L(ij,:),RB_L(ij,:,:),FNC_L(ij,:),...
                    Nreserve_L(ij,:),Preserve_L(ij,:),Kreserve_L(ij,:),rNc_L(ij,:),rPc_L(ij,:),rKc_L(ij,:),ManIL(ij,:),...
                    TexC_H(ij,:),TexN_H(ij,:),TexP_H(ij,:),TexK_H(ij,:),TNIT_H(ij,:),TPHO_H(ij,:),TPOT_H(ij,:),...
                    SupN_H(ij,:),SupP_H(ij,:),SupK_H(ij,:),ISOIL_H(ij,:,:),...
                    TexC_L(ij,:),TexN_L(ij,:),TexP_L(ij,:),TexK_L(ij,:),TNIT_L(ij,:),TPHO_L(ij,:),TPOT_L(ij,:),...
                    SupN_L(ij,:),SupP_L(ij,:),SupK_L(ij,:),ISOIL_L(ij,:,:),...
                    BA_H(ij,:),Tden_H(ij,:),AgePl_H(ij,:),BA_L(ij,:),Tden_L(ij,:),AgePl_L(ij,:),Ccrown_t(ij,:)]=VEGETATION_MODULE_PAR(cc_max,Ccrown,ZR95_H,ZR95_L,B_Htm1(ij,:,:),...
                    PHE_S_Htm1(ij,:),dflo_Htm1(ij,:),AgeL_Htm1(ij,:),AgeDL_Htm1(ij,:),...
                    Ta_t(ij,:),PAR_t(ij,:),Tdp_H_t(ij,:,:),Psi_x_H_t(ij,:,:),Psi_l_H_t(ij,:,:),An_H_t(ij,:,:),Rdark_H_t(ij,:,:),NPP_Htm1(ij,:),jDay,Datam_S,...
                    NPPI_Htm1(ij,:),TdpI_Htm1(ij,:),Bfac_weekHtm1(ij,:),...
                    Stoich_H,aSE_H,VegH_Param_Dyn,...
                    Nreserve_Htm1(ij,:),Preserve_Htm1(ij,:),Kreserve_Htm1(ij,:),Nuptake_H(ij,:),Puptake_H(ij,:),Kuptake_H(ij,:),FNC_Htm1(ij,:),Tden_Htm1(ij,:),AgePl_Htm1(ij,:),...
                    fab_H,fbe_H,ParEx_H,Mpar_H,TBio_H(ij,:),SAI_Htm1(ij,:),hc_Htm1(ij,:),...
                    B_Ltm1(ij,:,:),PHE_S_Ltm1(ij,:),dflo_Ltm1(ij,:),AgeL_Ltm1(ij,:),AgeDL_Ltm1(ij,:),...
                    Tdp_L_t(ij,:,:),Psi_x_L_t(ij,:,:),Psi_l_L_t(ij,:,:),An_L_t(ij,:,:),Rdark_L_t(ij,:,:),NPP_Ltm1(ij,:),...
                    NPPI_Ltm1(ij,:),TdpI_Ltm1(ij,:),Bfac_weekLtm1(ij,:),...
                    NupI_Htm1(ij,:,:),NupI_Ltm1(ij,:,:),NuLit_Htm1(ij,:,:),NuLit_Ltm1(ij,:,:),NBLeaf_Htm1(ij,:),NBLeaf_Ltm1(ij,:),...
                    PARI_Htm1(ij,:,:),NBLI_Htm1(ij,:),PARI_Ltm1(ij,:,:),NBLI_Ltm1(ij,:),...
                    Stoich_L,aSE_L,VegL_Param_Dyn,...
                    NavlI(ij,:),Bam(ij),Bem(ij),Ccrown_t_tm1(ij,:),...
                    Nreserve_Ltm1(ij,:),Preserve_Ltm1(ij,:),Kreserve_Ltm1(ij,:),Nuptake_L(ij,:),Puptake_L(ij,:),Kuptake_L(ij,:),FNC_Ltm1(ij,:),Tden_Ltm1(ij,:),AgePl_Ltm1(ij,:),...
                    fab_L,fbe_L,ParEx_L,Mpar_L,TBio_L(ij,:),SAI_Ltm1(ij,:),hc_Ltm1(ij,:),...
                    ExEM,Lmax_day,L_day,Se_bio,Tdp_bio,OPT_EnvLimitGrowth,OPT_VD,OPT_VCA,OPT_ALLOME,OPT_SoilBiogeochemistry);
                
                BLit(ij,:)= 0.0 ; % %% %%[kg DM / m2]
            end
            %%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%% HYDROLOGY MODULE
            %%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [V(ij,:),O(ij,:),Vice(ij,:),Oice(ij,:),ZWT(ij),OF(ij),OS(ij),OH(ij,:),OL(ij,:),Psi_s_H(ij,:),Psi_s_L(ij,:),Rd(ij),Qi_out(ij,:),...
                Rh(ij),Lk(ij),f(ij),WIS(ij),Ts(ij),Csno(ij),Cice(ij),NDVI(ij),...
                Pr_sno(ij),Pr_liq(ij),rb_H(ij,:),rb_L(ij,:),rs_sunH(ij,:),rs_sunL(ij,:),rs_shdH(ij,:),rs_shdL(ij,:),...
                rap_H(ij,:),rap_L(ij,:),r_soil(ij),b_soil(ij),alp_soil(ij),ra(ij),r_litter(ij,:),...
                WR_SP(ij),U_SWE(ij),NIn_SWE(ij),dQ_S(ij),DQ_S(ij),DT_S(ij),...
                WAT(ij),ICE(ij),ICE_D(ij),IP_wc(ij),WR_IP(ij),NIce(ij),Cicew(ij),Csnow(ij),FROCK(ij),...
                Dr_H(ij,:),Dr_L(ij,:),SE_rock(ij),SE_urb(ij),Lk_wat(ij),Lk_rock(ij),...
                An_L(ij,:),An_H(ij,:),Rdark_L(ij,:),Rdark_H(ij,:),Ci_sunH(ij,:),Ci_sunL(ij,:),Ci_shdH(ij,:),Ci_shdL(ij,:),Rn(ij),...
                H(ij),QE(ij),Qv(ij),Lpho(ij),T_H(ij,:),T_L(ij,:),EIn_H(ij,:),EIn_L(ij,:),EG(ij),ELitter(ij),ESN(ij),ESN_In(ij),...
                EWAT(ij),EICE(ij),EIn_urb(ij),EIn_rock(ij),dw_SNO(ij),Imelt(ij),Smelt(ij),...
                G(ij),Gfin(ij),Tdp(ij,:),Tdeb(ij,:),Tice(ij,:),Tdamp(ij),Tdp_H(ij,:),Tdp_L(ij,:),SWE(ij),SND(ij),ros(ij),In_SWE(ij),SP_wc(ij),Qfm(ij),t_sls(ij),...
                In_H(ij,:),In_L(ij,:),In_Litter(ij),In_urb(ij),In_rock(ij),...
                gsr_H(ij,:),Psi_x_H(ij,:),Psi_l_H(ij,:),Jsx_H(ij,:),Jxl_H(ij,:),Kleaf_H(ij,:),Kx_H(ij,:),Vx_H(ij,:),Vl_H(ij,:),...
                gsr_L(ij,:),Psi_x_L(ij,:),Psi_l_L(ij,:),Jsx_L(ij,:),Jxl_L(ij,:),Kleaf_L(ij,:),Kx_L(ij,:),Vx_L(ij,:),Vl_L(ij,:),...
                fapar_H(ij,:),fapar_L(ij,:),SIF_H(ij,:),SIF_L(ij,:),...
                Ws_under(ij),er(ij),snow_albedo(ij,:),tau_sno(ij),e_sno(ij),...
                HV(ij),QEV(ij),dQVEG(ij),TsVEG(ij),EK(ij),POT(ij,:),CK1(ij)]=HYDROLOGY_MODULE_PAR(Vtm1(ij,:),Oicetm1(ij,:),...
                aR,Zs,...
                EvL_Zs,Inf_Zs,Zinf,RfH_Zs,RfL_Zs,dz,Dz,ms,Kbot,Pr_S(ij),Ta_S(ij),Ds_S(ij),Ws_S(ij),zatm,Tstm1(ij),dt,dth,ea_S(ij),N_S(ij),Pre_S(ij),Tstm0(ij),...
                LAI_H(ij,:),SAI_H(ij,:),LAI_L(ij,:),SAI_L(ij,:),LAIdead_H(ij,:),LAIdead_L(ij,:),...
                Rrootl_H(ij,:),Rrootl_L(ij,:),BLit(ij,:),Sllit,Kct,...
                Datam_S,DeltaGMT,Lon,Lat,Elev,t_bef,t_aft,...
                Ccrown,Cbare,Crock,Curb,Cwat,...
                SAB1_S(ij),SAB2_S(ij),SAD1_S(ij),SAD2_S(ij),PARB_S(ij),PARD_S(ij),SvF(ij),SNDtm1(ij),...
                snow_albedotm1(ij,:),Asno_meas, Asno_meas_on,...
                Color_Class,OM_H,OM_L,PFT_opt_H,PFT_opt_L,hc_H(ij,:),hc_L(ij,:),d_leaf_H,d_leaf_L,...
                Soil_Param,Interc_Param,SnowIce_Param,VegH_Param,VegL_Param,...
                Ca_S(ij),Oa,Citm1_sunH(ij,:),Citm1_shdH(ij,:),Citm1_sunL(ij,:),Citm1_shdL(ij,:),...
                e_rel_H(ij,:),e_relN_H(ij,:),e_rel_L(ij,:),e_relN_L(ij,:),...
                e_snotm1(ij),In_Htm1(ij,:),In_Ltm1(ij,:),In_Littertm1(ij),In_urbtm1(ij),In_rocktm1(ij),SWEtm1(ij),In_SWEtm1(ij),....
                Tdebtm1(ij,:),Ticetm1(ij),Tdptm1(ij,:),Tdamptm1(ij),...
                WATtm1(ij),ICEtm1(ij),IP_wctm1(ij),ICE_Dtm1(ij),Cicewtm1(ij),...
                Vx_Htm1(ij,:),Vl_Htm1(ij,:),Vx_Ltm1(ij,:),Vl_Ltm1(ij,:),Psi_x_Htm1(ij,:),Psi_l_Htm1(ij,:),Psi_x_Ltm1(ij,:),Psi_l_Ltm1(ij,:),...
                ZR95_H,ZR95_L,...
                FROCKtm1(ij),Krock,...
                Deb_Par,Zs_deb,...
                Tdew_S(ij),t_slstm1(ij),rostm1(ij),SP_wctm1(ij),fpr,IrD_S(ij),...
                In_max_urb,In_max_rock,K_usle,tau_snotm1(ij),Ta_t(ij,:),...
                Slo_top(ij),Slo_head(ij,:),Asur(ij),Ared(ij),aTop(ij),EKtm1(ij),q_runon(ij),Qi_in(ij,:),...
                Ws_undertm1(ij),Pr_sno_t(ij,:),...
                pow_dis,a_dis,...
                SPAR,SNn(ij),OPT_VegSnow,OPT_SoilTemp,OPT_PlantHydr,Opt_CR,Opt_ST,OPT_SM,OPT_STh,OPT_FR_SOIL,OPT_PH);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii=1:ms_max
        Qi_out_Rout(:,:,ii) =  reshape(Qi_out(:,ii),m_cell,n_cell); %%[mm/h]
    end
    Rd = reshape(Rd,m_cell,n_cell); %%[mm]
    Rh = reshape(Rh,m_cell,n_cell); %%[mm]
    %%%%%%%%%%% ROUTING MODULE
    %%%%%%%%%%%%%%%%
    [q_runon,Q_channel,Qi_in_Rout,Slo_pot,Q_exit,Qsub_exit,T_pot,...
        QpointH,QpointC,UpointH,UpointC]= ROUTING_MODULE(dt,dth,Rd,Rh,Qi_out_Rout,Q_channel,...
        cellsize,Area,DTM,NMAN_H,NMAN_C,MRough,WC,SN,T_flow,T_pot,Slo_top,ms_max,POT,ZWT,OPT_HEAD,Xout,Yout);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii=1:ms_max
        Qi_in(:,ii) =  reshape(Qi_in_Rout(:,:,ii),num_cell,1); %%[mm]
        Slo_head(:,ii) = reshape(Slo_pot(:,:,ii),num_cell,1);
    end
    Rd = reshape(Rd,num_cell,1); %%[mm]
    Rh = reshape(Rh,num_cell,1); %%[mm]
    q_runon = reshape(q_runon,num_cell,1); %%[mm]
    %%%%%%%%%%%%%
    Qi_in=Qi_in/dth;%%% [mm/h]
    q_runon = q_runon/dth; %%% [mm/h]
    %%% Q_exit Qsub_exit [mm] over the entire domain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if not(isreal(sum(sum(q_runon))))
        disp('The Program fails because of runon numerical instability')
        break
    end
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% Glacier volume redistribution %%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if (Datam_S(2)==1) && (Datam_S(3)==1) && (Datam_S(4)==1) && (Idyn > 0) % e.g. first of January
%         % sum up total water equivalent in [m w.e.], store ratios
%         WEbef = ICE + SWE; %total w.e. before redistribution
%         WEym1 = ICEym1 + SWEym1; %total w.e. from the year before
%         raICE = ICE./WEbef; raICE(isnan(raICE))=0; %ice ratio of total
%         raSWE = SWE./WEbef; raSWE(isnan(raSWE))=0; %snow ratio of total mm w.e.
%         THbef = reshape(WEbef,m_cell,n_cell)./916; %conversion to thickness in m, rhi parameter needed
%         THym1 = reshape(WEym1,m_cell,n_cell)./916; %conversion to thickness in m, rhi parameter needed
%         
%         % adjust surface for average mass balance for each glacier
%         [THnew,DTM,GLA_ID] = apply_Huss_TC(THbef,THym1,DTM,GLA_ID,Gla_nEl_nDH,Idyn,resol);%mass supplied as m w.e.
%         
%         % reshape total water equivalent back to array map
%         WEnew = reshape(THnew,num_cell,1).*916; % conversion back to mm w.e., rhi parameter needed
%         % recalculate snow/ice w.e. maps based on old ratios
%         ICE = raICE.*WEnew; % re-calculate ice w.e.
%         SWE(WEnew>0) = raSWE(WEnew>0).*WEnew(WEnew>0); % re-calculate SWE, only where glacier
%         % store for next year redistribution
%         ICEym1 = ICE;     
%         SWEym1 = SWE;
%         % calculate new snow/ice depths
%         ICE_D = ICE/916;
%         SND = 0.001*SWE./ros;
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% AVALANCHES COMPONENT %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SND=reshape(SND,m_cell,n_cell);
    SWE=reshape(SWE,m_cell,n_cell);
    ros=reshape(ros,m_cell,n_cell);
    [SND,SWE,ros,Swe_exit]= AVALANCHES(DTM,cellsize,Area,reshape(Asur,m_cell,n_cell),Slo_top,SND,SWE,ros);
    SND= reshape(SND,num_cell,1);
    SWE = reshape(SWE,num_cell,1);
    ros = reshape(ros,num_cell,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% FRACTURED ROCK COMPONENT %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Q_channel,FROCK,Qflow_rock]= FRACTURED_ROCK(Q_channel,FROCK,SPRINGn,dth,m_cell,n_cell,num_cell,Kres_Rock);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if t==2
        V_tgtm1 =sum(Ared.*Asur.*sum(Vtm1,2))*(cellsize^2)/Area;
        Vice_tgtm1 =sum(Ared.*Asur.*sum(Vicetm1,2))*(cellsize^2)/Area;
        SWE_tgtm1 = sum(SWEtm1)*(cellsize^2)/Area; %%
        In_tgtm1 = (sum(sum(In_Htm1)) + sum(sum(In_Ltm1)) + sum(sum(In_Littertm1)) + sum(SP_wctm1) + ...
            sum(In_SWEtm1) + sum(In_urbtm1) + sum(In_rocktm1) +sum(IP_wctm1) )*(cellsize^2)/Area ;
        ICE_tgtm1 =  sum(ICEtm1)*(cellsize^2)/Area; %%
        WAT_tgtm1 =  sum(WATtm1)*(cellsize^2)/Area; %%
        FROCK_tgtm1 =  sum(FROCKtm1)*(cellsize^2)/Area; %%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% INITIAL CONDITION FOR THE NEXT STEP        %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Tstm0 =2*Ts-Tstm1;
    Vtm1=V;
    Tstm1 =Ts;
    SNDtm1 = SND ;
    snow_albedotm1 = snow_albedo ;
    Citm1_sunH = Ci_sunH;
    Citm1_sunL = Ci_sunL;
    Citm1_shdH = Ci_shdH;
    Citm1_shdL = Ci_shdL;
    e_snotm1 = e_sno;
    In_Htm1 = In_H ;
    In_Ltm1 = In_L;
    In_Littertm1 = In_Litter;
    In_urbtm1 = In_urb;
    In_rocktm1 = In_rock;
    SWEtm1 = SWE ;
    In_SWEtm1 = In_SWE;
    Tdamptm1 = Tdamp;
    Tdptm1 = Tdp;
    Vx_Htm1 =  Vx_H;
    Vl_Htm1 = Vl_H;
    Vx_Ltm1 =  Vx_L;
    Vl_Ltm1 =  Vl_L;
    Psi_x_Htm1 = Psi_x_H;
    Psi_l_Htm1 = Psi_l_H;
    Psi_x_Ltm1 = Psi_x_L;
    Psi_l_Ltm1 = Psi_l_L;
    WATtm1=WAT;
    ICEtm1=ICE;
    IP_wctm1=IP_wc;
    ICE_Dtm1=ICE_D;
    Cicewtm1=Cicew;
    FROCKtm1=FROCK;
    t_slstm1= t_sls;
    rostm1=ros;
    SP_wctm1=SP_wc ;
    tau_snotm1 = tau_sno;
    EKtm1 = EK;
    Ws_undertm1 = Ws_under;
    Tdebtm1= Tdeb;
    Ticetm1= Tice;
    Oicetm1 = Oice;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Datam_S(4)==1)
        B_Htm1= B_H;
        PHE_S_Htm1=PHE_S_H;
        dflo_Htm1= dflo_H;
        AgeL_Htm1= AgeL_H;
        AgeDL_Htm1 = AgeDL_H;
        NPP_Htm1 =  NPP_H;
        NPPI_Htm1 = NPPI_H;
        TdpI_Htm1 = TdpI_H;
        Bfac_weekHtm1 = Bfac_weekH ;
        Nreserve_Htm1 = Nreserve_H;
        Preserve_Htm1 = Preserve_H;
        Kreserve_Htm1 = Kreserve_H;
        FNC_Htm1 = FNC_H;
        SAI_Htm1= SAI_H;
        hc_Htm1= hc_H;
        Tden_Htm1 =Tden_H;
        AgePl_Htm1 = AgePl_H;
        NBLeaf_Htm1=NBLeaf_H;
        PARI_Htm1 = PARI_H;
        NBLI_Htm1 = NBLI_H;
        NuLit_Htm1=NuLit_H;
        %%%
        B_Ltm1= B_L;
        PHE_S_Ltm1= PHE_S_L;
        dflo_Ltm1= dflo_L ;
        AgeL_Ltm1= AgeL_L ;
        AgeDL_Ltm1 = AgeDL_L;
        NPP_Ltm1 =  NPP_L;
        NPPI_Ltm1 = NPPI_L;
        TdpI_Ltm1 = TdpI_L;
        Bfac_weekLtm1 = Bfac_weekL ;
        Nreserve_Ltm1 = Nreserve_L;
        Preserve_Ltm1 = Preserve_L;
        Kreserve_Ltm1 = Kreserve_L;
        FNC_Ltm1 = FNC_L;
        SAI_Ltm1= SAI_L;
        hc_Ltm1= hc_L;
        Tden_Ltm1=  Tden_L;
        AgePl_Ltm1=  AgePl_L;
        NBLeaf_Ltm1=NBLeaf_L;
        PARI_Ltm1 = PARI_L;
        NBLI_Ltm1 = NBLI_L;
        NuLit_Ltm1=NuLit_L;
        %%%
        Ccrown_t_tm1 = Ccrown_t;
        NupI_Htm1 =  NupI_H;
        NupI_Ltm1 =  NupI_L;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% MEMORY CONDITION FOR VEGETATION  MODEL    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% 1 day
    if t > 24
        Ta_t(:,1:23)=Ta_t(:,2:24);
        Ta_t(:,24)=Ta_S;
        PAR_t(:,1:23)= PAR_t(:,2:24);
        PAR_t(:,24)=PARB_S + PARD_S;
        Tdp_t(:,:,1:23)=Tdp_t(:,:,2:24);
        Tdp_t(:,:,24)=Tdp;
        O_t(:,:,1:23)=O_t(:,:,2:24);
        O_t(:,:,24)=O;
        V_t(:,:,1:23)=V_t(:,:,2:24);
        V_t(:,:,24)=V;
        Pr_sno_t(:,1:23)=Pr_sno_t(:,2:24);
        Pr_sno_t(:,24)=Pr_sno;
        %%%
        An_H_t(:,:,1:23)=An_H_t(:,:,2:24);
        An_H_t(:,:,24)=An_H;
        An_L_t(:,:,1:23)=An_L_t(:,:,2:24);
        An_L_t(:,:,24)=An_L;
        Rdark_H_t(:,:,1:23)=Rdark_H_t(:,:,2:24);
        Rdark_H_t(:,:,24)=Rdark_H;
        Rdark_L_t(:,:,1:23)=Rdark_L_t(:,:,2:24);
        Rdark_L_t(:,:,24)=Rdark_L;
        Psi_x_H_t(:,:,1:23)=Psi_x_H_t(:,:,2:24);
        Psi_x_H_t(:,:,24)=Psi_x_H;
        Psi_l_H_t(:,:,1:23)=Psi_l_H_t(:,:,2:24);
        Psi_l_H_t(:,:,24)=Psi_l_H;
        Tdp_H_t(:,:,1:23)=Tdp_H_t(:,:,2:24);
        Tdp_H_t(:,:,24)=Tdp_H;
        Psi_x_L_t(:,:,1:23)=Psi_x_L_t(:,:,2:24);
        Psi_x_L_t(:,:,24)=Psi_x_L;
        Psi_l_L_t(:,:,1:23)=Psi_l_L_t(:,:,2:24);
        Psi_l_L_t(:,:,24)=Psi_l_L;
        Tdp_L_t(:,:,1:23)=Tdp_L_t(:,:,2:24);
        Tdp_L_t(:,:,24)=Tdp_L;
    else
        Ta_t(:,t)=Ta_S;
        PAR_t(:,24)=PARB_S + PARD_S;
        Tdp_t(:,:,t)=Tdp;
        O_t(:,:,t)=O;
        V_t(:,:,t)=V;
        An_H_t(:,:,t)=An_H;
        An_L_t(:,:,t)=An_L;
        Rdark_H_t(:,:,t)=Rdark_H;
        Rdark_L_t(:,:,t)=Rdark_L;
        Psi_x_H_t(:,:,t)=Psi_x_H;
        Psi_l_H_t(:,:,t)=Psi_l_H;
        Tdp_H_t(:,:,t)=Tdp_H;
        Psi_x_L_t(:,:,t)=Psi_x_L;
        Psi_l_L_t(:,:,t)=Psi_l_L;
        Tdp_L_t(:,:,t)=Tdp_L;
        Pr_sno_t(:,t)=Pr_sno;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%% OUTPUT WRITING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %run('OUTPUT_MANAGER_PAR');
    %%%%%%%%%
    if  mod(t,25)==0
        save(Fstep);
    end
    if  mod(t,8760)==0  ||  t==N_time_step
        Fstep2= strcat(Fstep,'_',num2str(t));
        save(Fstep2);
    end
end
%close(bau)
Computational_Time =toc;
%profile off
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MASS CHECK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('COMPUTATIONAL TIME [h] ')
disp(Computational_Time/3600)
disp(' COMPUTATIONAL TIME [s/cycle] ')
disp(Computational_Time/N_time_step)

%matlabpool close;