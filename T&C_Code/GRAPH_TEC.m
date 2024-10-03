%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  GRAPHICAL OUTPUT FUNCTION     %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C)ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('D:\T&C_setups\DATA_FOLDER\dtm_rietholzbach.mat')
ksv=reshape(VEG_CODE,numel(DTM),1);
%%%
load('D:\TeCbeta\RESULT_GRASS_RIETHOLZBACH\Long_7cm\OUTPUT_SIM_Rietholzbach_SPATIAL_262992.mat')
%%%%
MASK=ones(m_cell,n_cell); 
MASK(isnan(DTM))=0;
MASKn=reshape(MASK,numel(MASK),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstep=262992; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_cell = cellsize:cellsize:cellsize*length(x_cell);
y_cell = cellsize:cellsize:cellsize*length(y_cell);
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
Ccrown =[ 1 0 ; 0.5 0.5 ; 1  0 ; 1 0];  %% Ccrown fraction for PFT
cc_max = 2 ; %% max number 
Vcode = 4;  %% number of codes 
EVcode = [ 1 2 3 4 ];  %% code of each PFTs 

PSAN = 5*MASK;
PCLA = 48*MASK;
PORG = 10*MASK;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%% ENERGY 
figure(1001)
Rsw_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(Rsw_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Shortwave Radiation {n} [W m^{-2}]')

figure(1002)
QE_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(QE_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Latent Heat {n} [W m^{-2}]')

figure(1003)
H_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(H_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Sensible Heat {n}  [W m^{-2}]')

figure(1004)
Rn_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(Rn_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Net Radiation {n} [W m^{-2}]')

figure(1005)
Ts_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(Ts_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Surface Temperature [°C]')


%%%%%%%%% EVAPOTRANSPIRATION 

figure(1011)
T_H_spatial(MASKn==0)= NaN;
T_L_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*T_L_spatial,m_cell,n_cell)+reshape(8760*T_H_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Transpiration  {v} [mm yr^{-1}]')

figure(1012)
EG_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*EG_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Soil Evaporation {v} [mm yr^{-1}]')

figure(1013)
EWAT_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*EWAT_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Water Evaporation {v} [mm yr^{-1}]')

figure(1014)
EIn_H_spatial(MASKn==0)= NaN;
EIn_L_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*EIn_H_spatial,m_cell,n_cell)+reshape(8760*EIn_L_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Interception Evaporation {v} [mm yr^{-1}]')

figure(1015)
ih=imagesc(x_cell,y_cell,reshape(8760*(T_H_spatial+T_L_spatial+EG_spatial+EIn_L_spatial+EIn_H_spatial+ESN_spatial+EWAT_spatial),m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('TOTAL ET {v} [mm yr^{-1}]')

%%%%%%%%% SNOW 

figure(1021)
ESN_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*ESN_spatial ,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Snow evaporation/sublimation {v} [mm yr^{-1}]')

figure(1022)
SWE_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(SWE_spatial ,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Snow Water Equivalent {v} [mm]')

figure(1023)
In_SWE_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(In_SWE_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Intercepted Snow Water Equivalent {v} [mm]')

figure(1024)
Csno_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(Csno_spatial/tstep,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Fraction of time with snow cover [-]')

figure(1020)
WR_SP_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*WR_SP_spatial ,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Water released from Snowpack {n} [mm yr^{-1}]')

figure(1021)
SND_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(SND_spatial ,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Snow depth {n}  [m]')

%%%%%%%%% SURFACE HYDROLOGY 

figure(1031)
f_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*f_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Infiltration {n} [mm yr^{-1}]')

figure(1032)
Rd_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*Rd_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Saturation Excess Runoff {v} [mm yr^{-1}]')

figure(1033)
Rh_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*Rh_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Infiltration Excess Runoff {v} [mm yr^{-1}]')

figure(1034)
Q_channel_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,Q_channel_spatial); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Channel Discharge [mm]')

figure(1035)
q_runon_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*q_runon_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Runon {v} [mm yr^{-1}]')

%%%%%%%%% SUBSURFACE HYDROLOGY 

figure(1041)
Qlat_out_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*Qlat_out_spatial ,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Lateral flow Out {v} [mm yr^{-1}]')

figure(1042)
Qlat_out_spatial(MASKn==0)= NaN;
Qlat_in_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*(Qlat_in_spatial-Qlat_out_spatial),m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Lateral flow Net {v} [mm yr^{-1}]')

figure(1043)
Lk_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(8760*Lk_spatial ,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Bedrock Leakage {v} [mm yr^{-1}]')

Pss = [100]; Pwp = [1000]; %%% [kPa]
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
[Osat,L,Pe,Ks,O33]=Soil_parameters_spatial(PSAN/100,PCLA/100,PORG/100);
[Ofc,Oss,Owp,Ohy]=Soil_parametersII_spatial(Osat,L,Pe,Ks,O33,Kfc,Pss,Pwp,Phy);
%%%%%%%%%%%%%%%%%%%%
figure(1044)
O_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,(-Ohy+reshape(O_spatial,m_cell,n_cell))./(Osat-Ohy)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Effective Saturation Soil Profile [-]')

figure(1045)
SAT_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(SAT_spatial/tstep ,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Saturation Time Fraction [-]')

figure(1046)
ZWT_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(-ZWT_spatial ,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Water table depth {n} [mm]')


%%%%%%%%% VEGETATION 
S_LAI=MASKn*0; 
S_ANPP =MASKn*0;
S_NPP= MASKn*0;
S_GPP = MASKn*0;
S_An= MASKn*0;
%%%%%%%%%%%%
for j=1:cc_max;
    for i=1:Vcode
        r=EVcode(i);
        S_LAI = S_LAI + Ccrown(i,j)*(ksv==r).*LAI_H_spatial(:,j)+ Ccrown(i,j)*(ksv==r).*LAI_L_spatial(:,j);
        S_ANPP = S_ANPP  + Ccrown(i,j)*(ksv==r).*ANPP_H_spatial(:,j)+ Ccrown(i,j)*(ksv==r).*ANPP_L_spatial(:,j);
        S_NPP = S_NPP + Ccrown(i,j)*(ksv==r).*NPP_H_spatial(:,j)+ Ccrown(i,j)*(ksv==r).*NPP_L_spatial(:,j);
        S_GPP = S_GPP + Ccrown(i,j)*(ksv==r).*GPP_H_spatial(:,j)+ Ccrown(i,j)*(ksv==r).*GPP_L_spatial(:,j);
        S_An = S_An + Ccrown(i,j)*(ksv==r).*An_H_spatial(:,j)+ Ccrown(i,j)*(ksv==r).*An_L_spatial(:,j);
    end
end


figure(1051)
S_LAI(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(S_LAI,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Leaf Area Index {n}  [-]')


figure(1052)
S_ANPP(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(365*S_ANPP,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('ANPP {n}  [gC m^{-2} yr^{-1}]')


figure(1053)
S_GPP(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(365*S_GPP,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('GPP {n}  [gC m^{-2} yr^{-1}]')

figure(1054)
S_NPP(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(365*S_NPP,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('NPP {n}  [gC m^{-2} yr^{-1}]')

figure(1055)
S_An(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(12*31.536*S_An,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Net Assimilation {n}  [gC m^{-2} yr^{-1}]')


%%%%%%%%% METEO INPUT 

figure(1061)
Ta_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(Ta_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Air Temperature  [°C]')

figure(1062)
Ws_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(Ws_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Wind speed  [m s^{-1}]')

figure(1063)
ea_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(ea_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Vapor Pressure [Pa]')

figure(1064)
Ca_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(Ca_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('CO_2 concentration [ppm]')

figure(1065)
Pr_spatial(MASKn==0)= NaN;
ih=imagesc(x_cell,y_cell,reshape(Pr_spatial,m_cell,n_cell)); axis xy ; axis equal;
set(ih,'alphadata',MASK)
colorbar ; xlabel('[m]'); ylabel('[m]')
title('Precipitation [mm h^{-1}]')



