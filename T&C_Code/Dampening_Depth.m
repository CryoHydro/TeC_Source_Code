%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Soil Heat Dampening Depth     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Zdep]=Dampening_Depth(lan_dry,lan_s,cv_s,Osat,Ohy,Cwat,Curb,Crock)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Evaluated for a dry soil %%%%%%%%%%%%%
%%%INPUTS
%%% lan_dry=  Thermal conductivity dry soil [W/m K] 
%%% cv_s =  Volumetric heat capacity soil solid [J/m^3 K]
%Osat 
%Owp; %% Residual/Wilting Point Water Content  
%%REFERENCES 
%%% Boone et al., 2004; Hu and Islam 1995; 
%%%%%%%%%%%%%
%%%%%%%%%% THERMAL PROPERTIES SOIL 
tau= 86400; %% [s] time constant 
lan_wat = 0.6; %%% [W/m K ] Thermal conductivity water 
cv_w =  4186000;  % [J/m^3 K] Volumetric heat capcity water 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
if Cwat == 1
    KS = lan_wat/cv_w;   %% Thermal heat  Diffusivity Water [m^2/s]
else
    %%% Hip Ts>0 
    lan_sat = (lan_wat^Osat)*(lan_s^(1-Osat)); %% Saturated Conductivity [W/m K]
    Ke = log((Ohy)/Osat)+ 1; %% Kersten number
    Ke= Ke*(Ke>=0);
    lanS = Ke*lan_sat + (1-Ke)*lan_dry ; 
    %%%%%%%%%%%%%%%%%%%%%%
    cv_Soil = cv_s*(1-Osat) + Ohy*cv_w;  %  Volumetric heat capacity Soil  [J/m^3 K]
    %cs_Soil = cv_Soil/rsoil ; %%% [J/kg K]  %% Specific Heat Soil
    KS = lanS/cv_Soil; %%% [m^2/s]  %% Thermal heat  Diffusivity  Soil
end
%%%%%%%%%%%%%% COMPUTATION
w1= 2*pi/tau; %%% daily frequency of oscillation [1/s]
Zdep= sqrt(2*KS/w1); %%%  soil heat wave damping depth [m]
Zdep = Zdep*1000; %% [mm] 
return 
