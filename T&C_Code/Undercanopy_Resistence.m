%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Undercanopy_Resistence     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%   - Mahat et al 2013
function[rap,ums]=Undercanopy_Resistence(Ws,Ta,Ts,hc,LAI,zatm,disp_h,zom,zom_under,SND)
%%%INPUTS
%Ta = %% air temperature [°C] --
%Ts =  ; %% surface temperature [°C] --
%z= %  reference height [m] % ---
%zom
%zoh_under = %% roughness  eddy diffusivities for heat  [m]  undercanopy land use
%zom_under = roughness eddy diffusivities for momentum [m] undercanopy land use
%hc = canpy height [m]
%%% Undercanopy resistence
%%% z = measurement height
%%% d = displacement height
%%% hc = canopy height
u = Ws; %% [m/s] wind speed  ---
% LAI = [Leaf Area Index ]
%%% OUTPUTS
%%% rap % [s/m] Undercanopy  Resistence
OPT=1 ;
if OPT == 1
    %%% Mahat et al 2013
    g=9.81; %%[m/s2]
    k= 0.41; %% Von Karman Constant
    d = disp_h ; %% Zero plane displacement [m]
    z=zatm; %% Measurement Height [m]
    us =  k*u/log((z-d)/zom); %%% Friction Velocity Canopy above [m/s]
    u_hc = (us/k)*log((hc-d)/zom); %% Wind Speed top Canopy [m/s]
    %%% alpha = 0.6 -1.5 -- alpha = LAI/2;
    %%% Match exponential and logarithmic
    alpha = log(u/u_hc)/(z/hc -1); %% Attenuation Coefficient
    %%% Yi et al 2008
    alpha = alpha*LAI*0.5/2; 
    % zom_under
    zms = 2 + SND ; %% reference height above the surface
    %%%
    ums = u_hc*exp(-alpha*(1-zms/hc)); %% Wind at reference height below canopy
    %%%
    Kh = k^2*u*(hc-d)/log((z-d)/zom);
    rap_n = hc*exp(alpha)/(Kh*alpha)*(exp(-alpha*(zms/hc))-exp(-alpha*((d+zom)/hc))) + 1/(k^2*ums)*log(zms/zom_under)^2;
    %%% Stability correction
    Ri = (g*(Ta-Ts)*zms)/(ums^2*(0.5*(Ta+Ts)+273.15)); %% [-]
    if Ri < 0 %% unstable
        rap = rap_n/(1-5*Ri)^(3/4);
    else %% Stable
        rap = rap_n/(1-5*Ri)^2;
    end
else
    %%% PARAMETERS
    k= 0.4; %% Von Karman Constant
    bet = 0.7*LAI; %% Par Zeng 2005
    d = disp_h ; % 0.67*hc; %% Zero plane displacement [m]
    z=zatm; %% Measurement Height [m]
    zoh_under = 0.1*zom_under; 
    %%% Hypothesis Logaritmic distribution of wind speed
    us =  k*u/log((z-d)/zom); %%% Friction Velocity  [m/s]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rt1 = (hc/(d*(bet+0.1)))*(1-exp(-bet*d/hc))*exp(bet); %% []
    rt2 = (hc/(d*bet))*(1-exp(-bet*d/hc))*exp(bet); % []
    %%%
    r1= (d/(k*(hc-d)))*rt1; %%% % []
    r2= (rt2^0.45)*log(zom_under/zoh_under)/k; %%% []
    Cs = 1/(r1 + r2); %% []
    rap = 1/(Cs*us);%%% [s/m]
    %%%%%%%%%% Sagakuchi and Zeng 2009
    %%%% Correction for stability
    Ta = Ta+273.15; %% air temperature [K] --
    Ts = Ts+273.15; %% surface temperature [K] --
    g=9.81; %%[m/s^2]
    %%%%%%%%%%%%%%
    if Ta > Ts %% Stability
        S=g*hc*(Ta-Ts)/(Ta*(us^2)); %% Correction factor
    else
        S=0;
    end
    rap = rap*(1+0.5*min(S,10));%%% [s/m]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ums = NaN; 
end
return