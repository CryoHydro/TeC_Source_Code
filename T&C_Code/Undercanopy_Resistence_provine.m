%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Undercanopy_Resistence     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%   - Modified from Zeng et al., 2005
function[rap,ums]=Undercanopy_Resistence_new(Ws,Ta,Ts,hc,LAI,z,d,zom,zom_under,zoh_under,SND)
OPT = 2;
switch OPT
    case 1
        %%% MEANINGLESS FORMULATION 
        %%%INPUTS
        %Ta = %% air temperature [°C] --
        %Ts =  ; %% surface temperature [°C] --
        %z= %  reference height [m] % ---
        %zom
        %zoh_under = %% roughness  eddy diffusivities for heat  [m]  undercanopy land use
        %zom_under = roughness eddy diffusivities for momentum [m] undercanopy land use
        %hc = canpy height [m]
        u = Ws; %% [m/s] wind speed  ---
        % LAI = [Leaf Area Index ]
        %%% OUTPUTS
        %%% rap % [s/m] Undercanopy  Resistence
        %%% PARAMETERS
        k= 0.4; %% Von Karman Constant
        bet = 0.7*LAI; %% Par Zeng 2005
        %d = 0.67*hc; %% Zero plane displacement [m]
        %z=z+d; %% Measurement Height [m]
        %%% Hypothesis Logaritmic distribution of wind speed
        us =  k*u/log((z-d)/zom); %%% Friction Velocity Canopy above [m/s]
        u_hc = (us/k)*log((hc-d)/zom); %% Wind Speed top Canopy [m/s]
        alpha = log(u/u_hc)/(z/hc -1); %% Attenuation Coefficient
        %%%%%%
        zbel = 2+ SND; %% [m] Reference height below the canopy 2 meter above the snow
        u_zbel = u_hc*exp(alpha*(zbel/hc-1));
        us =  k*u_zbel/log((zbel)/zom_under); %%% Friction Velocity below canopy [m/s]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        rt1 = (hc/(d*(bet+0.1)))*(1-exp(-bet*d/hc))*exp(bet); %% []
        rt2 = (hc/(d*bet))*(1-exp(-bet*d/hc))*exp(bet); % []
        %%%
        r1= (d/(k*(hc-d)))*rt1; %%% % []
        r2= (rt2^0.45)*log(zom_under/zoh_under)/k; %%% []
        Cs = 1/(r1 + r2); %% []
        rap = 1/(Cs*us);%%% [s/m]
    case 2
        %%%%%%%%%%
        
        alpha=1.0; 
        %%% Yi et al 2008 
        alpha = LAI/2; 
        %%% Match exponential and logarithmic
        %alpha = log(u/u_hc)/(z/hc -1); %% Attenuation Coefficient
        
        %%% Mahat et al 2013
        g=9.81; %%[m/s2] 
        k= 0.4; %% Von Karman Constant
        %%y 1-2
        %%% alpha = 0.6 -1.5
        %d = hc*(0.05 + LAI^0.02/2 + (y-1)/20);
        %zom = hc*(0.23 - LAI^0.25/10 - (y-1)/67);
        % zom_under
        %%% Undercanopy resistence
        %%% z = measurement height
        %%% d = displacement height
        %%% hc = canopy height
        u = Ws; %% [m/s] wind speed  ---
        zms = 2+ SND ; %% reference height above the surface
        %%%
        us =  k*u/log((z-d)/zom); %%% Friction Velocity Canopy above [m/s]
        u_hc = (us/k)*log((hc-d)/zom); %% Wind Speed top Canopy [m/s]
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
        
    case 3
        % Ellis et al 2010 Essery et al 2003
         u = Ws; %% [m/s] wind speed  ---
         k= 0.4; %% Von Karman Constant
        us =  k*u/log((z-d)/zom); %%% Friction Velocity Canopy above [m/s]
        u_hc = (us/k)*log((hc-d)/zom); %% Wind Speed top Canopy [m/s]
        %%%%% Within Canopy wind speed
        alpha = 2.5; 
        zms = 2+ SND ; %% reference height above the surface
        %uwc = u_hc*exp(-gam*(1-z/hc)); %%% gam = k*LAI;  2.5
        ums = u_hc*exp(-alpha*(1-zms/hc)); %% Wind at reference height below canop
        rap=0;
end



