%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Undercanopy_Resistence     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%   - Mahat et al 2013
function[rap_H,rap_L,rb_H,rb_L,ums,ums2]=Undercanopy_Leaf_Resistence(Ws,Ta,Ts,hc_H,hc_L,LAI_H,LAI_L,d_leaf_H,d_leaf_L,...
    zatm,disp_h,zom,zom_under,SND,disp_h_L,zom_L)
%%%INPUTS
%Ta = %% air temperature [°C] --
%Ts =  ; %% surface temperature [°C] --
%z= %  reference height [m] % ---
%zom
%zoh_under = %% roughness  eddy diffusivities for heat  [m]  undercanopy land use
%zom_under = roughness eddy diffusivities for momentum [m] undercanopy land use
%hc = canpy height [m]
%%% Undercanopy resistence
%%% zatm = measurement height
%%% disp_h = displacement height
%%% hc = canopy height
% LAI = [Leaf Area Index ]
u = Ws; %% [m/s] wind speed  ---
g=9.81; %%[m/s2]
k= 0.41; %% Von Karman Constant
%%% OUTPUTS
%%% rap % [s/m] Undercanopy  Resistence
%%% rb [s/m] Leaf Boundary Layer Resistence
%%%%%%%%%
%%% Generic Plot --
%z = zatm;
%d = disp_h;
%zom
%zom_under
%%%%
%%%%%
d = disp_h ; %% Zero plane displacement [m]
z = zatm; %% Measurement Height [m]
%%%%%
if (hc_L>0) && (hc_H>0)  %%% Two vegetations
    if SND >= (d+zom)
        rap_H=0;
        rb_H=0;
        rap_L=0;
        rb_L=0;
        ums=Ws; ums2=0;
    else
        %%% Aerodynamic lower vegetation in case of two vegetations
        d2 = disp_h_L ;
        zom2  = zom_L;
        hc2 =max([hc_L,d2+zom2+0.01,0.05]);
        %%%%%
        hc = max([hc_H,d+zom+0.01,0.05]);
        %%%%%
        us =  k*u/log((z-d)/zom); %%% Friction Velocity Canopy above [m/s]
        u_hc = (us/k)*log((hc-d)/zom); %% Wind Speed top Canopy [m/s]
        alpha  = log(u/u_hc)/(z/hc -1); %% Attenuation Coefficient
        %%% Yi et al 2008
        alpha  = alpha*LAI_H*0.5/2;
        %%%%%
        zms = 2 + max(hc2,SND) ; %% reference height above the surface
        if zms > d+zom
            zms=d+zom;
        end
        %%%
        ums = u_hc*exp(-alpha*(1-zms/hc)); %% Wind at reference height for Low vegetation
        %%%
        Kh = k^2*u*(hc-d)/log((z-d)/zom);
        rap_n = hc*exp(alpha)/(Kh*alpha)*(exp(-alpha*(zms/hc))-exp(-alpha*((d+zom)/hc))) + 1/(k^2*ums)*log(zms/zom2)^2;
        %%% Stability correction
        Ri = (g*(Ta-Ts)*zms)/(ums^2*(0.5*(Ta+Ts)+273.15)); %% [-]
        Ri(Ri>0.16)=0.16; %% Max. Stability 
        if Ri < 0 %% unstable
            rap = rap_n/((1-5*Ri)^(3/4));
        else %% Stable
            rap = rap_n/((1-5*Ri)^2);
        end
        rap_H=rap;
        rb_H=Leaf_BR(u_hc,Ts,Ta,d_leaf_H,alpha);
        %%%%%%%%%%%%%%%%%%%%%%%%
        if SND > (d2+zom2)
            rap_L=0;
            rb_L=0;
            ums2=Ws;
        else
            us2 =  k*ums/log((zms-d2)/zom2); %%% Friction Velocity Canopy above [m/s]
            u_hc2 = (us2/k)*log((hc2-d2)/zom2); %% Wind Speed top Canopy [m/s]
            alpha2  = log(ums/u_hc2)/(zms/hc2 -1); %% Attenuation Coefficient
            %%% Yi et al 2008
            alpha2  = alpha2*LAI_L*0.5/2;
            %%%%%
            zms2 = 2 +SND ; %% reference height above the surface
            if zms2 > d2+zom2
                zms2=d2+zom2;
            end
            %%%
            ums2 = u_hc2*exp(-alpha2*(1-zms2/hc2)); %% Wind at reference height for Low vegetation
            %%%
            Kh2 = k^2*ums*(hc2-d2)/log((zms-d2)/zom2);
            rap_n2 = hc2*exp(alpha2)/(Kh2*alpha2)*(exp(-alpha2*(zms2/hc2))-exp(-alpha2*((d2+zom2)/hc2))) + 1/(k^2*ums2)*log(zms2/zom_under)^2;
            %%% Stability correction
            Ri = (g*(Ta-Ts)*zms2)/(ums2^2*(0.5*(Ta+Ts)+273.15)); %% [-]
            Ri(Ri>0.16)=0.16; %% Max. Stability 
            if Ri < 0 %% unstable
                rap2 = rap_n2/((1-5*Ri)^(3/4));
            else %% Stable
                rap2 = rap_n2/((1-5*Ri)^2);
            end
            rap_L=rap2;
            rb_L=Leaf_BR(u_hc2,Ts,Ta,d_leaf_L,alpha2);
        end
    end
else
    if hc_H > 0 %% High Vegetation
        hc =max([hc_H,d+zom+0.01,0.05]);
        LAI = LAI_H;
        d_leaf = d_leaf_H;
    else %%% Low vegetation
        hc = max([hc_L,d+zom+0.01,0.05]);
        LAI = LAI_L;
        d_leaf = d_leaf_L;
    end
    if SND >= (d+zom)
        rb=0;
        rap=0;
        ums=Ws;
    else
        us =  k*u/log((z-d)/zom); %%% Friction Velocity Canopy above [m/s]
        u_hc = (us/k)*log((hc-d)/zom); %% Wind Speed top Canopy [m/s]
        %%% alpha = 0.6 -1.5 -- alpha = LAI/2;
        %%% Match exponential and logarithmic
        alpha = log(u/u_hc)/(z/hc -1); %% Attenuation Coefficient
        %%% Yi et al 2008
        alpha = alpha*LAI*0.5/2;
        % zom_under
        zms = 2 + SND ; %% reference height above the surface
        if zms > d+zom
            zms=d+zom;
        end
        %%%
        ums = u_hc*exp(-alpha*(1-zms/hc)); %% Wind at reference height below canopy
        %%%
        Kh = k^2*u*(hc-d)/log((z-d)/zom);
        rap_n = hc*exp(alpha)/(Kh*alpha)*(exp(-alpha*(zms/hc))-exp(-alpha*((d+zom)/hc))) + 1/(k^2*ums)*log(zms/zom_under)^2;
        %%% Stability correction
        Ri = (g*(Ta-Ts)*zms)/(ums^2*(0.5*(Ta+Ts)+273.15)); %% [-]
        Ri(Ri>0.16)=0.16; %% Max. Stability 
        if Ri < 0 %% unstable
            rap = rap_n/((1-5*Ri)^(3/4));
        else %% Stable
            rap = rap_n/((1-5*Ri)^2);
        end
        [rb]=Leaf_BR(u_hc,Ts,Ta,d_leaf,alpha);
    end
    if hc_H > 0 %% High Vegetation
        rap_H=rap;
        rb_H=rb;
        rap_L=0;
        rb_L=0;
        ums2=0;
    else %%% Low vegetation
        rap_L=rap;
        rb_L=rb;
        rap_H=0;
        rb_H=0;
        ums2=0;
    end
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Leaf_Boundary_Resistence   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%   - Vesala 1998 -- Chodhury and Monteith 1988 - Ivanov 2008
%%% Leigh et al 2011 ; Leuning et al 1995
function[rb]=Leaf_BR(u_hc,Ts,Ta,d_leaf,alpha)
%%%INPUTS
% Ta = %% air temperature [°C] --
% Ts =  ; %% surface temperature [°C] --
% Pre = %pressure [Pa]--
% z= %  reference height [m] % ---
% zoh = %% roughness  eddy diffusivities for heat  [m]
% zom = roughness eddy diffusivities for momentum [m]
% u = Ws; %% [m/s] wind speed  ---
% hc = canpy height [m]
% d_leaf = [cm] Leaf Dimension
% LAI = [Leaf Area Index ]
%%% OUTPUTS
%%% rb % [s/m] Leaf Boundary Resistence
%%% PARAMETERS
d_leaf = d_leaf/100; %% [m]
%k= 0.4; %% Von Karman Constant
a = 0.01; %% [m/s^0.5] --  Chodhury and Monteith 1988
%d = disp_h; %% Zero plane displacement [m]
%z= zatm ; %% Measurement Height [m]
%%% Hypothesis Logaritmic distribution of wind speed
%us =  k*u/log((z-d)/zom); %%% Friction Velocity  [m/s]
%u_hc = (us/k)*log((hc-d)/zom); %% Wind Speed top Canopy [m/s]
%%%
%alpha = log(u/u_hc)/(z/hc -1); %% Attenuation Coefficient
%alpha = 0.5*alpha*LAI/2;
%%% u(z) = u_hc*exp(alpha*((z/hc-1));
%%%%%%%%% Expression of Leaf Boundary Layer Resistence
%%% gb(z)= a*(u(z)/d_leaf)^0.5; %% [Jones 1992]  Integral between 0-LAI_TOT with LAI(z) = LAI_TOT*z/Hc
gb = (2*a/alpha)*((u_hc/d_leaf)^0.5)*(1-exp(-alpha/2)); %% [m/s]
%rb= (1/(gb*LAI)); %%   Leaf Boundary Layer Resistnce [s/m] Plant
%rb=1/gb; %%   Leaf Boundary Layer Resistnce [s/m] one-sided for unit leaf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Expression for free convection  (Leuning 1995, Monteith 1973)
Dh= 1.9*1e-5; %% [m2/s]
Gr= 1.6*1e+8*(Ts-Ta)*(d_leaf^3).*(Ts>Ta); %% [-]
gb_free = 0.5*Dh*Gr^(0.25)/d_leaf; %%[m/s]
%%%%%
gb=gb+gb_free;
rb=1/gb; %%   Leaf Boundary Layer Resistnce [s/m] one-sided for unit leaf
return