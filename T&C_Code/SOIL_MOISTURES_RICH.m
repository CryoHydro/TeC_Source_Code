%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Richards Model Integration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dV]=SOIL_MOISTURES_RICH(t,V,...
Lk,f,EG,T_H,T_L,...
Qi_in,Slo_pot,IS,SPAR,Param) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Osat =Param.Osat; Ohy= Param.Ohy; O33 =Param.O33 ; 
dz =Param.dz; Ks_Zs = Param.Ks_Zs; Dz= Param.Dz; numn = Param.ms; 
L = Param.L ; Pe = Param.Pe; aR= Param.aR; aT = Param.aTop ; 
alpVG=Param.alpVG; nVG=Param.nVG; Zs=Param.Zs; 
cosalp = Param.cosalp; sinalp=Param.sinalp; SN=Param.SN;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Need only - Oi 
%%% --> numn=length(V);
O = V'./dz + Ohy ; 
I1 = O >= Osat -1e-5 ; I2 = O <= Ohy + 1e-5 ; 
O(I1)=Osat(I1); O(I2)=Ohy(I2)+ 1e-5; 
%O(O>Osat)=Osat; O(O<Ohy)=Ohy;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Soil Water
dV=zeros(numn,1);
K= zeros(1,numn);
P= zeros(1,numn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch SPAR
    case 1
        mVG= 1-1./nVG; Se = (O-Ohy)./(Osat-Ohy);
        %%% Van Genuchten  1980
        P = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        K=  Ks_Zs.*((Se).^(0.5)).*(1-(1-(Se).^(1./mVG)).^mVG).^2; %%% [mm/h]
    case 2
        B=1./L;
        A= exp(log(33)+B.*log(O33)); % Coefficient of moisture tension
        %%%%%%%%%%%%%%%%%%%%%%
        for i = 1:numn
            %[K(i),P(i)]=Conductivity_Suction(Ks_Zs(i),Osat(i),L(i),Pe(i),O33(i),O(i));
            %%%%%%%%%%%%%%%%%%%%
            K(i) = Ks_Zs(i)*(O(i)/Osat(i))^(3+(2/L(i))); %%% [mm/h]
            if O(i) < O33(i)
                P(i) = A(i).*(O(i)^-B(i)); %% [kPa]
            else
                P(i) =33-((O(i)-O33(i))*(33-Pe(i))/(Osat(i)-O33(i)));%% [kPa]
            end
            %%%%%%%%%%%%%%%
        end
        P =-101.9368*P;%%[mm]%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
To = K*aR.*dz;  %%% Trasmissivity  unsaturated/saturated [mm^2/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Slo_pot: Slope Total Hydraulic Head 
if SN==1
    %Qi_out = IS.*(To/aT).*cosalp.*(Slo_pot); %%%% [mm/h]
    Qi_out = IS.*(To/aT).*(sinalp); %%%% [mm/h]
else
    %Qi_out = (To/aT).*cosalp.*(Slo_pot); %%%% [mm/h]
    Qi_out = (To/aT).*(sinalp); %%%% [mm/h]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Qi_in [mm/h] Lateral Flow incoming 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%     Richards Model      %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Khalf(1:numn-1)= 0.5*(K(1:numn-1) + K(2:numn));  
%Khalf(1:numn-1)=exp(0.5*(log(K(1:numn-1)) + log(K(2:numn))));  
q=Khalf.*(1*cosalp - (P(2:numn)-P(1:numn-1))./Dz(2:end)); %% [mm/h] %%% Flux positive downward from layer i+1 (above) to i (below) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qf=(1-I1(2:numn)).*q; 
q(q>0)=qf(q>0); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% SOIL WATER BALANCE without other terms 
dV(1) =  f - q(1)  - T_H(1) - T_L(1) - EG(1) + Qi_in(1) - Qi_out(1);
for i =2:numn-1
    dV(i)= q(i-1) - q(i) - T_H(i) - T_L(i) - EG(i) + Qi_in(i) - Qi_out(i) ; %
end
dV(numn) = q(numn-1) - Lk  - T_H(numn) - T_L(numn) - EG(numn) + Qi_in(numn) - Qi_out(numn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dV=dV'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 