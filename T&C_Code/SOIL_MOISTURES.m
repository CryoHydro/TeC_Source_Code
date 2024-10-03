%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Richards Model Integration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dV]=SOIL_MOISTURES(t,V,...
  Osat,Ohy,O33,dz,Ks_Zs,Dz,n,...
L,Pe,Lk,f,EG,T_H,T_L,...
Qi_in,aR,Slo_pot,aT) 
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Soil Water 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Need only - Oi 
%%% --> n=length(V);
O = V'./dz + Ohy ; 
O(O>Osat)=Osat; O(O<Ohy)=Ohy;  
%[Qi_out]=Lateral_Subsurface_Flow(V,n,dz,Ks_Zs,aR,L,Osat,Ohy,ZWT,Zs,Slo_top,Slo_acq,aT); 
%Qi_in [mm/h] Lateral Flow incoming 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Simplify Richards Model %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Zs(i) dz(i) Ko(i) Pe(i)  
dV=zeros(1,n);
%Ko= zeros(1,n);  q=zeros(1,n-1);
%Ko_half=zeros(1,n-1); Po=zeros(1,n); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for i = 1:n  
%   [Ko(i),Po(i)]=Conductivity_Suction(Ks_Zs(i),Osat,L,Pe,O33,O(i));  
%end 
Ko = Ks_Zs.*(O./Osat).^(3+(2/L)); %%% [mm/h] %%% Estimation Conductivty outside loop 
%%%%%%%%%%%%%%%%%%%%%%%%
[Qi_out]=Lateral_Subsurface_Flow(Ko*aR,dz,Slo_pot,aT); 
%%%%%%%%%%%% Input--> Po; Ko; Dz  
%for i=2:n 
        %Ko_half(i-1)= exp((log(Ko(i)) + log(Ko(i-1)))/2); %% Geometric Mean of Ks %%[mm/h] 
        %Ko_half(i-1)= (Ko(i)+ Ko(i-1))/2; %% Normal Mean of Ks %%[mm/h] 
        %Ko_half(i-1)= Ko(i-1); %%% Upper Layer Ks %%[mm/h] 
        %q(i-1) = Ko_half(i-1)*(1 + (Po(i)-Po(i-1))/Dz(i)); %%% [mm/h] Flux Between i and i+1 directed downward 
        %q(i-1) = Ko_half(i-1); %%% [mm/h] Flux Between i and i+1 direct downward only gravity effects 
%end
q=Ko; %%% Estimation flux outside loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Ko_half(i)  % hydraulic conductivty between O[i i-1] [mm/h]
%%%%%% Po(i)  % Tension at O(i) [mm] - Head (Positive) 
%%%%%% f(i) % infiltration [mm/h] 
%%%%%% Lk [mm/h] leakage between layer n and bedrock  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% SOIL WATER BALANCE
dV(1) =  f - q(1) - T_H(1) - T_L(1) - EG(1) + Qi_in(1) - Qi_out(1);
for i =2:(n-1)
    dV(i)= q(i-1) - q(i) - T_H(i) - T_L(i) - EG(i) + Qi_in(i) - Qi_out(i) ; %
end
dV(n) = q(n-1) - Lk - T_H(n) - T_L(n) - EG(n) + Qi_in(n) - Qi_out(n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dV=dV'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 