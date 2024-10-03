%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Thermal Mod %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dS]=SOIL_THERMAL_COND(t,S,lan,cv,dz,nz,Dz,G0,Gn,Tup,Tdown,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%
% S = % [J °C/m^2 K]
%%% Method of lines
%%% Soil Water
%dS=zeros(nz,1); % [J °C/m^2 s K]
dS = S*0; 
T=S./(0.001*dz)./cv; %% [°C]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opt == 1 
        G0=-lan(1).*((T(1)-Tup)./(0.001*Dz(1)));  %  [W °C/m^2 K]
else
        G0=G0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=-lan(2:nz).*((T(2:nz)-T(1:nz-1))./(0.001*Dz(2:nz))); %% [W °C/m^2 K] %%% Flux positive downward from layer i (above) to i+1 (below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% SOIL Thermal BALANCE
dS(1) =  G0 - q(1)  ;
for i =2:nz-1
    dS(i)= q(i-1) - q(i) ;
end
dS(nz) = q(nz-1) - Gn;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 