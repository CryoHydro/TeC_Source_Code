%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Temperature Humidity Biogeochemistry Zone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Se_bio,Se_fc,Psi_bio,Tdp_bio,Vtot,VT]=Biogeo_environment_old(Tdp,O,V,Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,Ofc,SPAR,Bio_Zs)
O=mean(O,1); Tdp=mean(Tdp,1);  V=mean(V,1); 
%%%% Biogeochemistry Layer Water Content
Obio = sum(Bio_Zs.*O);
Ohy = sum(Bio_Zs.*Ohy); 
Osat = sum(Bio_Zs.*Osat); 
Ofc = sum(Bio_Zs.*Ofc); 
Obio(Obio<=Ohy)= Ohy + 1e-5; 
Obio(Obio>=Osat)= Osat - 1e-5; 
Se_bio = (Obio-Ohy)./(Osat-Ohy);%%
Se_fc = (Ofc-Ohy)./(Osat-Ohy);%%
[Ko,Psi_bio]=Conductivity_Suction_old(SPAR,sum(Bio_Zs.*Ks_Zs),Osat,Ohy,...
    sum(Bio_Zs.*L),sum(Bio_Zs.*Pe),sum(Bio_Zs.*O33),sum(Bio_Zs.*alpVG),sum(Bio_Zs.*nVG),Obio);
Tdp_bio = sum(Bio_Zs.*Tdp); %% [°C]
Psi_bio=-(Psi_bio/1000)*1000*9.81/1e+6; %%[MPa]
Vtot = sum(V); %% [mm] Volume of entire depth 
VT = sum(Bio_Zs.*V); %% [mm] Volume in the Biogeochemistry layer
end
