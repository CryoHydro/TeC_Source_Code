% ========================================================================
%% INITIAL CONDITIONS just for multipoint
%==========================================================================

L_day=zeros(NNd,1);
for j=2:24:NN
    [~,~,~,~,~,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('L_day')

%% HYDROLOGY
%==========================================================================
% Initial conditions Hydrology/non-Vegetation
%==========================================================================
SWE(1)=SWEtm1(ij); %% [mm]
SND(1)=SNDtm1(ij); %% [m]
Ts(1)=Ta(1)+2;
Tdamp(1)=Ta(1);
Tdp(1,:)= Ta(1)*ones(1,ms);

%%% Snow_alb = soil_alb initial
snow_alb.dir_vis = snow_albedotm1(ij,1); %Was 0.6
snow_alb.dif_vis = snow_albedotm1(ij,2);
snow_alb.dir_nir = snow_albedotm1(ij,3);
snow_alb.dif_nir = snow_albedotm1(ij,4);

In_L(1,:)=In_Ltm1(ij); In_H(1,:)=In_Htm1(ij);
In_urb(1)=In_urbtm1(ij); In_rock(1)= In_rocktm1(ij);
In_Litter(1)=In_Littertm1(ij);
SP_wc(1)=SP_wctm1(ij) ; %%[mm]
In_SWE(1)= In_SWEtm1(ij);
ros(1)= rostm1(ij);
t_sls(1)= t_slstm1(ij);
e_sno(1) = e_snotm1(ij);
tau_sno(1) = tau_snotm1(ij);
EK(1)=EKtm1(ij);
WAT(1) = WATtm1(ij);% 
ICE(1) = ICEtm1(ij);% 
IP_wc(1)= IP_wctm1(ij);
ICE_D(1)= ICE_Dtm1(ij);%   ; 
FROCK(1)=FROCKtm1(ij);
Ws_under(1)=Ws_undertm1(ij); 
%%%%%%%%%%%%%% Volume [mm]
O(1,:)= Soil_Param.Ofc;
%%%%%%%%%%%%%%%%%%%
V(1,:) = Vtm1(ij,:); %(O(1,:) - Soil_Param.Ohy).*dz.*0; % TRYING WITH 0 initial water content in soil
%cd(cur_dir)
%%%%%%%%%%%%%%%%%

%% VEGETATION - Carbon
%==========================================================================
% Initial conditions Vegetation 
%==========================================================================

Ci_sunL(1,:) = [Ca(1)]; % [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; % [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; % [umolCO2/mol]
%%%%%%%%%%%%%%%%%%
%%% B1 Leaves - Grass  %%% B2 Sapwood  %%% B3 Fine Root  %%% B4 Carbohydrate Reserve
%%% B5 Fruit and Flower %%% B6 Heartwood - Dead Sapwood %%% B7 Leaves - Grass -- Standing Dead
%%%%%%%%%%%%%%%%%%
LAI_H(1,:) = LAI_Htm1(ij,:);  
Rrootl_H(1,:)= Rrootl_Htm1(ij,:); 
PHE_S_H(1,:)= PHE_S_Htm1(ij,:); 
dflo_H(1,:)=dflo_Htm1(ij,:); 
AgeL_H(1,:)=AgeL_Htm1(ij,:);
e_rel_H(1,:)=e_rel_Htm1(ij,:); 
hc_H(1,:) =hc_Htm1(ij,:); 
SAI_H(1,:) = SAI_Htm1(ij,:);

B_H(1,:,1)= full(B_Htm1{:,:,1}(ij,:));
B_H(1,:,2)= full(B_Htm1{:,:,2}(ij,:));
B_H(1,:,3)= full(B_Htm1{:,:,3}(ij,:));
B_H(1,:,4)= full(B_Htm1{:,:,4}(ij,:));
B_H(1,:,5)= full(B_Htm1{:,:,5}(ij,:));
B_H(1,:,6)= full(B_Htm1{:,:,6}(ij,:));
B_H(1,:,7)= full(B_Htm1{:,:,7}(ij,:));
B_H(1,:,8)= full(B_Htm1{:,:,8}(ij,:));

B_L(1,:,1)= full(B_Ltm1{:,:,1}(ij,:));
B_L(1,:,2)= full(B_Ltm1{:,:,2}(ij,:));
B_L(1,:,3)= full(B_Ltm1{:,:,3}(ij,:));
B_L(1,:,4)= full(B_Ltm1{:,:,4}(ij,:));
B_L(1,:,5)= full(B_Ltm1{:,:,5}(ij,:));
B_L(1,:,6)= full(B_Ltm1{:,:,6}(ij,:));
B_L(1,:,7)= full(B_Ltm1{:,:,7}(ij,:));
B_L(1,:,8)= full(B_Ltm1{:,:,8}(ij,:));

%%%%%%%%%%%%%%%%%%
LAI_L(1,:) = LAI_Ltm1(ij,:);  Rrootl_L(1,:)= Rrootl_Ltm1(ij,:);
PHE_S_L(1,:)=PHE_S_Ltm1(ij,:); dflo_L(1,:)=dflo_Ltm1(ij,:); AgeL_L(1,:)=AgeL_Ltm1(ij,:);
e_rel_L(1,:)=e_rel_Ltm1(ij,:); hc_L(1,:) =hc_Ltm1(ij,:); SAI_L(1,:) =SAI_Ltm1(ij,:);
%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%%%%%%%%%
RexmyI(1,:)= [0 0 0];
%%%%%%%%%%%%%%%%%%
Nreserve_H(1,:) = Nreserve_Htm1(ij,:); 
Preserve_H(1,:) = Preserve_Htm1(ij,:); 
Kreserve_H(1,:) = Kreserve_Htm1(ij,:);

NupI_H(1,:,1)= full(NupI_Htm1{:,:,1}(ij,:)); 
NupI_H(1,:,2)= full(NupI_Htm1{:,:,2}(ij,:)); 
NupI_H(1,:,3)= full(NupI_Htm1{:,:,3}(ij,:)); 

NupI_L(1,:,1)= full(NupI_Ltm1{:,:,1}(ij,:));
NupI_L(1,:,2)= full(NupI_Ltm1{:,:,2}(ij,:));
NupI_L(1,:,3)= full(NupI_Ltm1{:,:,3}(ij,:));

FNC_H(1,:)=FNC_Htm1(ij,:); 

Nreserve_L(1,:)= Nreserve_Ltm1(ij,:);
Preserve_L(1,:)= Preserve_Ltm1(ij,:); 
Kreserve_L(1,:)= Kreserve_Ltm1(ij,:);

FNC_L(1,:)=FNC_Ltm1(ij,:);
%%%%%%%%%%%%%%%%%%
TdpI_H(1,:)=TdpI_Htm1(ij,:); 
TdpI_L(1,:)=TdpI_Ltm1(ij,:);