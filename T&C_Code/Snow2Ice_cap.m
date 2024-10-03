%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Convert snow thicker than than a given threshold (SWE_max in mm) into ice  %
% Achille Jouberton, 21/07/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ICE,ICE_D,SWE,SND]=Snow2Ice_cap(Asur, SWE_max, SWEtm1, ICEtm1,ros)

%%%INPUTS

% SWE_max: maximum snow height at the end of a hydrological year in [mm w.e.]. The
% extra will be converted to ice
% ros: snow density
% Asur: 

row = 1000; % water density [kg/m^3]
roi = 916.2; %% ice density [kg/m^3]

if isscalar(SWEtm1)  % For the point-scale version

    SWE_toice = max(SWEtm1 - SWE_max,0); % Snow to be converted into ice [mm]
    newSWE = SWEtm1 - SWE_toice; % New snowpack SWE [mm w.e.]
    newSND = 0.001*(newSWE./Asur).*(row./ros); % New snow depth after removal of snow coverted to ice [m]

    ICE = ICEtm1 + SWE_toice;  %%% [mm] Icepack addition of snow converted to ice [mm]
    ICE_D= Asur.*ICE*(row/roi).*0.001; %% New Icepack Depth [m]

    SWE = newSWE;
    SND = newSND;

else  % for fully distributed distributed version

    SWE_toice = max(SWEtm1 - SWE_max,0); % Snow to be converted into ice [mm]
    newSWE = SWEtm1 - SWE_toice; % New snowpack SWE [mm w.e.]
    newSND = newSWE.*0;
    newSND(newSWE>0) = 0.001*(newSWE(newSWE>0)./Asur(newSWE>0)).*(row./ros(newSWE>0)); % New snow depth after removal of snow coverted to ice [m]

    ICE = ICEtm1 + SWE_toice;  %%% [mm] Icepack addition of snow converted to ice [mm]
    ICE_D= Asur.*ICE*(row/roi).*0.001; %% New Icepack Depth [m]

    SWE = newSWE;
    SND = newSND;
end 

