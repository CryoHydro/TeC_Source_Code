%==========================================================
% INITIALISING OUTPUTS
%==========================================================

%Script run on first day of the month to initialise the outputs

% function  Date_generator requires yy and mth as numbers
date_forMonth = Date_generator(str2double(yy),str2double(mth));
t_store = date_forMonth(end); % Date to store data

xdim = size(MASK, 1);
ydim = size(MASK, 2);
zdim = length(date_forMonth);

VAR_OUT = [... % Outputs for vegetation
           %---------------------------------------------------------------
           "ANPP_H";  % Above ground net primary production (for each Ccrown) daily
           "ANPP_L";
           "LAI_H";   % Leaf Area Index  (for each Ccrown)
           "LAI_L";
           "NDVI";    %Normalised difference vegetation index
           "NPP_H";   % Net Primary Production (for each Ccrown)
           "NPP_L";
           "hc_H";    % Vegetation Height (for each Ccrown)
           "hc_L"; 
           "An_H";    % Net assimilation Vegetation (for each Ccrown)
           "An_L";          
           %"B_H";     % Carbon pool biomass (for each Ccrown) Note *8 for each carbon pool
           %"B_L";     
           %"Bfac_dayL"; %Plant stress factor (daily) Low (for each Ccrown)
           %"Bfac_dayH"; %Plant stress factor (daily) High (for each Ccrown)
           %"Ca";      % CO2 atmospheric concentration
           %"Ci_shdH"; % CO2 shaded leaf internal concentration (for each Ccrown)
           %"Ci_shdL";      
           %"Ci_sunH"; % CO2 sunlit leaf internal concentration (for each Ccrown) 
           %"Ci_sunL";     
           %"LAIdead_H"; % Dead Leaf Area Index (for each Ccrown) 
           %"LAIdead_L";  
%            "GPP_H";      %Don't see in variable list?? need to know if
%            for each Ccrown
%            "GPP_L";        
           "RA_H";      % Autotrophic Respiration (for each Ccrown) 
           "RA_L";      
           %"Rdark_H";   % Leaf Dark Respiration (for each Ccrown) 
           %"Rdark_L";   
           "Rg_H";      % Growth Respiration (for each Ccrown) 
           "Rg_L";         
           %"Rmc_H";     % Maintenance Respiration Carbohydrate reserve (for each Ccrown) 
           %"Rmc_L";       
           %"Rmr_H";     % Maintenance Respiration roots (for each Ccrown) 
           %"Rmr_L";     
           %"Rms_H";     % Maintenance Respiration sapwood (for each Ccrown) 
           %"Rms_L";         
           "SAI_H";     % Stem Area Index (for each Ccrown) 
           "SAI_L";        
           %"SAT";   %Not in variable list
           "PHE_S_H";    % Phenology state (for each Ccrown) 
           "PHE_S_L";    
           %"TsV";        % Vegetation temperature for the energy balance -
           %supposedly doesn't exist?

          % Outputs for ET
          %--------------------------------------------------------------------------
          "EG";         % Evaporation from Bare soil (mm/h)
          "EICE";       % Evaporation/sublimation from Ice (mm/h)
          "EIn_H";      % Evaporation from intercepted water (mm/h) per Ccrown
          "EIn_L";        
          "EIn_rock";   % ET from rock  (mm/h)
          "EIn_urb";    % ET litter (mm/h)
          "EWAT";       % Evaporation from water and ponds (mm/h)
          "T_H";        % Transpiration (mm/h) per Ccrown
          "T_L";                  
          "ESN";        % Evap from snowpack from the ground (mm/h)
          "ELitter";    % ET litter (mm/h)
          "ESN_In";     % Evap from intercepted snow
          %"ET";         % Doesn't exist??
          %"ETen";       % Doesn't exist??

          % Outputs for snow
          %--------------------------------------------------------------------------
          %"ALB";             % Albedo - supposedly no snow albedo??
          "SND";             % Snow depth
          "SWE";             % Snow water equivalent
          "ros";             % Snow density
          "snow_albedo";
          "Smelt";           % Snow melt in w.e./h 
          "Imelt";           % Ice melt in w.e./h 
          "SP_wc";           % Snowpack water content
          "SWE_avalanched"; 
          "Cicew";           % Boolean operator for presence or absence of frozen water
          "Cice";            % Boolean operator for presence or absence of ice water
          "IP_wc";          % Ice pack water content
          %"NICe";           % New formed ice - doesn't seem to exist
          "ICE_D";          % Ice thickness
          "ICE";            % Ice water equivalent                   
          "U_SWE";          % Unloaded snow water equivalent from intercepted snow

          % Outputs for groundwater
        %--------------------------------------------------------------------------
        "FROCK";          % Storage in fractured rocks
        "Gfin";           % Ground Heat Flux heat diffusion
        "Lk_rock";        % Leakage rock surface to bedrock (recharge)
        "Lk";             % Bottom Leakage soil to bedrock (recharge)
        "Lk_wat";         % Leakage water pond to bedrock (recharge)
        "ZWT";            % Water table depth (mm)

        % Outputs for soil moisture
        %----------------------------------------------------------------------
        "O";            % Soil Moisture – Soil Water Content
        "OF";           % Soil Moisture First Soil Layer
        "OH";           % Soil Moisture available to roots (high)
        "OL";           % Soil Moisture available to roots (low)
        "OS";           % Soil Moisture for Bare Evaporation Layers
        
        % Outputs of runoff
        %------------------------------------------------------------------
        "QpointC";        % Discharge
        "q_runon";        % Runon
        "Rd";             % Saturation excess runoff
        "Rh";             % Infiltration excess runoff
        "f";              % Infiltration
        "In_rock";        % Intercepted water storage (rock)
        %"In";            % Intercepted water (storage) 
        "In_SWE";         % Intercepted snow water equivalent (storage)
        "In_L";           % Intercepted water storage (low) (per Ccrown)
        "In_H";           % Intercepted water storage (high)  (per Ccrown)
        %"Vice";         
        %"V";             % Volume of water stored in the soil layer
        %"WAT";           % Volume of water in the lakes/ponds
        %"WIS";           % Water flux incoming to the soil
        %"WR_IP";         % Water released from the ice pack
        %"WR_SP";         % Water released from the snow pack    

        % Meteorological and energy variables
        %--------------------------------------------------------------------------
        "Pr_liq";         % Liquid Precipitation
        "Pr_sno";         % Solid (snow)
        "Pr_St";          % Precipitation - input variable
        "Ta_St";          % Air temperature - input variable
        "Tdew_St";        % Dew Point temperature - input variable
        "Ws_St";          % Wind speed - input variable
        "U_St";           % Relative humidity - input variable
        "N_St";           % Cloudiness - input variable
        "Ts";             % Soil/snow Prognostic Temperature for the energy balance
        "ea_St";          % Vapor Pressure - input variable
        "Rn";             % Net radiation 
        "Pre_St";         % Atmospheric Pressure - input variable
        "SAD1_St";        % First band diffuse radiation - input variable
        "SAD2_St";        % Second band diffuse radiation - input variable
        "SAB1_St";        % First band direct radiation - input variable
        "SAB2_St";        % Second band direct radiation - input variable
        "QE";             % Latent Heat
        "H";              % Sensible Heat Flux
        "HV";             % Sensible heat flux from vegetation in presence of snow
        "G";              % Ground Heat Flux force restore method
        "Qfm";            % Heat for freezing or melting
        "Qv";             % Heat advected by Precipitation
        "dQ_S";             % Residual from the energy balance
        "dQVEG";          % Residual from the energy balance from snow free vegetation

        % Other potential variables
        %--------------------------------------------------------------------------
        %"CK1";            % Check on Mass Balance
        %"Csnow";          % Boolean operator for presence or absence of snow above frozen water 
        "Csno";           % Boolean operator for presence or absence of snow 
        "DQ_S";           % Residual of the energy budget
        %"Dr_H";           % Total Drainage from intercepted water
        %"Dr_L";           %
        %"Ds";             % Vapor Pressure Deficit
        %"DT_S";           % Residual temperature difference in the energy budget
        %"dw_SNO";         % Fraction of leaf covered by snow
        %"er";             % Splash erosion
        %"Inveg";         
        %"PAR";            
        %"Qlat_in";        
        %"Qlat_out";       
        %"Q_channel";      %This is output as a grid?  
        %"Tdamp";         % Soil/snow Temperature at Dampening depth
        "Tdp";            % Soil Temperature of the layer
        "Tdp_L";          % Soil temperature of the root zone (low)
        "Tdp_H";          % Soil temperature of the root zone (high)
        %"Tice";           %I don't think this exists
        %"TsVEG";         
          ];
% Names of variables
%--------------------------------------------------------------------------
        VAR_OUT_HOUR = VAR_OUT+'_hourly'; 
        VAR_OUT_DAY = VAR_OUT+'_daily'; 
        
        VAR_OUT_HOUR = [VAR_OUT_HOUR;
           "snow_albedo1_hourly"; "snow_albedo2_hourly"; "snow_albedo3_hourly"; "snow_albedo4_hourly"; 
           "O1_hourly"; "O2_hourly"; "O3_hourly"; "O4_hourly"; "O5_hourly";
           "O6_hourly"; "O7_hourly"; "O8_hourly"; "O9_hourly"; "O10_hourly"];
    
        VAR_OUT_DAY = [VAR_OUT_DAY;
       "snow_albedo1_daily"; "snow_albedo2_daily"; "snow_albedo3_daily"; "snow_albedo4_daily"; 
       "O1_daily"; "O2_daily"; "O3_daily"; "O4_daily"; "O5_daily";
       "O6_daily"; "O7_daily"; "O8_daily"; "O9_daily"; "O10_daily"];

% Creation of variables to store results per hour
%--------------------------------------------------------------------------
for i=1:length(VAR_OUT_HOUR)     
    if VAR_OUT_HOUR(i) == "QpointC_hourly" 
        QpointC_hourly = single(zeros(zdim, length(Xout))); % Discharge
    else
        commandStr = [char(VAR_OUT_HOUR(i)) ' = single(zeros(xdim , ydim, zdim));']; %Other matrices
        eval(commandStr);
    end
end

% Creation of variables to store results per day
%--------------------------------------------------------------------------
numdays = eomday(str2num(yy), str2num(mth));

for i=1:length(VAR_OUT_DAY)     
    if VAR_OUT_DAY(i) == "QpointC_daily"
        QpointC_daily = single(zeros(numdays, length(Xout))); % Discharge
    else
        commandStr = [char(VAR_OUT_DAY(i)) ' = single(zeros(xdim,ydim,numdays));']; %Other matrices
        eval(commandStr); 
    end
end
