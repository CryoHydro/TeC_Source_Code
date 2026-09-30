%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% CHECK OF VARIABLE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FUNCTION
function check_var_3D(Ta,Tdew,Pr,Latm,WS,Pre,es,ea,RH,SAD1,SAD2,SAB1,SAB2,PARB,PARD)

%{
AUTHOR: MAXIMILIANO RODRIGUEZ
DATE: JULY/2025

EDITED CAT FYFFE
DATE 30/9/2026

Names of Variables must be inserted in this order:
    Ta = Temperature
    Tdew = Dew Point temperature
    Pr = Precipitation
    Latm = Downdward Long wave radiationIn_Littertm1
    WS = Wind speed
    Pre = Air pressure
    es = saturation vapor pressure
    ea = actual vapor pressure
    RH = Relative humidity
    SAD1 = SAD1
    SAD2 = SAD2
    SAB1 = SAB1
    SAB2 = SAB2
    PARB = PARB
    PARD = PARD
%}

%% CHECK FOR NANS
errors = 0;
mistakes = cell(1, 1); 
warnings_code = cell(1, 1); 
z = 1; 
z2 = 1; 

if sum(isnan(Ta), "all") > 0   
errors = errors + 1;    
mistakes{z,1} = "Error: There are NaNs in your forcings for Temperature";
z =z+1;
end

if sum(isnan(Tdew), "all") > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for Dew Point Temperature";
z =z+1;
end

if sum(isnan(Pr), "all") > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for Precipitation";
z =z+1;
end


if sum(isnan(Latm), "all") > 0  
mistakes{z,1} = "Error: There are NaNs in your forcings for Long wave radiation";
z =z+1;   
end

if sum(isnan(WS), "all") > 0 
mistakes{z,1} = "Error: There are NaNs in your forcings for Wind speed";
z =z+1;  
end

if sum(isnan(Pre), "all") > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for Air pressure";
z =z+1;      
end

if sum(isnan(es)) > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for saturation vapor pressure";
z =z+1; 
end


if sum(isnan(ea)) > 0 
mistakes{z,1} = "Error: There are NaNs in your forcings for actual vapor pressure";
z =z+1;   
end
%}

if sum(isnan(RH), "all") > 0  
mistakes{z,1} = "Error: There are NaNs in your forcings for relative humidity";
z =z+1; 
end

if sum(isnan(SAD1), "all") > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for SAD1";
z =z+1; 
end

if sum(isnan(SAD2), "all") > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for SAD2";
z =z+1; 
end

if sum(isnan(SAB1), "all") > 0   
mistakes{z,1} = "Error: There are NaNs in your forcings for SAB1";
z =z+1; 
end

if sum(isnan(SAB2), "all") > 0  
mistakes{z,1} = "Error: There are NaNs in your forcings for SAB2";
z =z+1; 
end

if sum(isnan(PARB), "all") > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for PARB";
z =z+1; 
end

if sum(isnan(PARD), "all") > 0 
mistakes{z,1} = "Error: There are NaNs in your forcings for PARD";
z =z+1; 
end


%% CHECK TEMPERATURES

Tmin = min(Ta, [], "all");
Tmax = max(Ta, [], "all");

TDPmin = min(Tdew, [], "all");
TDPmax = max(Tdew, [], "all");

if Tmin < -100
mistakes{z,1} = "Error: Temperature seems to not be in celsius degrees. There values of temperature lower than 100 celsius";
z =z+1; 
end

if Tmax > 100
mistakes{z,1} = "Error: Temperature seems to not be in celsius degrees. There values of temperature higher than 100 celsius";
z =z+1; 
end

if TDPmin < -100 
mistakes{z,1} = "Error: Dew Point Temperature seems to not be in celsius degrees. There values of Dew Point temperature lower than 100 celsius";
z =z+1; 
end

if TDPmax > 100
mistakes{z,1} = "Error: Dew Point Temperature seems to not be in celsius degrees. There values of Dew Point temperature higher than 100 celsius";
z =z+1; 
end

if any(Tdew > Ta, "all")
mistakes{z,1} = "Error: Some values of Dew Point Temperature are higher than Air Temperature. Check.";
z =z+1; 
end


%% CHECK PRECIPITATION

if any(Pr > 200, "all")
mistakes{z,1} = "Error: Precipitation seems to be very high. More than 200 mm in a single hour or Inf";
z =z+1;
end

if any(Pr < 0, "all")
mistakes{z,1} = "Error: Precipitation is negative";
z =z+1;
end

%% CHECK SOLAR COMPONENTS

if any(SAD1 > 1000, "all")
mistakes{z,1} = "Error: SAD1 from the radiation partition has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end

if any(SAD2 > 1000, "all")
mistakes{z,1} = "Error: SAD2 from the radiation partition has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end

if any(SAB1 > 1000, "all")
mistakes{z,1} = "Error: SAB1 from the radiation partition has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end

if any(SAB2 > 1000, "all")
mistakes{z,1} = "Error: SAB2 from the radiation partition has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end

if any(PARB > 1000, "all")
mistakes{z,1} = "Error: PARB from the radiation partition has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end

if any(PARD > 1000, "all")
mistakes{z,1} = "Error: PARD from the radiation partition has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end

if any(Latm > 1000, "all")
mistakes{z,1} = "Error: Latm (incoming longwave) has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end


%Less than zero
if any(SAD1 < 0, "all")
warnings_code{z,1} = "Warning: SAD1 from the radiation partition has negative values";
z =z+1;
end

if any(SAD2 < 0, "all")
warnings_code{z,1} = "Warning: SAD2 from the radiation partition has negative values";
z =z+1;
end

if any(SAB1 < 0, "all")
warnings_code{z,1} = "Warning: SAB1 from the radiation partition has negative values";
z =z+1;
end

if any(SAB2 < 0, "all")
warnings_code{z,1} = "Warning: SAB2 from the radiation partition has negative values";
z =z+1;
end

if any(PARB < 0, "all")
warnings_code{z,1} = "Warning: PARB from the radiation partition has negative values";
z =z+1;
end

if any(PARD < 0, "all")
warnings_code{z,1} = "Warning: PARD from the radiation partition has negative values";
z =z+1;
end

if any(Latm < 0, "all")
warnings_code{z,1} = "Warning: Latm (incoming longwave) has negative values";
z =z+1;
end


%% Vapor pressure
if any(ea> 10000, "all")
mistakes{z,1} = "Error: Actual vapor pressure must be in Pa. There are values over 10000 Pa in es." + ...
    " A Temperature over 50 C is needed to get ea>10000 Pa. Typical values are within 600 Pa and 5000 Pa";
z =z+1;
end

if any(es > 10000)
mistakes{z,1} = "Error: Saturated vapor pressure must be in Pa. There are values over 10000 Pa in ea. " + ...
    "A Temperature over 50 C is needed to get es>10000 Pa. Typical values are within 600 Pa and 5000 Pa";
z =z+1;
end

%% Wind speed
if any(WS > 100, "all")
mistakes{z,1} = "Error: There are wind speeds over 100 m/s. Wind speed must be in m/s";
z =z+1;
end

if any(WS < 0, "all")
mistakes{z,1} = "Error: There are negative values. Wind speed must be in m/s";
z =z+1;
end

%% Relative humidity
if any(RH > 1)
mistakes{z,1} = "Error: Relative humidity must be between 0 and 1 (as a fraction). RH is greater than 1.";
z =z+1;
end

if any(RH < 0) 
mistakes{z,1} = "Error: Relative humidity must be between 0 and 1 (as a fraction). RH is lower than 0.";
z =z+1;
end

%% FINAL WARNINGS
if length(warnings_code{1}) > 0
cellfun(@disp, warnings_code);
warning('check_var:Var', ['Warnings: Some data is out of range. Please check details']);
end

%% FINAL ERROR
if length(mistakes{1}) > 0
cellfun(@disp, mistakes);
error('check_var:Var', ['Error: Errors found in the forcings . Please check details']);
else
disp(['No errors found in your forcings'])    
end


end % end function

