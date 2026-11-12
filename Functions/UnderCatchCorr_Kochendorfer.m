
function [PrecipOUT,CorrFactor] = UnderCatchCorr_Kochendorfer(PrecipIN,WindSpeed,AirTempC,PrecipH,DisplacementH,Z0,Shield)

% PrecipIN = un-corrected precipitation [mm]
% Wind Speed - wind speed [m s^1]
% AirTempC = air temperature [°C]
% PrecipH =  gauge height [m]
% WindH =  anemometer height [m]
% Z0 roughness length [m]
% DisplacementH = displacement height [m]


%% Wind speed at gauge height

% Logarithmic wind profile (Kochendorfer et al., 2017a)
WindCorr = WindSpeed .* log((PrecipH-DisplacementH)/Z0) ./ log((PrecipH-DisplacementH)/Z0);

%% Undercatch correction of Kochendorfer et al., 2017b

% Catch efficiency, eq. (3): CE = exp(-a*U*(1 - atan(b*Tair) + c)), with U the wind speed at gauge
% height [m/s] and Tair the air temperature [degC]. Coefficients [a b c]:
if strcmp(Shield,'none')
    coef = [0.0785 0.729 0.407];
elseif strcmp(Shield,'Single')
    coef      = [0.0348 1.366 0.779]; % single Alter shield
else
    coef      = [0.0348 1.366 0.779]; % single Alter shield
end

% The transfer function is only validated for U < 7.2 m/s at gauge height; U is not capped here.
CE_fun = @(U,Tair,c) exp(-c(1).*U.*(1 - atan(c(2).*Tair) + c(3)));

CorrFactor = NaN(size(PrecipIN));
CorrFactor = CE_fun(WindCorr, AirTempC, coef);
CorrFactor(CorrFactor>1) = 1; % no "overcatch" correction

PrecipOUT = PrecipIN ./ CorrFactor;