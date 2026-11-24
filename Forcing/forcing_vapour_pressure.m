function [esat, ea, Ds] = forcing_vapour_pressure(Ta, U)
%FORCING_VAPOUR_PRESSURE Saturation vapour pressure, vapour pressure and
% vapour pressure deficit from air temperature and relative humidity.
%
% Inputs (any size, e.g. time x cells):
%   Ta    air temperature [degC]
%   U     relative humidity [fraction]
% Outputs (same size) [Pa]:
%   esat  saturation vapour pressure
%   ea    vapour pressure
%   Ds    vapour pressure deficit (>= 0)

a = 17.27; b = 237.3;
esat = 611.*exp(a.*Ta./(b+Ta)); % Vapour pressure at saturation (Pa)
ea = U.*esat;                   % Vapour pressure (Pa)
Ds = esat - ea;                 % Vapour pressure deficit (Pa)
Ds(Ds<0) = 0;
end
