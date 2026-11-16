function [Ta, Tdew, U] = forcing_glacier_cooling(Ta, Tdew, U, Pre, Ws, SW, DTM, FPL)
%FORCING_GLACIER_COOLING Glacier cooling effect on air and dew point temperature
% (Shaw et al., 2025), applied to the forcing of the launcher (OPT_Decouple = 1).
%
% Inputs (time x cells, one row per hour, one column per cell):
%   Ta    air temperature [degC]
%   Tdew  dew point temperature [degC]
%   U     relative humidity [fraction]
%   Pre   air pressure [hPa]
%   Ws    wind speed [m s-1]
%   SW    incoming shortwave on a flat surface, SAD1 + SAB1 + SAD2 + SAB2 [W m-2]
% Grids (m_cell x n_cell, as in the launcher; cells in the same order as the columns):
%   DTM   elevation [m a.s.l.]
%   FPL   glacier flow path length (NaN outside glaciers)
%
% Outputs: Ta and Tdew multiplied by the cooling factor K (0.2 to 1) in the
% hours when all cells receive more than 30 W m-2 of shortwave, in the
% glacier cells with Ta > 0 degC; U recomputed from Ta and Tdew [fraction].
%
% The model coefficients are read from ModelCoefficients_Decoupling.mat
% (in Functions/, on the path).

C = load('ModelCoefficients_Decoupling.mat', 'DECOUPLE_COEFF');

SpecHum = convert_humidity(Pre*100, Ta+273.15, U*100, 'relative humidity', 'specific humidity');

D = DTM(:)'; F = FPL(:)';

for tt = 1:size(Ta,1)
    if min(SW(tt,:), [], 'all') > 30 % day time only (without very overcast conditions)

        T = Ta(tt,:); S = SpecHum(tt,:); FF = Ws(tt,:); Td = Tdew(tt,:);
        K = predict(C.DECOUPLE_COEFF, [T', S', FF', F', D'])'; % one factor per cell
        K(K<0.2) = 0.2; K(K>1) = 1; % limits of the observations

        ValidCooling = ~isnan(F) & ~isnan(D) & T > 0;

        Td(ValidCooling) = Td(ValidCooling) .* K(ValidCooling);
        T(ValidCooling) = T(ValidCooling) .* K(ValidCooling);

        Ta(tt,:) = T;
        Tdew(tt,:) = Td;
    end
end

% Relative humidity from the cooled Ta and Tdew [fraction]
U = exp((237.3*17.227.*(Tdew - Ta))./((237.3+Ta).*(237.3+Tdew)));
U(U>1) = 1; % Cannot be more than 100%
end
