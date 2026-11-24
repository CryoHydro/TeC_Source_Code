function [Ta, Tdew, U] = forcing_glacier_cooling(Ta, Tdew, U, Pre, Ws, SW, DTM, FPL, OPT_Decouple, Tmod)
%FORCING_GLACIER_COOLING Glacier cooling effect on air and dew point temperature,
% applied to the forcing of the launcher (OPT_Decouple > 0).
%
% Methods (OPT_Decouple):
%   1 = Shaw et al. (2025): Ta and Tdew multiplied by the cooling factor K (0.2 to 1)
%       in the hours when all cells receive more than 30 W m-2 of shortwave, over
%       clean ice (FPL not NaN) with Ta > 0 degC; no correction over debris-covered ice
%   2 = bias (former OPT_Tmod = 1): TmodB removed from Ta and Tdew over clean ice
%       where Ta > 0 degC
%   3 = bias and multiplier (former OPT_Tmod = 2): TmodB removed from Ta and Tdew over
%       clean ice, then Ta and Tdew multiplied by TmodM over clean ice where Ta > 0 degC
% In all methods U is recomputed from the corrected Ta and Tdew.
%
% Inputs (time x cells, one row per hour, one column per cell):
%   Ta    air temperature [degC]
%   Tdew  dew point temperature [degC]
%   U     relative humidity [fraction]
%   Pre   air pressure [hPa]                                         (method 1)
%   Ws    wind speed [m s-1]                                         (method 1)
%   SW    incoming shortwave on a flat surface, SAD1+SAB1+SAD2+SAB2 [W m-2] (method 1)
% Grids (m_cell x n_cell, as in the launcher; cells in the same order as the columns):
%   DTM   elevation [m a.s.l.]                                       (method 1)
%   FPL   glacier flow path length (NaN outside glaciers)            (method 1)
% OPT_Decouple  method, 1, 2 or 3
% Tmod  structure:
%   clean_ice  cells of clean ice (glacier without debris), one value per cell (all methods)
%   TmodB      bias [degC] (positive reduces the temperature over ice) (methods 2 and 3)
%   TmodM      multiplier [-] (method 3)
%
% Outputs: Ta, Tdew [degC] and U [fraction], time x cells.
%
% The coefficients of method 1 are read from ModelCoefficients_Decoupling.mat
% (in Functions/, on the path).

if nargin < 10 || ~isstruct(Tmod) || ~isfield(Tmod, 'clean_ice')
    error('forcing_glacier_cooling: the clean ice mask Tmod.clean_ice is needed')
end
clean_ice = logical(Tmod.clean_ice(:))'; % one value per cell, as a row

switch OPT_Decouple

    case 1 % Shaw et al. (2025)
        C = load('ModelCoefficients_Decoupling.mat', 'DECOUPLE_COEFF');

        SpecHum = convert_humidity(Pre*100, Ta+273.15, U*100, 'relative humidity', 'specific humidity');

        D = DTM(:)'; F = FPL(:)';

        for tt = 1:size(Ta,1)
            if min(SW(tt,:), [], 'all') > 30 % day time only (without very overcast conditions)

                T = Ta(tt,:); S = SpecHum(tt,:); FF = Ws(tt,:); Td = Tdew(tt,:);
                K = predict(C.DECOUPLE_COEFF, [T', S', FF', F', D'])'; % one factor per cell
                K(K<0.2) = 0.2; K(K>1) = 1; % limits of the observations

                ValidCooling = ~isnan(F) & ~isnan(D) & T > 0 & clean_ice; % clean ice only (not debris-covered)

                Td(ValidCooling) = Td(ValidCooling) .* K(ValidCooling);
                T(ValidCooling) = T(ValidCooling) .* K(ValidCooling);

                Ta(tt,:) = T;
                Tdew(tt,:) = Td;
            end
        end

        % Relative humidity from the cooled Ta and Tdew [fraction]
        U = exp((237.3*17.227.*(Tdew - Ta))./((237.3+Ta).*(237.3+Tdew)));
        U(U>1) = 1; % Cannot be more than 100%

    case {2, 3} % Bias (2), bias and multiplier (3) over clean ice
        if ~isfield(Tmod, 'TmodB')
            error('forcing_glacier_cooling: OPT_Decouple = %d needs Tmod.TmodB', OPT_Decouple)
        end
        idcli_r = repmat(clean_ice, size(Ta,1), 1); % clean ice, time x cells

        if OPT_Decouple == 2
            idapply = (Ta > 0) & idcli_r; % where the temperature threshold is exceeded, over clean ice
            Ta(idapply) = Ta(idapply) - Tmod.TmodB; % only remove bias
            Tdew(idapply) = Tdew(idapply) - Tmod.TmodB; % also applied to Tdew
        else
            if ~isfield(Tmod, 'TmodM')
                error('forcing_glacier_cooling: OPT_Decouple = 3 needs Tmod.TmodM')
            end
            Ta(idcli_r) = Ta(idcli_r) - Tmod.TmodB; % remove bias first (over clean ice)
            Tdew(idcli_r) = Tdew(idcli_r) - Tmod.TmodB; % also applied to Tdew
            idapply = (Ta > 0) & idcli_r; % multiplier where Ta > 0 after the bias, over clean ice
            Ta(idapply) = Ta(idapply) .* Tmod.TmodM;
            Tdew(idapply) = Tdew(idapply) .* Tmod.TmodM;
        end

        % Relative humidity from the corrected Ta and Tdew [fraction]
        c = 237.3; b = 17.27;
        U = exp((c*b.*(Tdew - Ta))./((c+Ta).*(c+Tdew)));
        U(U>1) = 1; % Cannot be more than 100%

    otherwise
        error('forcing_glacier_cooling: unknown OPT_Decouple = %g (1, 2 or 3)', OPT_Decouple)
end
end
