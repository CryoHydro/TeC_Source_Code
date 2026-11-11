%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% KYZYLSU - USER FORCING FUNCTION (OPT_Forcing = 4) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Distributed forcing as in the former launcher
% (Case_study/Kyzylsu_distributed/Launcher_Kyzylsu_Distributed.m):
%   - Air temperature: Pluvio + hourly lapse rate between Pluvio (3371 m)
%     and AWS_3900m (3900 m)
%   - Relative humidity, longwave, wind speed, precipitation, air pressure:
%     Pluvio, same value in all cells
%   - Shortwave (SAB1, SAB2, SAD1, SAD2, PARB, PARD): Pluvio, already
%     partitioned, flat surface (the terrain correction is applied by
%     Launcher_Dev_Distributed)
%   - Dew point temperature: from distributed Ta and uniform RH, same
%     formula as the former launcher
% Temperature modification over clean ice (Tmods = 0 in the former launcher)
% is left to the launcher (OPT_Tmod in Options_and_NoVegParams.xlsx).
%
% Called by Launcher_Dev_Distributed.m:
%   Forc_data = Kyzylsu_forcing('load', forc_dir)
%   F = Kyzylsu_forcing('maps', Forc_data, Date_t, Forc_grid)
% with Date_t the datetime of the time step and Forc_grid a structure with
% the cell vectors DTM, MASK, GLH, DEB_MAP and num_cell.
%==========================================================================

function out = Kyzylsu_forcing(mode, varargin)

switch mode

    case 'load'
        forc_dir = varargin{1};
        S1 = load(fullfile(forc_dir,'Pluvio.mat'));    S1 = S1.forcing_all;
        S2 = load(fullfile(forc_dir,'AWS_3900m.mat')); S2 = S2.forcing_all;
        if ~isequal(S1.Time, S2.Time)
            error('Kyzylsu_forcing: Pluvio and AWS_3900m do not have the same time steps')
        end
        out.Station_1 = S1;
        out.Station_2 = S2;
        out.Elev = [3371, 3900]; % Elevation of Pluvio and AWS_3900m [m a.s.l.]
        fprintf('Kyzylsu_forcing: station data from %s to %s\n', string(S1.Time(1)), string(S1.Time(end)));

    case 'maps'
        Forc_data = varargin{1}; Date_t = varargin{2}; G = varargin{3};
        S1 = Forc_data.Station_1; S2 = Forc_data.Station_2; Elev = Forc_data.Elev;

        % Index of the time step in the station data
        ind = find(S1.Time == Date_t);
        if numel(ind) ~= 1
            error('Kyzylsu_forcing: no station data for %s', string(Date_t))
        end

        % Air temperature: lapse rate between the two stations
        Ta_1 = S1.TA(ind); Ta_2 = S2.TA(ind);
        LR = (Ta_2 - Ta_1)/(Elev(2) - Elev(1)); % [C m-1]
        out.Ta = Ta_1 + LR.*(G.DTM - Elev(1));

        % Relative humidity [fraction]
        U = S1.RH(ind)/100;
        U(U<0) = 0; U(U>1) = 1;
        out.U = U;

        % Longwave radiation [W m-2]
        out.N = S1.LWIN(ind);

        % Wind speed [m s-1]
        Ws = S1.FF(ind);
        if isnan(Ws) || Ws < 0.01; Ws = 0.01; end
        out.Ws = Ws;

        % Precipitation [mm h-1]
        Pr = S1.PP(ind);
        if isnan(Pr) || Pr < 0.001; Pr = 0; end
        out.Pr = Pr;

        % Air pressure [hPa]
        out.Pre = S1.PRESS(ind);

        % Shortwave radiation [W m-2], flat surface
        out.SAB1 = S1.SAB1(ind); out.SAB2 = S1.SAB2(ind);
        out.SAD1 = S1.SAD1(ind); out.SAD2 = S1.SAD2(ind);
        out.PARB = S1.PARB(ind); out.PARD = S1.PARD(ind);

        % Dew point temperature [C], as in the former launcher
        a = 17.27; b = 237.3;
        xr = a.*out.Ta./(b+out.Ta) + log10(U);
        out.Tdew = b.*xr./(a-xr);

    otherwise
        error('Kyzylsu_forcing: unknown mode %s', mode)
end
end
