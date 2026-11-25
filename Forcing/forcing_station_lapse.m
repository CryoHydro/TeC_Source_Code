%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% FORCING FROM ONE STATION AND LAPSE RATES %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Selected in the site config with CONFIG_vals.forc_function = 'forcing_station_lapse'.
% Distributes the data of one meteorological station over the catchment:
%   - Air temperature and dew point: lapse rate per hour and month (file)
%   - Precipitation: ratio = p1*elevation + p2 per month (file)
%   - Longwave: lapse rate of Marty et al. (2002)
%   - Direct shortwave (SAB1, SAB2, PARB): lapse rates (Alps values)
%   - Air pressure: from elevation (barometric formula)
%   - Relative humidity: from the lapsed Ta and Tdew
%   - Wind speed and diffuse shortwave: same value in all cells
% The shortwave is partitioned once at load (Automatic_Radiation_Partition_I),
% which also gives t_bef and t_aft.
%
% Called by Launcher_Dev_Distributed.m (same interface as a user forcing function):
%   Forc_data = forcing_station_lapse('load', forc_dir, Info)
%   F = forcing_station_lapse('maps', Forc_data, Dates, Forc_grid)
% Info: structure with CONFIG (CONFIG_vals), Date (run period, not used here), Lat,
% Lon, DeltaGMT.
% Settings read from CONFIG: forc_meteo_file, forc_meteo_name, ta_lapse_file,
% ta_lapse_name, pr_lapse_file, pr_lapse_name, Meteo_el.
% Outputs of 'maps': every variable as a grid of hours (rows) x cells (columns).
% Names: _P = station values (one per hour), _Sm = grids hours x cells.
%==========================================================================

function out = forcing_station_lapse(mode, varargin)

% Constants of the dew point / relative humidity formula (Magnus), as in the launcher
a = 17.27; b = 237.3;

switch mode

    case 'load'
        forc_dir = varargin{1}; Info = varargin{2}; C = Info.CONFIG;

        % Station data
        Meteo_data = load([forc_dir '/' C.forc_meteo_file], C.forc_meteo_name);
        Meteo_data = Meteo_data.(C.forc_meteo_name);
        out.Forcing_Date = Meteo_data.DateTime;
        out.Ta = Meteo_data.Ta;     % Air temperature, degrees C
        out.U = Meteo_data.RH;      % Relative humidity, %
        out.Ws = Meteo_data.WS;     % Wind speed, m s-1
        out.Pr = Meteo_data.Pr;     % Precipitation, mm
        SWin = Meteo_data.SWin;     % Incoming shortwave, W m-2
        out.Nin = Meteo_data.LWin;  % Incoming longwave, W m-2
        out.Meteo_el = C.Meteo_el;  % Elevation of the station, m a.s.l.

        % Lapse rates
        L = load(C.ta_lapse_file);
        out.hm_Ta_lapse = L.(C.ta_lapse_name); % Ta lapse rates per hour and month (table with columns hour, month, Ta_lapse)
        
        L = load(C.pr_lapse_file);
        out.Pr_lapse = L.(C.pr_lapse_name);    % Pr ratio per month: ratio = p1*elevation + p2, multiplied by Pr

        % Dew point temperature at the station, degrees C
        xr = a.*out.Ta./(b+out.Ta) + log10(out.U/100);
        out.Tdew = b.*xr./(a-xr);

        % Radiation partition, once: shortwave components for the distribution, and t_bef, t_aft
        GRAPH_VAR = 0;
        Forcing_Date_num = datenum(out.Forcing_Date);
        [~,~,out.SAD1,out.SAD2,out.SAB1,out.SAB2,out.PARB,out.PARD,~,~,out.t_bef,out.t_aft] = ...
            Automatic_Radiation_Partition_I(Forcing_Date_num,Info.Lat,Info.Lon,C.Meteo_el,Info.DeltaGMT,out.Pr,out.Tdew,SWin,GRAPH_VAR);

        % All station series as columns (one row per hour), for the distribution in 'maps'
        series = {'Forcing_Date','Ta','U','Ws','Pr','Nin','Tdew','SAD1','SAD2','SAB1','SAB2','PARB','PARD'};
        for k = 1:numel(series)
            out.(series{k}) = out.(series{k})(:);
        end

    case 'maps'
        Data = varargin{1}; Dates = varargin{2}; Grid = varargin{3};
        nt = numel(Dates); num_cell = Grid.num_cell;
        z = Grid.DTM(:)';           % elevation of the cells, as a row
        dZ = z - Data.Meteo_el;     % elevation difference from the station, as a row

        % Fixed lapse rates
        Ldown_lapse = -0.031; % W m-2 per m, from Marty et al. (2002) p145
        SAB1_lapse = 0.0049;  % W m-2 per m (from Simone, for the Alps, may be a little high for Andes), radiation increases with elevation
        SAB2_lapse = 0.0080;  % W m-2 per m (from Simone, for the Alps)
        PARB_lapse = 0.0047;  % W m-2 per m (from Simone, for the Alps)

        % Station data of the requested hours (columns, as stored by 'load')
        [found, id] = ismember(Dates(:), Data.Forcing_Date);
        if ~all(found)
            error('forcing_station_lapse: no station data for %s', string(Dates(find(~found,1))))
        end
        Ta_P = Data.Ta(id); Tdew_P = Data.Tdew(id); Pr_P = Data.Pr(id); N_P = Data.Nin(id); Ws_P = Data.Ws(id);
        SAB1_P = Data.SAB1(id); SAB2_P = Data.SAB2(id); PARB_P = Data.PARB(id);
        SAD1_P = Data.SAD1(id); SAD2_P = Data.SAD2(id); PARD_P = Data.PARD(id);

        % Ta lapse rate of each hour (per hour and month), made negative
        T = Data.hm_Ta_lapse;
        [found, row] = ismember([hour(Dates(:)) month(Dates(:))], [T.hour T.month], 'rows');
        if ~all(found)
            error('forcing_station_lapse: no Ta lapse rate for %s', string(Dates(find(~found,1))))
        end
        Ta_lapse = -T.Ta_lapse(row);

        % Precipitation ratio coefficients of each hour (per month)
        p1_month = arrayfun(@(k) Data.Pr_lapse(k).month.p1, (1:12)'); % one slope per month (column)
        p2_month = arrayfun(@(k) Data.Pr_lapse(k).month.p2, (1:12)'); % one intercept per month (column)
        mth = month(Dates(:));
        p1 = p1_month(mth); p2 = p2_month(mth); % coefficients of each hour (column)

        % Distribution over the cells: hours (column) combined with cells (row) gives hours x cells
        Ta_Sm = Ta_P + Ta_lapse.*dZ;         % Air temperature, degrees C
        Tdew_Sm = Tdew_P + Ta_lapse.*dZ;     % Dew point temperature, degrees C (same lapse rate as Ta)
        Pr_Sm = (p1.*z + p2).*Pr_P;          % Precipitation, mm (ratio from the WRF ratio v elevation relationship)
        N_Sm = N_P + Ldown_lapse.*dZ;        % Longwave radiation, W m-2

        % Direct radiation, lapsed only in the hours with radiation at the station, never negative
        SAB1_Sm = SAB1_P + SAB1_lapse.*dZ; SAB1_Sm(~(SAB1_P > 0),:) = 0; SAB1_Sm(SAB1_Sm < 0) = 0;
        SAB2_Sm = SAB2_P + SAB2_lapse.*dZ; SAB2_Sm(~(SAB2_P > 0),:) = 0; SAB2_Sm(SAB2_Sm < 0) = 0;
        PARB_Sm = PARB_P + PARB_lapse.*dZ; PARB_Sm(~(PARB_P > 0),:) = 0; PARB_Sm(PARB_Sm < 0) = 0;

        % Air pressure from elevation (barometric formula), hPa, same at all hours
        P0 = 101325; T0 = 288.15; Lr = 0.0065; g = 9.81; M = 0.0289644; R = 8.31447;
        Pre_Sm = repmat(P0*(1 - Lr.*z./T0).^(g*M/(R*Lr))./100, nt, 1);

        % Relative humidity from the lapsed Ta and Tdew, % (turned into a fraction by the launcher)
        U_Sm = 100*exp((b*a.*(Tdew_Sm - Ta_Sm))./((b+Ta_Sm).*(b+Tdew_Sm)));

        % Outputs (the launcher sets Ws >= 0.01)
        out.Ta = Ta_Sm; out.Tdew = Tdew_Sm; out.U = U_Sm; out.Pr = Pr_Sm; out.N = N_Sm; out.Pre = Pre_Sm;
        out.SAB1 = SAB1_Sm; out.SAB2 = SAB2_Sm; out.PARB = PARB_Sm;
        
        % Same value in all cells: wind speed and diffuse radiation
        out.Ws = repmat(Ws_P,1,num_cell);
        out.SAD1 = repmat(SAD1_P,1,num_cell); 
        out.SAD2 = repmat(SAD2_P,1,num_cell);
        out.PARD = repmat(PARD_P,1,num_cell);

    otherwise
        error('forcing_station_lapse: unknown mode %s', mode)
end
end
