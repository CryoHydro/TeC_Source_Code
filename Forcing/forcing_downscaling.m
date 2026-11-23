%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% ON-THE-FLY DOWNSCALING OF ERA5-LAND (OPT_Forcing = 3) %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Downscales ERA5-Land to the model grid once per month with downscale_4_T_C
% (repository TC_Downscaling: downscale_4_T_C.m and src/, to be added to the
% MATLAB path), with optional bias correction and radiation partition.
% Files (on the path, e.g. in <site>/Forcing), with <meteo_name> = forc_meteo_name:
%   <meteo_name>_ERAdata.mat          ERA5-Land data (dateTimes, TLR, 3D fields)
%   <meteo_name>_GRIDdata.mat         DEM grid for the downscaling (demZs, demLats, ...)
%   <meteo_name>_BC_Coefficients.mat  monthly bias correction coefficients (COEFFS)
%
% Called by Launcher_Dev_Distributed.m (same interface as a user forcing function):
%   Forc_data = forcing_downscaling('load', forc_dir, Info)
%   F = forcing_downscaling('maps', Forc_data, Dates, Forc_grid)
% Info: structure with CONFIG (CONFIG_vals), Date (run period), Lat, Lon, DeltaGMT.
% Settings read from CONFIG: forc_meteo_name, OPT_RunBC.
% t_bef and t_aft are read by the launcher from CONFIG_vals.
% Outputs of 'maps': hours (rows) x cells (columns); U in %.
%==========================================================================

function out = forcing_downscaling(mode, varargin)

switch mode

    case 'load'
        forc_dir = varargin{1}; Info = varargin{2}; C = Info.CONFIG; %#ok<NASGU>
        Date = Info.Date;
        meteo_name = C.forc_meteo_name; % e.g. 'Hybas_Lvl7_Catchment_2070497830_250m'
        if exist('downscale_4_T_C','file') ~= 2
            error('forcing_downscaling: downscale_4_T_C not found; add the TC_Downscaling repository (and its src folder) to the MATLAB path')
        end

        % -------------------------------------------------------------------------
        % LOAD AND ADJUST ERA5-LAND INPUTS
        % -------------------------------------------------------------------------
        OPT_Run_BC = C.OPT_RunBC;
        dsERA = load([meteo_name '_ERAdata.mat']);
        % Ensure alignment of ERA data to model run period
        [~,ia,ib] = intersect(Date,dsERA.dateTimes);
        FN = fieldnames(dsERA);
        for vv = 1:length(FN)
            if length(size(dsERA.(FN{vv})))==3
                dsERA.(FN{vv}) = dsERA.(FN{vv})(:,:,ib);
            end
        end
        dsERA.dateTimes = dsERA.dateTimes(ib);
        dsERA.TLR = dsERA.TLR(ib);

        % -------------------------------------------------------------------------
        % LOAD AND ADJUST DEM GRID DATA FOR DOWNSCALING USE
        % -------------------------------------------------------------------------
        dsGRID = load([meteo_name '_GRIDdata.mat']);
        dsGRID = rmfield(dsGRID,{'dsLats','dsLons','dsZs','dsNames'}); % will remove later

        nYr = years(Date(end)-Date(1));
        [~,ia,ib] = intersect(dsERA.dateTimes,Date);
        dsERA.tpMean = nansum(dsERA.tp(:,:,ia),3)/nYr;

        dsERA.tpLr = 0.0003; % NEWER VERSIONS HAVE THIS ALREADY INSIDE dsERA.   SO COULD REMOVE IN FUTURE

        n = 1500;
        dsGRID.demCurv = getcurvature(dsGRID.demZs,n);
        dsGRID.pixelSize = abs(dsGRID.demLats(1,1)-dsGRID.demLats(2,1))*111.1e3; % (m)
        [dsGRID.demSlope,dsGRID.demAspect] = slopeaspect(dsGRID.demZs,dsGRID.pixelSize);

        BC = load([meteo_name '_BC_Coefficients.mat']);
        COEFFS = BC.COEFFS;

        % -------------------------------------------------------------------------
        % CONFIGURATION FIELDS FOR ON-THE-FLY GRID DOWNSCALING
        % -------------------------------------------------------------------------
        cfg.runMode = 'OTF'; %  'OTF' or 'save'
        cfg.dim = 'grid'; %  'grid', 'cluster', or 'points'
        cfg.DoubleOutput = 1; % output both downscaled AND bias-corrected data (if relevant)
        cfg.RunRadiationPartition = 1; % If desire to run radiation partition. 1 = yes
        cfg.Date = Date; % datetime vector at hourly resolution (replaced by the hours of the month in 'maps')
        cfg.UTC = Info.DeltaGMT;
        cfg.VarNames = {'TA','RH','PP','FF','SWIN','LWIN','PRESS'}; % cell array, e.g. {'TA','RH','PP','SWIN','FF'}
        cfg.VarDS = {'dsT2m','dsRH','dsTp','dsWs10','dsSsrd','dsStrd','dsSp'}; % output field names corresponding to VarNames
        cfg.VarOut = {'Ta','RH','Pr','Ws','Swin','N','Pre'}; % output field names corresponding to VarNames
        cfg.FullNames = {'Air Temperature','Relative Humidity','Precipitation','Wind Speed','Shortwave Radiation','Longwave Radiation','Air Pressure'}; % display names corresponding to VarNames

        cfg.dsLats = []; %: grid latitudes
        cfg.dsLons = []; %: grid longitudes
        cfg.dsZs = []; %: grid elevations

        cfg.RunBC = OPT_Run_BC; % Run Bias Correction? 1 = yes
        cfg.BC = COEFFS; % structure containing monthly EQM information.

        cfg.bcLats = COEFFS.bcLats; %: POI latitudes
        cfg.bcLons = COEFFS.bcLons; %: POI longitudes
        cfg.bcZs = COEFFS.bcZs; %: POI elevations

        cfg.foOut = []; % output directory (required if cfg.runMode = 'save')

        cfg.dsNames = []; % point names, required for dim = 'points' and save
        cfg.clusterNumber = []; % cluster IDs, required for dim = 'cluster' and save

        cfg.precipLimit = 0.1; % precipLimit;    % precipitation threshold (mm hr); default 0.1   - Provide as parameter OPT?
        cfg.windMinimum = 0.1; % windMinimum;   % wind-speed minimum (m s-1); default 0.3;                - Provide as parameter OPT?
        cfg.writeBiasCorr = true; % optional logical; default true

        out.dsERA = dsERA; out.dsGRID = dsGRID; out.cfg = cfg;

    case 'maps'
        Forc_data = varargin{1}; Dates = varargin{2};
        dsERA = Forc_data.dsERA; dsGRID = Forc_data.dsGRID; cfg = Forc_data.cfg;

        dsERAsub = dsERA; % Make a copy of the ERA5 data
        [~,ia,ib] = intersect(Dates,dsERA.dateTimes); % Find the overlapping timesteps
        if numel(ib) ~= numel(Dates)
            error('forcing_downscaling: no ERA5-Land data for some hours from %s to %s', string(Dates(1)), string(Dates(end)))
        end
        FN = fieldnames(dsERAsub); % get field names from the ERA5 data structure
        for vv = 1:length(FN) % Adjust the timesteps on the 3D matrices
            if length(size(dsERAsub.(FN{vv})))==3
                dsERAsub.(FN{vv}) = dsERAsub.(FN{vv})(:,:,ib);
            end
        end
        % adjust the timing of the date itself and the lapse rates
        dsERAsub.dateTimes = dsERAsub.dateTimes(ib);
        dsERAsub.TLR = dsERAsub.TLR(ib);
        cfg.Date = dsERAsub.dateTimes; % Hours of this month (radiation partition and monthly bias correction)

        % MAIN DOWNSCALING CODE
        %--------------------
        OTF = downscale_4_T_C(cfg,dsERAsub,dsGRID);
        %--------------------

        % Bias-corrected fields if available, otherwise the downscaled fields
        if isfield(OTF,'biasCorrected')
            P = OTF.biasCorrected; names_in = cfg.VarOut;
        else
            P = OTF.downscaled; names_in = cfg.VarDS;
        end
        rad = {'SAB1','SAB2','SAD1','SAD2','PARB','PARD'}; % from the radiation partition
        names_in = [names_in rad]; names_out = [cfg.VarOut rad];
        for ff = 1:numel(names_in)
            [m,n,q] = size(P.(names_in{ff}));
            V.(names_out{ff}) = transpose(reshape(flipud(P.(names_in{ff})),[m*n,q])); % time (rows) x cells (columns)
        end

        out.Ta = V.Ta; out.Pr = V.Pr; out.Ws = V.Ws; out.N = V.N; out.Pre = V.Pre;
        out.SAB1 = V.SAB1; out.SAB2 = V.SAB2; out.SAD1 = V.SAD1; out.SAD2 = V.SAD2;
        out.PARB = V.PARB; out.PARD = V.PARD;
        out.U = V.RH; % Relative humidity, % (turned into a fraction by the launcher)

        U = V.RH / 100;
        out.Tdew = (243.12 .* (log(U) + (17.62 .* V.Ta)./(243.12 + V.Ta))) ./ ...
            (17.62 - (log(U) + (17.62 .* V.Ta)./(243.12 + V.Ta)));

    otherwise
        error('forcing_downscaling: unknown mode %s', mode)
end
end
