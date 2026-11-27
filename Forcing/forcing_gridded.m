%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% FULLY GRIDDED FORCING (OPT_Forcing = 2) %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reads monthly gridded forcing files (NOT compressed),
% one file per variable and month:
%   <forc_dir>/Bias_corrected/<forc_meteo_folder>/<yyyy>/<var>_<yyyy>_<m>.mat
%       for t2m, d2m, tp, N, PARB, PARD, SAB1, SAB2, SAD1, SAD2
%   <forc_dir>/Output_Downscaling/<forc_meteo_folder>/<yyyy>/<var>_<yyyy>_<m>.mat
%       for RH, sp, ws10, strd
% with <m> the month without leading zero. Each file holds a grid
% rows x columns x hours (all hours of the month); the name of the variable
% in each file is given by forc_meteo_name (14 names separated by spaces,
% same order as the variables above: t2m d2m tp RH sp ws10 strd N PARB PARD
% SAB1 SAB2 SAD1 SAD2).
%
% Called by Launcher_Dev_Distributed.m (same interface as a user forcing function):
%   Forc_data = forcing_gridded('load', forc_dir, Info)
%   F = forcing_gridded('maps', Forc_data, Dates, Forc_grid)
% Info: structure with CONFIG (CONFIG_vals), Lat, Lon, DeltaGMT.
% Settings read from CONFIG: forc_meteo_folder, forc_meteo_name.
% t_bef and t_aft are read by the launcher from CONFIG_vals.
% Outputs of 'maps': hours (rows) x cells (columns); U in %, Pre in Pa
% (converted by the launcher).
%==========================================================================

function out = forcing_gridded(mode, varargin)

switch mode

    case 'load'
        forc_dir = varargin{1}; Info = varargin{2}; C = Info.CONFIG;
        out.forc_meteo = [forc_dir '/Gridded_Downscaling/biasCorrected/'];
        out.forc_meteoOD = [forc_dir '/Gridded_Downscaling/downscaled/'];
        out.forc_names = split(C.forc_meteo_name);

    case 'maps'

        Forc_data = varargin{1}; Dates = varargin{2}; G = varargin{3};
        yy = char(num2str(year(Dates(1))));
        mth = char(num2str(month(Dates(1)),'%02d'));
        %So Dates is the date required by the run
        
        %Load variables
        forc_vars = ["Ta" "RH" "Pre" "Pr" "Ws" "N" "PARB" "PARD" "SAB1" "SAB2" "SAD1" "SAD2" "Swin"];
        num_forc_vars = size(forc_vars,2);
        forc_date = strcat(yy,'_',mth);
        File_names = strcat(forc_vars,'_',forc_date,'.nc')';
    
         for v = 1:num_forc_vars
             var_file = strcat(Forc_data.forc_meteo,forc_vars{v},'/',File_names(v));
             Forc_S.(forc_vars{v}) = ncread(var_file,forc_vars{v});
             [m,n,q] = size(Forc_S.(forc_vars{v}));
             Forc_S.(forc_vars{v}) = transpose(reshape(flipud(Forc_S.(forc_vars{v})),[m*n,q])); %flipud here to match T&C grids. Transpose so each row is a grid for one time step
         end

        % Sort time - presume time is the same for all vars
        Forcing_Date = ncread(var_file{1},'time');
        t0 = datetime(1970, 1, 1, 0, 0, 0);
        Forcing_Date = t0 + hours(Forcing_Date);
        id_date_for_month = ismember(Forcing_Date,Dates); % Index of the requested hours in the forcing data of this month

        % Bring into correct names; only the requested hours
        out.Ta = double(Forc_S.Ta(id_date_for_month,:)); % Air temperature, degrees C
        out.Pr = double(Forc_S.Pr(id_date_for_month,:)); % Precipitation, mm
        out.U = double(Forc_S.RH(id_date_for_month,:)); % Relative humidity, %
        out.N = double(Forc_S.N(id_date_for_month,:)); % Incoming Longwave, Wm^-2
        out.Pre = double(Forc_S.Pre(id_date_for_month,:)); % Surface pressure
        out.Ws = double(Forc_S.Ws(id_date_for_month,:)); % Wind speed, ms-1
        out.SAB1 = double(Forc_S.SAB1(id_date_for_month,:)); % First band direct radiation, Wm^-2
        out.SAB2 = double(Forc_S.SAB2(id_date_for_month,:)); % Second band direct radiation, Wm^-2
        out.SAD1 = double(Forc_S.SAD1(id_date_for_month,:)); % First band diffuse radiation, Wm^-2
        out.SAD2 = double(Forc_S.SAD2(id_date_for_month,:)); % Second band diffuse radiation, Wm^-2
        out.PARB = double(Forc_S.PARB(id_date_for_month,:)); % PAR radiation direct, Wm^-2
        out.PARD = double(Forc_S.PARD(id_date_for_month,:)); % PAR radiation diffuse Wm^-2
        out.Swin = double(Forc_S.Swin(id_date_for_month,:)); % Incoming Shortwave radiation Wm^-2

    otherwise
        error('forcing_gridded: unknown mode %s', mode)
end
end
