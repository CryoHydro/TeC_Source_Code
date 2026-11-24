%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% FULLY GRIDDED FORCING (OPT_Forcing = 2) %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reads monthly gridded forcing files (compressed with Data_compressor_grid),
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
        out.forc_meteo = [forc_dir '/Bias_corrected/' C.forc_meteo_folder '/'];
        out.forc_meteoOD = [forc_dir '/Output_Downscaling/' C.forc_meteo_folder '/'];
        out.forc_names = split(C.forc_meteo_name);

    case 'maps'
        Forc_data = varargin{1}; Dates = varargin{2}; G = varargin{3};
        num_cell = G.num_cell;
        yy = char(num2str(year(Dates(1))));
        mth = char(num2str(month(Dates(1))));

        % Sort time
        StartMonth = datetime(str2double(yy),str2double(mth),1,0,0,0,0);
        EndMonth = dateshift(StartMonth,'end','month'); EndMonth = EndMonth + hours(23);
        Forcing_Date = StartMonth:hours(1):EndMonth; Forcing_Date = Forcing_Date';
        id_date_for_month = ismember(Forcing_Date,Dates); % Index of the requested hours in the forcing data of this month
        % Load variables
        forc_vars = ["t2m" "d2m" "tp" "RH" "sp" "ws10" "strd" "N" "PARB" "PARD" "SAB1" "SAB2" "SAD1" "SAD2"];
        forc_names = Forc_data.forc_names;
        OD_vars = contains(forc_vars,["RH","strd","sp","ws10"]);
        num_forc_vars = size(forc_vars,2);
        for v = 1:num_forc_vars
            if OD_vars(v)==0
                fn = char(strcat(Forc_data.forc_meteo,yy,'/',forc_vars(v), '_',yy,'_',mth,'.mat'));
                Forc_S.forc_vars{v} = load(fn);
            else
                fn = char(strcat(Forc_data.forc_meteoOD,yy,'/',forc_vars(v), '_',yy,'_',mth,'.mat'));
                Forc_S.forc_vars{v} = load(fn);
            end
        end

        % Then decompress, flip and reshape
        for v = 1:num_forc_vars
            %** Warning the decompress code changes the magnitude of the
            % variables so they should have been compressed before. **
            Forc_S_dc.(forc_vars{v}) = Data_compressor_grid(Forc_S.forc_vars{v}.(forc_names{v}),forc_vars(v),"decompress");
            Forc_S_dc.(forc_vars{v}) = flipud(Forc_S_dc.(forc_vars{v}));
            for ht = 1:size(Forcing_Date,1) % Probably this could be faster
                Forc_S_dc_r.(forc_vars{v})(ht,:) = reshape(Forc_S_dc.(forc_vars{v})(:,:,ht),[1,num_cell]); % Reshape to be time (rows) by grid cells (columns)
            end
        end
        % Bring into correct names; only the requested hours
        out.Ta = double(Forc_S_dc_r.t2m(id_date_for_month,:)); % Air temperature, degrees C
        out.Pr = double(Forc_S_dc_r.tp(id_date_for_month,:)); % Precipitation, mm
        out.U = double(Forc_S_dc_r.RH(id_date_for_month,:)); % Relative humidity, %
        out.N = double(Forc_S_dc_r.strd(id_date_for_month,:)); % Incoming Longwave (note N in the ERA data is cloudiness), Wm^-2
        out.Tdew = double(Forc_S_dc_r.d2m(id_date_for_month,:)); % Dew point temperature, degrees C
        out.Pre = double(Forc_S_dc_r.sp(id_date_for_month,:)); % Surface pressure
        out.Ws = double(Forc_S_dc_r.ws10(id_date_for_month,:)); % Wind speed, ms-1
        out.SAB1 = double(Forc_S_dc_r.SAB1(id_date_for_month,:)); % First band direct radiation, Wm^-2
        out.SAB2 = double(Forc_S_dc_r.SAB2(id_date_for_month,:)); % Second band direct radiation, Wm^-2
        out.SAD1 = double(Forc_S_dc_r.SAD1(id_date_for_month,:)); % First band diffuse radiation, Wm^-2
        out.SAD2 = double(Forc_S_dc_r.SAD2(id_date_for_month,:)); % Second band diffuse radiation, Wm^-2
        out.PARB = double(Forc_S_dc_r.PARB(id_date_for_month,:)); % PAR radiation direct, Wm^-2
        out.PARD = double(Forc_S_dc_r.PARD(id_date_for_month,:)); % PAR radiation diffuse Wm^-2

    otherwise
        error('forcing_gridded: unknown mode %s', mode)
end
end
