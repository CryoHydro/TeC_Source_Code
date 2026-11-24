%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% FORCING FROM ONE STATION AND LAPSE RATES (OPT_Forcing = 1) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Distributes the data of one meteorological station over the catchment:
%   - Air temperature and dew point: lapse rate per hour and month (file)
%   - Precipitation: ratio = p1*elevation + p2 per month (file)
%   - Longwave: lapse rate of Marty et al. (2002)
%   - Direct shortwave (SAB1, SAB2, PARB): lapse rates (Alps values)
%   - Air pressure: from elevation
%   - Relative humidity: from the lapsed Ta and Tdew
%   - Wind speed and diffuse shortwave: same value in all cells
% The shortwave is partitioned once at load (Automatic_Radiation_Partition_I),
% which also gives t_bef and t_aft.
%
% Called by Launcher_Dev_Distributed.m (same interface as a user forcing function):
%   Forc_data = forcing_station_lapse('load', forc_dir, Info)
%   F = forcing_station_lapse('maps', Forc_data, Dates, Forc_grid)
% Info: structure with CONFIG (CONFIG_vals), Lat, Lon, DeltaGMT.
% Settings read from CONFIG: forc_meteo_file, forc_meteo_name, ta_lapse_file,
% ta_lapse_name, pr_lapse_file, pr_lapse_name, Meteo_el.
% Outputs of 'maps': every variable as a grid of hours (rows) x cells (columns).
%==========================================================================

function out = forcing_station_lapse(mode, varargin)

switch mode

    case 'load'
        forc_dir = varargin{1}; Info = varargin{2}; C = Info.CONFIG;

        Meteo_data = load([forc_dir '/' C.forc_meteo_file], C.forc_meteo_name);
        Meteo_data = Meteo_data.(C.forc_meteo_name);
        out.Forcing_Date = Meteo_data.DateTime;
        out.Ta = Meteo_data.Ta;     %Air temperature, degrees C
        out.U = Meteo_data.RH;      %Relative humidity, %
        out.Ws = Meteo_data.WS;     %Wind speed, ms-1
        out.Pr = Meteo_data.Pr;     %Precipitation, mm
        SWin = Meteo_data.SWin;     %Incoming shortwave, Wm-2
        out.Nin = Meteo_data.LWin;  %Incoming longwave, Wm-2
        out.Meteo_el = C.Meteo_el;  %Elevation of the station

        %Lapse rates
        hm_Ta_lapse = load(C.ta_lapse_file); %Ta lapse rates per hour and month
        out.hm_Ta_lapse = hm_Ta_lapse.(C.ta_lapse_name); %The lapse rate column should be named 'Ta_lapse'
        Pr_lapse = load(C.pr_lapse_file); %Pr lapse - based on ratio per month
        out.Pr_lapse = Pr_lapse.(C.pr_lapse_name); %Fit of ratio = (p1 * elevation) + p2. Just multiply ratio by Pr.

        %Dew point temperature at the station
        a=17.27; b=237.3;
        Uf = out.U/100; %Change to fraction
        xr=a.*out.Ta./(b+out.Ta)+log10(Uf);
        out.Tdew=b.*xr./(a-xr);               %Presumed dewpoint temperature (C)

        %Run radiation partition once to get t_bef and t_aft, and the radiation for distribution
        GRAPH_VAR=0;
        Forcing_Date_num = datenum(out.Forcing_Date);
        [~,~,out.SAD1,out.SAD2,out.SAB1,out.SAB2,out.PARB,out.PARD,~,~,out.t_bef,out.t_aft] = ...
            Automatic_Radiation_Partition_I(Forcing_Date_num,Info.Lat,Info.Lon,C.Meteo_el,Info.DeltaGMT,out.Pr,out.Tdew,SWin,GRAPH_VAR);

    case 'maps'
        D = varargin{1}; Dates = varargin{2}; G = varargin{3};
        DTMn = G.DTM; num_cell = G.num_cell; Meteo_el = D.Meteo_el;

        %Fixed lapse rates
        Ldown_lapse = -0.031; %from Marty et al. (2002) p145
        SAB1_lapse = 0.0049; %Wm-2 per m (from Simone, for the Alps, may be a little high for Andes), radiation increases with elevation
        SAB2_lapse = 0.0080; %Wm-2 per m (from Simone, for the Alps)
        PARB_lapse = 0.0047; %Wm-2 per m (from Simone, for the Alps)

        % Station data of the requested hours
        [found, id] = ismember(Dates(:), D.Forcing_Date);
        if ~all(found)
            error('forcing_station_lapse: no station data for %s', string(Dates(find(~found,1))))
        end
        Ta_P = D.Ta(id);
        N_P = D.Nin(id);
        Ws_P = D.Ws(id);
        Ws_P(Ws_P < 0.01) = 0.01; %Should not be zero, but double check
        Pr_P = D.Pr(id);
        Tdew_P = D.Tdew(id);
        SAD1_P = D.SAD1(id); SAD2_P = D.SAD2(id);
        SAB1_P = D.SAB1(id); SAB2_P = D.SAB2(id);
        PARB_P = D.PARB(id); PARD_P = D.PARD(id);

        % Distribution over the cells: time in rows, space in columns
        nt = numel(Dates);
        Ta_Sm = NaN(nt,num_cell); Pr_Sm = NaN(nt,num_cell); N_Sm = NaN(nt,num_cell);
        Tdew_Sm = NaN(nt,num_cell); SAB1_Sm = NaN(nt,num_cell); SAB2_Sm = NaN(nt,num_cell);
        PARB_Sm = NaN(nt,num_cell); Pr_ratio = NaN(nt,num_cell);

        for tinm = 1:nt % run over each hour
            m_temp = month(Dates(tinm));
            h_temp = hour(Dates(tinm));
            mh_id = D.hm_Ta_lapse.hour==h_temp & D.hm_Ta_lapse.month==m_temp;
            Ta_lapse = -D.hm_Ta_lapse.Ta_lapse(mh_id); % the lapse rate per hour and month - note converting to negative rate
            Tdew_lapse = Ta_lapse; %Its ok to use the same lapse rate
            %Air temperature, degrees C
            Ta_Sm(tinm,:) = Ta_P(tinm) + (Ta_lapse.*(DTMn-Meteo_el)); %Lapse rates should be negative
            %Precipitation, mm - This uses the WRF ratio v el relationship
            Pr_ratio(tinm,:) = (D.Pr_lapse(m_temp).month.p1.*DTMn) + D.Pr_lapse(m_temp).month.p2;  %val(x) = p1*x + p2
            Pr_Sm(tinm,:) = Pr_ratio(tinm,:).*Pr_P(tinm); %Apply ratio
            %Longwave radiation, W m-2
            N_Sm(tinm,:) = N_P(tinm) + (Ldown_lapse.*(DTMn-Meteo_el));
            %Dew point temperature, degrees C
            Tdew_Sm(tinm,:) = Tdew_P(tinm) + (Tdew_lapse.*(DTMn-Meteo_el)); %Lapse rates should be negative
            %Direct radiation with elevation (only if radiation values)
            if SAB1_P(tinm)>0
                SAB1_Sm(tinm,:) = SAB1_P(tinm) + (SAB1_lapse.*(DTMn-Meteo_el)); %Radiation should increase with elevation
            else
                SAB1_Sm(tinm,:) = 0;
            end
            if SAB2_P(tinm)>0
                SAB2_Sm(tinm,:) = SAB2_P(tinm) + (SAB2_lapse.*(DTMn-Meteo_el));
            else
                SAB2_Sm(tinm,:) = 0;
            end
            if PARB_P(tinm)>0
                PARB_Sm(tinm,:) = PARB_P(tinm) + (PARB_lapse.*(DTMn-Meteo_el));
            else
                PARB_Sm(tinm,:) = 0;
            end
        end
        %Prevent negative values
        SAB1_Sm(SAB1_Sm<0)=0; SAB2_Sm(SAB2_Sm<0)=0; PARB_Sm(PARB_Sm<0)=0;

        %Air pressure (calculated directly from elevation)
        Pre_temp = (101325*((1-(0.0065.*DTMn./288.15)).^(9.18*0.0289644./(8.31447*0.0065))))./100; %/100 converts to mbar
        out.Pre = repmat(Pre_temp(:)',nt,1); %Copy to all time steps

        %Use the lapsed dew point temperature to derive relative humidity [%]
        %(turned into a fraction by the launcher)
        c=237.3; b=17.27;
        out.U = 100*exp((c*b.*(Tdew_Sm - Ta_Sm))./((c+Ta_Sm).*(c+Tdew_Sm)));

        out.Ta = Ta_Sm; out.Pr = Pr_Sm; out.N = N_Sm; out.Tdew = Tdew_Sm;
        out.SAB1 = SAB1_Sm; out.SAB2 = SAB2_Sm; out.PARB = PARB_Sm;
        %Same value in all cells: wind speed and diffuse radiation
        out.Ws = repmat(Ws_P(:),1,num_cell); out.SAD1 = repmat(SAD1_P(:),1,num_cell);
        out.SAD2 = repmat(SAD2_P(:),1,num_cell); out.PARD = repmat(PARD_P(:),1,num_cell);

    otherwise
        error('forcing_station_lapse: unknown mode %s', mode)
end
end
