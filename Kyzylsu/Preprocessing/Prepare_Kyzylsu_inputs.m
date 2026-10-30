%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% KYZYLSU - INPUTS FOR Launcher_Dev_Distributed %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================================================================
% Converts the former Kyzylsu set-up (Case_study/Kyzylsu_distributed,
% Launcher_Kyzylsu_Distributed.m) into the files expected by
% Launcher_Dev_Distributed.m (same structure as the Shallap site folder):
%
%   Preprocessing/OUTPUTS/<PreProc_folder>/dtm_Kyzylsu_100m.mat
%   Multipoint/OUTPUTS/Kyzylsu_MultiPoints.txt
%   Forcing/Kyzylsu_Pluvio_meteo.mat        (station forcing, OPT_Forcing = 1)
%   Forcing/Kyzylsu_hm_Ta_lapse.mat         (hour x month Ta lapse rates)
%   Forcing/Kyzylsu_Pr_ratio_v_el_const.mat (precipitation ratio vs elevation)
%
% Run this script from the Kyzylsu/Preprocessing folder.
% Parameter tables are in Kyzylsu/Parameters (xlsx files).
%==========================================================================

clc; clear;

%% SETTINGS
%==========================================================================
PreProc_folder = 'Kyzylsu_100m_from_Kyzylsu_distributed'; % Must match PreProc_folder in Site_Config_TC.xlsx
dtm_in = 'INPUTS/dtm_Kyzylsu_100m.mat'; % dtm used by Launcher_Kyzylsu_Distributed
dtm_out = 'dtm_Kyzylsu_100m.mat'; % Must match dtm_file in Site_Config_TC.xlsx

% Initial snow depth and albedo (diff_IniSND = 1 in the former launcher)
use_ini_snow = 1; % 1 = replace SNOWD/SNOWALB from the dtm with the files below
fn_IniSnowDepth = 'INPUTS/Ini_SnowDepth/Kyzylsu_Init_SnowDepth_LSS2_1999-10-16.mat';
fn_IniSnowAlbedo = 'INPUTS/Ini_SnowDepth/Kyzylsu_Init_SnowAlbedo_LSS2_1999-10-16.mat';

% Stations
Elev = [3371, 3900]; % Elevation of Pluvio and AWS_3900m stations (Meteo_el = Elev(1))
fn_station_1 = '../Forcing/Pluvio.mat';
fn_station_2 = '../Forcing/AWS_3900m.mat';

% Soil thickness: in the former PARAMETERS_SOIL.m SOIL_TH was a percentage of
% Zs(end) = 1000 mm. Here it is converted to cm (OPT_soil_th = 2).
Zs_end_former = 1000; % [mm]
% Launcher_Dev_Distributed builds Zs from Zs_all up to max(SOIL_TH) and its
% point outputs (BC_O) assume 10 soil layers, which requires a maximum soil
% thickness between 500 and 999 mm (Zs = [0 10 20 50 100 150 200 250 300 350 500]).
% A thickness of 1000 mm would give 11 layers and the run would stop.
soil_th_max_cm = 99;

UTM_zone_text = '42 N'; % Kyzylsu, used for the multipoint coordinates

%% LAND COVER REMAPPING
%==========================================================================
% Former VEG_CODE (Launcher_Kyzylsu_Distributed / PARAMETERS_SOIL):
%   1 Fir, 2 Larch, 3 Grass C3, 4 Shrub, 5 BL evergreen, 6 BL deciduous,
%   7 Rock, 8 Ice (8 was treated as 7 = Rock)
% New classes (VPAR_T.xlsx) keep the former codes 1-7, and Ice (8) becomes
% Rock (7) as in the former launcher (ksv(ksv==8)=7).
old_codes = [1 2 3 4 5 6 7 8];
new_codes = [1 2 3 4 5 6 7 7];

%% LOAD FORMER DTM
%==========================================================================
S = load(dtm_in);
VEG_old = S.VEG_CODE;
VEG_CODE = zeros(size(VEG_old));
for i = 1:numel(old_codes)
    VEG_CODE(VEG_old == old_codes(i)) = new_codes(i);
end
missing = setdiff(unique(VEG_old(S.MASK==1)), old_codes);
if ~isempty(missing)
    error('VEG_CODE classes %s in the mask are not remapped', mat2str(missing))
end
VEG_CODE(S.MASK~=1) = 0;
S.VEG_CODE = VEG_CODE;

% Land cover fractions as grids (OPT_Veg_Param = 2), one value per class as
% in VPAR_T.xlsx / former PARAMETERS_SOIL.m (Ccrown, Cwat, Curb, Crock, Cbare)
%            Fir Larch Grass Shrub BL_ev BL_dec Rock
C_class = [  1   1     1     0.9   1     1      0;    % CCROWN
             0   0     0     0     0     0      0;    % CWATER
             0   0     0     0     0     0      0;    % CURB
             0   0     0     0     0     0      1;    % CROCK
             0   0     0     0.1   0     0      0];   % CBARE
C_names = {'CCROWN','CWATER','CURB','CROCK','CBARE'};
for k = 1:numel(C_names)
    G = zeros(size(VEG_CODE));
    for c = 1:size(C_class,2)
        G(VEG_CODE == c) = C_class(k,c);
    end
    S.(C_names{k}) = G;
end

% Soil thickness, % of 1000 mm -> cm
SOIL_TH = min(S.SOIL_TH,100)./100.*Zs_end_former./10;
SOIL_TH = min(SOIL_TH,soil_th_max_cm);
SOIL_TH(S.MASK~=1) = 0;
S.SOIL_TH = SOIL_TH;

% Initial snow
if use_ini_snow == 1
    tmp = load(fn_IniSnowDepth); S.SNOWD = tmp.SNOWD;
    tmp = load(fn_IniSnowAlbedo); S.SNOWALB = tmp.SNOWALB;
end

%% SAVE NEW DTM
%==========================================================================
out_dir = ['OUTPUTS/' PreProc_folder '/'];
if ~exist(out_dir,'dir'); mkdir(out_dir); end
save([out_dir dtm_out],'-struct','S');
disp(['Saved ' out_dir dtm_out])
for c = unique(VEG_CODE(S.MASK==1))'
    fprintf('  Class %d: %d cells\n',c,sum(VEG_CODE(S.MASK==1)==c));
end

%% MULTIPOINTS
%==========================================================================
% Same format as Shallap_MultiPoints.txt (COL = Xout, ROW = Yout)
[m_cell,~] = size(S.DTM);
ID = (1:numel(S.Xout))';
Name = strtrim(cellstr(S.POI_names));
COL = S.Xout(:); ROW = S.Yout(:);
LON_UTM = S.x(COL); LON_UTM = LON_UTM(:);
LAT_UTM = S.y(ROW); LAT_UTM = LAT_UTM(:);
[LAT,LON] = utm2deg(LON_UTM,LAT_UTM,repmat(UTM_zone_text,numel(ID),1));
idx = (COL-1).*m_cell + ROW;
DT = S.DEB_MAP(idx); DT = DT(:);
Points = table(ID,Name,COL,ROW,LON_UTM,LAT_UTM,LON,LAT,DT,idx);
writetable(Points,'../Multipoint/OUTPUTS/Kyzylsu_MultiPoints.txt','Delimiter','\t');
disp('Saved ../Multipoint/OUTPUTS/Kyzylsu_MultiPoints.txt')

%% METEOROLOGICAL FORCING (OPT_Forcing = 1)
%==========================================================================
St1 = load(fn_station_1); St1 = St1.forcing_all;
St2 = load(fn_station_2); St2 = St2.forcing_all;
if ~isequal(St1.Time,St2.Time); error('Station time steps differ'); end

% Former launcher: all variables from Pluvio, except air temperature which
% was lapsed using Pluvio and AWS_3900m. The shortwave was given already
% partitioned, here the total is passed as SWin and re-partitioned by the
% launcher (Automatic_Radiation_Partition_I).
DateTime = St1.Time;
Ta = St1.TA;
RH = St1.RH;
WS = St1.FF;
Pr = St1.PP;
SWin = St1.SAB1 + St1.SAB2 + St1.SAD1 + St1.SAD2;
LWin = St1.LWIN;
SWout = NaN(size(Ta)); % Not used by the launcher
LWout = NaN(size(Ta)); % Not used by the launcher

% Same NaN handling as the former launcher
fprintf('NaNs in forcing: Ta %d, RH %d, WS %d, Pr %d, SWin %d, LWin %d\n', ...
    sum(isnan(Ta)),sum(isnan(RH)),sum(isnan(WS)),sum(isnan(Pr)),sum(isnan(SWin)),sum(isnan(LWin)));
Pr(isnan(Pr)) = 0; Pr(Pr<0.001) = 0;
WS(isnan(WS)) = 0.01; WS(WS<0.01) = 0.01;
RH(RH<0) = 0; RH(RH>100) = 100;
SWin(isnan(SWin) | SWin<0) = 0;

Meteo_data = timetable(DateTime,Ta,RH,WS,SWout,LWin,LWout,Pr,SWin);
save('../Forcing/Kyzylsu_Pluvio_meteo.mat','Meteo_data');
disp('Saved ../Forcing/Kyzylsu_Pluvio_meteo.mat')

%% AIR TEMPERATURE LAPSE RATES PER HOUR AND MONTH
%==========================================================================
% Former launcher: hourly lapse rate from Pluvio and AWS_3900m.
% Here: mean of these hourly lapse rates per hour and month (positive
% values = decrease of temperature with elevation, in degC/m).
LR = (St2.TA - St1.TA)./(Elev(2) - Elev(1));
[hh,mm] = ndgrid(0:23,1:12);
hour = hh(:); month = mm(:);
Ta_lapse = NaN(size(hour)); GroupCount = NaN(size(hour));
for i = 1:numel(hour)
    id = St1.Time.Hour == hour(i) & St1.Time.Month == month(i) & ~isnan(LR);
    Ta_lapse(i) = -mean(LR(id));
    GroupCount(i) = sum(id);
end
if any(isnan(Ta_lapse)); error('Missing Ta lapse rate for some hour/month'); end
hm_Ta_lapse = table(hour,month,GroupCount,Ta_lapse);
save('../Forcing/Kyzylsu_hm_Ta_lapse.mat','hm_Ta_lapse');
disp('Saved ../Forcing/Kyzylsu_hm_Ta_lapse.mat')
fprintf('Ta lapse rate: min %.4f, mean %.4f, max %.4f degC/m\n',min(Ta_lapse),mean(Ta_lapse),max(Ta_lapse));

%% PRECIPITATION RATIO VERSUS ELEVATION
%==========================================================================
% Former launcher: Pmod = 0, i.e. the Pluvio precipitation everywhere.
% Ratio = p1*elevation + p2 for each month -> p1 = 0, p2 = 1.
for i = 1:12
    Pr_lapse(i).month = struct('p1',0,'p2',1);
end
save('../Forcing/Kyzylsu_Pr_ratio_v_el_const.mat','Pr_lapse');
disp('Saved ../Forcing/Kyzylsu_Pr_ratio_v_el_const.mat')

%% GLACIER ALBEDO VERSUS ELEVATION
%==========================================================================
% Former launcher: piecewise profile from Kyzylsu_Albedo_vs_elev.mat, capped
% at 0.4. Launcher_Dev_Distributed uses a linear relationship, fitted here
% (copy Alb_v_el_MM / Alb_v_el_CC into Options_and_NoVegParams.xlsx).
A = load('INPUTS/Kyzylsu_Albedo_vs_elev.mat');
pf = polyfit(A.Alb_el_tt.Elevation,min(A.Alb_el_tt.Glacie_albedo,0.4),1);
fprintf('Albedo vs elevation fit: Alb_v_el_MM = %.9g, Alb_v_el_CC = %.9g\n',pf(1),pf(2));
fprintf('Albedo profile min %.3f, max %.3f\n',min(A.Alb_el_tt.Glacie_albedo),max(A.Alb_el_tt.Glacie_albedo));
