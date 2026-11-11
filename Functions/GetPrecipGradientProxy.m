%% estimate_precip_gradient_from_vertical_velocity.m
%
% Estimate a precipitation-gradient proxy from a single-column ERA5
% pressure-level time series, using ERA5 pressure vertical velocity.
%
% Input NetCDF structure:
%   z, q, t, w : [pressureLevel x valid_time]
%   pressureLevel : [pressureLevel x 1], hPa
%   valid_time    : [valid_time x 1], Unix seconds
%
% In ERA5:
%   z = geopotential                      [m^2 s^-2]
%   q = specific humidity                 [kg kg^-1]
%   t = temperature                       [K]
%   w = pressure vertical velocity omega  [Pa s^-1]
%
% Convention:
%   w < 0 = ascent
%   w > 0 = descent
%
% Output:
%   meanGradient_mm_hr_m
%      Condensation / precipitation-generation proxy [mm hr^-1 m^-1].
%
% Draft research code: validate against precipitation data before using
% the estimate as an actual precipitation-elevation gradient.
%

function [Gradient] = GetPrecipGradientProxy(ncfile,elevationRange_m,startDate,endDate,plot)

%% -------------------------- USER SETTINGS ------------------------------

% ncfile = "ERA5_PL_timeseries.nc";

% Geopotential-height interval to analyse [m above mean sea level].
% elevationRange_m = [500, 3500];

% Optional UTC analysis period.
% Leave [] to retain every 6-hourly record in the file.
% startDate = [];
% endDate   = [];


% Minimum upward vertical velocity magnitude [Pa s^-1].
%
% ERA5 ascent is w < 0. This threshold retains w <= -minAscent_Pa_s.
% Use 0 to retain every ascending sample.
minAscent_Pa_s = 0.01;

% Physical constants
g = 9.80665;          % m s^-2
Rd = 287.05;          % J kg^-1 K^-1
Rv = 461.5;           % J kg^-1 K^-1
epsilon = Rd / Rv;    % dimensionless
rhoWater = 1000;      % kg m^-3
secondsPerHour = 3600;

%% ---------------------------- READ DATA --------------------------------

p_hPa = double(ncread(ncfile, "pressureLevel"));
tUnix = double(ncread(ncfile, "valid_time"));

z = double(ncread(ncfile, "z"));  % m^2 s^-2
q = double(ncread(ncfile, "q"));  % kg kg^-1
T = double(ncread(ncfile, "t"));  % K
omega = double(ncread(ncfile, "w")); % Pa s^-1; negative = ascent

latitude = double(ncread(ncfile, "latitude"));
longitude = double(ncread(ncfile, "longitude"));

% Verify expected dimensions: [pressureLevel x valid_time].
assert(size(q,1) == numel(p_hPa), ...
    "Expected q to have pressureLevel in its first dimension.");

assert(size(q,2) == numel(tUnix), ...
    "Expected q to have valid_time in its second dimension.");

assert(isequal(size(q), size(z), size(T), size(omega)), ...
    "q, z, t, and w must have identical dimensions.");

p_Pa = p_hPa(:) .* 100;
time = datetime(tUnix, "ConvertFrom", "posixtime", "TimeZone", "UTC");

% ERA5 geopotential -> geopotential height [m].
Z_m = z ./ g;

%% --------------------------- TIME FILTER -------------------------------

keepTime = true(size(time));

if ~isempty(startDate)
    keepTime = keepTime & time >= startDate;
end

if ~isempty(endDate)
    keepTime = keepTime & time <= endDate;
end

if ~any(keepTime)
    error("No records remain after applying the time filter.");
end

time  = time(keepTime);
q     = q(:,keepTime);
T     = T(:,keepTime);
Z_m   = Z_m(:,keepTime);
omega = omega(:,keepTime);

nLevel = numel(p_Pa);
nTime = numel(time);

%% ---------------------- ELEVATION-RANGE FILTER -------------------------

inElevationRange = Z_m >= elevationRange_m(1) & ...
                   Z_m <= elevationRange_m(2);

if ~any(inElevationRange, "all")
    error("No pressure-level samples fall within elevationRange_m.");
end

%% --------------------- THERMODYNAMIC QUANTITIES ------------------------
%
% Estimate saturation specific humidity q_s(T,p), then calculate its
% vertical gradient dq_s/dz.
%
% For upward motion, cooling reduces saturation humidity. The approximate
% condensate generation rate is:
%
% C = rho_air * w_z * max(-dq_s/dz, 0) * RH_factor
%
% where:
%   C     = kg m^-3 s^-1
%   w_z   = upward geometric vertical velocity, m s^-1
%   RH_factor = q/q_s, clipped to [0,1].
%
% Dividing C by water density produces a water-equivalent tendency per
% metre: [m^-1 s^-1], converted below to [mm hr^-1 m^-1].

P_Pa = repmat(p_Pa, 1, nTime);

% Saturation vapour pressure over liquid water [Pa].
Tc = T - 273.15;
es_Pa = 611.2 .* exp((17.67 .* Tc) ./ (Tc + 243.5));

% Avoid invalid denominator if pressure is close to saturation pressure.
es_Pa = min(es_Pa, 0.99 .* P_Pa);

% Saturation mixing ratio -> saturation specific humidity.
rs = epsilon .* es_Pa ./ (P_Pa - es_Pa);
qs = rs ./ (1 + rs);

% Virtual temperature and air density.
Tv = T .* (1 + 0.61 .* q);
rhoAir = P_Pa ./ (Rd .* Tv);

% Relative saturation, bounded physically for this diagnostic.
relativeSaturation = min(max(q ./ qs, 0), 1);

%% ------------- CONVERT PRESSURE VELOCITY TO HEIGHT VELOCITY ------------
%
% Hydrostatic approximation:
%
%   dp/dz = -rho_air * g
%
% omega = dp/dt, so:
%
%   wz = dz/dt = -omega / (rho_air * g)
%
% Therefore omega < 0 gives wz > 0 (upward motion).

wz_m_s = -omega ./ (rhoAir .* g);

% Keep upward motion exceeding user-selected threshold.
ascentMask = omega <= -minAscent_Pa_s;

%% ------------------- VERTICAL SATURATION GRADIENT ----------------------

dqsdz = nan(nLevel, nTime);

for it = 1:nTime

    good = isfinite(Z_m(:,it)) & isfinite(qs(:,it));

    if nnz(good) < 3
        continue
    end

    zCol = Z_m(good,it);
    qsCol = qs(good,it);

    % gradient() requires monotonic coordinates.
    [zSorted, order] = sort(zCol, "ascend");
    qsSorted = qsCol(order);

    dqsdzSorted = gradient(qsSorted, zSorted);

    dqsdzCol = nan(size(zCol));
    dqsdzCol(order) = dqsdzSorted;

    dqsdz(good,it) = dqsdzCol;
end

% Positive value means saturation specific humidity decreases with height.
saturationHumidityDecrease = max(-dqsdz, 0);

%% ---------------- CONDENSATION / PRECIPITATION PROXY -------------------
%
% Condensation per unit volume:
%
% C = rho_air * vertical_velocity * (-dq_s/dz) * RH_fraction
%
% Units:
% rho_air [kg m^-3]
% wz      [m s^-1]
% dq_s/dz [m^-1]
% C       [kg m^-3 s^-1]
%
% Water-equivalent rate per vertical metre:
%
% C / rho_water * 1000 * 3600
% = mm hr^-1 m^-1

condensation_kg_m3_s = rhoAir .* wz_m_s .* ...
                       saturationHumidityDecrease .* ...
                       relativeSaturation;

gradientEachSample_mm_hr_m = ...
    (condensation_kg_m3_s ./ rhoWater) .* ...
    1000 .* secondsPerHour;

% The physical construction should already give non-negative values during
% ascent, but explicitly prevent numerical negative values.
gradientEachSample_mm_hr_m = max(gradientEachSample_mm_hr_m, 0);

%% ----------------------- COMBINE QUALITY FILTERS -----------------------

validMask = inElevationRange & ascentMask & ...
            isfinite(gradientEachSample_mm_hr_m) & ...
            isfinite(q) & isfinite(T) & isfinite(omega);

selectedGradient = gradientEachSample_mm_hr_m(validMask);

if isempty(selectedGradient)
    error(["No valid ascending samples remain. Try a lower ascent ", ...
           "threshold or inspect the specified elevation range."]);
end

%% ----------------------------- RESULTS ---------------------------------

meanGradient_mm_hr_m = mean(selectedGradient, "omitnan");
medianGradient_mm_hr_m = median(selectedGradient, "omitnan");
stdGradient_mm_hr_m = std(selectedGradient, "omitnan");

% Sampling fractions are useful diagnostics.
nPotential = nnz(inElevationRange);
nAscending = nnz(inElevationRange & ascentMask);
nValid = nnz(validMask);

fprintf("\n========================================================\n");
fprintf("ERA5 single-column ascent-based condensation proxy\n");
fprintf("Location: %.4f deg N, %.4f deg E\n", latitude, longitude);
fprintf("Period: %s to %s\n", string(time(1)), string(time(end)));
fprintf("Elevation interval: %.0f to %.0f m\n", ...
    elevationRange_m(1), elevationRange_m(2));
fprintf("Ascent criterion: omega <= -%.4g Pa s^-1\n", minAscent_Pa_s);
fprintf("--------------------------------------------------------\n");
fprintf("Elevation-range samples: %d\n", nPotential);
fprintf("Ascending samples:       %d\n", nAscending);
fprintf("Valid retained samples:  %d\n", nValid);
fprintf("--------------------------------------------------------\n");
fprintf("Mean gradient:   %.6g mm hr^-1 m^-1\n", ...
    meanGradient_mm_hr_m);
fprintf("Median gradient: %.6g mm hr^-1 m^-1\n", ...
    medianGradient_mm_hr_m);
fprintf("Std. deviation:  %.6g mm hr^-1 m^-1\n", ...
    stdGradient_mm_hr_m);
fprintf("Per 100 m:       %.6g mm hr^-1 per 100 m\n", ...
    meanGradient_mm_hr_m * 100);
fprintf("Per kilometre:   %.6g mm hr^-1 per km\n", ...
    meanGradient_mm_hr_m * 1000);
fprintf("========================================================\n");

%% ------------------ MEAN PROFILE WITHIN ANALYSIS RANGE -----------------

meanGradientByLevel = nan(nLevel, 1);
meanHeightByLevel_m = nan(nLevel, 1);
ascentFractionByLevel = nan(nLevel, 1);

for il = 1:nLevel

    levelRangeMask = inElevationRange(il,:);
    levelValidMask = validMask(il,:);

    meanHeightByLevel_m(il) = mean(Z_m(il,levelRangeMask), "omitnan");

    meanGradientByLevel(il) = mean( ...
        gradientEachSample_mm_hr_m(il,levelValidMask), "omitnan");

    if nnz(levelRangeMask) > 0
        ascentFractionByLevel(il) = ...
            nnz(ascentMask(il,:) & levelRangeMask) / nnz(levelRangeMask);
    end
end

%% ------------------------------- PLOTS ---------------------------------

if plot == 1
figure("Color", "w", ...
       "Name", "ERA5 ascent-based precipitation-gradient proxy");

subplot(1,3,1)
plot(meanGradientByLevel, meanHeightByLevel_m, "-o", ...
    "LineWidth", 1.5, "MarkerFaceColor", [0.1 0.4 0.8]);
xlabel("Gradient proxy (mm hr^{-1} m^{-1})");
ylabel("Mean geopotential height (m)");
yline(elevationRange_m(1), "--k");
yline(elevationRange_m(2), "--k");
grid on;
title("Mean gradient by level");

subplot(1,3,2)
plot(100 .* ascentFractionByLevel, meanHeightByLevel_m, "-o", ...
    "LineWidth", 1.5, "MarkerFaceColor", [0.8 0.3 0.1]);
xlabel("Ascending fraction (%)");
ylabel("Mean geopotential height (m)");
yline(elevationRange_m(1), "--k");
yline(elevationRange_m(2), "--k");
grid on;
title("Frequency of ERA5 ascent");

subplot(1,3,3)
histogram(selectedGradient, 50, ...
    "FaceColor", [0.1 0.4 0.8], "EdgeColor", "none");
xlabel("Gradient proxy (mm hr^{-1} m^{-1})");
ylabel("Number of samples");
grid on;
title("Distribution during ascent");

end

%% ------------------------------- SAVE ----------------------------------

results = struct();

results.ncfile = ncfile;
results.latitude_degN = latitude;
results.longitude_degE = longitude;
results.time = time;

results.elevationRange_m = elevationRange_m;
results.minAscent_Pa_s = minAscent_Pa_s;

results.meanGradient_mm_hr_m = meanGradient_mm_hr_m;
results.medianGradient_mm_hr_m = medianGradient_mm_hr_m;
results.stdGradient_mm_hr_m = stdGradient_mm_hr_m;

results.pressureLevel_hPa = p_hPa;
results.meanGeopotentialHeight_m = meanHeightByLevel_m;
results.meanGradientByLevel_mm_hr_m = meanGradientByLevel;
results.ascentFractionByLevel = ascentFractionByLevel;

results.validMask = validMask;

% save("ERA5_ascent_precipitation_gradient_proxy.mat", "results");

Gradient = results.meanGradient_mm_hr_m;