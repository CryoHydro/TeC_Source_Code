%% Initialize
clear 
close all
addpath('C:\Users\kneibm\Documents\UsefulCodes')
addpath('C:\Users\kneibm\Documents\UsefulCodes\ImGRAFT')
addpath('C:\Users\kneibm\Documents\UsefulCodes\topotoolbox')

homedir = 'C:\Users\kneibm\Documents\Projects\PI\2024_CAIRN-GLOBAL\SMB_inversions\Rhone\code\';
cd(homedir)

%% 

glacier_shp = shaperead('../data/gis/inventory_sgi2016_r2020/SGI_2016_glaciers.shp');

[dem,xdem,ydem] = geoimread_upd('..\output\dh_results\meanDEM-2020_02_29.tif');
[vx,xv,yv] = geoimread_upd('..\output\velocity\vx_2015-2021_clip_raw_epsg2056.tif');
[vy,~,~] = geoimread_upd('..\output\velocity\vy_2015-2021_clip_raw_epsg2056.tif');
[thx,xt,yt,~] = geoimread_upd('..\output\thx\thk_distr_grab2021_shifted.tif');
[~,R] = readgeoraster('..\output\thx\thk_distr_grab2021_shifted.tif');
dgeot = geotiffinfo('..\output\thx\thk_distr_grab2021_shifted.tif');

thx_obs = readgeotable('..\output\thx\GPR_points_Rhone_shifted.shp');
%dgeot = geotiffinfo([homedir,'output\dh_results\meanDEM-2017_02_15.tif']);

thx_obs_thx = thx_obs.thk;
thx_obs_x = thx_obs.Shape.X;
thx_obs_y = thx_obs.Shape.Y;

thx_obs_x(thx_obs_thx<0) = [];
thx_obs_y(thx_obs_thx<0) = [];
% thx_obs_thx(thx_obs_thx<0) = [];

thx(thx<0) = 0;
 
% resample to resolution of Millan thx (50 m)
[Xdem,Ydem] = meshgrid(xdem,ydem);
[Xv,Yv] = meshgrid(xv,yv);
[Xt,Yt] = meshgrid(xt,yt);

dem_r = interp2(Xdem,Ydem,dem,Xt,Yt,'bilinear');
vx_r = interp2(Xv,Yv,vx,Xt,Yt,'bilinear');
vy_r = interp2(Xv,Yv,vy,Xt,Yt,'bilinear');

[Glacier_r,~] = geotiffcrop_shp(glacier_shp,thx,R);
Glacier_r(isnan(Glacier_r)) = 0;

% Define the size and resolution of the DEM grid
demSizeX = length(xt); % Number of cells in the x direction
demSizeY = length(yt); % Number of cells in the y direction

% Create an empty raster grid filled with NaN values
thx_obs_grid = NaN(demSizeY, demSizeX);

resol = xt(2)-xt(1);

% Loop through each cell in the raster
for y = 1:demSizeY-1
    for x = 1:demSizeX-1
        % Define the boundaries of the current cell
        cellXMin = xt(x);
        cellXMax = xt(x + 1);
        cellYMin = yt(y+1);
        cellYMax = yt(y);
        
        pointMeasurements = [thx_obs_x thx_obs_y thx_obs_thx];
        
        % Find point measurements that intersect with the current cell
        intersectingPoints = pointMeasurements(...
            pointMeasurements(:, 1) >= cellXMin & pointMeasurements(:, 1) < cellXMax & ...
            pointMeasurements(:, 2) >= cellYMin & pointMeasurements(:, 2) < cellYMax, :);
        
        % Calculate the average value of intersecting point measurements
        if ~isempty(intersectingPoints)
            averageValue = mean(intersectingPoints(:, 3)); 
            thx_obs_grid(y, x) = averageValue;
        end
    end
end


% prepare data
YY = Yt';
XX = Xt';
ZZ = dem_r;
ZZ(isnan(ZZ)) = 0;
THKINIT = thx;
THKINIT(~Glacier_r) = 0;
ICEMASKOBS = Glacier_r;
UVELSURFOBS = vx_r;
UVELSURFOBS(ICEMASKOBS==0) = NaN;
VVELSURFOBS = vy_r;
VVELSURFOBS(ICEMASKOBS==0) = NaN;
THKOBS = thx_obs_grid;

geotiffwrite('igm/run_2016-2023_2024_10_01/usurf.tif', ZZ, R, 'GeoKeyDirectoryTag',dgeot.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('igm/run_2016-2023_2024_10_01/usurfobs.tif', ZZ, R, 'GeoKeyDirectoryTag',dgeot.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('igm/run_2016-2023_2024_10_01/thkinit.tif', THKINIT, R, 'GeoKeyDirectoryTag',dgeot.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('igm/run_2016-2023_2024_10_01/thkobs.tif', THKOBS, R, 'GeoKeyDirectoryTag',dgeot.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('igm/run_2016-2023_2024_10_01/icemaskobs.tif', ICEMASKOBS, R, 'GeoKeyDirectoryTag',dgeot.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('igm/run_2016-2023_2024_10_01/icemask.tif', ICEMASKOBS, R, 'GeoKeyDirectoryTag',dgeot.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('igm/run_2016-2023_2024_10_01/uvelsurfobs.tif', UVELSURFOBS, R, 'GeoKeyDirectoryTag',dgeot.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('igm/run_2016-2023_2024_10_01/vvelsurfobs.tif', VVELSURFOBS, R, 'GeoKeyDirectoryTag',dgeot.GeoTIFFTags.GeoKeyDirectoryTag);




