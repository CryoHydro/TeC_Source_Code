%Create table with info on multipoints for Tom's data

% NOTE THIS IS A MESS AS WAS MADE WHEN THE TABLE INFO WAS FROM UNFLIPPED
% DEM

%NO LONGER NEEDED

%Requires Cluster_Points.mat

%To be in form
%ID,Name,COL,ROW,LON_UTM,LAT_UTM,LON,LAT,DT,idx
%Don't absolutelyneed debris thickness

Catch = 'Pasterze';

load('Pasterze_CLUSTERpoints.mat');
mm = load('B:\group\pelligrp\Project_folders\MountainWater\Tethys_Chloris\Workflow_dev\Scripts\Cat\TeC_Source_Code\Pasterze\Grids\dtm_Pasterze_250m.mat');

%So here the DEM and probablyLat and Lon are not flipped but need to be! 
%So row and col are wrong

NAME = clusterTable.Name;
Num_cluster = size(NAME,1);
ID = clusterTable.ClusterID;
LON = clusterTable.Lon;
LAT = clusterTable.Lat;
X = clusterTable.X_UTM;
Y = clusterTable.Y_UTM;
ELEV = clusterTable.Elev_m;
% COL = clusterTable.GridCol;
% ROW = clusterTable.GridRow;

% IDX = zeros(Num_cluster,1);
% for c=1:Num_cluster
%     GridT = zeros(279,445);
%     GridT(ROW(c),COL(c))=1;
%     IDX(c) = find(GridT==1);
% end

demLatsF = flipud(demLats);
demLonsF = flipud(demLons);
demZsF = flipud(demZs);

zone = utmzone(LAT,LON);
zone=zone(1:end-1);

%Easting 72.... Northing 52....
[LAT_UTM,LON_UTM,utmZ1] = ll2utm(LAT,LON,'wgs84',str2double(zone)); %Need to specify zone as the same for all points

%% Find row and column
ROW = zeros(Num_cluster,1);
COL = zeros(Num_cluster,1);
idx = zeros(Num_cluster,1);

for c=1:Num_cluster
    diffLat_a = abs(demLatsF - LAT(c));
    min_Lat = min(min(diffLat_a));
    loc_Lat = min_Lat == diffLat_a;
    
    diffLon_a = abs(demLonsF - LON(c));
    min_Lon = min(min(diffLon_a));
    loc_Lon = min_Lon == diffLon_a;
   
    idx(c) = find(loc_Lat & loc_Lon);
    [ROW(c), COL(c)] = find(loc_Lat & loc_Lon);
end

ELEV2 = zeros(Num_cluster,1);
for c=1:Num_cluster
ELEV2(c) = demZsF(ROW(c), COL(c));
end

%% 
cellsize = mm.cellsize;
xllcorner =  mm.xllcorner;
yllcorner = mm.yllcorner;
DTM = mm.DTM_orig;

for k = 1:Num_cluster
id_location = char(string(NAME(k))); %id

y_coord = Y(k); %Y - check as Shallap might be wrong!
x_coord = X(k); %X

pos_col = floor((x_coord - xllcorner) / cellsize) + 1;
pos_row = floor((y_coord - yllcorner) / cellsize) + 1;

%ij = POI.idx(loc);
ij = sub2ind(size(DTM),pos_row, pos_col); % Location
[j, i] = ind2sub(size(DTM), ij); % Location

POI.ij(k) = ij;
POI.i(k) = i; %Same as column
POI.j(k) = j; %Same as row

POI.Zbas(k) = DTM(j,i); % Altitude
end

%% In original code

%==========================================================================
% PREPARE MAP-COORDINATES
%==========================================================================
[~,x_v,y_v] = GRIDobj2mat(demZs);    % x y coordinates
y_v = flip(y_v);    %y has to go from smallest to largest coordinate (S-N)

%==========================================================================
% PREPARE MAP-CELL INDICES FOR POI
%==========================================================================
%Y-coordinate has to be flipped to be in T&C format
% correct cell coordinate can be checked here:
fldem = flipud(demZs);
imagesc(fldem); 
axis xy

%x & y indices to all POI cells ("tracked points")- so these are cell
%indices
X_v = table2array(poi_utm(:,'LON_UTM'));
Y_v = table2array(poi_utm(:,'LAT_UTM'));
[npoi,~] = size(poi_utm);
X_POI = zeros(npoi,1); Y_POI = zeros(npoi,1);
for i = 1:npoi
    [~,X_POI(i)] = min(abs(x_v - X_v(i)));
    [~,Y_POI(i)] = min(abs(y_v - Y_v(i)));
end
clear('X_v','Y_v');

%==========================================================================
% PREPARE MAP-CELL INDICES FOR CLUSTERS
%==========================================================================
%Y-coordinate has to be flipped to be in T&C format
% correct cell coordinate can be checked here:
% fldem = flipud(DEM.Z);
% imagesc(fldem);

figure; imagesc(flipud(DEM.Z))
hold on; scatter(X_POI(1), Y_POI(1),'+','r'); %So can try different POIs here and check make sense
set(gca,'YDir','normal')
% axis xy

%% 
Multipoints_T = table(ID,NAME,LAT,LON,X,Y,ELEV,COL,ROW,idx,POI.i',POI.j',POI.ij');
Multipoints_T.Properties.VariableNames = {'ID','NAME','LAT','LON','X','Y','ELEV','COL','ROW','IDX','i','j','ij'};

writetable(Multipoints_T,strcat(Catch,"_MultiPoints.txt"));




