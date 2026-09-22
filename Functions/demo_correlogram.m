%% Prepare the data
% Load the example dataset and compute correlation matrix
D = load('accidents.mat');
C = corr(D.hwydata);
axislabels = D.hwyheaders;

%% Draw the correlogram
% Prepare the figure and set the size
figure(1); clf();
set(gcf, 'Position', [0 0 1080 640]); movegui('center');
correlogram(C, 'AxisLabels', axislabels); 

%% Nodes as sorted by default, set 'Sorting' false to disable it
figure(2); clf();
set(gcf, 'Position', [0 0 1080 640]); movegui('center');
correlogram(C, 'AxisLabels', axislabels, 'Sorting', false);

