%% Prepare the data
% Load the example dataset and compute correlation matrix
D = load('patients.mat');
D.IsFemale = ismember(D.Gender, 'Female');
[~, D.HealthStatus] = ismember(D.SelfAssessedHealthStatus, {'Poor', 'Fair', 'Good', 'Excellent'});

% Prepare the matrix
X = [D.IsFemale, D.Age, D.Weight, D.Height, D.Smoker, ...
     D.Systolic, D.Diastolic, D.HealthStatus];
axislabels = {'GenderFemale', 'Age', 'Weight', 'Height', ...
    'IsSmoker', 'Systolic', 'Diastolic', 'HealthStatus'};
C = corr(X);
%% Draw the correlogram
% Prepare the figure and set the size
figure(1); clf();
set(gcf, 'Position', [0 0 640 480]); movegui('center');
% Draw the figure
correlogram(C, 'AxisLabels', axislabels); 
set(gca, 'FontSize', 12)