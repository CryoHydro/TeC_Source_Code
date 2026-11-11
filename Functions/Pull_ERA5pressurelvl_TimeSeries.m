function [TAB] = Pull_ERA5pressurelvl_TimeSeries(templateFile,OutName,latitude,longitude,startDate,endDate,LoadExisting)

% OutName   = 'ERA5pressurelvl_TimeSeries_NAME';
% latitude  = 47.1;
% longitude = 12.7;
% 
% startDate = datetime(1950,1,2);
% endDate   = datetime(2026,8,22);


% Python template and output filenames
% templateFile = 'Request_ERA5Land_Timeseries_template.py';
pythonFile   = 'ERA5pressurelvls_download_run.py';

NewFolder = fullfile('./Reanalysis',strrep(OutName,'nc',''));
if ~exist(NewFolder,'dir')
    mkdir(NewFolder)


end

cd(NewFolder) 
FILES = dir('*.nc');
cd ..
cd ..

if isempty(FILES) | LoadExisting == 0

%% PYTHON EXECUTABLE

% Get the Python executable configured in MATLAB
pyExe = pyenv;

% if pyExe.Status == "NotLoaded"
%     error(['No Python environment is configured in MATLAB. ' ...
%            'Set one using pyenv(Version="path_to_python.exe").']);
% end

pythonExe = pyExe.Executable;

fprintf('Using Python:\n%s\n\n', pythonExe);

%% Read Python template

txt = fileread(templateFile);


%% Replace variables

txt = regexprep(txt, ...
    "OutName\s*=\s*'.*?'", ...
    sprintf("OutName = '%s'", OutName));

txt = regexprep(txt, ...
    '"longitude":\s*[-+]?\d*\.?\d+', ...
    sprintf('"longitude": %.8f', longitude));

txt = regexprep(txt, ...
    '"latitude":\s*[-+]?\d*\.?\d+', ...
    sprintf('"latitude": %.8f', latitude));

txt = regexprep(txt, ...
    '"date":\s*\[".*?"\]', ...
    sprintf('"date": ["%s/%s"]', ...
    datestr(startDate,'yyyy-mm-dd'), ...
    datestr(endDate,'yyyy-mm-dd')));


%% Save modified Python file

fid = fopen(pythonFile,'w');
fwrite(fid,txt);
fclose(fid);

pause(2)

%% RUN PYTHON

command = sprintf('"%s" "%s"', pythonExe, pythonFile);

fprintf('Sending Request to ECMWF...\n\n');

[status,result] = system(command);

if status ~= 0
    error('Python script failed:\n%s',result);
end


fprintf('\n Data Downloaded....\n');

%% UNZIP AND PROCESS

FILES = dir('*.nc');

else
    disp('Data Already Exist... reprocessing..')
end


for ii = 1:length(FILES)
   DAT = load_era5_data_region_v2([FILES(ii).folder '/' FILES(ii).name]);
   if ii == 1
    TAB.Date = DAT.time(:);
   end
   F = fieldnames(DAT.variables);
    for vv = 1:length(F)
        TAB.(F{vv}) = ncread([FILES(ii).folder '/' FILES(ii).name],F{vv});  
    end
end

%% MOVE FILE INTO FOLDER
try
    for ii = 1:length(FILES)
        movefile(FILES(ii).name,fullfile(NewFolder,FILES(ii).name));
    end
catch
    
end
try
        movefile(pythonFile,fullfile(NewFolder,pythonFile));
catch
    
end



