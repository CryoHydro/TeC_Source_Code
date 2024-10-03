function checkInputsValidity(variables,variable_names)

invalidVariables = {};

for i = 1:numel(variables)

    % Don't check variables which won't make T&C crash if NaNs or weird values
    if any(strcmp(variable_names{i},{'OM_H', 'pow_dis', 'a_dis'}))  % Variable which won't make T&C crash if NaNs
        continue

    % Check specific T&C variables
    elseif strcmp(variable_names{i},'Ws') && ~and(variables{i}> 0, variables{i}< 300)
        disp('Problem with wind speed')
    elseif strcmp(variable_names{i},'Ta') && ~and(variables{i}> 0, variables{i}< 300)
        disp('Problem with air temperature')

    % Check all other T&C variables
    elseif ~isstruct(variables{i})
        if any(isnan(variables{i}(:))) || any(isinf(variables{i}(:)))
         invalidVariables{end+1} = variable_names{i};
        end
    end 
end

if ~isempty(invalidVariables)
    disp('The following variables have NaN or inf values:');
    disp(invalidVariables);
end 
end