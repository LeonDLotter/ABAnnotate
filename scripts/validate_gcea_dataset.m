function [valid] = validate_gcea_dataset(dataset_path)
    
% is file?
if ~isfile(dataset_path)
    error(['GCEA dataset file not found! (%s). \nType <abannotate_get_datasets> ',...
           'for a list of available datasets.'], dataset_path);
end

% contains cTable?
try
    load(dataset_path, 'cTable');
catch
    error('GCEA dataset *.mat must contain a table named <cTable>!');
end

% has correct columns?    
if ~any(ismember(cTable.Properties.VariableNames, {'cLabel', 'cSize', 'cGenes'}))
    error('GCEA dataset must contain at least 3 columns named <cLabel>, <cSize>, <cGenes>!');
end

% label column is cellstring?
if ~iscellstr(cTable.cLabel)
    error('"cLabel" column of dataset table must be cellstring not %s!', class(cTable.cLabel));
end

% size column is numeric?
if ~isnumeric(cTable.cSize(1))
    error('"cSize" column of dataset table must be numeric not %s!', class(cTable.cSize));
end

% has >= 1 row?
if height(cTable) < 1
    error('GCEA dataset must contain at least 1 row (= category)!');
end

valid = true;

end