function [dataset_names, dataset_files] = get_gcea_datasets(disp_names)
%function [dataset_names, dataset_files] = get_gcea_datasets()
%
% Outputs two cell arrays of ABAnnotate GCEA datasets, first only dataset
% names, second full file paths.

wd = fileparts(which('ABAnnotate'));
dataset_dir = dir(fullfile(wd, 'datasets'));

if nargin < 1
    disp_names = true;
end

j = 1;
for i=1:length(dataset_dir)
    if dataset_dir(i).isdir == 0
        file = dataset_dir(i).name;
        folder = dataset_dir(i).folder;
        [~, dataset_names{j}, ~] = fileparts(file);
        dataset_files{j} = fullfile(folder, file) ;
        
        if disp_names
            disp(dataset_names{j});
        end
        
        j = j+1;
    end
end


    