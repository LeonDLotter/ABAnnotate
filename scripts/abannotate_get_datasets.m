function [datasets] = abannotate_get_datasets(dataset_type, disp_names)
%function [dataset_names, dataset_files] = abannotate_get_gcea_datasets(disp_names)
%
% Outputs a table of ABAnnotate dataset names and OSF ids.
% dataset_type = one of {expression_data, 'parcellation', 'gcea_dataset'}
% disp_names = bool

if ~exist('dataset_type', 'var')
    dataset_type = 'gcea_dataset';
end
if ~exist('disp_names', 'var')
    disp_names = true;
end

% paths
wd = fileparts(which('ABAnnotate'));
sources_path = fullfile(wd, 'dataset_sources.csv');

% get sources (download if necessary)
if exist(sources_path, 'file')
    sources = readtable(sources_path);
else
    sources = abannotate_get_sources();
end

% get sub dataset
sources = sources(strcmp(sources.type, dataset_type), :);
datasets = sources(:, {'name', 'osf'});

% print

if disp_names
    switch dataset_type
    case 'gcea_dataset'
        disp('Available GCEA datasets:')
    case 'expression_data'
        disp('Source of gene expression datasets:')
    case 'parcellation'
        disp('Available parcellations with processed gene expression data:')
    end
    for i=1:height(datasets)
        disp(['- ' datasets.name{i}]);   
    end
end

end


    