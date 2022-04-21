function [dataset_path] = get_gcea_dataset(dataset_name)
%function [dataset_path] = get_gcea_dataset(dataset_name)
%
% Gets path to ABAnnotate GCEA dataset or downloads it from OSF

wd = fileparts(which('Abannotate'));
dataset_name = char(dataset_name);

% get list of datasets
dataset_list = abannotate_get_datasets('gcea_dataset', false);

% check if exist and get osf id
try
    osf_id = dataset_list.osf{strcmp(dataset_list.name, dataset_name)};
catch
    error('GCEA dataset %s not found!', dataset_name);
end

% get path or download
dataset_path = fullfile(wd, 'datasets', [dataset_name '.mat']);
if ~exist(dataset_path, 'file')
    disp('Dataset not found. Downloading from OSF...')
    abannotate_download_osf(osf_id, dataset_path);
    fprintf('Saved to: %s\n', dataset_path);
end

end

