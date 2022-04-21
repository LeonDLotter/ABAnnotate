function [] = abannotate_delete_datasets(which_datasets)
%function [] = abannotate_delete_datasets(which_datasets)
% 
% Function to clear all ABAnnotate datasets (atlas and gcea datasets)
% Input:
%   which_datasets = one of {'all', 'atlas', 'gcea'}, default to 'all'

if nargin < 1
    which_datasets = 'all';
end

wd = fileparts(which('ABAnnotate'));
atlas_dir = fullfile(wd, 'atlas');
gcea_dir = fullfile(wd, 'datasets');

atlas_dirs = {dir(fullfile(atlas_dir, '*')).name};
atlas_dirs = atlas_dirs(~startsWith(atlas_dirs, '.'));
gcea_files = {dir(fullfile(gcea_dir, '*.mat')).name};

fprintf('%u atlas and %u GCEA datasets found. Deleting %s:\n', ...
        length(atlas_dirs), length(gcea_files), which_datasets);

switch which_datasets
    case 'all'
        files_to_delete = [fullfile(atlas_dir, atlas_dirs)
                           fullfile(gcea_dir, gcea_files)]; 
    case 'atlas'
        files_to_delete = fullfile(atlas_dir, atlas_dirs); 
    case 'gcea'
        files_to_delete = fullfile(gcea_dir, gcea_files); 
end

for i=1:length(files_to_delete)
    disp(files_to_delete{i}); 
end

del = input(sprintf('Deleting %s null data, proceed? Y/N [N]', which_datasets), 's');
if isempty(del)
    del = 'N'; 
end

if ismember(del, {'Y', 'y', 'yes', 'Yes', 'YES'})
    for i=1:length(files_to_delete)
        file_to_delete = files_to_delete{i};
        fprintf('Deleting file: %s\n', file_to_delete);
        if endsWith(file_to_delete, '.mat')
            delete(file_to_delete); 
        else
            rmdir(file_to_delete, 's')        
        end
    end
else
    disp('No files deleted.');
end

end