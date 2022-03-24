function [] = abannotate_delete_null_data(which_data)
%function [] = abannotate_delete_null_data(which_data)
% 
% Function to clear all ABAnnotate null data (phenotype and category
% nulls) stored in the '/nulls' directory.
% Input:
%   which_data = one of {'all', 'categories', 'phenotypes'}, default to 'all'

if nargin < 1
    which_data = 'all';
end

wd = fileparts(which('ABAnnotate'));
null_dir = fullfile(wd, 'nulls');

cat_nulls = {dir(fullfile(null_dir, 'categonulls_*.mat')).name};
pheno_nulls = {dir(fullfile(null_dir, 'phenonulls_*.mat')).name};

fprintf('%u category and %u phenotype null datasets found. Deleting %s:\n', ...
        length(cat_nulls), length(pheno_nulls), which_data);

switch which_data
    case 'all'
        files_to_delete = [cat_nulls, pheno_nulls]; 
    case 'categories'
        files_to_delete = cat_nulls; 
    case 'phenotypes'
        files_to_delete = pheno_nulls; 
end

for i=1:length(files_to_delete)
    disp(files_to_delete{i}); 
end

del = input(sprintf('Deleting %s null data, proceed? Y/N [N]', which_data), 's');
if isempty(del)
    del = 'N'; 
end

if ismember(del, {'Y', 'y', 'yes', 'Yes', 'YES'})
    for i=1:length(files_to_delete)
        fullfile_to_delete = fullfile(null_dir, files_to_delete{i});
        disp(['Deleting file: ' fullfile_to_delete]);
        delete(fullfile_to_delete); 
    end
else
    disp('No files deleted.');
end

end