function import_aba_data(aba_csv, atlas, save_path, prefix)
% function import_aba_data(aba_csv, atlas, save_path, prefix)
% convert ABA csv exported from python into matlab expression matrix and 
% gene symbol vector

% read aba data
disp('Reading ABA data.')
aba = readtable(aba_csv, 'Format', 'auto');
try
    aba.label = [];
end

% get gene symbols
gene_symbols = aba.Properties.VariableNames;

% get expression matrix
expression_matrix = table2array(aba);

% save expression data
if ~exist('save_path', 'var')
    save_path = fullfile(fileparts(which('ABAnnotate')), 'atlas');
end
if ~exist('prefix', 'var')
    [~, prefix, ~] = fileparts(aba_csv);
end
save(fullfile(save_path, sprintf('%s_expression.mat', prefix)), 'expression_matrix', 'gene_symbols');

% save atlas
[~, ~, ext] = fileparts(atlas);
if strcmp(ext, '.gz')
    atlas_unzip = gunzip(atlas, save_path);
    movefile(char(atlas_unzip), fullfile(save_path, sprintf('%s_atlas.nii', prefix)));
elseif strcmp(ext, '.nii')
    copyfile(atlas, fullfile(save_path, [char(prefix) '_atlas.nii']))
else
    error('Unsupported atlas volume file format!');
end
fprintf('Imported data saved to: %s\n', save_path);

end