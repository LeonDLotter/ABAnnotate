function import_aba_data(aba_csv, atlas, save_path, prefix)
% function import_aba_data(aba_csv, save_path) 
% convert ABA csv exported from python into matlab expression matrix and 
% gene symbol vector

% read aba data
disp('Reading ABA data.')
aba = readtable(aba_csv);
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
save(fullfile(save_path, [prefix '_expression.mat']), 'expression_matrix', 'gene_symbols');

% save atlas
[~, ~, ext] = fileparts(atlas);
if strcmp(ext, '.gz')
    gunzip(atlas, fullfile(save_path, [prefix '_atlas.nii']))
elseif strcmp(ext, '.nii')
    copyfile(atlas, fullfile(save_path, [prefix '_atlas.nii']))
else
    error('Unsupported atlas volume file format!');
end
disp(['Imported data saved to: ', save_path]);

end