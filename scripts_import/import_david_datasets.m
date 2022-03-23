function [cTable] = import_david_datasets(david_file, save_path)
%function [cTable] = import_david_datasets(david_file, save_path)
%
% Function to import gene-symbol-to-category mapping files as provided by
% DAVID (https://david-d.ncifcrf.gov/). This function works for all files
% with the structure: "gene-symbol \tab ID~description" or 
% "gene-symbol \tab ID"
% 
% Examples: 
%   ABCA4	601718~Retinitis pigmentosa 19
%   APOA2	GO:0006641~triglyceride metabolic process
% Source (requires registration with email, name, and institution):
%   https://david.ncifcrf.gov/knowledgebase/login.html

% read table
tab = readtable(david_file, 'ReadVariableNames', false, 'Format', '%s%s');
tab = tab(:,1:2);
tab.Properties.VariableNames = {'genes', 'categories'};

% get categories
categories = unique(tab.categories);

% get genes
genes = cellfun(@(x) tab.genes(strcmp(x, tab.categories)), ...
                categories, 'UniformOutput', false);
sizes = cellfun(@length, genes);
            
% ID and description are saved as "XXXXX~Description"
categories = split(categories, '~');

% make cTable
% case 1: two category labels separated by ~
if width(categories)==2
    cTable = table(categories(:,1), categories(:,2), sizes, genes, ...
                   'VariableNames', {'cLabel', 'cDesc', 'cSize', 'cGenes'});
% case 2: one category label
elseif width(categories)==1
    cTable = table(categories, sizes, genes, ...
                   'VariableNames', {'cLabel', 'cSize', 'cGenes'});
else
    error('Does the category column contain multiple occurences of "~"!?');
end
           
% save
save(save_path, 'cTable');

end
