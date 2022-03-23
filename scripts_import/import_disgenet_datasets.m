function [cTable] = import_disgenet_datasets(disgenet_file, save_path, ...
    select_diseaseType, select_diseaseSemanticType)
%function [cTable] = import_disgenet_datasets(disgenet_file, select_type, save_path)
%
% Function to import disGeNET datasets retrieved from
% https://www.disgenet.org/downloads
%
% select_diseaseType can be empty or >=1 of:
%   {'disease', 'group', 'phenotype'}
% select_diseaseSemanticType can be empty or >=1 of:
%   all dataset:
%   {'Body Substance';'Cell Component';'Cell Function';'Cell or Molecular Dysfunction';'Classification';'Clinical Attribute';'Congenital Abnormality';'Diagnostic Procedure';'Disease or Syndrome';'Disease or Syndrome; Anatomical Abnormality';'Disease or Syndrome; Congenital Abnormality';'Experimental Model of Disease';'Finding';'Individual Behavior';'Injury or Poisoning';'Laboratory Procedure';'Laboratory or Test Result';'Mental Process';'Mental or Behavioral Dysfunction';'Molecular Function';'Neoplastic Process';'Neoplastic Process; Experimental Model of Disease';'Organ or Tissue Function';'Organism Attribute';'Organism Function';'Pathologic Function';'Physiologic Function';'Sign or Symptom';'Social Behavior';'Temporal Concept'}
%   curated dataset:
%   {'Anatomical Abnormality';'Body Substance';'Cell Component';'Cell Function';'Cell or Molecular Dysfunction';'Clinical Attribute';'Congenital Abnormality';'Diagnostic Procedure';'Disease or Syndrome';'Disease or Syndrome; Congenital Abnormality';'Experimental Model of Disease';'Finding';'Individual Behavior';'Injury or Poisoning';'Laboratory Procedure';'Laboratory or Test Result';'Mental Process';'Mental or Behavioral Dysfunction';'Molecular Function';'Neoplastic Process';'Neoplastic Process; Experimental Model of Disease';'Organ or Tissue Function';'Organism Attribute';'Organism Function';'Pathologic Function';'Sign or Symptom'}
%   befree dataset:
%   {'Acquired Abnormality';'Anatomical Abnormality';'Congenital Abnormality';'Disease or Syndrome';'Disease or Syndrome; Anatomical Abnormality';'Disease or Syndrome; Congenital Abnormality';'Mental or Behavioral Dysfunction';'Neoplastic Process';'Neoplastic Process; Experimental Model of Disease';'Sign or Symptom'}

if ~exist('select_diseaseType', 'var')
    select_diseaseType = [];
end
if ~exist('select_diseaseSemanticType', 'var')
    select_diseaseSemanticType = [];
end

% load data
tab = readtable(disgenet_file, 'FileType', 'delimitedtext');

% filter categories
if ~isempty(select_diseaseType)
    tab = tab(ismember(tab.diseaseType, select_diseaseType), :);
end
if ~isempty(select_diseaseSemanticType)
    tab = tab(ismember(tab.diseaseSemanticType, select_diseaseSemanticType), :);
end

% get categories
[cat_id, ia, ~] = unique(tab.diseaseId);
cat_lab1 = tab.diseaseName(ia);
cat_lab2 = tab.diseaseType(ia);
cat_lab3 = tab.diseaseSemanticType(ia);

% get genes
genes = cellfun(@(x) tab.geneSymbol(strcmp(x, tab.diseaseId)), ...
                cat_id, 'UniformOutput', false);
sizes = cellfun(@length, genes);

% combine categories with exactly the same genes 
%...

% make cTable
cTable = table(cat_id, cat_lab1, cat_lab2, cat_lab3, sizes, genes, ...
               'VariableNames', {'cLabel', 'cDesc1', 'cDesc2', 'cDesc3', 'cSize', 'cGenes'});

% save
save(save_path, 'cTable');

end