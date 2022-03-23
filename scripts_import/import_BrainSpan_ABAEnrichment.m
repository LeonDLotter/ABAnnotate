function [cTable] = import_BrainSpan_ABAEnrichment(ABAEnrichment_file, save_path)
%function [cTable] = import_aba_BrainSpan_ABAEnrichment(ABAEnrichment_file, save_path)
% 
% Function to import AHBA BrainSpan data derived in processed form from
% the ABAEnrichment R-package 
% (https://doi.org/doi:10.1038/nature13185)
% (https://doi.org/10.1093/bioinformatics/btw392)
% (http://bioconductor.org/packages/release/data/experiment/html/ABAData.html)

% get dataset
brainspan = readtable(ABAEnrichment_file);

% get categories
tab = brainspan(:, {'age_category', 'structure'});
tab = unique(tab, 'rows');

% make table
for i=1:height(tab)
    % age labels
    switch tab.age_category(i)
        case 1; tab.age_label{i} = 'prenatal';
        case 2; tab.age_label{i} = 'infant';
        case 3; tab.age_label{i} = 'child';
        case 4; tab.age_label{i} = 'adolescent';
        case 5; tab.age_label{i} = 'adult';
    end
    % brain region labels
    switch tab.structure(i)
        case 10163; tab.struct_label{i} = 'primary motor cortex (area M1, area 4)';
                    tab.struct_abbr{i} = 'M1C';
        case 10173; tab.struct_label{i} = 'dorsolateral prefrontal cortex';
                    tab.struct_abbr{i} = 'DFC';
        case 10185; tab.struct_label{i} = 'ventrolateral prefrontal cortex';
                    tab.struct_abbr{i} = 'VFC';
        case 10194; tab.struct_label{i} = 'orbital frontal cortex';
                    tab.struct_abbr{i} = 'OFC';            
        case 10209; tab.struct_label{i} = 'primary somatosensory cortex (area S1, areas 3,1,2)';
                    tab.struct_abbr{i} = 'S1C';            
        case 10225; tab.struct_label{i} = 'posteroventral (inferior) parietal cortex';
                    tab.struct_abbr{i} = 'IPC';            
        case 10236; tab.struct_label{i} = 'primary auditory cortex (core)';
                    tab.struct_abbr{i} = 'A1C';            
        case 10243; tab.struct_label{i} = 'posterior (caudal) superior temporal cortex (area 22c)';
                    tab.struct_abbr{i} = 'STC';
        case 10252; tab.struct_label{i} = 'inferolateral temporal cortex (area TEv, area 20)';
                    tab.struct_abbr{i} = 'ITC';            
        case 10269; tab.struct_label{i} = 'primary visual cortex (striate cortex, area V1/17)';
                    tab.struct_abbr{i} = 'V1C';            
        case 10278; tab.struct_label{i} = 'anterior (rostral) cingulate (medial prefrontal) cortex';
                    tab.struct_abbr{i} = 'MFC';            
        case 10294; tab.struct_label{i} = 'hippocampus (hippocampal formation)';
                    tab.struct_abbr{i} = 'HIP';
        case 10333; tab.struct_label{i} = 'striatum';
                    tab.struct_abbr{i} = 'STR';            
        case 10361; tab.struct_label{i} = 'amygdaloid complex';
                    tab.struct_abbr{i} = 'AMY';            
        case 10398; tab.struct_label{i} = 'mediodorsal nucleus of thalamus';
                    tab.struct_abbr{i} = 'MD';            
        case 10657; tab.struct_label{i} = 'cerebellar cortex';
                    tab.struct_abbr{i} = 'CBC';
    end
    % new category label
    tab.label{i} = [tab.age_label{i} '_' tab.struct_abbr{i}];
    % new category desc
    tab.desc{i} = [tab.age_label{i} ' - ' tab.struct_label{i}];
    
    % genes
    tab.genes{i} = brainspan.hgnc_symbol(...
        brainspan.age_category==tab.age_category(i) & ...
        brainspan.structure==tab.structure(i));
    
    % gene weights
    tab.weights{i} = brainspan.signal(...
        brainspan.age_category==tab.age_category(i) & ...
        brainspan.structure==tab.structure(i));
    tab.size(i) = length(tab.weights{i});
end

% category name
cTable = table(tab.label, tab.desc, tab.size, tab.genes, tab.weights, ...
               'VariableNames', {'cLabel', 'cDesc', 'cSize', 'cGenes', 'cWeights'});
           
% save
save(save_path, 'cTable');

