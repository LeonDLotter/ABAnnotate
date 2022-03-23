function cTablePheno = ensemble_enrichment(...
    dataset_mat, aba_mat, category_nulls, phenotype_data)
%function cTablePheno = ensemble_enrichment(...
%    dataset_mat, aba_mat, category_nulls, phenotype_data)
%
% Compute enrichment across GO categories for a given null model.
% Assumes that nulls have been precomputed using ComputeAllCategoryNulls.
% The parameters are taken from the null ensemble enrichment file (params)
%
% INPUTS:
% - 
%
% OUTPUT:
% - cTablePheno: a table with p-values estimated from the null ensemble

%--------------------------------------------------------------------------
% Process inputs and set defaults:
%--------------------------------------------------------------------------
if nargin < 1
    error('You must specify the dataset file!');
end
if nargin < 2
    error('You must specify a file containing the gene expression data!');
end
if nargin < 3
    error('You must specify a file containing the precomputed ensemble nulls!');
end
if nargin < 4
    error('You must provide a parcellated phenotype vector!');
end

%-------------------------------------------------------------------------------
% Load null distributions into GOTableNull
%-------------------------------------------------------------------------------
load(category_nulls, 'cTable', 'gcea_opt');
cTableNull = cTable;

%-------------------------------------------------------------------------------
% Now compute scores for the input phenotype using the same parameter settings
% as for the null distribution computation:
%-------------------------------------------------------------------------------
cTablePheno = generate_category_nulls(dataset_mat, ... 
                                      aba_mat, ...
                                      gcea_opt, ...
                                      [], ...
                                      phenotype_data, ...
                                      false, false);  

% Check that we have the same category Labels in both cases:
n_categories = height(cTablePheno);
if ~(height(cTableNull)==n_categories) && ~all(strcmp(cTableNull.cLabel, cTablePheno.cLabel))
    error('Error while matching categories to precomputed null data!');
end

%--------------------------------------------------------------------------
% Estimate p-values
%--------------------------------------------------------------------------
if isfield(gcea_opt, 'p_tail')
    tail = gcea_opt.p_tail;
else
    tail = 'right';
    fprintf(1,['Right-tailed test by default: larger values of %s %s ', ...
               'correlation are interesting.\n'], ...
            gcea_opt.aggregation_method, gcea_opt.correlation_method);
end

cTablePheno = estimate_pvals(cTableNull.cScores, ...
                             [cTablePheno.cScores{:}], ...
                             tail, cTablePheno);

% make nice output table & sort by gaussean z p-values and 
cTablePheno.cScoresNull = cTableNull.cScores;
cTablePheno = renamevars(cTablePheno, 'cScores', 'cScorePheno');
cTablePheno.cScorePheno = cell2mat(cTablePheno.cScorePheno);
cTablePheno = movevars(cTablePheno, 'cScoresNull', 'Before', 'cScorePheno');
cTablePheno = sortrows(cTablePheno, 'pValZ', 'ascend');

% Give a basic output about significance 
fprintf(1,'%u significant categories at p < %.3f derived from gaussian distributions.\n', ...
    sum(cTablePheno.pValZCorr < gcea_opt.p_thresh), gcea_opt.p_thresh);
fprintf(1,'%u significant categories at p < %.3f derived from null distributions. \n', ...
    sum(cTablePheno.pValPermCorr < gcea_opt.p_thresh), gcea_opt.p_thresh);

end
