function [cTable, save_path] = generate_category_nulls(...
    dataset_mat, aba_mat, gcea_opt, phenotype_nulls_mat, ...
    phenotype_data, save_result, verbose)
%function [cTable, save_path] = generate_category_nulls(...
%    dataset_mat, aba_mat, phenotype_nulls_mat, gcea_opt, ...
%    phenotype_parc, save_result, verbose)
%
% Compute an ensemble-based null distribution for all given gene categories.
% Can also compute for a single phenotype by specifying a non-empty 
% phenotype_parc vector.
% All parameters of the calculation are set in the gcea_opt structure.
%
% INPUTS:
% -
  
%%-------------------------------------------------------------------------
% Check inputs,  set defaults, get data
%--------------------------------------------------------------------------

if nargin < 3
    error(['You must provide the dataset file, processed gene-expression ', ...
           'data & GCEA parameters!']);
end
if ~exist('phenotype_data', 'var')
    phenotype_data = [];
end
if ~exist('save_result', 'var')
    save_result = true;
end
if ~exist('verbose', 'var')
    verbose = true;
end

%--------------------------------------------------------------------------
% get gene-expression data
load(aba_mat, 'gene_symbols', 'expression_matrix');
%n_rois = size(expression_matrix,1);

%--------------------------------------------------------------------------
% filter dataset (exclude categories not matching size filter and missing genes
cTable = filter_dataset(dataset_mat, gcea_opt, gene_symbols);
n_categories = height(cTable);

%--------------------------------------------------------------------------
corr_method = gcea_opt.correlation_method;
n_null = gcea_opt.n_category_nulls;

%--------------------------------------------------------------------------
% choose if calculation based on null maps

% run calculation only on phenotype vector (no null maps)
if ~isempty(phenotype_data) && isempty(phenotype_nulls_mat)
    nullMaps = phenotype_data;
    n_null = 1;
    fprintf(1,'Computing category scores for the spatial phenotype provided.\n');

% run calculation for provided null maps
elseif ~isempty(phenotype_nulls_mat) && isempty(phenotype_data)
    load(phenotype_nulls_mat, 'nullMaps');
    fprintf(1,'Computing category scores for %u generated %u-region null-phenotype-maps.\n',...
            size(nullMaps, 2), size(nullMaps, 1));
    fprintf(1,'(Spatial null maps loaded from %s)\n', phenotype_nulls_mat);
else
    error('You must provide either phenotype data or phenotype nulls!');
end

%--------------------------------------------------------------------------
% Prepare for category-wise agglomeration by first computing the correlation
% of each gene with a given spatial map (or a null ensemble of spatial maps)

% get all unique annoted genes
all_annotated_genes = unique(vertcat(cTable.cGenes{:}));
fprintf(1, '%u of %u annotated genes are unique.\n', ...
        length(all_annotated_genes), sum(cTable.cSize));
all_annotated_genes = intersect(all_annotated_genes, gene_symbols);
n_annotated_genes = length(all_annotated_genes);
fprintf(1,'Check again: We should have data for all %u of %u unique annotated genes.\n',...
        length(all_annotated_genes), n_annotated_genes);

% get gene scores
gene_scores = nan(n_annotated_genes, n_null);
fprintf(1,['Computing category null distributions corresponding to %u null ', ...
           'phenotypes for all %u genes annotated to categories.\n'], ...
           n_null, n_annotated_genes);
fprintf(['Progress (this might take some time...):\n', repmat('.', 1, 100) '\n']);
progress_step = round(n_annotated_genes/100,0);
for g = 1:n_annotated_genes
    % Get this gene's expression vector:
    match_me = (strcmp(gene_symbols, all_annotated_genes(g)));
    if sum(match_me)~=1
        fprintf(1,'How the heck did I get %u matches for gene %s?!\n',...
                sum(match_me), string(all_annotated_genes(g)));
    end
    expression_vector = expression_matrix(:, match_me);

    % The correlation to compute for this gene: (r-to-z transformed)
    theCorr_fun = @(x) fishers_r_to_z(corr(x, expression_vector, 'type', corr_method, 'rows', 'pairwise'));
    if n_null==1
        gene_scores(g) = theCorr_fun(nullMaps(:, 1));
    else
        parfor n = 1:n_null
            gene_scores(g, n) = theCorr_fun(nullMaps(:, n));
        end
    end
    
    % print progress
    if mod(g, progress_step)==0; fprintf('.'); end
end
fprintf('\n');

%--------------------------------------------------------------------------
% Agglomerate gene scores by gene category
%--------------------------------------------------------------------------
category_scores = cell(n_categories, 1);
for i = 1:n_categories
    if verbose
        fprintf(1,'Looking in at Category %u/%u: %s (%u genes).\n',...
            i, n_categories, cTable.cLabel{i}, cTable.cSize(i));
    end

    % Match genes in this category to all genes
    category_genes = cTable.cGenes{i};
    match_me_genes = find(ismember(all_annotated_genes, category_genes));
    category_genes_matched = all_annotated_genes(match_me_genes);
    % Get distribution of gene scores across phenotypes
    scores_here = gene_scores(match_me_genes, :);

    % Get gene weights corresponding to genes
    category_weights = cTable.cWeights{i};
    match_me_scores = find(ismember(category_genes_matched, category_genes));
    weights_here = category_weights(match_me_scores);

    if verbose
        fprintf(1,['%u/%u genes from this category have matching records ', ...
                   'in the expression data.\n'],...
                length(category_genes_matched), length(category_genes));
    end   

    %----------------------------------------------------------------------
    % Aggregate gene-wise scores into an overall category score (for each phenotype)
    category_scores{i} = aggregate_scores(scores_here, weights_here, ...
                                          gcea_opt.aggregation_method);
                                      
end

%--------------------------------------------------------------------------
% Assign to the table
cTable.cScores = category_scores;

%--------------------------------------------------------------------------
% Save results to .mat file
save_path = gcea_opt.category_nulls;
if save_result
    fprintf(1,'Saving nulls from %u iterations to "%s"\n', n_null, save_path);
    save(save_path, 'cTable', 'gcea_opt', '-v7');
end

end
