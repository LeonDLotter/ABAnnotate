function cTable = filter_dataset(dataset_mat, gcea_opt, gene_symbols)
%function cTable = filter_dataset(dataset_mat, gcea_opt.category_sizes, gene_symbols)
% 
% Returns a table of gene annotations after filtering according to given
% category sizes, threshold, and removing genes not in provided expression 
% dataset.

%--------------------------------------------------------------------------
% Check Inputs:
if nargin < 3
    error('You must provide dataset file, a vector of category sizes, & a vector of gene symbols!')
end

cat_select = gcea_opt.categories_select;
%weights_rescale = gcea_opt.weights_rescale;
weights_thresh = gcea_opt.weights_thresh;
weights_quant = gcea_opt.weights_quant;
weights_cutoff = gcea_opt.weights_cutoff;
gene_coocc_thresh = gcea_opt.gene_coocc_thresh;
size_filter = gcea_opt.size_filter;

%--------------------------------------------------------------------------
% Load dataset
load(dataset_mat, 'cTable');
fprintf(1,'Loaded %u annotated categories from %s.\n', height(cTable), dataset_mat);

%--------------------------------------------------------------------------
% Check that all cGenes in GOTable are column vectors:
is_row_vector = find(cellfun(@(x)size(x, 2), cTable.cGenes)~=1);
for i = 1:length(is_row_vector)
    cTable.cGenes{is_row_vector(i)} = cTable.cGenes{is_row_vector(i)}';
end
fprintf(1,'Transposed %u category annotation vectors from row to column vector.\n', ...
        length(is_row_vector));

%--------------------------------------------------------------------------
% Remove unselected categories 
%--------------------------------------------------------------------------
if isempty(cat_select)
    disp('Using all categories in dataset.')
else
    fprintf(1, 'Keeping %u/%u selected categories in dataset.\n', ...
            length(cat_select), height(cTable))
    % search selected cagegories in dataset
    match_me = ismember(cTable.cLabel, cat_select);
    % check if at least one category selected
    if sum(match_me) > 0
        % filter
        cTable = cTable(match_me, :);
    else
        disp('No selected categories match those in  dataset. Defaulting to all categories!')        
    end
end
    
%--------------------------------------------------------------------------
% Match genes 
%--------------------------------------------------------------------------
fprintf(1,'Matching annotated genes to %u ABA-genes.\n', length(gene_symbols));

% get category sizes
category_sizes_old = cellfun(@length, cTable.cGenes);
    
% Make sure all cGenes (& cWeights) match those in our set of gene_symbols
match_me = cellfun(@(x) ismember(x, gene_symbols), cTable.cGenes, 'UniformOutput', false);
cTable.cGenes = cellfun(@(x,y) x(y), cTable.cGenes, match_me, 'UniformOutput', false);
try
    cTable.cWeights = cellfun(@(x,y) x(y), cTable.cWeights, match_me, 'UniformOutput', false);
end

% Renew category size column
category_sizes = cellfun(@length, cTable.cGenes);
fprintf(1, '%u/%u annotated genes are in gene expression dataset.\n', ...
        sum(category_sizes), sum(category_sizes_old));
cTable.cSize = category_sizes;

fprintf(1, '%u categories have no annotated genes matching our %u genes.\n',...
        sum(category_sizes==0), length(gene_symbols));

%--------------------------------------------------------------------------
% Filter genes by weights and adjust weights
%--------------------------------------------------------------------------
% ADD RESCALE AND STANDARDIZE? (FOR WHOLE DATASET OR PER CATEGORY?)

% Modify genes annotated to categories by weights
if any(ismember(cTable.Properties.VariableNames, 'cWeights'))
    
    % Treshold genes at quantile if quantile given
    if ~isempty(weights_quant)
        % get weights distribution of every unique gene
        all_genes = vertcat(cTable.cGenes{:});
        all_weights = cell2mat(cTable.cWeights(:));
        [~, match_me, ~] = unique(all_genes);
        weights_distr = all_weights(match_me);
        % get threshold value
        quant = quantile(weights_distr, weights_quant);
        % apply threshold
        fprintf('Thresholding genes at %.2fth quantile (weight > %2.2f).\n', weights_quant, quant);
        match_me = cellfun(@(x) x > quant, cTable.cWeights, 'UniformOutput', false);
        cTable.cGenes = cellfun(@(x,y) x(y), cTable.cGenes, match_me, 'UniformOutput', false);  
        cTable.cWeights = cellfun(@(x,y) x(y), cTable.cWeights, match_me, 'UniformOutput', false);
    end
    % Treshold genes at value if threshold given
    if ~isempty(weights_thresh) 
        fprintf('Thresholding genes at weight > %.2f.\n', weights_thresh);
        match_me = cellfun(@(x) x > weights_thresh, cTable.cWeights, 'UniformOutput', false);
        cTable.cGenes = cellfun(@(x,y) x(y), cTable.cGenes, match_me, 'UniformOutput', false);  
        cTable.cWeights = cellfun(@(x,y) x(y), cTable.cWeights, match_me, 'UniformOutput', false);
    end
end

% Set all weights to 1 if no weights in dataset or if cutoff == true 
%   -> 'weightedmean' will also calculate mean
if ~any(ismember(cTable.Properties.VariableNames, 'cWeights')) || weights_cutoff 
    disp('Using binarized category-gene-annotations.');
    cTable.cWeights(:) = {[]};
    cTable.cWeights = cellfun(@(x) ones(length(x),1), cTable.cGenes, 'UniformOutput', false); 
else
    disp('Using weighted category-scores (if aggregation_method is <weightedmean>).');
end

%--------------------------------------------------------------------------
% Filter genes by co-occurence across categories
%--------------------------------------------------------------------------

% remove genes that are annotated for >= gene_coocc_thresh of categories
if gene_coocc_thresh
    % get proportion of genes 
    genes_unique = unique(vertcat(cTable.cGenes{:}));
    genes_unique_categ = cellfun(@(x) ismember(genes_unique, x), cTable.cGenes, 'UniformOutput', false);
    is_genes_unique_categ = horzcat(genes_unique_categ{:});
    genes_unique_prop = sum(is_genes_unique_categ, 2) / height(cTable);
    % get over-threshold genes
    genes_to_remove = genes_unique(genes_unique_prop >= gene_coocc_thresh);
    % remove these genes
    match_me = cellfun(@(x) ~ismember(x, genes_to_remove), cTable.cGenes, 'UniformOutput', false);
    cTable.cGenes = cellfun(@(x,y) x(y), cTable.cGenes, match_me, 'UniformOutput', false);  
    cTable.cWeights = cellfun(@(x,y) x(y), cTable.cWeights, match_me, 'UniformOutput', false);
    fprintf(1, 'Removed %u genes annotated in >= %u percent of categories.\n', ...
            length(genes_to_remove), gene_coocc_thresh*100);
end

% Renew category size column
category_sizes = cellfun(@length, cTable.cGenes);
cTable.cSize = category_sizes;

%--------------------------------------------------------------------------
% Filter categories by size
%--------------------------------------------------------------------------
fprintf(1, 'Categories have between %u and %u annotated genes.\n',...
        min(category_sizes), max(category_sizes));
    
% filter categories by size
is_good_size = (category_sizes >= size_filter(1)) & (category_sizes <= size_filter(2));
cTable = cTable(is_good_size, :);
fprintf(1,'Filtered to %u categories with between %u and %u annotated genes.\n',...
        height(cTable), size_filter(1), size_filter(2));

%--------------------------------------------------------------------------
% Stop if no categories left empty
%--------------------------------------------------------------------------
if isempty(cTable)
    error('Dataset is empty after filtering. Check size filter and thresholding settings!')
end

end
