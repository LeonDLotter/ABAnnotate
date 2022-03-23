function [default_opt] = gcea_default_settings()
%function [] = gcea_default_settings()


%% default settings -------------------------------------------------------

% choose categories
default_opt.categories_select = [];

% filter categories by size (number of genes; [min, max])
default_opt.size_filter = [1, inf];

% nr of category null samples to compute
default_opt.n_category_nulls = 1000;

% what type of correlation to use for phenotype x gene correlations 
%   ('Pearson', 'Spearman')
default_opt.correlation_method = 'Spearman';

% rescale weights? (default to no)
%if ~isfield(gcea_opt, 'weights_rescale'); gcea_opt.weights_rescale = false; end

% threshold gene weights at quantile? (default to no)
default_opt.weights_quant = []; 

% threshold gene weights at value? (default to no)
default_opt.weights_thresh = []; 

% binarize gene weights? (default to no)
default_opt.weights_cutoff = false; 

% remove genes that co-occurre in >= gene_coocc_thresh * 100% of categories? 
%   (default to 1 -> will remove genes existing in every category)
default_opt.gene_coocc_thresh = 1; 

% how to aggregate scores within a category ('mean', 'absmean', 'median', 
%   'absmedian', 'weightedmean', 'absweightedmean')
default_opt.aggregation_method = 'weightedmean'; 

% higher or lower scores are 'better' (compute p-value from right or left tail)
default_opt.p_tail = 'right'; 

% significance threshold for categories (only for display)
default_opt.p_thresh = 0.05; 


end