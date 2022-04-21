function category_score = aggregate_scores(...
    gene_scores, gene_weights, aggregation_method)
% Aggregates a set of scores of genes in a category in a category-wide score
%
% geneScores are aggregated down columns (assuming each phenotype is a column)
%-------------------------------------------------------------------------------

if nargin < 2
    aggregation_method = 'weightedmean';
end

%---------------------------------------------------------------------------
% Aggregate gene-wise scores into an overall GO category score (for each phenotype)
switch aggregation_method
    case 'mean'
        category_score = nanmean(gene_scores,1);
    case 'absmean'
        category_score = nanmean(abs(gene_scores),1);
    case 'median'
        category_score = nanmedian(gene_scores,1);
    case 'absmedian'
        category_score = nanmedian(abs(gene_scores),1);
    case 'weightedmean'
        category_score = nansum(gene_scores .* gene_weights, 1) ./ nansum(gene_weights, 1);
    case 'absweightedmean'
        category_score = nansum(abs(gene_scores) .* gene_weights, 1) ./ nansum(gene_weights, 1);
    otherwise
        error('Unknown aggregation option: "%s"', aggregation_method);
end

end
