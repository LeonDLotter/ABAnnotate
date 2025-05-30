function cTable = estimate_pvals(null_scores, real_scores, tail, cTable)
%function cTable = estimate_pvals(null_scores, real_scores, tail, cTable)
%
% Estimates p-values for each gene category from a given null distribution.
% Computes using both a permutation-based and fitted Gaussian approximation.
%
% INPUTS:
% - nullScores (cell): a set of null scores for each category
% - realScores (vector): scores for the same set of categories obtained 
%       from real data
% - tail ('right' or 'left'): whether to use a right- or left-tailed test.
%       'right': tests whether each realScores is greater than the null.
%       'left': tests whether each realScores is less than the null.
% - cTable (table): the table to save the computations back to.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Checks on matching null_scores, real_scores, cTable:
%--------------------------------------------------------------------------
n_categories = height(cTable);
assert(n_categories==length(null_scores))
assert(n_categories==length(real_scores))

%--------------------------------------------------------------------------
% Estimate category p-values, looping over categories
%--------------------------------------------------------------------------
pValPerm = zeros(n_categories,1);
pValZ = zeros(n_categories,1);
for i = 1:n_categories
    score_here = real_scores(i);
    null_here = null_scores{i};
    switch tail
    case 'right'
        pValPerm(i) = mean(null_here >= score_here);
        pValZ(i) = 1 - normcdf(score_here, mean(null_here), std(null_here));
    case 'left'
        pValPerm(i) = mean(null_here <= score_here);
        pValZ(i) = normcdf(score_here, mean(null_here), std(null_here));
    otherwise
        error('Unknown tail setting, ''%s''',tail)
    end
end

%--------------------------------------------------------------------------
% Multiple hypothesis correction using Benjamini-Hochberg (FDR):
%--------------------------------------------------------------------------
[~, ~, ~, pValPermCorr] = fdr_bh(pValPerm);
[~, ~, ~, pValZCorr] = fdr_bh(pValZ);

%--------------------------------------------------------------------------
% Assign values to categories of cTable:
%--------------------------------------------------------------------------
cTable.pValZ = pValZ;
cTable.pValZCorr = pValZCorr;
cTable.pValPerm = pValPerm;
cTable.pValPermCorr = pValPermCorr;

end
