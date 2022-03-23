function go_desc = get_go_desc(go_ids)
%function go_desc = get_go_desc(go_terms)
%
% returns a list of go descriptors given a list of go IDs

go_dataset = fullfile(fileparts(which('ABAnnotate')), 'datasets', 'GO-biologicalProcessProp-discrete.mat');
load(go_dataset, 'cTable');

go_desc = cell([length(go_ids), 1]);
for i=1:length(go_ids)
    id = go_ids{i};
    match_me = ismember(cTable.cLabel, id);
    go_desc(i) = cTable.cDesc(match_me);
end

end
