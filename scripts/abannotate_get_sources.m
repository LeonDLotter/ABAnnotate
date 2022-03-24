function [sources] = abannotate_get_sources()
%function [sources] = abannotate_get_sources()
%
% Download an overview of available ABAnnotate datasets from OSF.
    
% file paths
wd = fileparts(which('ABAnnotate'));
sources_path = fullfile(wd, 'dataset_sources.csv');

% download
disp(['(Re-)downloading available ABAnnotate datasets to: ' sources_path])
try
    abannotate_download_osf('623c2c81e6b58b0b8ad6b79a', sources_path);
catch
    error('Download failed. Check your internet connection?!');
end

% read table
sources = readtable(sources_path);

end
    