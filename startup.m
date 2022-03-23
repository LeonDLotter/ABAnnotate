% Add paths required for the project (ignoring hidden, including version control)
files = dir;
mainpath = files.folder;
directories = files([files.isdir]);
directories(strmatch('.',{files([files.isdir]).name})) = []; % remove hidden
paths = arrayfun(@(x)fullfile(directories(x).folder,directories(x).name),1:length(directories),'UniformOutput',false);

fprintf('Adding %s to matlab path\n', mainpath);
addpath(genpath(mainpath))
for j = 1:length(paths)
    fprintf('Adding %s to matlab path\n', paths{j});
    addpath(genpath(paths{j}))
end

clear;
