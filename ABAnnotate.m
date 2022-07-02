function [cTablePheno, opt] = ABAnnotate(opt)
%[cTablePheno, opt] = ABAnnotate(opt)
%
% Performs Category Enrichment Analysis on gene expression data from the
% Allen Human Brain Atlas. Analysis settings are defined in a settings 
% structure. 'phenotype' OR 'phenotype_data' are the only mandatory input 
% variables.

%% ------------------------------------------------------------------------
% Start
% -------------------------------------------------------------------------

% stop if no settings structure
if nargin < 1
    error('You must provide settings structure as input!')
end
    
% stop if no pheno volume or data, or if only empty pheno datafields 
if ~isfield(opt, 'phenotype') && ~isfield(opt, 'phenotype_data')
    error(['You must provide either phenotype volume (opt.phenotype) ' ...
           'or data (opt.phenotype_data)!']);
end

% get directories & start logging
wd = fileparts(which('ABAnnotate'));
tempdir = fullfile(wd, 'temp');
if ~exist(tempdir, 'dir')
    mkdir(tempdir);
end
logfile = fullfile(tempdir, 'logfile_temp.txt');
diary(logfile);

% update sources table
try
    abannotate_get_sources();
catch
    disp('Could not update dataset_sources.csv from OSF, using existing file.');
end

% get analysis name and print section
if ~isfield(opt, 'analysis_name') && isfield(opt, 'phenotype')
    [~, opt.analysis_name, ~] = fileparts(opt.phenotype);
elseif ~isfield(opt, 'analysis_name')
    opt.analysis_name = sprintf('enrich_%s', date);
else
    opt.analysis_name = char(opt.analysis_name);
end
print_section(sprintf('Starting ABAnnotate: %s', opt.analysis_name));


%% ------------------------------------------------------------------------
% Prepare imaging files & null maps
% -------------------------------------------------------------------------
 
% set atlas & aba data ----------------------------------------------------

% default: SchaeferTian-116
if ~isfield(opt, 'atlas'); opt.atlas = 'Schaefer100-7_TianS1'; end
% choose atlas - currently only SchaeferTian116, Schaefer100, and Neuromorphometrics defined
opt.atlas = char(opt.atlas);
switch opt.atlas
    case {'SchaeferTian', 'Schaefer100-7_TianS1'}
        atlas = 'Schaefer100-7_TianS1';
        atlas_type = 'integrated';
    case {'Schaefer', 'Schaefer100-7', 'Schaefer100'}
        atlas = 'Schaefer100-7';
        atlas_type = 'integrated';
    case {'Neuromorphometrics', 'neuromorphometrics'}
        atlas = 'Neuromorphometrics';
        atlas_type = 'integrated';
    otherwise
        atlas_type = 'user';
end

% get atlas
switch atlas_type 
    
    case 'integrated'
        atlas_dir = fullfile(wd, 'atlas', atlas);
        if ~exist(atlas_dir, 'dir')
            % download zipped atlas dir from OSF
            fprintf('Atlas %s not found. Downloading...\n', atlas);
            atlas_list = abannotate_get_datasets('parcellation', false);
            osf_id = atlas_list.osf{strcmp(atlas_list.name, atlas)};
            atlas_zip = abannotate_download_osf(osf_id, [atlas_dir '.zip'], true);
            % unzip & delete
            unzip(atlas_zip, atlas_dir);
            delete(atlas_zip);
        end
        % set atlas
        opt.atlas = fullfile(wd, 'atlas', atlas, sprintf('%s_atlas.nii', atlas));
        opt.aba_mat = fullfile(wd, 'atlas', atlas, sprintf('%s_expression.mat', atlas));
        load(opt.aba_mat, 'expression_matrix')
        opt.n_rois = size(expression_matrix,1);
        clear('expression_matrix');
        fprintf('Setting atlas & ABA data to %s - %u parcels.\n', opt.atlas, opt.n_rois);
        
    case 'user'
        % check if custom volume. if yes, unzip if ending *.nii.gz
        try
            % load & unzip atlas
            [~, ~, atlas_ext] = fileparts(opt.atlas);
            if strcmp(atlas_ext, '.gz')
                atlas_unzip = gunzip(opt.atlas, fullfile(wd, 'temp'));
                opt.atlas = atlas_unzip{1};
            elseif strcmp(atlas_ext, '.nii')
            else
                error('Filetype of %s is not supported!', opt.phenotype);
            end
        catch
            error(['Atlas not valid. Choose one of the integrated atlases or ', ...
                   'provide a path to a custom atlas volume! ', ...
                   'Run <abannotate_get_datasets("parcellation");> for a list of available atlases.']);
        end
        % get roi number
        opt.n_rois = get_number_of_rois(opt.atlas);
        fprintf('Using custom provided atlas: %s - %u parcels.\n', opt.atlas, opt.n_rois);
        % check associated ABA expression matrix
        if ~isfield(opt, 'aba_mat')
            error('If you provide a custom atlas, you must also provide matching ABA data!')
        end
        try
            load(opt.aba_mat, 'expression_matrix', 'gene_symbols');
        catch
            error('Problems loading variables "expression_matrix" and "gene_symbols" from opt.aba_mat!');
        end
        % compare numbers
        if opt.n_rois~=size(expression_matrix,1)
            error('Number of ROIs in atlas (%u) does not match number of values in expression matrix (%u)!', opt.n_rois, height(expression_matrix))
        end
        if length(gene_symbols)~=size(expression_matrix,2)
            error('Length of gene_symbols width of expression_matrix do not match!');
        end
        clear('expression_matrix', 'gene_symbols');
end

% parcellate pheno volume -------------------------------------------------

% if pheno data provided, will use data and ignore volume!
if ~isfield(opt, 'phenotype_data')
    fprintf('Applying parcellation to %s.\n', opt.phenotype);
    % unzip and check file type
    [~, ~, phenotype_ext] = fileparts(opt.phenotype);
    if strcmp(phenotype_ext, '.gz')
        phenotype_unzip = gunzip(opt.phenotype, fullfile(wd, 'temp'));
    elseif strcmp(phenotype_ext, '.nii')
        phenotype_unzip = cellstr(opt.phenotype);
    else
        error('Filetype of %s is not supported!', opt.phenotype);
    end
    % parcellate
    [phenotype_parc, ~] = mean_time_course(phenotype_unzip, opt.atlas, 1:opt.n_rois);
    opt.phenotype_data = phenotype_parc';
    
% else check if input data is provided
elseif isfield(opt, 'phenotype_data')
    disp('Phenotype data array provided, overriding opt.phenotype with "data".');
    opt.phenotype = 'data';
    if ~isequal(size(opt.phenotype_data), [opt.n_rois 1])
        error('Phenotype data has to be vector of size [%d, 1]!', opt.n_rois);
    end

% else something did not work
else
    error('It seems you did not provide opt.phenotype or opt.phenotype_data!'); 
    end

% check for NaNs
if any(isnan(opt.phenotype_data))
    warning('Parcellated phenotype data contains NaNs! Respective regions will be ignored.');
end
    
% get phenotype nulls -----------------------------------------------------

% number of phenotype nulls
if ~isfield(opt, 'n_nulls'); opt.n_nulls = 1000; end

% if nulls provided, check if correct array name and size
if isfield(opt, 'phenotype_nulls')
    disp('Using provided null maps.');
    % check array name
    try
        load(opt.phenotype_nulls, 'nullMaps')
    catch
        error(['Error while loading null map. Either the file does not exist ', ...
               'or does not contain an array named <nullMaps>!']);
    end
    % check if correct size
    phenotype_nulls_size = size(nullMaps);
    if phenotype_nulls_size(1)~=opt.n_rois || phenotype_nulls_size(2)<opt.n_nulls
        error(['Provided null maps have wrong size [%u %u], ',...
               'must be n_rois x n_nulls array [%u %u]!'], ...
              phenotype_nulls_size(1), phenotype_nulls_size(2), opt.n_rois, opt.n_nulls)
    end
    
% if not provided, load or generate new phenotype nulls
else
    disp('Checking if phenotype null maps exist...')
    % path to pheno nulls (n_rois x n_nulls array)
    opt.phenotype_nulls = ...
        fullfile(wd, 'nulls', sprintf('phenonulls_%s.mat', opt.analysis_name));
    
    % if exist, check n_nulls and load
    if isfile(opt.phenotype_nulls)
        load(opt.phenotype_nulls, 'nullMaps')
        
        % check if size fits n_nulls
        phenotype_nulls_size = size(nullMaps);
        if phenotype_nulls_size(1)==opt.n_rois && phenotype_nulls_size(2)>=opt.n_nulls
            fprintf('Using existing phenotype null maps.\n')
            
        % if not enough maps, generate new
        else
            fprintf(['Existing null maps file has < %u maps. ', ...
                     'Generating new, may take a while...\n'], opt.n_nulls);
            % run null generation
            generate_phenotype_nulls(...
                opt.phenotype_data', opt.atlas, opt.n_nulls, opt.phenotype_nulls);
        end
        
    % if not exist, generate new
    else
        fprintf(['No null maps file found. Generating %u null maps, ', ...
                 'may take a while...\n'], opt.n_nulls);
        % run null generation
        generate_phenotype_nulls(...
            opt.phenotype_data', opt.atlas, opt.n_nulls, opt.phenotype_nulls);
    end
    fprintf('Phenotype null maps in: %s\n', opt.phenotype_nulls);
end


% -------------------------------------------------------------------------
% Gene Category Enrichment Analysis
% -------------------------------------------------------------------------

% if GCEA settings structure field does not exist, make one and set n_nulls
if ~isfield(opt, 'GCEA'); opt.GCEA = struct(); end
if ~isfield(opt.GCEA, 'n_category_nulls'); opt.GCEA.n_category_nulls = opt.n_nulls; end

% Dataset -----------------------------------------------------------------
% default to Gene Ontology - Biological Process
% case 1: none set, default to GO-BP
if ~isfield(opt.GCEA, 'dataset') && ~isfield(opt.GCEA, 'dataset_mat')
    opt.GCEA.dataset = 'GO-biologicalProcessDirect-discrete';
    fprintf('No dataset provided, defaulting to: %s\n', opt.GCEA.dataset);
    % get path or download
    opt.GCEA.dataset_mat = get_gcea_dataset(opt.GCEA.dataset);
% case 2: one of integrated datasets
elseif isfield(opt.GCEA, 'dataset') && ~strcmp(opt.GCEA.dataset, 'custom')
    % get path or download
    opt.GCEA.dataset_mat = get_gcea_dataset(opt.GCEA.dataset);
    fprintf('GCEA dataset: %s\n', opt.GCEA.dataset);
    validate_gcea_dataset(opt.GCEA.dataset_mat);
% case 3: a custom mat file
elseif strcmp(opt.GCEA.dataset, 'custom') && isfield(opt.GCEA, 'dataset_mat')
    fprintf('Custom GCEA dataset: %s\n', opt.GCEA.dataset_mat);
    validate_gcea_dataset(opt.GCEA.dataset_mat);
% case 4: other
else 
    disp(['You must provide ether a valid dataset name or set <opt.gcea.dataset> to ', ...
          '<custom> and provide an own dataset *.mat file. Possible datasets are:\n']);
    get_gcea_datasets();
    error('Invalid dataset.');      
end

print_section('Running GCEA');

% Generate Category null samples ------------------------------------------

if ~isfield(opt.GCEA, 'category_nulls')
    % get default settings for each field if not provided by user
    opt.GCEA = merge_struct(gcea_default_settings, opt.GCEA);
    % file path to save results to
    opt.GCEA.category_nulls = ...
        fullfile(wd, 'nulls', sprintf('categonulls_%s_%s.mat', opt.GCEA.dataset, opt.analysis_name));
    % print complete settings
    print_ABAnnotate_input(opt);
    % run
    disp('Generating category null samples...')
    generate_category_nulls(opt.GCEA.dataset_mat, ... 
                            opt.aba_mat, ...
                            opt.GCEA, ...
                            opt.phenotype_nulls, ...
                            [], ...
                            true);  
else
    fprintf('Loading category null samples from %s\n', opt.GCEA.category_nulls);
    % check if correct dataset
    load(opt.GCEA.category_nulls, 'gcea_opt');
    if ~strcmp(gcea_opt.dataset, opt.GCEA.dataset)
        error(['Provided dataset nulls seem to be for wrong dataset, define ', ...
               'correct nulls, change dataset or rerun null sample generation!']);
    end
    % print complete settings
    print_ABAnnotate_input(opt);
end

% Compute Enrichment ------------------------------------------------------

cTablePheno = ensemble_enrichment(opt.GCEA.dataset_mat, ...
                                  opt.aba_mat, ...
                                  opt.GCEA.category_nulls, ...
                                  opt.phenotype_data);
                            
                               
%% ------------------------------------------------------------------------
% Save
% -------------------------------------------------------------------------
if isfield(opt, 'dir_result')
 
    % write csv
    cTablePheno_csv = cTablePheno;
    cTablePheno_csv = removevars(cTablePheno_csv, {'cWeights', 'cScoresNull'});
    cTablePheno_csv.cGenes = cellfun(@(x) strjoin(x, ', '), cTablePheno.cGenes, 'UniformOutput', false);
    writetable(cTablePheno_csv, fullfile(opt.dir_result, sprintf('GCEA_%s.csv', opt.analysis_name)));
        
    % write mat
    save(fullfile(opt.dir_result, sprintf('GCEA_%s.mat', opt.analysis_name)), 'cTablePheno', 'opt');
     
    % write setting xml
    writestruct(opt, fullfile(opt.dir_result, sprintf('GCEA_%s_settings.xml', opt.analysis_name)))
    
    % copy logfile
    copyfile(logfile, fullfile(opt.dir_result, sprintf('GCEA_%s_log.txt', opt.analysis_name)));

    % print
    fprintf('Saved results to %s.\n', fullfile(opt.dir_result, sprintf('GCEA_%s*.*', opt.analysis_name)));
end

% delete temp data
diary('off');
delete(fullfile(tempdir, '*'));

print_section(sprintf('Finished ABAnnotate: %s', opt.analysis_name));
     
end % function


    
    
    
