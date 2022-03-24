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

% get directory & start logging
wd = fileparts(which('ABAnnotate'));
logfile = fullfile(wd, 'temp', 'logfile_temp.txt');
diary(logfile);

% get analysis name and start section
if ~isfield(opt, 'analysis_name') && isfield(opt, 'phenotype')
    [~, opt.analysis_name, ~] = fileparts(opt.phenotype);
elseif ~isfield(opt, 'analysis_name')
    opt.analysis_name = ['enrich_' date];
end
print_section(['Starting ABAnnotate: ' opt.analysis_name]);


%% ------------------------------------------------------------------------
% Prepare imaging files & null maps
% -------------------------------------------------------------------------
 
% set atlas & aba data ---------------------------------------------------
% default: SchaeferTian-116
if ~isfield(opt, 'atlas'); opt.atlas = 'SchaeferTian'; end
% choose atlas - currently only SchaeferTian116 and Neuromorphometrics defined
switch opt.atlas
    case 'SchaeferTian'
        atlas = 'Schaefer100-7_TianS1';
        atlas_type = 'integrated';
    case 'Neuromorphometrics'
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
            disp('Atlas not found. Downloading...');
            atlas_list = abannotate_get_datasets('parcellation', false);
            osf_id = atlas_list.osf{strcmp(atlas_list.name, atlas)};
            atlas_zip = abannotate_download_osf(osf_id, [atlas_dir '.zip'], true);
            % unzip & delete
            unzip(atlas_zip, atlas_dir);
            delete(atlas_zip);
        end
        % set atlas
        opt.atlas = fullfile(wd, 'atlas', [atlas '_atlas.nii']);
        opt.aba_mat = fullfile(wd, 'atlas', [atlas '_expression.mat']);
        load(opt.aba_mat, 'expression_matrix')
        opt.n_rois = height(expression_matrix);
        clear('expression_matrix');
        fprintf('Setting atlas & ABA data to %s - %u parcels.\n', opt.atlas, opt.n_rois)
    case 'user'
        % check if custom volume. if yes, unzip if ending *.nii.gz
        try
            [~, ~, atlas_ext] = fileparts(opt.atlas);
            if strcmp(atlas_ext, '.gz')
                atlas_unzip = gunzip(opt.atlas, fullfile(wd, 'temp'));
                opt.atlas = atlas_unzip{1};
            elseif strcmp(atlas_ext, '.nii')
            else
                error('Filetype of %s is not supported!', opt.phenotype);
            end
            disp(['Using custom provided atlas: ' opt.atlas]);
            if ~isfield(opt, 'aba_mat')
                error('If you provide a custom atlas, you must also provide matching ABA data!')
            end
        catch
            error(['Current included options for opt.atlas are: <SchaeferTian>, ', ...
                   '<Neuromorphometrics> or a path to a custom atlas volume!']);
        end
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
        phenotype_unzip = {opt.phenotype};
    else
        error('Filetype of %s is not supported!', opt.phenotype);
    end
    % parcellate
    [phenotype_parc, ~] = mean_time_course(phenotype_unzip, opt.atlas, 1:opt.n_rois);
    opt.phenotype_data = phenotype_parc';
% else check if input data is provided
elseif isfield(opt, 'phenotype_data')
    if size(opt.phenotype_data) ~= [opt.n_rois 1]
        error('Phenotype data has to be opt.n_rois x 1 vector!');
    end
    disp('Phenotype data array provided, ignoring opt.phenotype');
% else something did not work
else
    error('It seems you did not provide opt.phenotype or opt.phenotype_data!'); 
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
    if phenotype_nulls_size(1) ~= opt.n_rois
        error('Null maps have wrong size, must be n_rois x n_nulls array!')
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
% case 1: one of integrated datasets
if isfield(opt.GCEA, 'dataset') && ~strcmp(opt.GCEA.dataset, 'custom')
    % get path or download
    opt.GCEA.dataset_mat = get_gcea_dataset(opt.GCEA.dataset);
    disp(['GCEA dataset: ' opt.GCEA.dataset]);
    validate_gcea_dataset(opt.GCEA.dataset_mat);
% case 2: a custom mat file
elseif strcmp(opt.GCEA.dataset, 'custom') && isfield(opt.GCEA, 'dataset_mat')
    disp(['Custom GCEA dataset: ' opt.GCEA.dataset_mat]);
    validate_gcea_dataset(opt.GCEA.dataset_mat);
% case 3: none set, default to GO-BP
elseif ~isfield(opt.GCEA, 'dataset') && ~isfield(opt.GCEA, 'dataset_mat')
    opt.GCEA.dataset = 'GO-biologicalProcessDirect-discrete';
    % get path or download
    opt.GCEA.dataset_mat = get_gcea_dataset(opt.GCEA.dataset);
    disp(['No dataset provided, defaulting to: ' opt.GCEA.dataset]);
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
                            true, true);  
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
    writetable(cTablePheno_csv, fullfile(opt.dir_result, ['GCEA_' opt.analysis_name '.csv']));
        
    % write mat
    save(fullfile(opt.dir_result, ['GCEA_' opt.analysis_name '.mat']), 'cTablePheno', 'opt');
     
    % write setting xml
    writestruct(opt, fullfile(opt.dir_result, ['GCEA_' opt.analysis_name '_settings.xml']))
    
    % copy logfile
    copyfile(logfile, fullfile(opt.dir_result, ['GCEA_' opt.analysis_name '_log.txt']));

    % print
    fprintf('Saved results to %s.\n', fullfile(opt.dir_result, ['GCEA_' opt.analysis_name '*.*']));
end

% delete temp data
diary('off');
delete(fullfile(wd, 'temp', '*'));

print_section(['Finished ABAnnotate: ' opt.analysis_name]);
     
end % function


    
    
    
