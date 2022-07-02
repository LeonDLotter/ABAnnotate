# ABAnnotate Customization Options: 

## Default

```matlab
opt.analysis_name = 'pain_celltypes';
opt.n_nulls = 1000;
opt.phenotype = 'neuroquery_pain.nii';
opt.dir_result = 'gcea';
opt.GCEA.dataset = 'PsychEncode-cellTypesTPM-discrete';
% pass 'opt.GCEA.verbose_cat = true' to print process for each category
% opt.GCEA.verbose_cat = false

cTable = ABAnnotate(opt);
```

## Other integrated parcellation

```matlab
opt.analysis_name = 'pain_celltypes';
opt.n_nulls = 1000;
opt.atlas = 'Neuromorphometrics' % one of SchaeferTian, Neuromorpometrics, Schaefer
opt.phenotype = 'neuroquery_pain.nii';
opt.dir_result = 'gcea';
opt.GCEA.dataset = 'PsychEncode-cellTypesTPM-discrete';

cTable = ABAnnotate(opt);
```

## Pass parcellated phenotype data directly

```matlab
opt.analysis_name = 'pain_celltypes';
opt.n_nulls = 1000;
opt.phenotype_data = parcellated_volume; % a vector with size (number of rois, x 1), 
    % must be parcellated with one of the integrated atlases
opt.dir_result = 'gcea';
opt.GCEA.dataset = 'PsychEncode-cellTypesTPM-discrete';

cTable = ABAnnotate(opt);
```

## Pass null phenotype data directly

```matlab
opt.analysis_name = 'pain_celltypes';
opt.n_nulls = 1000;
opt.phenotype = 'neuroquery_pain.nii';
opt.phenotype_nulls = 'path/to/pheno_nulls.mat' % usefull to save the null map generation step, 
    % or to use custom null maps. ABAnnotate saves the null volumes in
    % 'ABAnnotate/nulls/phenonulls_analysis_name.mat'. If self-generated,
    % must be a matrix of size (number of rois, number of null maps)
opt.dir_result = 'gcea';
opt.GCEA.dataset = 'PsychEncode-cellTypesTPM-discrete';

cTable = ABAnnotate(opt);
```

## Use other integrated GCEA datasets

```matlab
% print a list of integrated datasets and replace 'opt.GCEA.dataset' with their name
abannotate_get_datasets;
```

```matlab
opt.analysis_name = 'pain_celltypes';
opt.n_nulls = 1000;
opt.phenotype = 'neuroquery_pain.nii';
opt.dir_result = 'gcea';
opt.GCEA.dataset = 'one_of_the_printed_datasets';

cTable = ABAnnotate(opt);
```

## Use custom GCEA datasets

```matlab
opt.analysis_name = 'pain_celltypes';
opt.n_nulls = 1000;
opt.phenotype = 'neuroquery_pain.nii';
opt.dir_result = 'gcea';
opt.GCEA.dataset = 'custom'; % this indicates a custom dataset
opt.GCEA.dataset_mat = 'path/to/custom_dataset.mat'; % the path to your dataset, must be a *.mat 
    % file with a table named 'cTable' and three mandatory columns: 
    % 'cLabel': category labels, 
    % 'cSize': the number of genes annotated to each category, and 
    % 'cGenes': for each row i a cell array with the official gene symbols annotated to each 
        % category with size (cSize(i), 1). Every additional column (e.g., category descriptions
        % or IDs) will be passed unchanged to the output table.
    % ('cWeights': if you want to include a 'weighted' dataset with abitray high numbers of genes 
        % annotated to each category and additional expression values associated to each gene in 
        % each category: a numeric vector with size (cSize(i), 1). See the BrainSpan dataset for an 
        % example (~12,000 expression values for each category))
    % see the scripts in 'ABAnnotate/scripts_import/' for info how integrated datasets were created

cTable = ABAnnotate(opt);
```

## GCEA options: general, and 'discrete' datasets

```matlab
opt.analysis_name = 'pain_celltypes';
opt.n_nulls = 1000;
opt.phenotype = 'neuroquery_pain.nii';
opt.dir_result = 'gcea';
opt.GCEA.dataset = 'PsychEncode-cellTypesTPM-discrete'; % 'discrete' datasets are those with a 
    % relatively small number of genes annotated to each category without additional weights 
    % associated with each gene in each category
    
% GCEA options:
opt.GCEA.categories_select = []; % choose categories by cLabel
opt.GCEA.size_filter = [1, inf]; % filter categories by cSize (number of genes; [min, max])
opt.GCEA.correlation_method = 'Spearman'; % what type of correlation to use for 
    % phenotype x gene correlations ('Pearson', 'Spearman')
opt.GCEA.aggregation_method = 'mean'; % how to aggregate scores within a category 
    % for discrete datasets one of ('mean', 'absmean', 'median', 'absmedian'), 
    % 'abs...' will ignore the sign of the phenotype-gene correlation 
    % if providing 'weightedmean' but a discrete datasets, will effectively use 'mean'.
opt.GCEA.p_tail = 'right'; % higher or lower scores are 'better' (compute p-value from right or left tail)
opt.GCEA.p_thresh = 0.05; % significance threshold for categories (only for display)

cTable = ABAnnotate(opt);
```

## GCEA options: 'weighted' datasets

```matlab
opt.analysis_name = 'pain_celltypes';
opt.n_nulls = 1000;
opt.phenotype = 'neuroquery_pain.nii';
opt.dir_result = 'gcea';
opt.GCEA.dataset = 'ABA-brainSpan-weights'; % Currently, the only 'non-discrete' dataset
    % included in ABAnnotate is the BrainSpan dataset. See "Use custom dataset" above
    
% GCEA options:
opt.GCEA.size_filter = [1, inf]; % filter categories by cSize (number of genes; [min, max])
    % with weighted datasets, the size filter is applied after category thresholding, see below
opt.GCEA.aggregation_method = 'weightedmean'; % how to aggregate scores within a category 
    % for weighted datasets one of ('mean', 'absmean', 'median', 'absmedian', 'weightedmean', 
    % 'absweightedmean'), 'weightedmean' will calculate category scores by taking the category-wise 
    %  mean of phenotype-gene associations weighted by the category-specific gene expression value ("weight")
default_opt.weights_quant = 0.9; % this would threshold gene weights at the 0.9th quantile of the whole dataset
default_opt.weights_thresh = 10; % this would threshold gene weights at the value 10
default_opt.weights_cutoff = true; % this would binarize gene weights after thresholding
    % if 'aggregation_method' is 'weightedmean', now effectively the 'mean' would be calculated
    % (cWeights are set to 1 for a each gene in each category)
default_opt.gene_coocc_thresh = 0.2; % this removes genes that co-occurre in >= 0.2 * 100% of categories
    % (defaults to 1 -> will only remove genes existing in every category)

cTable = ABAnnotate(opt);
```

## Use a custom parcellation and generate an ABA dataset to use for GCEA

Use `ABAnnotate/scripts_import/get_aba_data.py` in python to load and process ABA data and get parcel-wise expression values. This requires the [abagen](https://abagen.readthedocs.io/) toolbox and will download ~4GB of ABA data if not done before. Outputs a csv file with parcel data and a markdown wile with a report on the processing pipeline used.

```python
from get_aba_data import get_aba_data

get_aba_data(
    atlas='path/to/parcellation_volume.nii.gz',
    save_dir='path/to/save_dir',
    save_prefix='file_prefix'
    # any addtional arguments are passed to abagen.get_expression_data, otherwise, will use default settings
    # see https://abagen.readthedocs.io/en/stable/generated/abagen.get_expression_data.html#abagen.get_expression_data
)
```

Now, import the generated expression value csv file for use with ABAnnotate:

```matlab
import_aba_data( ...
    'path/to/aba_csv_from_python_script.csv', ...
    'path/to/parcellation_volume.nii', ...
    'path/to/save_dir', ...
    'file_prefix')
```

Use the new ABA dataset:

```matlab
opt.analysis_name = 'pain_celltypes';
opt.n_nulls = 1000;
opt.phenotype = 'neuroquery_pain.nii';
opt.dir_result = 'gcea';
opt.atlas = 'path/to/custom_parcellation.nii' % the parcellation volume you used to generate the dataset
opt.aba_mat = 'path/to/custom_aba_dataset.mat' % the output dataset generated by import_aba_data.m
opt.GCEA.dataset = 'PsychEncode-cellTypesTPM-discrete';

cTable = ABAnnotate(opt);
```

For a full working example, see [`example_pain.md`](example_pain.md).   
Feel free to [contact me](mailto:leondlotter@gmail.com), if you have questions or run into issues!