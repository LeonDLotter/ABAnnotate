function import_psychEncode_cellmarkers(tpm_file, umi_file, save_path) 
%function import_psychEncode_cellmarkers(tpm_file, umi_file, save_path) 
%
% Function to import PsychEncode single cell markers from 
% Wang et al., 2018, Science (https://doi.org/10.1126/science.aat8464)
% Sources:
% PsychEncode Resrouces: http://resource.psychencode.org/
% TPM file: http://resource.psychencode.org/Datasets/Derived/SC_Decomp/DER-20_Single_cell_expression_processed_TPM.tsv
% UMI file: http://resource.psychencode.org/Datasets/Derived/SC_Decomp/DER-21_Single_cell_markergenes_UMI.xlsx

% read PsychEncode data
psyEnc_tpm = readtable(tpm_file);
psyEnc_umi = readtable(umi_file);

% get clusters
tpm_unique_clusters = unique(psyEnc_tpm.cluster, 'stable');
umi_unique_clusters = unique(psyEnc_umi.cluster, 'stable');

% TPM
tpmTable = table();
for i=1:length(tpm_unique_clusters)
    tpmTable.cLabel(i) = tpm_unique_clusters(i);
    tpmTable.cID(i) = i;
    tpmTable.cGenes(i) = {psyEnc_tpm.gene(strcmp(psyEnc_tpm.cluster, tpm_unique_clusters(i)))};
    tpmTable.cSize(i) = length(tpmTable.cGenes{i});
end
tpmTable = movevars(tpmTable, 'cSize', 'Before', 'cGenes');

% UMI
umiTable = table();
for i=1:length(umi_unique_clusters)
    umiTable.cLabel(i) = umi_unique_clusters(i);
    umiTable.cID(i) = i;
    umiTable.cGenes(i) = {psyEnc_umi.gene(strcmp(psyEnc_umi.cluster, umi_unique_clusters(i)))};
    umiTable.cSize(i) = length(umiTable.cGenes{i});
end
umiTable = movevars(umiTable, 'cSize', 'Before', 'cGenes');

cTable = tpmTable;
save(fullfile(save_path, 'PsychEncode-cellTypesTPM-discrete.mat'), 'cTable');
cTable = umiTable;
save(fullfile(save_path, 'PsychEncode-cellTypesUMI-discrete.mat'), 'cTable');

end
