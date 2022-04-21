function [n_rois] = get_number_of_rois(vol)
%function [n_rois] = get_number_of_rois(vol)

vol_hdr = spm_vol(vol);
vol_data = spm_read_vols(vol_hdr);
roi_idc = unique(vol_data);
roi_idc(roi_idc==0) = [];

if ~isequal(roi_idc, round(roi_idc))
    warning('Parcellation volume seems to contain non-integers!')
end

n_rois = length(roi_idc);

end