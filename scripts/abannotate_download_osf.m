function [osf_file] = abannotate_download_osf(osf_id, save_path, zip_file, osf_resource)

if ~exist('zip_file', 'var')
    zip_file = false;
end

if ~exist('osf_resource', 'var')
    osf_resource = 'nvcmf';
end

OSF = 'https://files.osf.io/v1/resources/%s/providers/osfstorage/%s%s';

if zip_file
    osf_file = websave(save_path, sprintf(OSF, osf_resource, osf_id, '/?zip='));
else
    osf_file = websave(save_path, sprintf(OSF, osf_resource, osf_id, ''));
end

end