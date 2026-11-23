function folder = new_output_dir(folder)
% SPM jobs must write to a fresh directory, even when the caller changes cwd.
folder = spm_select('CPath', folder);
if ~isfolder(folder)
    mkdir(folder);
end
listing = dir(folder);
assert(all(ismember({listing.name}, {'.', '..'})), 'Choose an empty output directory: %s', folder);
end
