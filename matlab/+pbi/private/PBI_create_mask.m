function [mask, mapValues, geometry] = PBI_create_mask(boldMask, grayMatterMask, bioMap, boldVolume)
% Extract the template within coverage shared by BOLD, gray matter, and the biological map.

geometry = spm_vol(boldMask);
gm = spm_vol(grayMatterMask);
bio = spm_vol(bioMap);
[aligned, message] = spm_check_orientations([boldVolume; geometry; gm; bio], false);
assert(aligned, '%s', message);

mapVolume = spm_read_vols(bio);
mask = spm_read_vols(geometry) == 1 & spm_read_vols(gm) == 1 & ...
    isfinite(mapVolume) & mapVolume > 1e-8;
mapValues = mapVolume(mask);
assert(numel(mapValues) >= 2, 'The biological map and BOLD masks must contain at least two common voxels.');
end
