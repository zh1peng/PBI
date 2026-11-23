function [mask, pet, geometry] = PBI_create_mask(boldMask, grayMatterMask, petMap, boldVolume)
% Extract the template within coverage shared by BOLD, gray matter, and PET.

geometry = spm_vol(boldMask);
gm = spm_vol(grayMatterMask);
bio = spm_vol(petMap);
[aligned, message] = spm_check_orientations([boldVolume; geometry; gm; bio], false);
assert(aligned, '%s', message);

petVolume = spm_read_vols(bio);
mask = spm_read_vols(geometry) == 1 & spm_read_vols(gm) == 1 & ...
    isfinite(petVolume) & petVolume > 1e-8;
pet = petVolume(mask);
assert(numel(pet) >= 2, 'The PET and BOLD masks must contain at least two common voxels.');
end
