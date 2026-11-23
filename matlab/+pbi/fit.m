function PBI = fit(spmFile, petMap, grayMatterMask, outputDir, options)
% Run single-session volume PBI from an estimated SPM12 task GLM.
% Contrasts is an optional struct array with name and weights fields.
% Weights follow the condition order in the source SPM design.
arguments
    spmFile (1,:) char
    petMap (1,:) char
    grayMatterMask (1,:) char
    outputDir (1,:) char
    options.Contrasts struct = struct('name', {}, 'weights', {})
end
source = load(spmFile, 'SPM');
SPM = source.SPM;
assert(isscalar(SPM.nscan), 'PBI requires a single-session task GLM.');
assert(isfield(SPM, 'xCon') && ~isempty(SPM.xCon) && strcmp(SPM.xCon(1).STAT, 'F'), ...
    'The first contrast in the source GLM must be the task effects-of-interest F contrast.');
assert(all(arrayfun(@(u) numel(u.name) == 1, SPM.Sess.U)) && size(SPM.xBF.bf, 2) == 1, ...
    'PBI requires one canonical HRF basis function and no parametric or temporal modulators.');

outputDir = new_output_dir(outputDir);
[~, sourceFile] = fileattrib(spmFile);

spm('Defaults', 'fMRI');
spm_jobman('initcfg');

% Read the caller's mask/template paths before changing to the source GLM folder.
sourceFolder = fileparts(sourceFile.Name);
[mask, pet, geometry] = PBI_create_mask(fullfile(sourceFolder, 'mask.nii'), ...
    grayMatterMask, petMap, SPM.xY.VY(1));

caller = pwd;
restoreFolder = onCleanup(@() cd(caller));
cd(sourceFolder);
[neural, bold] = PBI_PEB_decon(SPM, mask);
[PBI, valid] = PBI_build_interaction(neural, bold, pet, SPM);

% Save the effective mask after excluding locations with nonfinite data.
mask(mask) = valid;
geometry.fname = fullfile(outputDir, 'seed_mask.nii');
geometry.dt = [spm_type('uint8') 0];
geometry.pinfo = [1; 0; 0];
spm_write_vol(geometry, double(mask));

save(fullfile(outputDir, 'PBI.mat'), 'PBI');
PBI_1st_level(SPM, PBI, outputDir, options.Contrasts);
end
