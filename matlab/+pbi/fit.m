function PBI = fit(spmFile, bioMap, grayMatterMask, outputDir, options)
% Run PBI with a volume-defined component and optional surface/grayordinate targets.
% TargetBold accepts paired GIFTI time series (left, right) or one CIFTI dtseries.
% TargetMeshes supplies the corresponding left and right surface geometries.
% Each contrast names its positive and negative conditions. For a weighted
% contrast, supply conditions and their corresponding weights instead.
arguments
    spmFile (1,:) char
    bioMap (1,:) char
    grayMatterMask (1,:) char
    outputDir (1,:) char
    options.Contrasts struct = struct([])
    options.TargetBold string = strings(1,0)
    options.TargetMeshes string = strings(1,0)
end
source = load(spmFile, 'SPM');
SPM = source.SPM;
assert(isscalar(SPM.nscan), 'PBI requires a single-session task GLM.');
assert(isfield(SPM, 'xCon') && ~isempty(SPM.xCon) && strcmp(SPM.xCon(1).STAT, 'F'), ...
    'The first contrast in the source GLM must be the task effects-of-interest F contrast.');
assert(all(arrayfun(@(u) numel(u.name) == 1, SPM.Sess.U)) && size(SPM.xBF.bf, 2) == 1, ...
    'PBI requires one canonical HRF basis function and no parametric or temporal modulators.');

conditionNames = arrayfun(@(u) string(u.name{1}), SPM.Sess.U);
contrasts = resolve_contrasts(options.Contrasts, conditionNames);
outputDir = new_output_dir(outputDir);
[~, sourceFile] = fileattrib(spmFile);

spm('Defaults', 'fMRI');
spm_jobman('initcfg');

% Read the caller's mask/template paths before changing to the source GLM folder.
sourceFolder = fileparts(sourceFile.Name);
[mask, mapValues, geometry] = PBI_create_mask(fullfile(sourceFolder, 'mask.nii'), ...
    grayMatterMask, bioMap, SPM.xY.VY(1));

caller = pwd;
restoreFolder = onCleanup(@() cd(caller));
cd(sourceFolder);
[PBI, valid] = PBI_PEB_decon(SPM, mask, mapValues);
PBI = PBI_build_interaction(PBI, SPM);

% Save the effective mask after excluding locations with nonfinite data.
mask(mask) = valid;
geometry.fname = fullfile(outputDir, 'bio_mask.nii');
geometry.dt = [spm_type('uint8') 0];
geometry.pinfo = [1; 0; 0];
spm_write_vol(geometry, double(mask));

save(fullfile(outputDir, 'PBI.mat'), 'PBI');
if isempty(options.TargetBold)
    assert(isempty(options.TargetMeshes), 'TargetMeshes requires TargetBold.');
    PBI_1st_level(SPM, PBI, outputDir, contrasts, cellstr(SPM.xY.P));
else
    cd(caller); % Resolve target paths relative to the caller, like the biological map.
    fit_surface_targets(SPM, PBI, outputDir, contrasts, options.TargetBold, options.TargetMeshes);
end
end
