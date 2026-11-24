function PBI_1st_level(SPM, PBI, outputDir, contrasts, scans)
% Estimate the PBI GLM with both main effects and condition interactions.

assert(isfield(PBI, 'Bio') && isequal(size(PBI.Bio), [SPM.nscan 1]) && ...
    all(isfinite(PBI.Bio)) && max(PBI.Bio) > min(PBI.Bio), ...
    'PBI requires a nonconstant physiological main effect at scan resolution.');

nInteraction = size(PBI.int, 2);
nuisance = SPM.Sess.C;
names = [PBI.intname, PBI.psyname, {'Bio'}, nuisance.name(:)'];
design = [PBI.int, PBI.psyhrf, PBI.Bio, nuisance.C];

spec.dir = {outputDir};
spec.timing.units = 'secs';
spec.timing.RT = SPM.xY.RT;
spec.timing.fmri_t = SPM.xBF.T;
spec.timing.fmri_t0 = SPM.xBF.T0;
spec.sess.scans = scans;
spec.sess.hpf = SPM.xX.K(1).HParam;
spec.global = 'None';
spec.mthresh = 0.8;
spec.cvi = 'AR(1)';
for j = 1:size(design, 2)
    spec.sess.regress(j).name = names{j};
    spec.sess.regress(j).val = design(:, j);
end

matlabbatch{1}.spm.stats.fmri_spec = spec;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', ...
    substruct('.','val','{}',{1},'.','val','{}',{1},'.','val','{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

for c = 1:numel(contrasts)
    weights = contrasts(c).weights;
    % Interactions come first. All main effects, confounds, and constant get zero weight.
    contrast.consess{c}.tcon.name = contrasts(c).name;
    contrast.consess{c}.tcon.weights = [weights, zeros(1, size(design, 2) + 1 - nInteraction)];
    contrast.consess{c}.tcon.sessrep = 'none';
end

if ~isempty(contrasts)
    contrast.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', ...
        substruct('.','val','{}',{2},'.','val','{}',{1},'.','val','{}',{1}), substruct('.','spmmat'));
    contrast.delete = 0;
    matlabbatch{3}.spm.stats.con = contrast;
end
spm_jobman('run', matlabbatch);
end
