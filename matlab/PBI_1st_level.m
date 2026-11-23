function PBI_1st_level(SPM, PBI, outputDir, contrasts)
% Estimate the PBI GLM with both main effects and condition interactions.

assert(isfield(PBI, 'phys') && isequal(size(PBI.phys), [SPM.nscan 1]) && ...
    all(isfinite(PBI.phys)) && max(PBI.phys) > min(PBI.phys), ...
    'PBI requires a nonconstant physiological main effect at scan resolution.');

nInteraction = size(PBI.int, 2);
nuisance = SPM.Sess.C;
names = [PBI.intname, PBI.psyname, {'seed_phys_bold'}, nuisance.name(:)'];
design = [PBI.int, PBI.psyhrf, PBI.phys, nuisance.C];

spec.dir = {outputDir};
spec.timing.units = 'secs';
spec.timing.RT = SPM.xY.RT;
spec.timing.fmri_t = SPM.xBF.T;
spec.timing.fmri_t0 = SPM.xBF.T0;
spec.sess.scans = cellstr(SPM.xY.P);
spec.sess.hpf = SPM.xX.K(1).HParam;
spec.global = 'None';
spec.mthresh = 0.8;
spec.cvi = 'AR(1)';
for j = 1:size(design, 2)
    spec.sess.regress(j).name = names{j};
    spec.sess.regress(j).val = design(:, j);
end

matlabbatch{1}.spm.stats.fmri_spec = spec;
for c = 1:numel(contrasts)
    weights = contrasts(c).weights(:)';
    assert(numel(weights) == nInteraction && all(isfinite(weights)) && any(weights ~= 0), ...
        'Contrast "%s" requires one weight per condition and at least one nonzero weight.', ...
        contrasts(c).name);
    % Interactions come first. All main effects, confounds, and constant get zero weight.
    contrast.consess{c}.tcon.name = contrasts(c).name;
    contrast.consess{c}.tcon.weights = [weights, zeros(1, size(design, 2) + 1 - nInteraction)];
    contrast.consess{c}.tcon.sessrep = 'none';
end

spm_jobman('run', matlabbatch);
spmFile = {fullfile(outputDir, 'SPM.mat')};
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = spmFile;
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', matlabbatch);

if ~isempty(contrasts)
    contrast.spmmat = spmFile;
    contrast.delete = 0;
    clear matlabbatch
    matlabbatch{1}.spm.stats.con = contrast;
    spm_jobman('run', matlabbatch);
end
end
