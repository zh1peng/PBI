function spmFile = task_glm(boldFile, conditions, tr, outputDir, options)
% Estimate a single-session SPM task GLM for pbi.fit.
% conditions has name, onset (seconds) and duration (seconds) fields.
% Confounds is a numeric matrix with one row per scan.
arguments
    boldFile (1,:) char
    conditions struct
    tr (1,1) double {mustBePositive}
    outputDir (1,:) char
    options.Confounds (:,:) double = []
end

spm('Defaults', 'fMRI');
spm_jobman('initcfg');
outputDir = new_output_dir(outputDir);
boldFile = spm_select('CPath', boldFile);
if endsWith(boldFile, '.nii.gz')
    files = gunzip(boldFile, outputDir);
    boldFile = files{1};
end

spec.dir = {outputDir};
spec.timing = struct('units', 'secs', 'RT', tr, 'fmri_t', 16, 'fmri_t0', 1);
spec.sess.scans = cellstr(spm_select('expand', boldFile));
spec.sess.cond = conditions;
spec.sess.hpf = 128;
spec.bases.hrf.derivs = [0 0];
spec.cvi = 'AR(1)';
for j = 1:size(options.Confounds, 2)
    spec.sess.regress(j).name = sprintf('nuisance_%d', j);
    spec.sess.regress(j).val = options.Confounds(:, j);
end
batch{1}.spm.stats.fmri_spec = spec;
spmFile = fullfile(outputDir, 'SPM.mat');
batch{2}.spm.stats.fmri_est.spmmat = {spmFile};
batch{2}.spm.stats.fmri_est.method.Classical = 1;
batch{3}.spm.stats.con.spmmat = {spmFile};
batch{3}.spm.stats.con.consess{1}.fcon = struct('name', 'Task effects', ...
    'weights', eye(numel(conditions)), 'sessrep', 'none');
spm_jobman('run', batch);
end
