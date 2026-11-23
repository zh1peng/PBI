function test_PBI_seed_effect
% Exercise the public entry point and verify the saved physiological column.
spm('Defaults', 'fMRI');
spm_jobman('initcfg');
rng(12);

root = tempname;
mkdir(root);
originalDir = pwd;
cleanup = onCleanup(@() cleanFixture(root, originalDir));
source = fullfile(root, 'source');
mkdir(source);

nScan = 72;
tr = 2;
dims = [4 4 4];
header = struct('fname', fullfile(root, 'bold.nii'), 'dim', dims, ...
    'dt', [spm_type('float32') 0], 'mat', diag([3 3 3 1]), 'pinfo', [1; 0; 0]);
pet = reshape(linspace(0.2, 1, prod(dims)), dims);
common = randn(nScan, 1);
motion = randn(nScan, 1);
scans = cell(nScan, 1);
for t = 1:nScan
    header.n = [t 1];
    spm_write_vol(header, 100 + 5 * motion(t) + pet * common(t) + 0.2 * randn(dims));
    scans{t} = sprintf('%s,%d', header.fname, t);
end
header.n = [1 1];
header.fname = fullfile(root, 'pet.nii');
spm_write_vol(header, pet);
header.fname = fullfile(root, 'gm.nii');
spm_write_vol(header, ones(dims));

spec.dir = {source};
spec.timing = struct('units', 'secs', 'RT', tr, 'fmri_t', 8, 'fmri_t0', 3);
spec.sess.scans = scans;
spec.sess.cond = struct('name', {'no_win', 'big_win'}, ...
    'onset', {10:24:130, 18:24:138}, 'duration', {3, 3});
spec.sess.regress = struct('name', 'motion', 'val', motion);
spec.cvi = 'AR(1)';
matlabbatch{1}.spm.stats.fmri_spec = spec;
spm_jobman('run', matlabbatch);

sourceSPM = fullfile(source, 'SPM.mat');
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_est.spmmat = {sourceSPM};
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run', matlabbatch);

clear matlabbatch
matlabbatch{1}.spm.stats.con.spmmat = {sourceSPM};
matlabbatch{1}.spm.stats.con.consess{1}.fcon = ...
    struct('name', 'Task effects', 'weights', eye(2), 'sessrep', 'none');
spm_jobman('run', matlabbatch);

out = fullfile(root, 'pbi');
PBI = PBI_run(sourceSPM, fullfile(root, 'gm.nii'), fullfile(root, 'pet.nii'), out, ...
    struct('name', 'big-no', 'weights', [-1 1]));
result = load(fullfile(out, 'SPM.mat'), 'SPM');
model = result.SPM;
assert(isequal(model.xBF.T, 8) && isequal(model.xBF.T0, 3));
seedCol = find(strcmp(model.xX.name, 'Sn(1) seed_phys_bold'));
assert(isscalar(seedCol), 'The design must contain exactly one physiological main effect.');

% Independently regress stored BOLD on the PET template over the effective mask.
mask = spm_read_vols(spm_vol(fullfile(out, 'seed_mask.nii'))) > 0;
input = load(sourceSPM, 'SPM');
% Use SPM's scan headers so the projection includes its intensity scaling.
bold = reshape(spm_read_vols(input.SPM.xY.VY), [], nScan);
petDisk = spm_read_vols(spm_vol(fullfile(root, 'pet.nii')));
coef = [ones(nnz(mask), 1), petDisk(mask)] \ bold(mask, :);
expectedSeed = coef(2, :)' - mean(coef(2, :));
assert(max(abs(PBI.phys - expectedSeed)) < 1e-10);
assert(max(abs(model.xX.X(:, seedCol) - expectedSeed)) < 1e-10);
assert(isequal(model.xCon(1).c(1:2)', [-1 1]) && all(model.xCon(1).c(3:end) == 0));

% Check the contrast can distinguish common coupling from task modulation.
X = model.xX.xKXs.X;
c = model.xCon(1).c;
y = 2 * X(:, seedCol);
interaction = 0.6 * X(:, 1) - 0.4 * X(:, 2);
assert(abs(c' * (X \ y)) < 1e-10);
assert(abs(c' * (X \ (y + interaction)) + 1) < 1e-10);
disp('PBI seed main effect regression passed.');
end

function cleanFixture(root, originalDir)
cd(originalDir);
rmdir(root, 's');
end
