function [PBI, valid] = PBI_build_interaction(neural, bold, pet, SPM)
% Project matched microtime and scan-level data onto the same PET template.

valid = all(isfinite(neural), 1) & all(isfinite(bold), 1);
pet = pet(valid);
assert(numel(pet) >= 2 && max(pet) > min(pet), ...
    'The PET template must vary across the retained voxels.');
centered = pet - mean(pet);
weights = centered / sum(centered.^2);
PBI.seed_microtime = neural(:, valid) * weights;
PBI.phys = spm_detrend(bold(:, valid) * weights);

nScan = SPM.nscan;
sample = ((0:nScan-1) * SPM.xBF.T + SPM.xBF.T0)';
hrf = spm_hrf(SPM.xBF.dt);
conditions = spm_get_ons(SPM, 1);
nCondition = numel(conditions);
PBI.int = zeros(nScan, nCondition);
PBI.psyhrf = zeros(nScan, nCondition);
PBI.psyname = cell(1, nCondition);
PBI.intname = cell(1, nCondition);

for c = 1:nCondition
    condition = conditions(c);
    PBI.psyname{c} = condition.name{1};
    PBI.intname{c} = [condition.name{1} '_int'];

    % SPM pads onset functions with 32 microtime bins.
    psy = spm_detrend(full(condition.u(33:end, 1)));
    interaction = conv(psy .* PBI.seed_microtime, hrf);
    PBI.int(:, c) = spm_detrend(interaction(sample));

    % Use the original task main-effect column at the scan time grid.
    name = ['Sn(1) ' condition.name{1} '*bf(1)'];
    column = strcmp(SPM.xX.name, name);
    assert(nnz(column) == 1, 'Expected one SPM task column named "%s".', name);
    PBI.psyhrf(:, c) = SPM.xX.X(:, column);
end
end
