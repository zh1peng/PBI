function [neural, bold] = PBI_PEB_decon(SPM, mask)
% Estimate microtime signals with voxelwise SPM PEB and retain observed BOLD.
% Run from the source GLM folder for relative SPM image paths.

nScan = SPM.nscan;
volumes = SPM.xY.VY;
[i, j, k] = ind2sub(size(mask), find(mask));
xyz = [i, j, k]';
bold = spm_data_read(volumes, 'xyz', xyz);
beta = spm_data_read(SPM.Vbeta, 'xyz', xyz);

% Adjust for effects outside the supplied F contrast after SPM filtering.
W = SPM.xX.W;
filtered = spm_filter(SPM.xX.K, W * bold);
adjusted = filtered - spm_FcUtil('Y0', SPM.xCon(1), SPM.xX.xKXs, beta);
nuisance = [SPM.xX.xKXs.X(:, [SPM.xX.iB SPM.xX.iG]), SPM.xX.K(1).X0];
nuisance = nuisance(:, any(nuisance, 1));

bins = SPM.xBF.T;
sample = (0:nScan-1) * bins + SPM.xBF.T0;
padding = 128;
basis = spm_dctmtx(nScan * bins + padding, nScan);
hrf = spm_hrf(SPM.xBF.dt);
response = zeros(nScan, nScan);
for j = 1:nScan
    convolved = conv(basis(:, j), hrf);
    response(:, j) = convolved(padding + sample);
end
basis = basis(padding+1:end, :);

% Preserve the prototype's PEB priors and SPM adjusted signal.
nNuisance = size(nuisance, 2);
model{1}.X = [W * response, nuisance];
model{1}.C = speye(nScan) / 4;
model{2}.X = sparse(nScan + nNuisance, 1);
model{2}.C = blkdiag(speye(nScan) * nScan / sum(response(:).^2), ...
    speye(nNuisance) * 1e6);

neural = nan(size(basis, 1), size(adjusted, 2));
for voxel = find(all(isfinite(adjusted), 1) & all(isfinite(bold), 1))
    estimate = spm_PEB(adjusted(:, voxel), model);
    neural(:, voxel) = spm_detrend(basis * estimate{2}.E(1:nScan));
end
end
