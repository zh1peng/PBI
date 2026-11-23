# MATLAB API

The first-level workflow uses MATLAB with SPM12. Add SPM12 and the repository's `matlab` directory to the MATLAB path.

```matlab
addpath('/path/to/spm12', '/path/to/PBI/matlab');
```

### Estimate the task GLM

`pbi.task_glm` takes one preprocessed 4D BOLD run (`.nii` or `.nii.gz`), an SPM condition struct array, the TR in seconds, and a fresh output directory. Supply one onset per event and either one duration for the condition or one per event, all in seconds. Optional numeric confounds have one row per scan. The function returns the path to the estimated `SPM.mat`.

```matlab
% Two-condition example. Replace the event times and TR with your run's values.
conditions = struct('name', {'no_win', 'big_win'}, ...
    'onset', {[10 42 74], [26 58 90]}, 'duration', {4, 4});
motion = readmatrix('motion.tsv');
spmFile = pbi.task_glm('bold.nii.gz', conditions, 2.2, 'task_model', Confounds=motion);
```

The model uses a canonical HRF, 16 microtime bins with reference bin 1, a 128 s high-pass filter and AR(1) errors. Its first contrast is an omnibus F contrast over the condition columns. SPM retains uncompressed scans at their original path or, for gzip inputs, in the task model directory.

### Fit PET-informed gPPI

`pbi.fit` takes an estimated task model, a PET template, a binary gray matter mask, and a fresh PBI output directory. Use aligned NIfTI volumes on the BOLD grid and a PET template with its reference regions excluded. Reuse the same task model for each PET template, choosing a separate output directory for each fit.

```matlab
bigVsNo = struct('name', 'anticipation_big-no', ...
    'positive', 'big_win', 'negative', 'no_win');
petMap = 'pet_template.nii';
gmMask = 'gm_mask.nii';
PBI = pbi.fit(spmFile, petMap, gmMask, 'pbi_model', Contrasts=bigVsNo);
```

You can also pass an existing single-session `SPM.mat` directly. It must have one canonical HRF basis function, one unmodulated regressor per condition, and the task effects of interest as its first F contrast. Keep the beta images, `mask.nii`, and original scan paths accessible.

Use the exact, case-sensitive condition names from `SPM.Sess.U`. The example estimates the `big_win` interaction coefficient minus the `no_win` interaction coefficient in any column order. Other conditions, main effects and nuisance columns receive zero weight. Supply several contrasts as a struct array, or omit `Contrasts` to estimate the design alone.

For a weighted combination, pair each weight with a named condition. For example, the mean of big and small win versus no win can be written as `struct('name', 'mean_win-no', 'conditions', ["big_win","small_win","no_win"], 'weights', [0.5 0.5 -1])`.

- The returned struct and `PBI.mat` contain the mean-centered observed BOLD projection (`phys`), the microtime seed (`seed_microtime`), condition columns, and condition-specific interactions. For the microtime seed, PBI projects the SPM-adjusted BOLD onto the same PET weights and deconvolves that series once. The shared design, fixed covariances, and zero prior mean make this equivalent to projecting the voxelwise deconvolutions.
- `seed_mask.nii` records the voxels used for both PET projections.
- The final `SPM.mat` includes task main effects, `seed_phys_bold`, interactions, and any source nuisance columns. SPM writes `con_*.nii` and `spmT_*.nii` for requested interaction contrasts.

See the [model specification](../docs/methods.md) for the interpretation of the interaction coefficients.

### Fit surface or grayordinate targets

Use the same source volume GLM, PET map and gray matter mask as above. `TargetBold` changes where PBI estimates coupling. The seed and interactions still come from the source volume run. Supply the same acquisition and retained frames in the same order, with compatible preprocessing. Use baseline-preserving BOLD suitable for SPM's intensity scaling and masking. Include the matching left and right anatomical meshes (for example midthickness); the mesh must have the same registration, vertex count and vertex order as its functional data.

**GIFTI** — pass the left and right time series and their meshes in that order. This example uses fsaverage5.

```matlab
PBI = pbi.fit(spmFile, petMap, gmMask, 'pbi_fsaverage5', Contrasts=bigVsNo, ...
    TargetBold=["lh.bold.func.gii", "rh.bold.func.gii"], ...
    TargetMeshes=["lh.midthickness.surf.gii", "rh.midthickness.surf.gii"]);
```

**CIFTI-2** — add the [Washington University cifti-matlab library](https://github.com/Washington-University/cifti-matlab) to the MATLAB path for CIFTI I/O. Its reader also invokes Connectome Workbench (`wb_command` on the system path) for files whose XML declares version `2.0` instead of `2`, including some NiBabel outputs. This example uses fsLR meshes for an HCP-style dense time series. The CIFTI brain models store vertex indices and counts, and the mesh files supply the triangles.

```matlab
addpath('/path/to/cifti-matlab');
PBI = pbi.fit(spmFile, petMap, gmMask, 'pbi_fslr', Contrasts=bigVsNo, ...
    TargetBold="bold.dtseries.nii", ...
    TargetMeshes=["L.32k_fs_LR.midthickness.surf.gii", "R.32k_fs_LR.midthickness.surf.gii"]);
```

Each hemisphere has a JobMan first-level model under `CORTEX_LEFT/` or `CORTEX_RIGHT/`, with its own `SPM.mat`, `con_0001.gii` and `spmT_0001.gii`. For a CIFTI input with voxels, all volume structures share a third model in `volume/`. The fits use the same task, seed, nuisance and interaction columns and contrast weights; SPM estimates the noise model and applies its intensity scaling within each fit. The prepared scans remain beside each model so that its `SPM.mat` can be reopened.

For each requested CIFTI contrast, `con_0001.dscalar.nii` and `spmT_0001.dscalar.nii` in the output root contain the effect and T statistic on the original dense axis. Locations excluded by SPM's analysis mask have value NaN. These are unthresholded first-level maps; the per-part models have separate search spaces for spatial inference.
