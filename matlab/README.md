# MATLAB API

The volume workflow uses MATLAB with SPM12. Add SPM12 and the repository's `matlab` directory to the MATLAB path.

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
bigVsNo = struct('name', 'anticipation_big-no', 'weights', [-1 1]);
PBI = pbi.fit(spmFile, 'pet_template.nii', 'gm_mask.nii', 'pbi_model', ...
    Contrasts=bigVsNo);
```

You can also pass an existing single-session `SPM.mat` directly. It must have one canonical HRF basis function, one unmodulated regressor per condition, and the task effects of interest as its first F contrast. Keep the beta images, `mask.nii`, and original scan paths accessible.

Each contrast has one weight per interaction in source condition order (`SPM.Sess.U`). For an eight-condition model starting with `no_win`, `small_win`, `big_win`, use `[-1 0 1 0 0 0 0 0]` for big versus no anticipation. Supply several contrasts as a struct array, or omit `Contrasts` to estimate the design alone.

- The returned struct and `PBI.mat` contain the mean-centered observed BOLD projection (`phys`), its deconvolved counterpart (`seed_microtime`), condition columns, and condition-specific interactions.
- `seed_mask.nii` records the voxels used for both PET projections.
- The final `SPM.mat` includes task main effects, `seed_phys_bold`, interactions, and any source nuisance columns. SPM writes `con_*.nii` and `spmT_*.nii` for requested interaction contrasts.

See the [model specification](../docs/methods.md) for the interpretation of the interaction coefficients.
