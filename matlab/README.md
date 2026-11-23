# MATLAB implementation

The single-session volume PBI workflow requires MATLAB and SPM12. Start with an estimated first-level task GLM using a single canonical HRF basis function and no parametric or temporal modulators. Set its first contrast to an omnibus **F contrast covering the task effects of interest**, so the deconvolution retains those task effects. The source folder must contain its `SPM.mat`, `mask.nii`, and beta images. The scan paths in the source model must resolve on the current machine.

Use a binary gray matter mask and a PET template already aligned to the BOLD voxel grid. Exclude PET reference regions and invalid template coverage during template preparation. Use a fresh output directory.

```matlab
addpath('/path/to/spm12');
addpath('/path/to/PBI/matlab');

% Example for an eight-condition MID design with no_win, small_win, big_win first.
contrast = struct('name', 'anticipation_big-no', 'weights', [-1 0 1 0 0 0 0 0]);
PBI = PBI_run('/path/to/source_glm/SPM.mat', '/path/to/gm_mask.nii', ...
    '/path/to/pet_template.nii', '/path/to/pbi_glm', contrast);
```

Each contrast needs one weight per source condition, in the order of `SPM.Sess.U`. Pass additional contrasts as a struct array. Omit the argument to estimate the model without generating contrasts.

- `PBI.mat` saves the scan-level physiological regressor (`phys`), microtime seed (`seed_microtime`), and task and interaction columns.
- `seed_mask.nii` records the voxels used when projecting the seed. The same PET weights are applied to observed BOLD and to the PEB-deconvolved signals.
- `SPM.mat` records the final design, including the `seed_phys_bold` main effect. SPM writes contrast and statistic images when contrasts are supplied.

The physiological main effect is the mean-centered projection of observed BOLD. SPM applies the target GLM's filtering and noise model to the design. Interaction terms use the same spatial projection of the microtime signals, multiplied by each centered condition before HRF convolution. See the [SPM PPI documentation](https://www.fil.ion.ucl.ac.uk/spm/docs/manual/ppi/ppi/) and the [shared method specification](../docs/methods.md).

Run the synthetic regression test from the repository root with SPM12 on the MATLAB path:

```matlab
addpath('matlab', 'matlab/tests');
test_PBI_seed_effect;
```

The test estimates an SPM task GLM and runs the PBI entry point on generated volume data. It verifies the saved seed column and checks that the interaction contrast distinguishes common coupling from a known slope difference.
